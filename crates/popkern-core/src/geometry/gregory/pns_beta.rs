//! PnS cap — β reparameterization for the n-gon cap.
//!
//! β is a change of variables applied to the **tensor-border** `b`
//! of the surrounding bi-2 spline complex. It reparameterises the
//! Hermite (position + transversal derivative) data along the cap's
//! outer boundary to:
//! 1. Speed up the transition from the tensor-structure to the
//!    arrangement of parameter lines at the cap's centre.
//! 2. Enable the G¹ construction (eq 2) to have a consistent
//!    solution via eq (3).
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` — Karčiauskas & Peters, "Smooth
//! multi-sided blending of biquadratic splines" (CAGD 2015),
//! §4 eq (3) + Appendix B for the **bi-4** construction (n ≥ 5);
//! §4.2 eq (3³) for the **bi-3** construction (n ∈ {3, 5}).
//!
//! # Shape of β
//!
//! Along one sector boundary parameterised by `u`, β acts on
//! `(u, v)` as:
//! ```text
//!   β(u, v) = (u + b(u) · v,  a(u) · v)
//! ```
//!
//! where `a(u)` and `b(u)` are polynomials in `u` given below. The
//! `v` direction is transversal to the boundary (pointing into the
//! cap).
//!
//! # Variants
//!
//! - **Bi-4** (general n, §4 eq 3): `a` is degree 2, `b` is degree 3.
//!   Paper defaults: `a_2^def = (n - 2)/(0.55n - 0.2)`, `b_2^def = 0`
//!   (Appendix B motivates these via σ's internal scaling heuristic).
//! - **Bi-3 n ∈ {3, 5}** (§4.2 eq 3³): closed-form
//!   `a(u) = (1 − c(1 − u))/(1 − c)`, `b(u) = c(1 − u)u / (1 − c)`
//!   where `c = cos(2π/n)`. No Appendix B knobs — fully specified.
//!
//! The bi-3 variant is strictly simpler and stable; the bi-4
//! variant's `a_2`, `b_2` parameters exist as empirical degrees of
//! freedom. We ship paper defaults, gated on `n`.

#![allow(dead_code)]

use popkern_math::PortableMath;

/// β reparameterization for one sector boundary. Evaluates `(β_x,
/// β_y)` given `(u, v)` ∈ `[0, 1]²`.
#[derive(Debug, Clone, Copy)]
pub struct BetaReparam {
    /// Polynomial `a(u) = Σ a_i B_i^2(u)` (degree 2 in u).
    pub a: [f64; 3],
    /// Polynomial `b(u) = Σ b_i B_i^3(u)` (degree 3 in u).
    pub b: [f64; 4],
}

impl BetaReparam {
    /// **Bi-4 variant** — for n ≥ 5 (and the `n ∈ {5, 6, 7, 8}`
    /// range where the bi-4 cap is the default). Uses paper's
    /// Appendix B defaults for the empirical `a_2`, `b_2`.
    ///
    /// `a_0 = 1,  a_1 = (2 − c)/2 · a_2,  a_2^def = (n − 2)/(0.55n − 0.2)`
    /// `b_0 = 0,  b_1 = −2/3 + (2 − c)/3 · a_2,  b_2^def = 0,  b_3 = 0`
    ///
    /// where `c := cos(2π/n)`. Panics (debug-only) for `n < 3` —
    /// the formula is ill-defined at `n = 2` (sector angle π).
    pub fn bi4(n: usize) -> Self {
        debug_assert!(n >= 3, "β bi-4: n must be ≥ 3, got {n}");
        let c = (2.0 * std::f64::consts::PI / (n as f64)).pcos();
        let a_2 = (n as f64 - 2.0) / (0.55 * (n as f64) - 0.2);
        let a_1 = (2.0 - c) / 2.0 * a_2;
        let a_0 = 1.0;
        let b_2 = 0.0;
        let b_3 = 0.0;
        let b_1 = -2.0 / 3.0 + (2.0 - c) / 3.0 * a_2;
        let b_0 = 0.0;
        BetaReparam {
            a: [a_0, a_1, a_2],
            b: [b_0, b_1, b_2, b_3],
        }
    }

    /// **Bi-3 variant** — for n ∈ {3, 5} where the cap is one bi-3
    /// patch per sector (§4.2). Closed-form β: no empirical knobs.
    ///
    /// `a(u) = (1 − c(1 − u))/(1 − c)`
    /// `b(u) = c(1 − u)u / (1 − c)`
    ///
    /// where `c := cos(2π/n)`. Sampled at Bézier control positions
    /// `u ∈ {0, 1/2, 1}` for `a` (degree 2) and `u ∈ {0, 1/3, 2/3,
    /// 1}` for `b` (degree 3) and converted to BB form.
    ///
    /// For higher valencies `(1 − c) → 0` producing very non-uniform
    /// parameter lines; the paper notes this variant is only safe
    /// for `n ∈ {3, 5}`.
    pub fn bi3_for_n3_n5(n: usize) -> Self {
        debug_assert!(
            n == 3 || n == 5,
            "β bi-3: only n ∈ {{3, 5}} supported, got {n}"
        );
        let c = (2.0 * std::f64::consts::PI / (n as f64)).pcos();
        let inv = 1.0 / (1.0 - c);

        // a(u) = (1 − c·(1 − u))·inv, a polynomial of degree 1 in u:
        //   a(u) = (1 − c)·inv + c·inv · u = 1 + c·inv · u − c·inv
        // At degree 2 in BB form (with monomial → Bernstein conversion):
        // a(0) = a_0 = (1 − c)·inv = 1
        // a(1) = a_2 = (1 − c·0)·inv = 1/(1 − c)  (note: not +1 because a is polynomial
        // of degree ≤ 2 and we sample at u=1)
        // Actually a(u) is linear so a_0 = a(0), a_2 = a(1), a_1 = ((a(0) + a(1))/2 is
        // NOT right for degree-2 BB — need degree elevation.
        //
        // Convert degree-1 polynomial `p(u) = α + β·u` to degree-2 BB by identity:
        //   [α, α + β/2, α + β]  (Bernstein CPs of a linear function at degree 2).
        // Here α = a(0) = 1, β = a'(u) = c·inv.
        let a_linear_alpha = 1.0;
        let a_linear_beta = c * inv;
        let a_0 = a_linear_alpha;
        let a_1 = a_linear_alpha + a_linear_beta / 2.0;
        let a_2 = a_linear_alpha + a_linear_beta;

        // b(u) = c·(1 − u)·u · inv, a polynomial of degree 2:
        //   b(u) = c·inv·u − c·inv·u². Monomial coeffs (0, c·inv, −c·inv).
        // Convert degree-2 monomial → degree-3 BB:
        //   p(u) = c_0 + c_1 u + c_2 u²
        // at degree 3: BB coefs b_i via
        //   b_0 = c_0
        //   b_1 = c_0 + c_1/3
        //   b_2 = c_0 + 2c_1/3 + c_2/3
        //   b_3 = c_0 + c_1 + c_2
        let c0 = 0.0;
        let c1 = c * inv;
        let c2 = -c * inv;
        let b_0 = c0;
        let b_1 = c0 + c1 / 3.0;
        let b_2 = c0 + 2.0 * c1 / 3.0 + c2 / 3.0;
        let b_3 = c0 + c1 + c2;
        BetaReparam {
            a: [a_0, a_1, a_2],
            b: [b_0, b_1, b_2, b_3],
        }
    }

    /// Evaluate `a(u) = Σ a_i B_i^2(u)` at parameter `u`.
    pub fn eval_a(&self, u: f64) -> f64 {
        let omu = 1.0 - u;
        self.a[0] * omu * omu + self.a[1] * 2.0 * u * omu + self.a[2] * u * u
    }

    /// Evaluate `b(u) = Σ b_i B_i^3(u)` at parameter `u`.
    pub fn eval_b(&self, u: f64) -> f64 {
        let omu = 1.0 - u;
        self.b[0] * omu.powi(3)
            + self.b[1] * 3.0 * u * omu * omu
            + self.b[2] * 3.0 * u * u * omu
            + self.b[3] * u.powi(3)
    }

    /// Evaluate β at `(u, v)`. Returns `(β_x, β_y)` = `(u + b(u)·v, a(u)·v)`.
    pub fn eval(&self, u: f64, v: f64) -> (f64, f64) {
        (u + self.eval_b(u) * v, self.eval_a(u) * v)
    }
}

#[cfg(test)]
#[path = "pns_beta_tests.rs"]
mod tests;
