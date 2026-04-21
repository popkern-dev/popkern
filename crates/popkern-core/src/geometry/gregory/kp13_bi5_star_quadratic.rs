//! Phase I.11c — Star-point quadratic identity per [KP09a §2.3(a)].
//!
//! # Paper statement
//!
//! > "The extraordinary point `g^i_{00}` and its five nearest
//! > coefficients `g^i_{10}, g^i_{01}, g^i_{20}, g^i_{11}, g^i_{02}`
//! > define the same quadratic polynomial (in degree-raised BB-
//! > form) for all `i`."
//!
//! That is, across **every** sector, the 6 near-star tri-5 BB-
//! coefficients are forced to be consistent with a single bivariate
//! quadratic `Q(x, y)` on the full n-gon domain Ω. So the whole
//! near-star block has 6 scalar DoFs globally, not `6n`.
//!
//! # Representation
//!
//! `Q(x, y) := α·x² + β·x·y + γ·y² + δ·x + ε·y + ζ`
//!
//! Each coefficient is `Point3`-valued — the quadratic is applied
//! independently per output axis.
//!
//! # Polar (blossom) form for degree-2 polynomials
//!
//! For two planar points `(x_1, y_1), (x_2, y_2)`:
//!
//! ```text
//!   B((x₁, y₁), (x₂, y₂)) := α·x₁·x₂
//!                          + β·(x₁·y₂ + y₁·x₂) / 2
//!                          + γ·y₁·y₂
//!                          + δ·(x₁ + x₂) / 2
//!                          + ε·(y₁ + y₂) / 2
//!                          + ζ
//! ```
//!
//! `B(P, P) = Q(P)`. The tri-2 BB-coef on a wedge with vertices
//! `(v_a, v_b, o)` at barycentric `(a, b, c)` with `a + b + c = 2`
//! is `B` evaluated at the multiset with `a` copies of `v_a`, `b`
//! copies of `v_b`, `c` copies of `o`.
//!
//! # Degree-raise to tri-5
//!
//! For a tri-2 polynomial `Q` represented by 6 coefficients
//! `q_{a,b,c}`, the standard degree-raise to tri-5 gives 21
//! tri-5 coefficients `b_{i,j,k}^{(5)}`. Only 6 of those sit in
//! the near-star region (`i, j ≤ 2, i + j ≤ 2`). The relevant
//! formulas (`C(m, d)` = binomial):
//!
//! ```text
//!   b^{(5)}_{i, j, k} = (1 / C(5, 2)) · Σ_{a+b+c=2, a≤i, b≤j, c≤k}
//!                          C(i, a)·C(j, b)·C(k, c) · q_{a, b, c}
//! ```
//!
//! Expanded for the 6 near-star coefs:
//!
//! | tri-5 coef    | formula (÷ 10)                               |
//! |---------------|----------------------------------------------|
//! | `b_{0,0,5}`   | `q_{0,0,2} · 10`                             |
//! | `b_{0,1,4}`   | `6·q_{0,0,2} + 4·q_{0,1,1}`                  |
//! | `b_{0,2,3}`   | `3·q_{0,0,2} + 6·q_{0,1,1} + q_{0,2,0}`      |
//! | `b_{1,0,4}`   | `6·q_{0,0,2} + 4·q_{1,0,1}`                  |
//! | `b_{1,1,3}`   | `3·q_{0,0,2} + 3·q_{1,0,1} + 3·q_{0,1,1} + q_{1,1,0}` |
//! | `b_{2,0,3}`   | `3·q_{0,0,2} + 6·q_{1,0,1} + q_{2,0,0}`      |
//!
//! The 6×6 map from `q_{a,b,c}` to near-star `b^{(5)}` is upper-
//! triangular in an appropriate ordering, hence invertible.
//!
//! # Near-star flat indices
//!
//! In the tri-5 flat layout, the 6 near-star CPs sit at flat
//! indices `[0, 1, 2, 6, 7, 11]`:
//!
//! | tri-5 coef | (i, j) | flat | paper label |
//! |------------|--------|-----:|-------------|
//! | `b_{0,0,5}` | `(0, 0)` |  0 | `500`      |
//! | `b_{0,1,4}` | `(0, 1)` |  1 | `401`      |
//! | `b_{0,2,3}` | `(0, 2)` |  2 | `302`      |
//! | `b_{1,0,4}` | `(1, 0)` |  6 | `410`      |
//! | `b_{1,1,3}` | `(1, 1)` |  7 | `311`      |
//! | `b_{2,0,3}` | `(2, 0)` | 11 | `320`      |
//!
//! (Paper labels reflect KP13's `(o-exp, m_k-exp, m_{k+1}-exp)`
//! triple — star is "500".)

#![allow(dead_code)]

use super::kp13_bi5_gtri5::{Gtri5, sector_corner};
use crate::geometry::Point3;

/// Flat indices of the 6 near-star tri-5 CPs.
/// Order aligns with `tri2_coefs`:
/// `[q_{0,0,2}, q_{0,1,1}, q_{0,2,0}, q_{1,0,1}, q_{1,1,0}, q_{2,0,0}]`.
pub const NEAR_STAR_FLAT: [usize; 6] = [0, 1, 2, 6, 7, 11];

/// Bivariate quadratic `Q(x, y) = α·x² + β·x·y + γ·y² + δ·x + ε·y + ζ`
/// with `Point3`-valued coefficients.
#[derive(Debug, Clone, Copy)]
pub struct StarQuadratic {
    pub alpha: Point3,
    pub beta: Point3,
    pub gamma: Point3,
    pub delta: Point3,
    pub epsilon: Point3,
    pub zeta: Point3,
}

impl StarQuadratic {
    /// All-zero quadratic (identity-valued Q).
    pub fn zero() -> Self {
        let z = Point3::new(0.0, 0.0, 0.0);
        Self {
            alpha: z,
            beta: z,
            gamma: z,
            delta: z,
            epsilon: z,
            zeta: z,
        }
    }

    /// Pure constant `Q(x, y) = c`.
    pub fn constant(c: Point3) -> Self {
        Self {
            zeta: c,
            ..Self::zero()
        }
    }

    /// Evaluate `Q(x, y)`.
    pub fn evaluate(&self, x: f64, y: f64) -> Point3 {
        let x2 = x * x;
        let xy = x * y;
        let y2 = y * y;
        Point3::new(
            self.alpha.x * x2
                + self.beta.x * xy
                + self.gamma.x * y2
                + self.delta.x * x
                + self.epsilon.x * y
                + self.zeta.x,
            self.alpha.y * x2
                + self.beta.y * xy
                + self.gamma.y * y2
                + self.delta.y * x
                + self.epsilon.y * y
                + self.zeta.y,
            self.alpha.z * x2
                + self.beta.z * xy
                + self.gamma.z * y2
                + self.delta.z * x
                + self.epsilon.z * y
                + self.zeta.z,
        )
    }

    /// Polar (blossom) form `B((x₁, y₁), (x₂, y₂))`. Invariant
    /// `B(P, P) = Q(P)`.
    pub fn polar(&self, p1: (f64, f64), p2: (f64, f64)) -> Point3 {
        let (x1, y1) = p1;
        let (x2, y2) = p2;
        Point3::new(
            self.alpha.x * x1 * x2
                + self.beta.x * (x1 * y2 + y1 * x2) * 0.5
                + self.gamma.x * y1 * y2
                + self.delta.x * (x1 + x2) * 0.5
                + self.epsilon.x * (y1 + y2) * 0.5
                + self.zeta.x,
            self.alpha.y * x1 * x2
                + self.beta.y * (x1 * y2 + y1 * x2) * 0.5
                + self.gamma.y * y1 * y2
                + self.delta.y * (x1 + x2) * 0.5
                + self.epsilon.y * (y1 + y2) * 0.5
                + self.zeta.y,
            self.alpha.z * x1 * x2
                + self.beta.z * (x1 * y2 + y1 * x2) * 0.5
                + self.gamma.z * y1 * y2
                + self.delta.z * (x1 + x2) * 0.5
                + self.epsilon.z * (y1 + y2) * 0.5
                + self.zeta.z,
        )
    }

    /// Tri-2 BB-coefs on the wedge `(o, v_a, v_b)` (o = origin).
    /// Returns `[q_{0,0,2}, q_{0,1,1}, q_{0,2,0}, q_{1,0,1}, q_{1,1,0}, q_{2,0,0}]`.
    pub fn tri2_bb_coefs_on_wedge(&self, v_a: (f64, f64), v_b: (f64, f64)) -> [Point3; 6] {
        let o = (0.0, 0.0);
        [
            self.polar(o, o),     // q_{0,0,2} = Q(o) = ζ
            self.polar(v_b, o),   // q_{0,1,1}
            self.polar(v_b, v_b), // q_{0,2,0} = Q(v_b)
            self.polar(v_a, o),   // q_{1,0,1}
            self.polar(v_a, v_b), // q_{1,1,0}
            self.polar(v_a, v_a), // q_{2,0,0} = Q(v_a)
        ]
    }

    /// Near-star tri-5 BB-coefs on the wedge `(o, v_a, v_b)`.
    /// Returns 6 coefs in the `NEAR_STAR_FLAT` order.
    pub fn tri5_near_star_on_wedge(&self, v_a: (f64, f64), v_b: (f64, f64)) -> [Point3; 6] {
        let q = self.tri2_bb_coefs_on_wedge(v_a, v_b);
        tri2_to_tri5_near_star(&q)
    }

    /// Populate sector `k`'s 6 near-star tri-5 coefs from this
    /// quadratic. Uses sector `k`'s wedge `(o, v_k, v_{k+1})`.
    ///
    /// **Note**: writes `Q(o) = ζ` into the sector's star slot
    /// (flat 0). If the caller pinned the star to the Catmull-Clark
    /// limit via [`Gtri5::set_star_point`], that value is silently
    /// overwritten — prefer [`Self::apply_to_all_sectors_pinned`]
    /// which enforces `Q.ζ = star`.
    pub fn apply_to_sector(&self, gtri5: &mut Gtri5, k: usize) {
        let v_a = sector_corner(gtri5.n, k);
        let v_b = sector_corner(gtri5.n, k + 1);
        let near = self.tri5_near_star_on_wedge(v_a, v_b);
        let sector = gtri5.sector_mut(k);
        for (slot, &flat) in NEAR_STAR_FLAT.iter().enumerate() {
            sector.cps[flat] = near[slot];
        }
    }

    /// Populate every sector's 6 near-star tri-5 coefs from this
    /// quadratic.
    ///
    /// See the caveat on [`Self::apply_to_sector`] — if the star
    /// was pre-pinned, prefer [`Self::apply_to_all_sectors_pinned`].
    pub fn apply_to_all_sectors(&self, gtri5: &mut Gtri5) {
        for k in 0..gtri5.n {
            self.apply_to_sector(gtri5, k);
        }
    }

    /// Populate every sector's 6 near-star tri-5 coefs, forcing
    /// `Q.ζ = pinned_star`.
    ///
    /// Under this path the star coef written is guaranteed to
    /// equal the pinned star — KP13 §3 "star = Catmull-Clark
    /// limit" and KP09a §2.3(a) "same quadratic Q across sectors"
    /// both hold simultaneously.
    pub fn apply_to_all_sectors_pinned(&self, gtri5: &mut Gtri5, pinned_star: Point3) {
        let q = Self {
            zeta: pinned_star,
            ..*self
        };
        q.apply_to_all_sectors(gtri5);
    }

    /// Recover the unique `StarQuadratic` whose sector-0 near-star
    /// tri-5 coefs equal the given 6 values. The 6 are passed in
    /// `NEAR_STAR_FLAT` order.
    ///
    /// Panics if `n < 3` (no valid wedge orientation) or if the
    /// recovery is singular (shouldn't happen for `n ≥ 3`).
    pub fn from_sector0_near_star(near_star: &[Point3; 6], n: usize) -> Self {
        assert!(n >= 3, "from_sector0_near_star needs n ≥ 3; got {n}");
        let q = tri5_near_star_to_tri2(near_star);
        tri2_coefs_to_quadratic(&q, n)
    }
}

/// Degree-raise tri-2 coefs `[q_{0,0,2}, q_{0,1,1}, q_{0,2,0},
/// q_{1,0,1}, q_{1,1,0}, q_{2,0,0}]` to the 6 near-star tri-5
/// coefs in `NEAR_STAR_FLAT` order.
pub fn tri2_to_tri5_near_star(q: &[Point3; 6]) -> [Point3; 6] {
    // Alias for readability.
    let q002 = q[0];
    let q011 = q[1];
    let q020 = q[2];
    let q101 = q[3];
    let q110 = q[4];
    let q200 = q[5];

    let b0_0_5 = q002;
    let b0_1_4 = lin2(6.0, q002, 4.0, q011, 10.0);
    let b0_2_3 = lin3(3.0, q002, 6.0, q011, 1.0, q020, 10.0);
    let b1_0_4 = lin2(6.0, q002, 4.0, q101, 10.0);
    let b1_1_3 = lin_n(&[(3.0, q002), (3.0, q101), (3.0, q011), (1.0, q110)], 10.0);
    let b2_0_3 = lin3(3.0, q002, 6.0, q101, 1.0, q200, 10.0);
    [b0_0_5, b0_1_4, b0_2_3, b1_0_4, b1_1_3, b2_0_3]
}

/// Inverse of [`tri2_to_tri5_near_star`]: given 6 tri-5 near-star
/// coefs, recover the 6 tri-2 coefs whose degree-raise produced them.
pub fn tri5_near_star_to_tri2(b: &[Point3; 6]) -> [Point3; 6] {
    let b005 = b[0];
    let b014 = b[1];
    let b023 = b[2];
    let b104 = b[3];
    let b113 = b[4];
    let b203 = b[5];

    // q_{0,0,2} = b_{0,0,5}
    let q002 = b005;
    // q_{1,0,1} = (10 b_{1,0,4} - 6 b_{0,0,5}) / 4
    let q101 = scale_sub(2.5, b104, 1.5, b005);
    // q_{0,1,1} = (10 b_{0,1,4} - 6 b_{0,0,5}) / 4
    let q011 = scale_sub(2.5, b014, 1.5, b005);
    // q_{2,0,0} = 10 b_{2,0,3} - 15 b_{1,0,4} + 6 b_{0,0,5}
    let q200 = lin3(10.0, b203, -15.0, b104, 6.0, b005, 1.0);
    // q_{0,2,0} = 10 b_{0,2,3} - 15 b_{0,1,4} + 6 b_{0,0,5}
    let q020 = lin3(10.0, b023, -15.0, b014, 6.0, b005, 1.0);
    // q_{1,1,0} = 10 b_{1,1,3} - 7.5 b_{1,0,4} - 7.5 b_{0,1,4} + 6 b_{0,0,5}
    let q110 = lin_n(
        &[(10.0, b113), (-7.5, b104), (-7.5, b014), (6.0, b005)],
        1.0,
    );
    [q002, q011, q020, q101, q110, q200]
}

/// Recover `Q`'s monomial coefs from sector-0 tri-2 BB-coefs.
/// `q` is in the same order as [`tri2_to_tri5_near_star`] inputs.
fn tri2_coefs_to_quadratic(q: &[Point3; 6], n: usize) -> StarQuadratic {
    let c = (2.0 * std::f64::consts::PI / n as f64).cos();
    let s = (2.0 * std::f64::consts::PI / n as f64).sin();
    debug_assert!(s.abs() > 1e-12, "sin(2π/n) ≈ 0 for n = {n}");

    let q002 = q[0];
    let q011 = q[1];
    let q020 = q[2];
    let q101 = q[3];
    let q110 = q[4];
    let q200 = q[5];

    // ζ = q_{0,0,2}
    let zeta = q002;
    // δ/2 + ζ = q_{1,0,1} → δ = 2·(q_{1,0,1} − ζ)
    let delta = scale_sub(2.0, q101, 2.0, zeta);
    // δ·c/2 + ε·s/2 + ζ = q_{0,1,1} → ε = (2·(q_{0,1,1} − ζ) − δ·c) / s
    let epsilon = {
        let rhs = scale_sub(2.0, q011, 2.0, zeta);
        // rhs = δ·c + ε·s → ε = (rhs − δ·c) / s
        scale_sub(1.0 / s, rhs, c / s, delta)
    };
    // α + δ + ζ = q_{2,0,0} → α = q_{2,0,0} − δ − ζ
    let alpha = lin3(1.0, q200, -1.0, delta, -1.0, zeta, 1.0);
    // β·s/2 = q_{1,1,0} − α·c − δ·(1+c)/2 − ε·s/2 − ζ
    // → β = (2/s)·(q_{1,1,0} − α·c − δ·(1+c)/2 − ε·s/2 − ζ)
    let beta = {
        let rhs = Point3::new(
            q110.x - alpha.x * c - delta.x * (1.0 + c) * 0.5 - epsilon.x * s * 0.5 - zeta.x,
            q110.y - alpha.y * c - delta.y * (1.0 + c) * 0.5 - epsilon.y * s * 0.5 - zeta.y,
            q110.z - alpha.z * c - delta.z * (1.0 + c) * 0.5 - epsilon.z * s * 0.5 - zeta.z,
        );
        scale(2.0 / s, rhs)
    };
    // γ·s² = q_{0,2,0} − α·c² − β·c·s − δ·c − ε·s − ζ
    let gamma = {
        let rhs = Point3::new(
            q020.x - alpha.x * c * c - beta.x * c * s - delta.x * c - epsilon.x * s - zeta.x,
            q020.y - alpha.y * c * c - beta.y * c * s - delta.y * c - epsilon.y * s - zeta.y,
            q020.z - alpha.z * c * c - beta.z * c * s - delta.z * c - epsilon.z * s - zeta.z,
        );
        scale(1.0 / (s * s), rhs)
    };

    StarQuadratic {
        alpha,
        beta,
        gamma,
        delta,
        epsilon,
        zeta,
    }
}

// ---- small scalar helpers ----

#[inline]
fn scale(s: f64, p: Point3) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

/// `(sa·pa − sb·pb)` component-wise.
#[inline]
fn scale_sub(sa: f64, pa: Point3, sb: f64, pb: Point3) -> Point3 {
    Point3::new(
        sa * pa.x - sb * pb.x,
        sa * pa.y - sb * pb.y,
        sa * pa.z - sb * pb.z,
    )
}

/// `(a·pa + b·pb) / denom`.
#[inline]
fn lin2(a: f64, pa: Point3, b: f64, pb: Point3, denom: f64) -> Point3 {
    Point3::new(
        (a * pa.x + b * pb.x) / denom,
        (a * pa.y + b * pb.y) / denom,
        (a * pa.z + b * pb.z) / denom,
    )
}

/// `(a·pa + b·pb + c·pc) / denom`.
#[inline]
fn lin3(a: f64, pa: Point3, b: f64, pb: Point3, c: f64, pc: Point3, denom: f64) -> Point3 {
    Point3::new(
        (a * pa.x + b * pb.x + c * pc.x) / denom,
        (a * pa.y + b * pb.y + c * pc.y) / denom,
        (a * pa.z + b * pb.z + c * pc.z) / denom,
    )
}

/// `(Σᵢ wᵢ·pᵢ) / denom`.
#[inline]
fn lin_n(pairs: &[(f64, Point3)], denom: f64) -> Point3 {
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for &(w, p) in pairs {
        x += w * p.x;
        y += w * p.y;
        z += w * p.z;
    }
    Point3::new(x / denom, y / denom, z / denom)
}

#[cfg(test)]
#[path = "kp13_bi5_star_quadratic_tests.rs"]
mod tests;
