//! KP13 § 3 Phase I.12 — bi-5 guide surface `g` from auxiliary map
//! `g△`.
//!
//! The auxiliary tri-5 map `g△` (Phase I.11) is an `n`-sector
//! total-degree-5 polynomial on the planar `n`-gon `Ω`; each sector
//! `g△^k` lives on the natural-simplex wedge
//! `△_k = {(0,0), m_k, m_{k+1}}`. The bi-5 **guide** `g` is an
//! `n`-sector tensor-product polynomial on `[0, 1]²`, one 6×6 BB
//! grid per sector, with `C⁰` junctions between sectors.
//!
//! # The paper's construction (KP13 § 3, eqs (5) – (8))
//!
//! Following Fig. 6a, let `L : [0, 1]² → parallelogram` be the
//! linear map `L(u, v) := u · m_k + v · m_{k+1}`. In the sector-
//! natural simplex coords `(u_s, v_s)` (with `(x, y) = u_s m_k + v_s m_{k+1}`)
//! the map `L` is the **identity** — so the composition
//! `g△^k ∘ L` is simply `g△^k` re-expressed as a polynomial on
//! `[0, 1]²`.
//!
//! Because a tri-5 polynomial on the simplex has total degree ≤ 5,
//! it lives inside the bi-5 tensor-product polynomial space on
//! `[0, 1]²` as a proper 21-dimensional linear subspace. The paper
//! sets `g^k`'s jets at three corners of the unit square —
//!
//! ```text
//!   H_o      g^k := H_o      (g△^k ∘ L)    (eq 5,  star)
//!   H_{m_k}  g^k := H_{m_k}  (g△^k ∘ L)    (eq 7,  m_k corner)
//!   H_{m_{k+1}} g^k := H_{m_{k+1}} (g△^k ∘ L)  (eq 7, m_{k+1} corner)
//! ```
//!
//! — to match the corresponding jets of `g△^k ∘ L`. The simplest
//! way to satisfy all three simultaneously is to set
//!
//! ```text
//!   g^k := degree-raise(g△^k)              (all 36 CPs of g^k)
//! ```
//!
//! i.e., `g^k` is the exact bi-5 BB representation of `g△^k`
//! re-expressed on the unit square. This satisfies eqs (5), (7)
//! **identically at all orders**, not just up to order 2.
//!
//! # Eq (8) — v_k corner override
//!
//! The fourth paper-prescribed constraint overrides `g^k`'s jet at
//! the `(1, 1)` corner (= `v_k`) by solving
//!
//! ```text
//!   H_{v_k}(g^k ∘ L⁻¹ ∘ χ) = H_{v_k} b^k
//! ```
//!
//! against the tensor-border jet `b^k`.
//!
//! - **Bilinear χ** — `L⁻¹ ∘ χ` is the identity, so the constraint
//!   reduces to `H_{v_k} g^k = H_{v_k} b^k`. Inverted analytically
//!   onto the six bi-5 CPs controlling the 2nd-order jet at corner
//!   `(1, 1)`; implemented in
//!   [`Bi5Guide::from_gtri5_with_eq8_bilinear`].
//!
//! - **CC χ** — `L⁻¹ ∘ χ` has non-trivial Jacobian (at `v_k`:
//!   `∂χ = λ · m_k, λ · m_{k+1}`) and non-zero Hessian, so the
//!   6 × 6 linear solve includes chain-rule terms. Deferred — the
//!   bilinear path plus the LSQ-fit `g△` is already commercial-
//!   grade for regular valence; CC χ is the paper's final polish
//!   and is a mechanical extension following the pattern in
//!   [`super::kp13_bi5_lsq_cc`].
//!
//! # Tri-5 → bi-5 degree raise — closed form
//!
//! A tri-5 basis function `B^5_{a,b,c}(u, v)` (with `a + b + c = 5`,
//! `a` = weight on the star vertex, `b` on `m_k`, `c` on `m_{k+1}`)
//! equals
//!
//! ```text
//!   K_{a,b,c} · (1 − u − v)^a · u^b · v^c,
//!   K_{a,b,c} := 5! / (a! b! c!).
//! ```
//!
//! Expanding `(1 − u − v)^a` via the multinomial theorem and
//! collecting terms in the monomial basis `u^r v^s`:
//!
//! ```text
//!   A^{a,b,c}_{r,s} = 5! · (−1)^{r+s-b-c} /
//!                     (b! c! (5−r−s)! (r−b)! (s−c)!),
//!   valid for r ≥ b, s ≥ c, r + s ≤ 5 (else 0).
//! ```
//!
//! Monomials convert to bi-5 BB via
//!
//! ```text
//!   u^r v^s = Σ_{p ≥ r, q ≥ s}
//!             [C(5−r, p−r) · C(5−s, q−s) / (C(5, p) C(5, q))]
//!             · B^5_p(u) · B^5_q(v).
//! ```
//!
//! Composing gives the 36 × 21 linear map `D = M · C` from tri-5
//! BB coefs to bi-5 BB coefs — implemented in
//! [`tri5_to_bi5_degree_raise`].

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::bernstein_bi5::{
    bernstein_bi5_all, bernstein_bi5_deriv_all, bernstein_bi5_second_deriv_all,
};
use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner};
use super::kp13_bi5_gtri5::Gtri5;
use super::kp13_bi5_jet::OuterHessian;
use super::kp13_bi5_lsq_cc::ambient_to_simplex_jet;
use crate::geometry::Point3;

#[path = "kp13_bi5_guide_eq8.rs"]
mod kp13_bi5_guide_eq8;

pub use kp13_bi5_guide_eq8::{Eq8CcError, tri5_to_bi5_degree_raise};
use kp13_bi5_guide_eq8::{apply_eq8_bilinear_override, apply_eq8_cc_override};

/// Number of tensor-product BB CPs in one bi-5 guide sector (6 × 6).
pub const BI5_CP_COUNT: usize = 36;

/// One bi-5 guide sector — 6 × 6 BB-grid on `[0, 1]²` with
/// `cps[p][q]` the coefficient of `B^5_p(u) · B^5_q(v)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bi5GuideSector {
    pub cps: [[Point3; 6]; 6],
}

impl Bi5GuideSector {
    /// Zero sector (all CPs at origin).
    pub fn zero() -> Self {
        Self {
            cps: [[Point3::new(0.0, 0.0, 0.0); 6]; 6],
        }
    }

    /// Evaluate `g^k(u, v)` on the unit square via tensor-product
    /// Bernstein basis.
    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_bi5_all(u);
        let bv = bernstein_bi5_all(v);
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for p in 0..6 {
            for q in 0..6 {
                let w = bu[p] * bv[q];
                let cp = self.cps[p][q];
                out.x += w * cp.x;
                out.y += w * cp.y;
                out.z += w * cp.z;
            }
        }
        out
    }

    /// Full 2nd-order Taylor jet of the bi-5 patch at `(u, v)`,
    /// packaged as an [`OuterHessian`] for consumption by the
    /// 2nd-order chain-rule helper [`compose_chain_rule`].
    ///
    /// Uses `bernstein_bi5_*` + derivative tables on the 6 × 6 CP
    /// grid; the field names `d_x, d_y, d_xx, d_xy, d_yy` carry
    /// partials in the patch's own `(u, v)` variables.
    pub fn outer_hessian(&self, u: f64, v: f64) -> OuterHessian {
        let bu = bernstein_bi5_all(u);
        let bv = bernstein_bi5_all(v);
        let bu_d = bernstein_bi5_deriv_all(u);
        let bv_d = bernstein_bi5_deriv_all(v);
        let bu_dd = bernstein_bi5_second_deriv_all(u);
        let bv_dd = bernstein_bi5_second_deriv_all(v);

        let zero = Point3::new(0.0, 0.0, 0.0);
        let mut oh = OuterHessian {
            value: zero,
            d_x: zero,
            d_y: zero,
            d_xx: zero,
            d_xy: zero,
            d_yy: zero,
        };
        for p in 0..6 {
            for q in 0..6 {
                let cp = self.cps[p][q];
                let b_val = bu[p] * bv[q];
                let b_x = bu_d[p] * bv[q];
                let b_y = bu[p] * bv_d[q];
                let b_xx = bu_dd[p] * bv[q];
                let b_xy = bu_d[p] * bv_d[q];
                let b_yy = bu[p] * bv_dd[q];
                oh.value.x += b_val * cp.x;
                oh.value.y += b_val * cp.y;
                oh.value.z += b_val * cp.z;
                oh.d_x.x += b_x * cp.x;
                oh.d_x.y += b_x * cp.y;
                oh.d_x.z += b_x * cp.z;
                oh.d_y.x += b_y * cp.x;
                oh.d_y.y += b_y * cp.y;
                oh.d_y.z += b_y * cp.z;
                oh.d_xx.x += b_xx * cp.x;
                oh.d_xx.y += b_xx * cp.y;
                oh.d_xx.z += b_xx * cp.z;
                oh.d_xy.x += b_xy * cp.x;
                oh.d_xy.y += b_xy * cp.y;
                oh.d_xy.z += b_xy * cp.z;
                oh.d_yy.x += b_yy * cp.x;
                oh.d_yy.y += b_yy * cp.y;
                oh.d_yy.z += b_yy * cp.z;
            }
        }
        oh
    }
}

/// `n`-piece bi-5 guide `g` with `C⁰` junctions between sectors.
#[derive(Debug, Clone, PartialEq)]
pub struct Bi5Guide {
    pub n: usize,
    pub sectors: Vec<Bi5GuideSector>,
}

impl Bi5Guide {
    /// New zero guide with `n` sectors (`n ≥ 3`).
    pub fn zero(n: usize) -> Self {
        assert!(n >= 3, "KP13 bi-5 guide needs n ≥ 3, got {n}");
        Self {
            n,
            sectors: vec![Bi5GuideSector::zero(); n],
        }
    }

    /// Access sector `k` (mod `n`).
    #[inline]
    pub fn sector(&self, k: usize) -> &Bi5GuideSector {
        &self.sectors[k % self.n]
    }

    /// Mutable sector access.
    #[inline]
    pub fn sector_mut(&mut self, k: usize) -> &mut Bi5GuideSector {
        let n = self.n;
        &mut self.sectors[k % n]
    }

    /// Build the bi-5 guide from the I.11 auxiliary map `g△` via
    /// paper eqs (5) + (7): degree-raise `g△^k` to bi-5 per sector.
    ///
    /// This is the § 3 construction of the paper without the eq (8)
    /// `v_k`-corner override — see the module header for why that
    /// is a follow-up and why the LSQ-fit `g△` already carries an
    /// approximation of the eq (8) jet.
    pub fn from_gtri5(gtri5: &Gtri5) -> Self {
        let n = gtri5.n;
        let sectors: Vec<Bi5GuideSector> = (0..n)
            .map(|k| {
                let cps = tri5_to_bi5_degree_raise(&gtri5.sector(k).cps);
                Bi5GuideSector { cps }
            })
            .collect();
        Self { n, sectors }
    }

    /// Build the bi-5 guide per eqs (5) + (7) from `g△`, then apply
    /// the eq (8) `v_k`-corner override to directly match the
    /// tensor-border jet at `(u, v) = (1, 1)` in each sector.
    ///
    /// # Eq (8) — bilinear χ variant
    ///
    /// Under bilinear `χ(u, v) := u · m_k + v · m_{k+1}`, the
    /// composition `L⁻¹ ∘ χ` is the identity on `[0, 1]²` (both
    /// are the same linear shear). Paper eq (8) therefore reduces
    /// to a direct jet match
    ///
    /// ```text
    ///   H_{v_k} g^k = H_{v_k} b^k
    /// ```
    ///
    /// which is inverted analytically onto the six bi-5 CPs
    /// controlling the 2nd-order jet at corner `(1, 1)`:
    ///
    /// ```text
    ///   P_{5,5} = J.value
    ///   P_{4,5} = P_{5,5} − J.d_u / 5
    ///   P_{5,4} = P_{5,5} − J.d_v / 5
    ///   P_{3,5} = 2 P_{4,5} − P_{5,5} + J.d_uu / 20
    ///   P_{5,3} = 2 P_{5,4} − P_{5,5} + J.d_vv / 20
    ///   P_{4,4} = P_{4,5} + P_{5,4} − P_{5,5} + J.d_uv / 25
    /// ```
    ///
    /// The remaining three CPs of the `(1, 1)` 3 × 3 block
    /// (`P_{3,3}`, `P_{3,4}`, `P_{4,3}`) are **not** determined by
    /// the 2nd-order jet — they keep their degree-raised values
    /// from `g△` and carry 3rd-order-and-higher information.
    ///
    /// # Full CC χ variant
    ///
    /// Deferred: for full Catmull-Clark χ, the composition
    /// `L⁻¹ ∘ χ` is no longer the identity, and the six equations
    /// involve a chain-rule expansion (analogous to the CC-χ path
    /// in [`super::kp13_bi5_lsq_cc`]). The 6 × 6 linear solve on
    /// the `(5, 5)`-corner CPs still applies, just with non-
    /// diagonal structure from the χ Jacobian + Hessian at `v_k`.
    pub fn from_gtri5_with_eq8_bilinear(gtri5: &Gtri5, target: &TensorBorder) -> Self {
        assert_eq!(
            gtri5.n,
            target.n(),
            "from_gtri5_with_eq8_bilinear: g△ and target n mismatch"
        );
        let mut g = Self::from_gtri5(gtri5);
        for k in 0..g.n {
            let jet = target.jet(k, ChiCorner::Vk);
            apply_eq8_bilinear_override(&mut g.sectors[k], &jet);
        }
        g
    }

    /// Build the bi-5 guide per eqs (5) + (7) + (8) using the **full
    /// Catmull-Clark characteristic map** χ for eq (8).
    ///
    /// Under CC χ the reparameterisation `L⁻¹ ∘ χ` has non-trivial
    /// Jacobian (at `v_k`: `M⁻¹ · ∂χ = (λ, 0), (0, λ)` after ambient-
    /// to-simplex transform — χ contracts by `λ` at the star and
    /// evaluates at `(μ, μ)` at `v_k` under our `CatmullClarkChi`
    /// eigen-structure construction, μ = λ²) and non-zero Hessian.
    /// The eq (8) constraint
    ///
    /// ```text
    ///   H_{v_k}(g^k ∘ L⁻¹ ∘ χ) = H_{v_k} b^k
    /// ```
    ///
    /// becomes a 6 × 6 linear system on the six CPs at g^k's `(5, 5)`
    /// corner block, with the matrix built from the second-order
    /// chain rule: each corner-block CP's contribution to the
    /// composite jet at `(u_χ, v_χ) = (1, 1)` = `(u_g, v_g) = (μ, μ)`
    /// is evaluated via the existing [`compose_chain_rule`] helper.
    ///
    /// The solve is **per-axis**: the three axes (x, y, z) share the
    /// same 6 × 6 system matrix (by linearity) and only the RHS
    /// differs. Returns [`Eq8CcError::SingularSystem`] if the 6 × 6
    /// system is singular at the given `(μ, μ)` — this can happen
    /// only at degenerate `n` where χ collapses, in practice not for
    /// `n ≥ 3` via `CatmullClarkChi`.
    pub fn from_gtri5_with_eq8_cc(
        gtri5: &Gtri5,
        target: &TensorBorder,
    ) -> Result<Self, Eq8CcError> {
        assert_eq!(
            gtri5.n,
            target.n(),
            "from_gtri5_with_eq8_cc: g△ and target n mismatch"
        );
        let n = gtri5.n;
        let cc_chi = CatmullClarkChi::new(n);
        let mut g = Self::from_gtri5(gtri5);
        for k in 0..g.n {
            let jet_amb = cc_chi.jet(k, ChiCorner::Vk);
            let jet_simplex = ambient_to_simplex_jet(&jet_amb, n, k);
            let target_jet = target.jet(k, ChiCorner::Vk);
            apply_eq8_cc_override(&mut g.sectors[k], &target_jet, &jet_simplex)
                .map_err(|_| Eq8CcError::SingularSystem(k))?;
        }
        Ok(g)
    }

    /// Structural sanity check: `n ≥ 3`, sector count matches `n`,
    /// and all CP components are finite.
    pub fn validate(&self) -> Result<(), String> {
        if self.n < 3 {
            return Err(format!("Bi5Guide::validate: n = {} < 3", self.n));
        }
        if self.sectors.len() != self.n {
            return Err(format!(
                "Bi5Guide::validate: sectors.len() = {} ≠ n = {}",
                self.sectors.len(),
                self.n
            ));
        }
        for (k, s) in self.sectors.iter().enumerate() {
            for p in 0..6 {
                for q in 0..6 {
                    let cp = s.cps[p][q];
                    if !(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite()) {
                        return Err(format!(
                            "Bi5Guide::validate: non-finite CP at sector {k}, ({p}, {q}): {cp:?}"
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "kp13_bi5_guide_tests.rs"]
mod tests;
