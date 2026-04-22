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

use super::bernstein_bi5::bernstein_bi5_all;
use super::bernstein_tri5::{TRI5_CP_COUNT, flat_index};
use super::kp13_bi5_boundary::TensorBorder;
use super::kp13_bi5_chi::ChiCorner;
use super::kp13_bi5_gtri5::Gtri5;
use super::kp13_bi5_jet::TaylorJet2;
use crate::geometry::Point3;

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

// ---- tri-5 → bi-5 degree-raise ----

/// Closed-form tri-5 → bi-5 degree raise.
///
/// Converts 21 tri-5 BB coefs (one sector of `g△^k`) into the 6 × 6
/// bi-5 BB coefs of the same polynomial expressed on `[0, 1]²`.
///
/// The tri-5 polynomial
/// `P(u, v) = Σ_{i+j ≤ 5} C_{ij} · B^5_{5−i−j, i, j}(u, v)`
/// (treating `(u, v)` as the natural-simplex coords `(u_s, v_s)`
/// via the linear shear `L`) equals the bi-5 polynomial
/// `Q(u, v) = Σ_{p, q ∈ 0..6} D_{pq} · B^5_p(u) · B^5_q(v)`
/// where each `D_{pq}` is a specific linear combination of the 21
/// `C_{ij}`. The coefficients are derived from the multinomial
/// expansion of `(1−u−v)^{5−i−j}` followed by the closed-form
/// monomial-to-Bernstein conversion.
pub fn tri5_to_bi5_degree_raise(cps_tri: &[Point3; TRI5_CP_COUNT]) -> [[Point3; 6]; 6] {
    let zero = Point3::new(0.0, 0.0, 0.0);

    // Step 1: tri-5 BB → bivariate monomials.
    // monomials[r][s] (for r + s ≤ 5) is the coefficient of u^r v^s.
    let mut monomials = [[zero; 6]; 6];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let cp = cps_tri[flat_index(5, i, j)];
            // B^5_{a,b,c} with (a, b, c) = (5 − i − j, i, j):
            //   K · (1 − u − v)^a · u^b · v^c
            // = K · Σ_{r ≥ b, s ≥ c, r+s ≤ 5}
            //       [a! / ((5 − r − s)! (r − b)! (s − c)!) · (−1)^{r+s-b-c}]
            //       · u^r · v^s
            // with K = 5! / (a! b! c!), so K · a! = 5! / (b! c!).
            for r in i..=5 {
                for s in j..=(5 - r) {
                    let r_b = r - i;
                    let s_c = s - j;
                    let i_exp = 5 - r - s;
                    let sign = if (r_b + s_c) % 2 == 0 { 1.0 } else { -1.0 };
                    let numer = factorial(5) as f64;
                    let denom = (factorial(i)
                        * factorial(j)
                        * factorial(i_exp)
                        * factorial(r_b)
                        * factorial(s_c)) as f64;
                    let coef = sign * numer / denom;
                    monomials[r][s].x += coef * cp.x;
                    monomials[r][s].y += coef * cp.y;
                    monomials[r][s].z += coef * cp.z;
                }
            }
        }
    }

    // Step 2: bivariate monomials → bi-5 BB.
    // u^r = Σ_{p=r}^{5} [C(5 − r, p − r) / C(5, p)] · B^5_p(u).
    // Bivariate: u^r v^s = Σ_{p ≥ r, q ≥ s} w(r, p) · w(s, q) · B^5_p(u) B^5_q(v).
    let mut bi5 = [[zero; 6]; 6];
    for p in 0..6 {
        for q in 0..6 {
            let denom = (binomial(5, p) * binomial(5, q)) as f64;
            for r in 0..=p {
                for s in 0..=q {
                    if r + s > 5 {
                        continue;
                    }
                    let numer = (binomial(5 - r, p - r) * binomial(5 - s, q - s)) as f64;
                    let w = numer / denom;
                    bi5[p][q].x += w * monomials[r][s].x;
                    bi5[p][q].y += w * monomials[r][s].y;
                    bi5[p][q].z += w * monomials[r][s].z;
                }
            }
        }
    }
    bi5
}

/// Apply the bilinear-χ eq (8) override to one sector's bi-5 CPs.
///
/// Writes `P_{5,5}, P_{4,5}, P_{5,4}, P_{3,5}, P_{5,3}, P_{4,4}` so
/// that the 2nd-order bi-5 jet at corner `(1, 1)` equals `jet`.
/// The other 30 CPs (including `P_{3,3}, P_{3,4}, P_{4,3}`) are
/// untouched.
fn apply_eq8_bilinear_override(sector: &mut Bi5GuideSector, jet: &TaylorJet2) {
    let value = jet.value;
    let d_u = jet.d_u;
    let d_v = jet.d_v;
    let d_uu = jet.d_uu;
    let d_uv = jet.d_uv;
    let d_vv = jet.d_vv;

    let p55 = value;
    let p45 = Point3::new(
        p55.x - d_u.x / 5.0,
        p55.y - d_u.y / 5.0,
        p55.z - d_u.z / 5.0,
    );
    let p54 = Point3::new(
        p55.x - d_v.x / 5.0,
        p55.y - d_v.y / 5.0,
        p55.z - d_v.z / 5.0,
    );
    let p35 = Point3::new(
        2.0 * p45.x - p55.x + d_uu.x / 20.0,
        2.0 * p45.y - p55.y + d_uu.y / 20.0,
        2.0 * p45.z - p55.z + d_uu.z / 20.0,
    );
    let p53 = Point3::new(
        2.0 * p54.x - p55.x + d_vv.x / 20.0,
        2.0 * p54.y - p55.y + d_vv.y / 20.0,
        2.0 * p54.z - p55.z + d_vv.z / 20.0,
    );
    let p44 = Point3::new(
        p45.x + p54.x - p55.x + d_uv.x / 25.0,
        p45.y + p54.y - p55.y + d_uv.y / 25.0,
        p45.z + p54.z - p55.z + d_uv.z / 25.0,
    );

    sector.cps[5][5] = p55;
    sector.cps[4][5] = p45;
    sector.cps[5][4] = p54;
    sector.cps[3][5] = p35;
    sector.cps[5][3] = p53;
    sector.cps[4][4] = p44;
}

// ---- small-integer helpers ----

#[inline]
const fn factorial(n: usize) -> u64 {
    let mut acc: u64 = 1;
    let mut k = 1;
    while k <= n {
        acc *= k as u64;
        k += 1;
    }
    acc
}

#[inline]
const fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    factorial(n) / (factorial(k) * factorial(n - k))
}

#[cfg(test)]
#[path = "kp13_bi5_guide_tests.rs"]
mod tests;
