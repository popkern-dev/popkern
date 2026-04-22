//! Phase I.13b — σ, the guide sampling reparameterisation.
//!
//! KP13 §4.3 defines σ as a symmetric planar bi-5 patch that relates
//! the guide `g`'s domain to the cap `p`'s pieces. It is **data-
//! independent** (depends only on valence `n`) and can be tabulated
//! once per `n`.
//!
//! # Paper specification
//!
//! > "σ is a symmetric planar 2 × 2 patch of degree bi-5 that matches
//! > the tensor-border of χ, reparameterized by β [...]. The map σ
//! > is chosen to be symmetric with respect to its diagonal `u = v`,
//! > to match the (equally symmetric) tensor-border of χ, to be
//! > internally C² and to enforce G² continuity with R_n β using
//! > the reparameterization ρ. Points along the boundary must lie on
//! > the sector partition lines, and symmetry across the diagonal
//! > remove half the degrees of freedom. After enforcing the G²
//! > constraints (see Fig. 8b) along the inter-sector boundaries, σ
//! > has 15 scalar parameters to be determined. These are fixed by
//! > minimizing F_3 σ."
//!
//! # Representation
//!
//! This module provides [`SigmaPatch`] — a single bi-5 tensor-
//! product patch with 6×6 = 36 control points, each a 2-D point
//! (σ maps `R² → R²`). Planar (hence the paper's "planar").
//!
//! Evaluation + 2nd-order jets at the four unit-square corners are
//! implemented via the tensor-product Bernstein basis; these are
//! what I.14 consumes for the cap assembly equations.
//!
//! # Scope
//!
//! This stage ships:
//!
//! * [`SigmaPatch`] data type + evaluation / jet extraction
//! * [`SigmaPatch::identity`] — the identity map `σ(u, v) = (u, v)`,
//!   useful as a bootstrapping default for I.14 development
//! * Symmetry predicate + enforcement helper
//! * [`f3_functional`] — the F_3 quadratic functional value for
//!   testing / as an objective scalar (ready to plug into a solver)
//!
//! Out of scope for this stage:
//!
//! * F_3 constrained minimisation solver (the proper KP13 σ values
//!   per `n`). A working bootstrap uses `SigmaPatch::identity`; the
//!   F_3 solve is a follow-on that reduces the jet-sampling error
//!   at sector boundaries.
//!
//! See `docs/journal/polar-kp-math.md` §3.5 for the math summary.

#![allow(dead_code)]

use super::bernstein_bi5::{
    bernstein_bi5_all, bernstein_bi5_deriv_all, bernstein_bi5_second_deriv_all,
};
use super::kp13_bi5_chi::ChiCorner;

/// A symmetric planar bi-5 patch `σ : [0, 1]² → R²`. 6×6 grid of
/// 2-D BB control points.
///
/// Symmetric about the diagonal `u = v` when `cps[i][j] ==
/// swap(cps[j][i])` where `swap((x, y)) = (y, x)` — see
/// [`SigmaPatch::is_symmetric`].
#[derive(Debug, Clone, Copy)]
pub struct SigmaPatch {
    /// `cps[i][j] = (σ_x, σ_y)` at BB-index `(i, j)`.
    pub cps: [[(f64, f64); 6]; 6],
}

/// 2-D Taylor-jet up to order 2 of a planar function at a point.
/// All six components are `(∂ to some order)` expressed as 2-D
/// vectors.
#[derive(Debug, Clone, Copy)]
pub struct SigmaJet {
    pub pos: (f64, f64),
    pub d_u: (f64, f64),
    pub d_v: (f64, f64),
    pub d_uu: (f64, f64),
    pub d_uv: (f64, f64),
    pub d_vv: (f64, f64),
}

impl SigmaPatch {
    /// σ(u, v) = (u, v) — the identity map. Fills the 6×6 CP grid
    /// with `cps[i][j] = (i/5, j/5)`.
    ///
    /// Useful as a bootstrap default for I.14 development. For
    /// paper-faithful σ values, the F_3 constrained minimisation
    /// (not implemented here) is required.
    pub fn identity() -> Self {
        let mut cps = [[(0.0, 0.0); 6]; 6];
        for (i, row) in cps.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (i as f64 / 5.0, j as f64 / 5.0);
            }
        }
        Self { cps }
    }

    /// Evaluate σ at `(u, v)` via tensor-product Bernstein basis.
    pub fn eval(&self, u: f64, v: f64) -> (f64, f64) {
        let bu = bernstein_bi5_all(u);
        let bv = bernstein_bi5_all(v);
        let mut x = 0.0;
        let mut y = 0.0;
        for (i, bu_i) in bu.iter().enumerate() {
            for (j, bv_j) in bv.iter().enumerate() {
                let w = bu_i * bv_j;
                x += w * self.cps[i][j].0;
                y += w * self.cps[i][j].1;
            }
        }
        (x, y)
    }

    /// Full second-order Taylor jet at the given unit-square corner.
    /// Each jet component is a 2-D vector (∂ direction in σ's image).
    pub fn jet_at_corner(&self, corner: ChiCorner) -> SigmaJet {
        let (u, v) = corner.uv();
        let bu = bernstein_bi5_all(u);
        let bv = bernstein_bi5_all(v);
        let bu_d = bernstein_bi5_deriv_all(u);
        let bv_d = bernstein_bi5_deriv_all(v);
        let bu_dd = bernstein_bi5_second_deriv_all(u);
        let bv_dd = bernstein_bi5_second_deriv_all(v);

        let combine = |wu: &[f64; 6], wv: &[f64; 6]| -> (f64, f64) {
            let mut x = 0.0;
            let mut y = 0.0;
            for (i, wu_i) in wu.iter().enumerate() {
                for (j, wv_j) in wv.iter().enumerate() {
                    let w = wu_i * wv_j;
                    x += w * self.cps[i][j].0;
                    y += w * self.cps[i][j].1;
                }
            }
            (x, y)
        };

        SigmaJet {
            pos: combine(&bu, &bv),
            d_u: combine(&bu_d, &bv),
            d_v: combine(&bu, &bv_d),
            d_uu: combine(&bu_dd, &bv),
            d_uv: combine(&bu_d, &bv_d),
            d_vv: combine(&bu, &bv_dd),
        }
    }

    /// Check KP13 §4.3's symmetry requirement: σ is symmetric about
    /// `u = v` iff `cps[i][j]` equals `swap(cps[j][i])` (componentwise)
    /// where `swap((x, y)) = (y, x)`.
    pub fn is_symmetric(&self) -> bool {
        const EPS: f64 = 1e-12;
        for i in 0..6 {
            for j in 0..6 {
                let p = self.cps[i][j];
                let q = self.cps[j][i];
                if (p.0 - q.1).abs() > EPS || (p.1 - q.0).abs() > EPS {
                    return false;
                }
            }
        }
        true
    }

    /// Enforce symmetry by averaging each pair `(cps[i][j], swap(cps[j][i]))`
    /// and writing the average back. After this, [`Self::is_symmetric`]
    /// returns `true`.
    pub fn enforce_symmetry(&mut self) {
        for i in 0..6 {
            for j in i..6 {
                let p = self.cps[i][j];
                let q = self.cps[j][i];
                let avg = ((p.0 + q.1) * 0.5, (p.1 + q.0) * 0.5);
                self.cps[i][j] = avg;
                self.cps[j][i] = (avg.1, avg.0);
            }
        }
    }
}

/// F_3 quadratic shape functional for a planar bi-5 patch over the
/// unit square (KP09b § 2.2 / KP13 § 4.3):
///
/// ```text
///   F_3(σ) := ∫∫_[0,1]² Σ_{i+j+k=3, i,j,k≥0} C(3; i, j, k)
///                (∂_u^i ∂_v^j (u - v)^k σ)² du dv
/// ```
///
/// The `(u - v)^k` factor is a *polynomial* multiplier, not a
/// derivative — it weighs higher-order partials by their distance
/// from the diagonal. This matches KP09b / KP13's definition.
///
/// Returns a non-negative scalar that the paper minimises over the
/// 15 free parameters of σ (once boundary + symmetry + G² are
/// fixed). Useful here as (1) a regression target to verify σ-
/// solver implementations and (2) a stand-alone quality measure.
///
/// We use numerical quadrature — 5-point Gauss-Legendre on each
/// axis (adequate for degree-5 polynomial integrands of degree
/// up to 2·5 = 10 post-differentiation, which 5-point Gauss-
/// Legendre integrates exactly).
pub fn f3_functional(sigma: &SigmaPatch) -> f64 {
    // 5-point Gauss-Legendre nodes + weights on [0, 1].
    let nodes = [
        (
            0.5 - 0.453_846_931_527_77 * 0.5,
            0.236_926_885_056_189 * 0.5,
        ),
        (
            0.5 - 0.269_356_648_453_358 * 0.5,
            0.478_628_670_499_366 * 0.5,
        ),
        (0.5, 0.568_888_888_888_889 * 0.5),
        (
            0.5 + 0.269_356_648_453_358 * 0.5,
            0.478_628_670_499_366 * 0.5,
        ),
        (
            0.5 + 0.453_846_931_527_77 * 0.5,
            0.236_926_885_056_189 * 0.5,
        ),
    ];

    let mut total = 0.0;
    for (u, wu) in nodes {
        for (v, wv) in nodes {
            let w = wu * wv;
            // Partials of σ at (u, v).
            let bu = bernstein_bi5_all(u);
            let bv = bernstein_bi5_all(v);
            let bu_d = bernstein_bi5_deriv_all(u);
            let bv_d = bernstein_bi5_deriv_all(v);
            let bu_dd = bernstein_bi5_second_deriv_all(u);
            let bv_dd = bernstein_bi5_second_deriv_all(v);
            let bu_ddd = super::bernstein_bi5::bernstein_bi5_third_deriv_all(u);
            let bv_ddd = super::bernstein_bi5::bernstein_bi5_third_deriv_all(v);

            let combine = |wu: &[f64; 6], wv: &[f64; 6]| -> (f64, f64) {
                let mut x = 0.0;
                let mut y = 0.0;
                for (i, wu_i) in wu.iter().enumerate() {
                    for (j, wv_j) in wv.iter().enumerate() {
                        let w = wu_i * wv_j;
                        x += w * sigma.cps[i][j].0;
                        y += w * sigma.cps[i][j].1;
                    }
                }
                (x, y)
            };

            // The three 3rd-order partials of σ.
            let uuu = combine(&bu_ddd, &bv);
            let uuv = combine(&bu_dd, &bv_d);
            let uvv = combine(&bu_d, &bv_dd);
            let vvv = combine(&bu, &bv_ddd);

            // F_3 integrand: Σ C(3; i, j, k) (∂^{i+j+k} σ · (u - v)^k)²
            // Trinomial coefs C(3; i, j, k) = 3!/(i!j!k!).
            //   i+j+k = 3 partitions:
            //     (3,0,0) → C = 1,  ∂_u^3 σ · (u-v)^0 = ∂³u
            //     (0,3,0) → 1,  ∂³v
            //     (0,0,3) → 1,  (u-v)^3 · σ  (no partials)
            //     (2,1,0) → 3,  ∂²u∂v
            //     (2,0,1) → 3,  ∂²u · (u-v)
            //     (1,2,0) → 3,  ∂u∂²v
            //     (1,0,2) → 3,  ∂u · (u-v)²
            //     (0,2,1) → 3,  ∂²v · (u-v)
            //     (0,1,2) → 3,  ∂v · (u-v)²
            //     (1,1,1) → 6,  ∂u∂v · (u-v)
            //
            // Terms involving `(u-v)^k` with k > 0 have no partials of σ
            // (the (u-v)^k is a polynomial multiplier, not a partial).
            // They represent σ's own value times (u-v)^k.
            //
            // Partials only: (3,0,0), (0,3,0), (2,1,0), (1,2,0)
            // Mixed (partial + (u-v)^k): (2,0,1), (0,2,1), (1,0,2),
            //                              (0,1,2), (1,1,1)
            // Pure (u-v)^3 σ: (0,0,3)
            //
            // Evaluate σ itself for the (u-v)^k-only and mixed terms.
            let sigma_val = combine(&bu, &bv);
            let dd_u = combine(&bu_d, &bv); // ∂u
            let dd_v = combine(&bu, &bv_d); // ∂v
            let dd_uu = combine(&bu_dd, &bv);
            let dd_uv = combine(&bu_d, &bv_d);
            let dd_vv = combine(&bu, &bv_dd);

            let uv_diff = u - v;

            let sq = |p: (f64, f64)| p.0 * p.0 + p.1 * p.1;
            let scale = |s: f64, p: (f64, f64)| (s * p.0, s * p.1);

            // Assemble F_3 integrand at this (u, v).
            let mut f = 0.0;
            f += 1.0 * sq(uuu);
            f += 1.0 * sq(vvv);
            f += 1.0 * sq(scale(uv_diff.powi(3), sigma_val));
            f += 3.0 * sq(uuv);
            f += 3.0 * sq(scale(uv_diff, dd_uu));
            f += 3.0 * sq(uvv);
            f += 3.0 * sq(scale(uv_diff * uv_diff, dd_u));
            f += 3.0 * sq(scale(uv_diff, dd_vv));
            f += 3.0 * sq(scale(uv_diff * uv_diff, dd_v));
            f += 6.0 * sq(scale(uv_diff, dd_uv));

            total += w * f;
        }
    }
    total
}

#[cfg(test)]
#[path = "kp13_bi5_sigma_tests.rs"]
mod tests;
