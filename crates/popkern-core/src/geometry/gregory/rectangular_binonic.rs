//! Phase H — `RectangularGregoryBiNonic`, a bidegree-(9, 9) pure
//! Bezier tensor-product patch with 10×10 = 100 control points.
//!
//! This is the substrate for the Karciauskas-Peters 2004
//! guided-surface construction (ADR
//! `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`).
//! The K-P construction derives:
//!
//! * Rows 0, 9 (boundary) from degree-elevated input spoke curves.
//! * Rows 1, 8 from the guide surface's first cross-boundary
//!   derivative — gives G1 continuity across every shared spoke
//!   by construction.
//! * Rows 2, 7 from the guide's second cross-boundary derivative
//!   — gives G2 continuity by construction.
//! * Rows 3..=6 from a least-squares approximation of the guide
//!   with the rows above pinned.
//!
//! The patch type itself is just a plain-Bezier tensor product;
//! the K-P construction machinery lives in `kp_*.rs` siblings
//! (pending H.4–H.6 — at H.3 none are landed yet, so this patch
//! type is only constructible via [`Self::zeros`] or direct grid
//! assignment in tests).
//!
//! # Representation
//!
//! `grid: [[Point3; 10]; 10]` — row-major, indexed `grid[i][j]`
//! where `i` is the u-parameter index and `j` is the v-parameter
//! index. This matches the other Gregory-family patch types in
//! this crate.

use crate::geometry::gregory::bernstein_binonic::{
    bernstein_binonic_all, bernstein_binonic_deriv_all, bernstein_binonic_second_deriv_all,
};
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Bidegree-(9, 9) Bezier tensor-product patch on `[0, 1]²`.
///
/// 100 plain-Bezier control points; no rational blending, no
/// inner-CP pair placeholders. The K-P guided-surface
/// construction sets every CP explicitly via the formulas in
/// `kp_g1_construction.rs` + `kp_g2_construction.rs`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularGregoryBiNonic {
    /// 10×10 control-point grid, row-major.
    pub grid: [[Point3; 10]; 10],
}

impl RectangularGregoryBiNonic {
    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`.
    ///
    /// `S(u, v) = Σ_{i,j} P[i][j] · B_i^9(u) · B_j^9(v)`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate(&self.grid, &bu, &bv)
    }

    /// `∂S/∂u` at `(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂S/∂v` at `(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²S/∂u²` at `(u, v)` — the u-direction second partial.
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²S/∂u∂v` at `(u, v)` — the mixed second partial.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²S/∂v²` at `(u, v)` — the v-direction second partial.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³S/∂u²∂v` at `(u, v)` — mixed third partial.
    pub fn third_partial_uuv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³S/∂u∂v²` at `(u, v)` — mixed third partial.
    pub fn third_partial_uvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂⁴S/∂u²∂v²` at `(u, v)` — mixed fourth partial.
    pub fn fourth_partial_uuvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// Flat zero patch with every CP at the origin. Useful as a
    /// sentinel / default state before the K-P construction
    /// writes the real CPs. **Not** a valid geometric surface
    /// on its own — every partial derivative is zero.
    pub fn zeros() -> Self {
        Self {
            grid: [[Point3::new(0.0, 0.0, 0.0); 10]; 10],
        }
    }
}

/// Scalar-basis tensor-product accumulate for point evaluation.
/// `out = Σ_{i, j} grid[i][j] · bu[i] · bv[j]` as a `Point3`.
#[inline]
fn accumulate(grid: &[[Point3; 10]; 10], bu: &[f64; 10], bv: &[f64; 10]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        // Inner v-combination first.
        let mut inner = Point3::new(0.0, 0.0, 0.0);
        for j in 0..10 {
            let p = grid[i][j];
            let w = bv[j];
            inner.x += p.x * w;
            inner.y += p.y * w;
            inner.z += p.z * w;
        }
        out.x += inner.x * wu;
        out.y += inner.y * wu;
        out.z += inner.z * wu;
    }
    out
}

/// Same shape as [`accumulate`] but returns a `Vec3` (for
/// derivative evaluations where the output is a direction,
/// not a position).
#[inline]
fn accumulate_vec(grid: &[[Point3; 10]; 10], bu: &[f64; 10], bv: &[f64; 10]) -> Vec3 {
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        let mut inner = Vec3::new(0.0, 0.0, 0.0);
        for j in 0..10 {
            let p = grid[i][j];
            let w = bv[j];
            inner.x += p.x * w;
            inner.y += p.y * w;
            inner.z += p.z * w;
        }
        out.x += inner.x * wu;
        out.y += inner.y * wu;
        out.z += inner.z * wu;
    }
    out
}

#[cfg(test)]
#[path = "rectangular_binonic_tests.rs"]
mod tests;
