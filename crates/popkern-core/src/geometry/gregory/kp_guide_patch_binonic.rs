//! Phase H.8 (Option C, Axis 1 ceiling) — bi-nonic Bezier guide
//! patch. 10×10 = 100 plain-Bezier CPs. Matches the output
//! bidegree-9 patch's native degree, giving the K-P construction
//! maximum expressiveness within the current output.
//!
//! # The equal-degree degeneracy
//!
//! When the guide has the same degree as the output, the K-P 2004
//! row-1 / row-2 formulas become identity operations under matching
//! boundaries:
//!
//! ```text
//!   row-1[j] = boundary[0][j] + (1/9) · [∂G/∂u(0, v)]_j^9
//! ```
//!
//! For a bi-nonic guide, the boundary cross-derivative at u=0 is
//! already in its native degree-9 form: `9·(P_G[1][j] − P_G[0][j])`.
//! Substituting:
//!
//! ```text
//!   row-1[j] = boundary[0][j] + (1/9) · 9 · (P_G[1][j] − P_G[0][j])
//!            = boundary[0][j] + (P_G[1][j] − P_G[0][j])
//! ```
//!
//! When the boundary matches the guide at the edge (`boundary[0][j]
//! = P_G[0][j]`), the output's row-1 equals the guide's row-1
//! exactly. The K-P guided-surface framework collapses to direct
//! construction of the output patch. The distinction between
//! "guide" and "output" becomes purely a code-organisation concern.
//!
//! This isn't a bug — it's the framework's expressive ceiling.
//! K-P 2004 formalises the guide-surface idea with the factor
//! `1/n` where `n` is the guide degree; making guide degree equal
//! output degree removes all information loss.
//!
//! # Why keep guide and output as distinct types
//!
//! Semantic clarity. `RectangularGregoryBiNonic` lives on a BRep
//! face; `KpGuidePatchBiNonic` is auxiliary scaffolding that
//! never reaches user-facing surfaces. Reading a function
//! signature `fn derive_row1<G: KpGuide>(guide: &G, ...)` where
//! `G` is the same concrete type as what we're building would be
//! confusing. The types encode intent; minor struct-storage
//! redundancy is a fair price.

use super::kp_guide_patch::KpGuidePatch;
use super::kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
use crate::geometry::gregory::bernstein_binonic::{
    bernstein_binonic_all, bernstein_binonic_deriv_all, bernstein_binonic_second_deriv_all,
};
use crate::geometry::{Point3, Vec3};

/// Bi-nonic tensor-product Bezier guide patch on `[0, 1]²`.
/// 10×10 = 100 plain-Bezier control points.
///
/// Equal-degree cousin of [`RectangularGregoryBiNonic`] with
/// distinct semantics — guide, not output. Never lives on a
/// face; only provides cross-derivative data to the K-P row-1 /
/// row-2 formulas where it collapses to direct CP construction
/// (see module doc).
#[derive(Debug, Clone)]
pub struct KpGuidePatchBiNonic {
    /// 10×10 control-point grid, row-major. `grid[i][j]` is the
    /// CP with `u`-index `i` and `v`-index `j`.
    pub grid: [[Point3; 10]; 10],
}

impl KpGuidePatchBiNonic {
    /// Flat zero guide — all 100 CPs at the origin.
    pub fn zeros() -> Self {
        Self {
            grid: [[Point3::new(0.0, 0.0, 0.0); 10]; 10],
        }
    }

    /// Evaluate `G(u, v)` at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate_point(&self.grid, &bu, &bv)
    }

    /// `∂G/∂u(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂G/∂v(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u²(u, v)`.
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u∂v(u, v)`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂v²(u, v)`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u²∂v(u, v)` — mixed third partial for G2 corner
    /// reconciliation.
    pub fn third_partial_uuv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u∂v²(u, v)` — mixed third partial.
    pub fn third_partial_uvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂⁴G/∂u²∂v²(u, v)` — fourth-order mixed partial.
    pub fn fourth_partial_uuvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_second_deriv_all(u);
        let bv = bernstein_binonic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// The `u = 0` boundary curve as a degree-9 Bezier in `v`.
    pub fn boundary_u0(&self) -> [Point3; 10] {
        self.grid[0]
    }

    /// The `u = 1` boundary curve as a degree-9 Bezier in `v`.
    pub fn boundary_u1(&self) -> [Point3; 10] {
        self.grid[9]
    }

    /// The `v = 0` boundary curve as a degree-9 Bezier in `u`.
    pub fn boundary_v0(&self) -> [Point3; 10] {
        core::array::from_fn(|i| self.grid[i][0])
    }

    /// The `v = 1` boundary curve as a degree-9 Bezier in `u`.
    pub fn boundary_v1(&self) -> [Point3; 10] {
        core::array::from_fn(|i| self.grid[i][9])
    }

    /// Degree-9 CPs of `∂G/∂u(0, v)`: `CP_j = 9 · (P[1][j] − P[0][j])`.
    /// Already at degree 9 — no elevation needed for the K-P row-1
    /// formula.
    pub fn cross_deriv_cps_u0(&self) -> [Vec3; 10] {
        core::array::from_fn(|j| first_diff9(self.grid[1][j], self.grid[0][j]))
    }

    /// `∂G/∂u(1, v)` CPs.
    pub fn cross_deriv_cps_u1(&self) -> [Vec3; 10] {
        core::array::from_fn(|j| first_diff9(self.grid[9][j], self.grid[8][j]))
    }

    /// `∂G/∂v(u, 0)` CPs.
    pub fn cross_deriv_cps_v0(&self) -> [Vec3; 10] {
        core::array::from_fn(|i| first_diff9(self.grid[i][1], self.grid[i][0]))
    }

    /// `∂G/∂v(u, 1)` CPs.
    pub fn cross_deriv_cps_v1(&self) -> [Vec3; 10] {
        core::array::from_fn(|i| first_diff9(self.grid[i][9], self.grid[i][8]))
    }

    /// Degree-9 CPs of `∂²G/∂u²(0, v)`: `CP_j = 72 · (P[2][j] −
    /// 2·P[1][j] + P[0][j])` (factor `n·(n-1) = 9·8 = 72`).
    pub fn cross_second_deriv_cps_uu_u0(&self) -> [Vec3; 10] {
        core::array::from_fn(|j| second_diff9(self.grid[2][j], self.grid[1][j], self.grid[0][j]))
    }

    /// `∂²G/∂u²(1, v)` CPs.
    pub fn cross_second_deriv_cps_uu_u1(&self) -> [Vec3; 10] {
        core::array::from_fn(|j| second_diff9(self.grid[9][j], self.grid[8][j], self.grid[7][j]))
    }

    /// `∂²G/∂v²(u, 0)` CPs.
    pub fn cross_second_deriv_cps_vv_v0(&self) -> [Vec3; 10] {
        core::array::from_fn(|i| second_diff9(self.grid[i][2], self.grid[i][1], self.grid[i][0]))
    }

    /// `∂²G/∂v²(u, 1)` CPs.
    pub fn cross_second_deriv_cps_vv_v1(&self) -> [Vec3; 10] {
        core::array::from_fn(|i| second_diff9(self.grid[i][9], self.grid[i][8], self.grid[i][7]))
    }

    /// Degree-elevate a bicubic guide (4×4 CPs) to bi-nonic. 3 → 9
    /// elevation steps in each parametric direction. Pointwise-
    /// identical to the source bicubic.
    #[allow(clippy::needless_range_loop)]
    pub fn from_bicubic_elevated(bicubic: &KpGuidePatch) -> Self {
        use crate::geometry::bernstein_poly::BernsteinPoly;
        // Helper: elevate a cubic (4 CPs) to degree 9 (10 CPs)
        // coordinate-wise.
        let elevate = |col: [Point3; 4]| -> [Point3; 10] {
            let xs =
                BernsteinPoly::new(col.iter().map(|p| p.x).collect::<Vec<_>>()).degree_elevate(6);
            let ys =
                BernsteinPoly::new(col.iter().map(|p| p.y).collect::<Vec<_>>()).degree_elevate(6);
            let zs =
                BernsteinPoly::new(col.iter().map(|p| p.z).collect::<Vec<_>>()).degree_elevate(6);
            core::array::from_fn(|i| Point3::new(xs.coeffs()[i], ys.coeffs()[i], zs.coeffs()[i]))
        };
        // u-direction first: 4 v-columns → each elevated 3→9.
        let mut ten_by_four = [[Point3::new(0.0, 0.0, 0.0); 4]; 10];
        for j in 0..4 {
            let col: [Point3; 4] = core::array::from_fn(|i| bicubic.grid[i][j]);
            let elevated = elevate(col);
            for (i, p) in elevated.iter().enumerate() {
                ten_by_four[i][j] = *p;
            }
        }
        // v-direction: each of 10 u-rows → elevated 3→9.
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for i in 0..10 {
            let row: [Point3; 4] = core::array::from_fn(|j| ten_by_four[i][j]);
            grid[i] = elevate(row);
        }
        Self { grid }
    }

    /// Degree-elevate a bi-quintic guide (6×6 CPs) to bi-nonic.
    /// 5 → 9 elevation in each direction. Pointwise-identical to
    /// the source bi-quintic.
    #[allow(clippy::needless_range_loop)]
    pub fn from_biquintic_elevated(bq: &KpGuidePatchBiQuintic) -> Self {
        use crate::geometry::bernstein_poly::BernsteinPoly;
        let elevate = |col: [Point3; 6]| -> [Point3; 10] {
            let xs =
                BernsteinPoly::new(col.iter().map(|p| p.x).collect::<Vec<_>>()).degree_elevate(4);
            let ys =
                BernsteinPoly::new(col.iter().map(|p| p.y).collect::<Vec<_>>()).degree_elevate(4);
            let zs =
                BernsteinPoly::new(col.iter().map(|p| p.z).collect::<Vec<_>>()).degree_elevate(4);
            core::array::from_fn(|i| Point3::new(xs.coeffs()[i], ys.coeffs()[i], zs.coeffs()[i]))
        };
        let mut ten_by_six = [[Point3::new(0.0, 0.0, 0.0); 6]; 10];
        for j in 0..6 {
            let col: [Point3; 6] = core::array::from_fn(|i| bq.grid[i][j]);
            let elevated = elevate(col);
            for (i, p) in elevated.iter().enumerate() {
                ten_by_six[i][j] = *p;
            }
        }
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for i in 0..10 {
            let row: [Point3; 6] = core::array::from_fn(|j| ten_by_six[i][j]);
            grid[i] = elevate(row);
        }
        Self { grid }
    }
}

/// `9 · (a − b)` — bi-nonic's first-derivative-at-boundary factor.
#[inline]
fn first_diff9(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(9.0 * (a.x - b.x), 9.0 * (a.y - b.y), 9.0 * (a.z - b.z))
}

/// `72 · (a − 2b + c)` — bi-nonic's second-derivative-at-boundary
/// factor (`n·(n-1) = 9·8`).
#[inline]
fn second_diff9(a: Point3, b: Point3, c: Point3) -> Vec3 {
    Vec3::new(
        72.0 * (a.x - 2.0 * b.x + c.x),
        72.0 * (a.y - 2.0 * b.y + c.y),
        72.0 * (a.z - 2.0 * b.z + c.z),
    )
}

#[inline]
fn accumulate_point(grid: &[[Point3; 10]; 10], bu: &[f64; 10], bv: &[f64; 10]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
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
#[path = "kp_guide_patch_binonic_tests.rs"]
mod tests;
