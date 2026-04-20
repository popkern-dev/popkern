//! Phase H.8 (Option C) — bi-quintic Bezier guide patch.
//!
//! 5×5 = 36 plain-Bezier CPs. Drop-in richer substrate for the
//! K-P 2004 guided-surface construction, closing the flexibility
//! gap the bicubic 16-CP guide (see `kp_guide_patch.rs`) hits on
//! extreme radius asymmetry (cube_extreme, 1:10 ratios). See
//! `docs/journal/wip.md` §Option C for the full rationale.
//!
//! # Flexibility over bicubic
//!
//! A bicubic guide's boundary cross-derivative is a cubic-in-v
//! (4 CPs per edge, degree-3 polynomial). A bi-quintic guide's
//! boundary cross-derivative is a quartic-in-v (5 CPs, degree-4
//! polynomial). One extra DoF per boundary per direction — enough
//! for the row-1 frame to represent asymmetric cross-curvature
//! variation the bicubic cannot.
//!
//! # Corner layout mirrors the bicubic
//!
//! ```text
//!   grid[0][0] = c0 (u=0, v=0)   grid[0][5] = c3 (u=0, v=1)
//!   grid[5][0] = c1 (u=1, v=0)   grid[5][5] = c2 (u=1, v=1)
//! ```
//!
//! # Derivative-CP identities
//!
//! For `G(u, v) = Σ_{i, j} P[i][j] · B_i^5(u) · B_j^5(v)` the
//! cross-derivatives at the `u = 0` boundary are
//!
//! ```text
//!   ∂G/∂u(0, v)   = Σ_j [5 · (P[1][j] − P[0][j])] · B_j^5(v)
//!   ∂²G/∂u²(0, v) = Σ_j [20 · (P[2][j] − 2·P[1][j] + P[0][j])]
//!                   · B_j^5(v)
//! ```
//!
//! i.e. each is a **quartic**-in-v Bezier curve (5 CPs). The
//! row-1 / row-2 frame construction degree-elevates these from
//! 4 to 9 (factor 5/9) — vs the bicubic's 3→9 (factor 3/9).

use super::kp_guide_patch::KpGuidePatch;
use crate::geometry::{Point3, Vec3};

/// Binomial coefficients `C(5, i)` for `i ∈ 0..=5`.
const C5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];

/// Binomial coefficients `C(4, i)` for `i ∈ 0..=4`.
const C4: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];

/// Binomial coefficients `C(3, i)` for `i ∈ 0..=3`.
const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];

/// Bi-quintic tensor-product Bezier patch on `[0, 1]²`. 6×6 = 36
/// plain-Bezier control points; no rational blend.
#[derive(Debug, Clone)]
pub struct KpGuidePatchBiQuintic {
    /// 6×6 control-point grid, row-major. `grid[i][j]` is the CP
    /// with `u`-index `i` and `v`-index `j`.
    pub grid: [[Point3; 6]; 6],
}

impl KpGuidePatchBiQuintic {
    /// Flat zero guide — all 36 CPs at the origin.
    pub fn zeros() -> Self {
        Self {
            grid: [[Point3::new(0.0, 0.0, 0.0); 6]; 6],
        }
    }

    /// Evaluate `G(u, v)` at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein5_all(u);
        let bv = bernstein5_all(v);
        accumulate_point(&self.grid, &bu, &bv)
    }

    /// `∂G/∂u(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_deriv_all(u);
        let bv = bernstein5_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂G/∂v(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_all(u);
        let bv = bernstein5_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u²(u, v)`.
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_second_deriv_all(u);
        let bv = bernstein5_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u∂v(u, v)`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_deriv_all(u);
        let bv = bernstein5_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂v²(u, v)`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_all(u);
        let bv = bernstein5_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u²∂v(u, v)` — mixed third partial for G2 corner-block
    /// reconciliation.
    pub fn third_partial_uuv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_second_deriv_all(u);
        let bv = bernstein5_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u∂v²(u, v)` — mixed third partial (mirror of uuv).
    pub fn third_partial_uvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_deriv_all(u);
        let bv = bernstein5_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂⁴G/∂u²∂v²(u, v)` — 4th-order mixed partial (top for G2).
    pub fn fourth_partial_uuvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein5_second_deriv_all(u);
        let bv = bernstein5_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// The `u = 0` boundary curve as a quintic Bezier in `v`. CPs
    /// are `grid[0][0..=5]`.
    pub fn boundary_u0(&self) -> [Point3; 6] {
        self.grid[0]
    }

    /// The `u = 1` boundary curve as a quintic Bezier in `v`. CPs
    /// are `grid[5][0..=5]`.
    pub fn boundary_u1(&self) -> [Point3; 6] {
        self.grid[5]
    }

    /// The `v = 0` boundary curve as a quintic Bezier in `u`. CPs
    /// are `grid[0..=5][0]`.
    pub fn boundary_v0(&self) -> [Point3; 6] {
        core::array::from_fn(|i| self.grid[i][0])
    }

    /// The `v = 1` boundary curve as a quintic Bezier in `u`. CPs
    /// are `grid[0..=5][5]`.
    pub fn boundary_v1(&self) -> [Point3; 6] {
        core::array::from_fn(|i| self.grid[i][5])
    }

    /// Quartic-Bernstein CPs of `∂G/∂u(0, v)`: `CP_j = 5 · (P[1][j]
    /// − P[0][j])` for `j ∈ 0..=5`.
    pub fn cross_deriv_cps_u0(&self) -> [Vec3; 6] {
        core::array::from_fn(|j| first_diff5(self.grid[1][j], self.grid[0][j]))
    }

    /// Quartic-Bernstein CPs of `∂G/∂u(1, v)`: `CP_j = 5 · (P[5][j]
    /// − P[4][j])`.
    pub fn cross_deriv_cps_u1(&self) -> [Vec3; 6] {
        core::array::from_fn(|j| first_diff5(self.grid[5][j], self.grid[4][j]))
    }

    /// Quartic-Bernstein CPs of `∂G/∂v(u, 0)`.
    pub fn cross_deriv_cps_v0(&self) -> [Vec3; 6] {
        core::array::from_fn(|i| first_diff5(self.grid[i][1], self.grid[i][0]))
    }

    /// Quartic-Bernstein CPs of `∂G/∂v(u, 1)`.
    pub fn cross_deriv_cps_v1(&self) -> [Vec3; 6] {
        core::array::from_fn(|i| first_diff5(self.grid[i][5], self.grid[i][4]))
    }

    /// Cubic-Bernstein CPs of `∂²G/∂u²(0, v)`: `CP_j = 20 · (P[2][j]
    /// − 2·P[1][j] + P[0][j])` for `j ∈ 0..=5`. Note: while the
    /// curve has 6 CPs (degree-5 like the patch), the second-
    /// derivative factor is `n·(n-1) = 5·4 = 20`.
    pub fn cross_second_deriv_cps_uu_u0(&self) -> [Vec3; 6] {
        core::array::from_fn(|j| second_diff5(self.grid[2][j], self.grid[1][j], self.grid[0][j]))
    }

    /// `∂²G/∂u²(1, v)` CPs.
    pub fn cross_second_deriv_cps_uu_u1(&self) -> [Vec3; 6] {
        core::array::from_fn(|j| second_diff5(self.grid[5][j], self.grid[4][j], self.grid[3][j]))
    }

    /// `∂²G/∂v²(u, 0)` CPs.
    pub fn cross_second_deriv_cps_vv_v0(&self) -> [Vec3; 6] {
        core::array::from_fn(|i| second_diff5(self.grid[i][2], self.grid[i][1], self.grid[i][0]))
    }

    /// `∂²G/∂v²(u, 1)` CPs.
    pub fn cross_second_deriv_cps_vv_v1(&self) -> [Vec3; 6] {
        core::array::from_fn(|i| second_diff5(self.grid[i][5], self.grid[i][4], self.grid[i][3]))
    }

    /// Build a bi-quintic guide by degree-elevating a bicubic guide.
    ///
    /// Bezier degree elevation: a bicubic (degree 3×3) patch with
    /// 4×4 CPs can be represented exactly as a degree-5×5 patch
    /// with 6×6 CPs by applying the standard elevation formula
    /// `Q_i = (i/(n+1)) · P_{i-1} + (1 − i/(n+1)) · P_i` in both
    /// parametric directions (twice: 3→4 then 4→5). The resulting
    /// bi-quintic evaluates **identically** — pointwise, per-partial
    /// — to the source bicubic. This is the "initial state" for the
    /// iterative K-P 2004 §4 optimizer: starts from the current
    /// bicubic output (no regression baseline) and the extra DoFs
    /// (CPs that were previously determined by elevation) become
    /// the optimizer's handles.
    #[allow(clippy::needless_range_loop)]
    pub fn from_bicubic_elevated(bicubic: &KpGuidePatch) -> Self {
        // u-direction elevation first (4 u-CPs → 6 u-CPs for each
        // of the 4 v-columns), then v-direction elevation on the
        // 6 resulting u-rows (4 v-CPs → 6 v-CPs).
        let mut six_by_four = [[Point3::new(0.0, 0.0, 0.0); 4]; 6];
        for j in 0..4 {
            let col: [Point3; 4] = core::array::from_fn(|i| bicubic.grid[i][j]);
            let elevated = elevate_cubic_column_to_degree_five(col);
            for (i, e) in elevated.iter().enumerate() {
                six_by_four[i][j] = *e;
            }
        }
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
        for i in 0..6 {
            let row: [Point3; 4] = core::array::from_fn(|j| six_by_four[i][j]);
            let elevated = elevate_cubic_column_to_degree_five(row);
            grid[i] = elevated;
        }
        Self { grid }
    }
}

/// Degree elevation on a single column of 4 CPs: 3 → 4 → 5 via
/// the standard Bezier formula applied twice.
fn elevate_cubic_column_to_degree_five(cps: [Point3; 4]) -> [Point3; 6] {
    // Step 3 → 4.
    let mut four = [Point3::new(0.0, 0.0, 0.0); 5];
    four[0] = cps[0];
    four[4] = cps[3];
    for i in 1..4 {
        let a = i as f64 / 4.0;
        four[i] = Point3::new(
            a * cps[i - 1].x + (1.0 - a) * cps[i].x,
            a * cps[i - 1].y + (1.0 - a) * cps[i].y,
            a * cps[i - 1].z + (1.0 - a) * cps[i].z,
        );
    }
    // Step 4 → 5.
    let mut five = [Point3::new(0.0, 0.0, 0.0); 6];
    five[0] = four[0];
    five[5] = four[4];
    for i in 1..5 {
        let a = i as f64 / 5.0;
        five[i] = Point3::new(
            a * four[i - 1].x + (1.0 - a) * four[i].x,
            a * four[i - 1].y + (1.0 - a) * four[i].y,
            a * four[i - 1].z + (1.0 - a) * four[i].z,
        );
    }
    five
}

/// `B_i^5(t) = C(5, i) · t^i · (1 − t)^{5 − i}` for `i ∈ 0..=5`.
#[inline]
fn bernstein5_all(t: f64) -> [f64; 6] {
    let omt = 1.0 - t;
    let omt2 = omt * omt;
    let omt3 = omt2 * omt;
    let omt4 = omt3 * omt;
    let omt5 = omt4 * omt;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    [
        C5[0] * omt5,
        C5[1] * t * omt4,
        C5[2] * t2 * omt3,
        C5[3] * t3 * omt2,
        C5[4] * t4 * omt,
        C5[5] * t5,
    ]
}

/// `B_i^5'(t) = 5 · (B_{i-1}^4(t) − B_i^4(t))` for `i ∈ 0..=5`.
#[inline]
fn bernstein5_deriv_all(t: f64) -> [f64; 6] {
    let b4 = bernstein4_all(t);
    [
        5.0 * (0.0 - b4[0]),
        5.0 * (b4[0] - b4[1]),
        5.0 * (b4[1] - b4[2]),
        5.0 * (b4[2] - b4[3]),
        5.0 * (b4[3] - b4[4]),
        5.0 * (b4[4] - 0.0),
    ]
}

/// `B_i^5''(t) = 20 · (B_{i-2}^3(t) − 2·B_{i-1}^3(t) + B_i^3(t))`.
#[inline]
fn bernstein5_second_deriv_all(t: f64) -> [f64; 6] {
    let b3 = bernstein3_all(t);
    let bm = |i: i32| -> f64 {
        if (0..=3).contains(&i) {
            b3[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 6];
    for i in 0..6_i32 {
        out[i as usize] = 20.0 * (bm(i - 2) - 2.0 * bm(i - 1) + bm(i));
    }
    out
}

#[inline]
fn bernstein4_all(t: f64) -> [f64; 5] {
    let omt = 1.0 - t;
    let omt2 = omt * omt;
    let omt3 = omt2 * omt;
    let omt4 = omt3 * omt;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    [
        C4[0] * omt4,
        C4[1] * t * omt3,
        C4[2] * t2 * omt2,
        C4[3] * t3 * omt,
        C4[4] * t4,
    ]
}

#[inline]
fn bernstein3_all(t: f64) -> [f64; 4] {
    let omt = 1.0 - t;
    let omt2 = omt * omt;
    let omt3 = omt2 * omt;
    let t2 = t * t;
    let t3 = t2 * t;
    [C3[0] * omt3, C3[1] * t * omt2, C3[2] * t2 * omt, C3[3] * t3]
}

/// `5 · (a − b)` as a [`Vec3`] — bi-quintic's first-derivative-at-
/// boundary factor.
#[inline]
fn first_diff5(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(5.0 * (a.x - b.x), 5.0 * (a.y - b.y), 5.0 * (a.z - b.z))
}

/// `20 · (a − 2b + c)` as a [`Vec3`] — bi-quintic's second-
/// derivative-at-boundary factor (`n · (n − 1) = 5 · 4`).
#[inline]
fn second_diff5(a: Point3, b: Point3, c: Point3) -> Vec3 {
    Vec3::new(
        20.0 * (a.x - 2.0 * b.x + c.x),
        20.0 * (a.y - 2.0 * b.y + c.y),
        20.0 * (a.z - 2.0 * b.z + c.z),
    )
}

/// Scalar-basis tensor-product accumulate for point evaluation.
#[inline]
fn accumulate_point(grid: &[[Point3; 6]; 6], bu: &[f64; 6], bv: &[f64; 6]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        let mut inner = Point3::new(0.0, 0.0, 0.0);
        for j in 0..6 {
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

/// Same as [`accumulate_point`] but returns `Vec3` (for derivatives).
#[inline]
fn accumulate_vec(grid: &[[Point3; 6]; 6], bu: &[f64; 6], bv: &[f64; 6]) -> Vec3 {
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        let mut inner = Vec3::new(0.0, 0.0, 0.0);
        for j in 0..6 {
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
#[path = "kp_guide_patch_biquintic_tests.rs"]
mod tests;
