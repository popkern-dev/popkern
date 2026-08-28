//! Bicubic Bezier "guide patch" used as the low-degree
//! auxiliary surface in the Karciauskas-Peters guided-surface
//! construction (ADR
//! `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`
//! §2.2, §2.3).
//!
//! The guide provides the cross-boundary derivative field (first
//! partials for G1, second partials for G2) that the bidegree-9
//! [`RectangularGregoryBiNonic`](super::rectangular_binonic::RectangularGregoryBiNonic)
//! inherits by the K-P row-1 / row-2 formulas. The guide is
//! **auxiliary** — it never lives on a face itself; the final
//! surface is always the bidegree-9 patch.
//!
//! # Design
//!
//! Standalone bicubic Bezier, 4×4 control points. Separate from
//! [`RectangularGregory`](super::RectangularGregory) because that
//! type has the 4 corner inner-CP pairs for rational Chiyokura-
//! Kimura blending — the guide is a **plain** Bezier with neither
//! the pairs nor the rational denominator. A plain bicubic is
//! both simpler to construct (no blending-denominator singularity
//! at corners) and simpler to differentiate — essential because
//! the row-2 G2 formula needs second partials of the guide
//! evaluated in closed form.
//!
//! # Derivative-CP identities
//!
//! For a bicubic `G(u, v) = Σ_{i, j} P[i][j] · B_i^3(u) · B_j^3(v)`
//! the cross-derivatives at the `u = 0` boundary factor as
//!
//! ```text
//!   ∂G/∂u(0, v)    = Σ_j [3 · (P[1][j] − P[0][j])] · B_j^3(v)
//!   ∂²G/∂u²(0, v)  = Σ_j [6 · (P[2][j] − 2·P[1][j] + P[0][j])]
//!                       · B_j^3(v)
//! ```
//!
//! i.e. each is a cubic-in-v Bezier curve whose 4 CPs come from a
//! closed-form finite-difference on the guide's `u`-columns. The
//! analogues on `u = 1`, `v = 0`, `v = 1` follow the same
//! template. These are the "`[∂G/∂u(0, v)]_j`" etc. quantities the
//! ADR's Notation paragraph refers to — at cubic degree, before
//! degree-elevation to 9.

use crate::geometry::{Point3, Vec3};

/// Binomial coefficients `C(3, i)` for `i ∈ 0..=3`.
const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];

/// Binomial coefficients `C(2, i)` for `i ∈ 0..=2`, for the first
/// derivative.
const C2: [f64; 3] = [1.0, 2.0, 1.0];

/// Binomial coefficients `C(1, i)` for `i ∈ 0..=1`, for the second
/// derivative.
const C1: [f64; 2] = [1.0, 1.0];

/// Bicubic tensor-product Bezier patch on `[0, 1]²`. 4×4 = 16
/// plain-Bezier control points; no rational blend.
///
/// Used exclusively as the K-P guide surface (H.4 onward); it is
/// not a `SurfaceGeometry` variant and never lives on a face.
#[derive(Debug, Clone)]
pub struct KpGuidePatch {
    /// 4×4 control-point grid, row-major. `grid[i][j]` is the
    /// CP with `u`-index `i` and `v`-index `j`.
    pub grid: [[Point3; 4]; 4],
}

impl KpGuidePatch {
    /// Flat zero guide — all 16 CPs at the origin. Useful as a
    /// sentinel before K-P construction fills it in.
    pub fn zeros() -> Self {
        Self {
            grid: [[Point3::new(0.0, 0.0, 0.0); 4]; 4],
        }
    }

    /// Evaluate `G(u, v)` at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_cubic_all(u);
        let bv = bernstein_cubic_all(v);
        accumulate_point(&self.grid, &bu, &bv)
    }

    /// `∂G/∂u(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_deriv_all(u);
        let bv = bernstein_cubic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂G/∂v(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_all(u);
        let bv = bernstein_cubic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u²(u, v)`.
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_second_deriv_all(u);
        let bv = bernstein_cubic_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂u∂v(u, v)`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_deriv_all(u);
        let bv = bernstein_cubic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂²G/∂v²(u, v)`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_all(u);
        let bv = bernstein_cubic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u²∂v(u, v)` — mixed third partial used by the
    /// corner 3 × 3 block reconciliation in H.4.f.
    pub fn third_partial_uuv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_second_deriv_all(u);
        let bv = bernstein_cubic_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂³G/∂u∂v²(u, v)` — mixed third partial (mirror of
    /// [`third_partial_uuv`](Self::third_partial_uuv)).
    pub fn third_partial_uvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_deriv_all(u);
        let bv = bernstein_cubic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// `∂⁴G/∂u²∂v²(u, v)` — the top-order mixed partial
    /// available on a bicubic guide. Used in the corner 3 × 3
    /// reconciliation.
    pub fn fourth_partial_uuvv(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_cubic_second_deriv_all(u);
        let bv = bernstein_cubic_second_deriv_all(v);
        accumulate_vec(&self.grid, &bu, &bv)
    }

    /// The `u = 0` boundary curve as a cubic Bezier in `v`. Its 4
    /// Bernstein CPs are `grid[0][0..=3]`.
    pub fn boundary_u0(&self) -> [Point3; 4] {
        self.grid[0]
    }

    /// The `u = 1` boundary curve as a cubic Bezier in `v`. Its 4
    /// Bernstein CPs are `grid[3][0..=3]`.
    pub fn boundary_u1(&self) -> [Point3; 4] {
        self.grid[3]
    }

    /// The `v = 0` boundary curve as a cubic Bezier in `u`. Its 4
    /// Bernstein CPs are `grid[0..=3][0]`.
    pub fn boundary_v0(&self) -> [Point3; 4] {
        [
            self.grid[0][0],
            self.grid[1][0],
            self.grid[2][0],
            self.grid[3][0],
        ]
    }

    /// The `v = 1` boundary curve as a cubic Bezier in `u`. Its 4
    /// Bernstein CPs are `grid[0..=3][3]`.
    pub fn boundary_v1(&self) -> [Point3; 4] {
        [
            self.grid[0][3],
            self.grid[1][3],
            self.grid[2][3],
            self.grid[3][3],
        ]
    }

    /// The 4 Bernstein CPs of the cross-derivative curve
    /// `∂G/∂u(0, v)` as a cubic Bezier in `v`. `CP_j = 3 ·
    /// (grid[1][j] − grid[0][j])`.
    pub fn cross_deriv_cps_u0(&self) -> [Vec3; 4] {
        core::array::from_fn(|j| first_diff(self.grid[1][j], self.grid[0][j]))
    }

    /// The 4 Bernstein CPs of `∂G/∂u(1, v)` as a cubic Bezier in
    /// `v`. `CP_j = 3 · (grid[3][j] − grid[2][j])`.
    pub fn cross_deriv_cps_u1(&self) -> [Vec3; 4] {
        core::array::from_fn(|j| first_diff(self.grid[3][j], self.grid[2][j]))
    }

    /// The 4 Bernstein CPs of `∂G/∂v(u, 0)` as a cubic Bezier in
    /// `u`. `CP_i = 3 · (grid[i][1] − grid[i][0])`.
    pub fn cross_deriv_cps_v0(&self) -> [Vec3; 4] {
        core::array::from_fn(|i| first_diff(self.grid[i][1], self.grid[i][0]))
    }

    /// The 4 Bernstein CPs of `∂G/∂v(u, 1)` as a cubic Bezier in
    /// `u`. `CP_i = 3 · (grid[i][3] − grid[i][2])`.
    pub fn cross_deriv_cps_v1(&self) -> [Vec3; 4] {
        core::array::from_fn(|i| first_diff(self.grid[i][3], self.grid[i][2]))
    }

    /// The 4 Bernstein CPs of `∂²G/∂u²(0, v)` as a cubic Bezier in
    /// `v`. `CP_j = 6 · (grid[2][j] − 2·grid[1][j] + grid[0][j])`.
    pub fn cross_second_deriv_cps_uu_u0(&self) -> [Vec3; 4] {
        core::array::from_fn(|j| second_diff(self.grid[2][j], self.grid[1][j], self.grid[0][j]))
    }

    /// The 4 Bernstein CPs of `∂²G/∂u²(1, v)` as a cubic Bezier in
    /// `v`. `CP_j = 6 · (grid[3][j] − 2·grid[2][j] + grid[1][j])`.
    pub fn cross_second_deriv_cps_uu_u1(&self) -> [Vec3; 4] {
        core::array::from_fn(|j| second_diff(self.grid[3][j], self.grid[2][j], self.grid[1][j]))
    }

    /// The 4 Bernstein CPs of `∂²G/∂v²(u, 0)` as a cubic Bezier in
    /// `u`. `CP_i = 6 · (grid[i][2] − 2·grid[i][1] + grid[i][0])`.
    pub fn cross_second_deriv_cps_vv_v0(&self) -> [Vec3; 4] {
        core::array::from_fn(|i| second_diff(self.grid[i][2], self.grid[i][1], self.grid[i][0]))
    }

    /// The 4 Bernstein CPs of `∂²G/∂v²(u, 1)` as a cubic Bezier in
    /// `u`. `CP_i = 6 · (grid[i][3] − 2·grid[i][2] + grid[i][1])`.
    pub fn cross_second_deriv_cps_vv_v1(&self) -> [Vec3; 4] {
        core::array::from_fn(|i| second_diff(self.grid[i][3], self.grid[i][2], self.grid[i][1]))
    }
}

/// `B_i^3(t) = C(3, i) · t^i · (1 − t)^{3 − i}` for `i ∈ 0..=3`.
#[inline]
fn bernstein_cubic_all(t: f64) -> [f64; 4] {
    let omt = 1.0 - t;
    [
        C3[0] * omt * omt * omt,
        C3[1] * t * omt * omt,
        C3[2] * t * t * omt,
        C3[3] * t * t * t,
    ]
}

/// `B_i^3'(t) = 3 · (B_{i - 1}^2(t) − B_i^2(t))` for `i ∈ 0..=3`.
#[inline]
fn bernstein_cubic_deriv_all(t: f64) -> [f64; 4] {
    let omt = 1.0 - t;
    let b2 = [C2[0] * omt * omt, C2[1] * t * omt, C2[2] * t * t];
    [
        3.0 * (0.0 - b2[0]),
        3.0 * (b2[0] - b2[1]),
        3.0 * (b2[1] - b2[2]),
        3.0 * (b2[2] - 0.0),
    ]
}

/// `B_i^3''(t) = 6 · (B_{i - 2}^1(t) − 2·B_{i - 1}^1(t) + B_i^1(t))`
/// for `i ∈ 0..=3`.
#[inline]
fn bernstein_cubic_second_deriv_all(t: f64) -> [f64; 4] {
    let omt = 1.0 - t;
    let b1 = [C1[0] * omt, C1[1] * t];
    let bm = |i: i32| -> f64 {
        if (0..=1).contains(&i) {
            b1[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 4];
    for i in 0..4_i32 {
        out[i as usize] = 6.0 * (bm(i - 2) - 2.0 * bm(i - 1) + bm(i));
    }
    out
}

/// `3 · (a − b)` as a [`Vec3`]. The `3` is the degree-3 Bezier
/// first-derivative-at-boundary factor.
#[inline]
fn first_diff(a: Point3, b: Point3) -> Vec3 {
    Vec3::new(3.0 * (a.x - b.x), 3.0 * (a.y - b.y), 3.0 * (a.z - b.z))
}

/// `6 · (a − 2b + c)` as a [`Vec3`]. The `6` is the degree-3
/// Bezier second-derivative-at-boundary factor `n · (n − 1) = 3 ·
/// 2`.
#[inline]
fn second_diff(a: Point3, b: Point3, c: Point3) -> Vec3 {
    Vec3::new(
        6.0 * (a.x - 2.0 * b.x + c.x),
        6.0 * (a.y - 2.0 * b.y + c.y),
        6.0 * (a.z - 2.0 * b.z + c.z),
    )
}

/// Scalar-basis tensor-product accumulate for point evaluation.
#[inline]
fn accumulate_point(grid: &[[Point3; 4]; 4], bu: &[f64; 4], bv: &[f64; 4]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        let mut inner = Point3::new(0.0, 0.0, 0.0);
        for j in 0..4 {
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

/// Same shape as [`accumulate_point`] but returns `Vec3` for
/// derivative-like values.
#[inline]
fn accumulate_vec(grid: &[[Point3; 4]; 4], bu: &[f64; 4], bv: &[f64; 4]) -> Vec3 {
    let mut out = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        let wu = bu[i];
        if wu == 0.0 {
            continue;
        }
        let mut inner = Vec3::new(0.0, 0.0, 0.0);
        for j in 0..4 {
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
#[path = "kp_guide_patch_tests.rs"]
mod tests;
