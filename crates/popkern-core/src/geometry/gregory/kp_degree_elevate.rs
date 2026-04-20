//! Phase H — degree-elevation helpers for the Karciauskas-Peters
//! bidegree-9 construction.
//!
//! The K-P row-1 / row-2 formulas (ADR
//! `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`
//! §2.3) consume **Bernstein control points of the guide's
//! cross-derivative curves expressed at degree 9**, not at the
//! guide's native degree 3. This module converts 4-CP cubic
//! curves to 10-CP degree-9 curves.
//!
//! # Why degree elevation matters (and direct point-sampling is
//! wrong)
//!
//! The row-1 formula is
//!
//! ```text
//!   P[1][j] = P[0][j] + (1/9) · [∂G/∂u(0, v)]_j^9
//! ```
//!
//! where `[F(v)]_j^9` denotes the j-th Bernstein CP of `F(v)` **as
//! a degree-9 Bezier curve**. `F(v) = ∂G/∂u(0, v)` is naturally
//! cubic in `v` (since the guide `G` is bicubic). To index the
//! j-th CP at degree 9, we must degree-elevate the cubic from 4
//! CPs to 10 CPs. Bernstein CPs are ≠ point samples except at the
//! two endpoints, so `[F(v)]_j^9 ≠ F(j / 9)` in general — every
//! K-P implementation must go through degree elevation, not point
//! evaluation.
//!
//! # Delegation
//!
//! For [`Vec3`] CPs we delegate to the generic
//! [`BernsteinPoly::degree_elevate`] which applies the standard
//! iterated-elevation formula
//!
//! ```text
//!   c'_0 = c_0
//!   c'_i = (i / (n + 1)) · c_{i - 1} + (1 − i / (n + 1)) · c_i
//!                                              for i ∈ 1..n
//!   c'_{n + 1} = c_n
//! ```
//!
//! six times to go from `n = 3` to `n = 9`. This is numerically
//! stable in `[0, 1]` and exact up to rounding. For [`Point3`]
//! CPs we perform coordinate-wise `f64` elevation, since `Point3`
//! is semantically a position and does not implement
//! [`PolyCoeff`](crate::geometry::bernstein_poly) (only [`Vec3`]
//! and `f64` do).

use crate::geometry::bernstein_poly::BernsteinPoly;
use crate::geometry::{Point3, Vec3};

/// Degree-elevate a cubic Bezier curve (4 [`Vec3`] CPs) to a
/// degree-9 Bezier curve (10 [`Vec3`] CPs).
///
/// This is the canonical helper for the K-P row-1 and row-2
/// formulas when the guide is bicubic, which consume the guide's
/// cross-derivative / cross-second-derivative CPs at degree 9.
pub fn elevate_cubic_vec_to_nine(cubic: [Vec3; 4]) -> [Vec3; 10] {
    let poly = BernsteinPoly::new(cubic.to_vec()).degree_elevate(6);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Degree-elevate a quintic Bezier curve (6 [`Vec3`] CPs) to a
/// degree-9 Bezier curve (10 [`Vec3`] CPs).
///
/// Used for the row-1 / row-2 formulas when the guide is
/// bi-quintic (Phase H.8 Option C). The bi-quintic guide's
/// boundary cross-derivative is a degree-5 Bezier (6 CPs),
/// so it only needs to be elevated by 4 degrees to reach 9.
pub fn elevate_quintic_vec_to_nine(quintic: [Vec3; 6]) -> [Vec3; 10] {
    let poly = BernsteinPoly::new(quintic.to_vec()).degree_elevate(4);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Degree-elevate a cubic Bezier curve (4 [`Point3`] CPs) to a
/// degree-9 Bezier curve (10 [`Point3`] CPs). Coordinate-wise
/// `f64` elevation.
pub fn elevate_cubic_point_to_nine(cubic: [Point3; 4]) -> [Point3; 10] {
    let xs = elevate_cubic_f64_to_nine([cubic[0].x, cubic[1].x, cubic[2].x, cubic[3].x]);
    let ys = elevate_cubic_f64_to_nine([cubic[0].y, cubic[1].y, cubic[2].y, cubic[3].y]);
    let zs = elevate_cubic_f64_to_nine([cubic[0].z, cubic[1].z, cubic[2].z, cubic[3].z]);
    core::array::from_fn(|i| Point3::new(xs[i], ys[i], zs[i]))
}

/// Degree-elevate an `N`-CP [`Point3`] Bezier curve to a 10-CP
/// degree-9 Bezier curve. Used for converting input boundary
/// curves (typically degree 7 from the biseptic promotion
/// pipeline) into the `row = 0` / `col = 0` CPs of the bidegree-9
/// patch.
///
/// Returns `None` if `cps.is_empty()` or `cps.len() > 10`
/// (requested degree exceeds target).
pub fn elevate_point_curve_to_nine(cps: &[Point3]) -> Option<[Point3; 10]> {
    let n = cps.len();
    if n == 0 || n > 10 {
        return None;
    }
    let by = 10 - n;
    let xs: Vec<f64> = cps.iter().map(|p| p.x).collect();
    let ys: Vec<f64> = cps.iter().map(|p| p.y).collect();
    let zs: Vec<f64> = cps.iter().map(|p| p.z).collect();
    let xs9 = BernsteinPoly::new(xs).degree_elevate(by);
    let ys9 = BernsteinPoly::new(ys).degree_elevate(by);
    let zs9 = BernsteinPoly::new(zs).degree_elevate(by);
    Some(core::array::from_fn(|i| {
        Point3::new(xs9.coeffs()[i], ys9.coeffs()[i], zs9.coeffs()[i])
    }))
}

/// Degree-elevate a cubic (4-CP) `f64` Bezier curve to degree-9
/// (10-CP). Private helper for [`elevate_cubic_point_to_nine`].
fn elevate_cubic_f64_to_nine(cubic: [f64; 4]) -> [f64; 10] {
    let poly = BernsteinPoly::new(cubic.to_vec()).degree_elevate(6);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Phase I — degree-elevate a degree-6 Bezier curve (7 [`Vec3`]
/// CPs) to a degree-9 Bezier curve (10 CPs). Used by the K-P 2015
/// bi-6 polar construction to produce bi-9 output uniform with the
/// Phase H rectangular path.
pub fn elevate_sextic_vec_to_nine(sextic: [Vec3; 7]) -> [Vec3; 10] {
    let poly = BernsteinPoly::new(sextic.to_vec()).degree_elevate(3);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Phase I — degree-elevate a degree-6 Bezier curve (7 [`Point3`]
/// CPs) to a degree-9 Bezier curve (10 CPs), coordinate-wise.
pub fn elevate_sextic_point_to_nine(sextic: [Point3; 7]) -> [Point3; 10] {
    let xs = elevate_sextic_f64_to_nine([
        sextic[0].x,
        sextic[1].x,
        sextic[2].x,
        sextic[3].x,
        sextic[4].x,
        sextic[5].x,
        sextic[6].x,
    ]);
    let ys = elevate_sextic_f64_to_nine([
        sextic[0].y,
        sextic[1].y,
        sextic[2].y,
        sextic[3].y,
        sextic[4].y,
        sextic[5].y,
        sextic[6].y,
    ]);
    let zs = elevate_sextic_f64_to_nine([
        sextic[0].z,
        sextic[1].z,
        sextic[2].z,
        sextic[3].z,
        sextic[4].z,
        sextic[5].z,
        sextic[6].z,
    ]);
    core::array::from_fn(|i| Point3::new(xs[i], ys[i], zs[i]))
}

fn elevate_sextic_f64_to_nine(sextic: [f64; 7]) -> [f64; 10] {
    let poly = BernsteinPoly::new(sextic.to_vec()).degree_elevate(3);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Phase I — degree-elevate a degree-5 Bezier curve (6 [`Point3`]
/// CPs) to a degree-9 Bezier curve (10 CPs), coordinate-wise.
/// Counterpart to the [`Vec3`] version
/// [`elevate_quintic_vec_to_nine`].
pub fn elevate_quintic_point_to_nine(quintic: [Point3; 6]) -> [Point3; 10] {
    let xs = elevate_quintic_f64_to_nine([
        quintic[0].x,
        quintic[1].x,
        quintic[2].x,
        quintic[3].x,
        quintic[4].x,
        quintic[5].x,
    ]);
    let ys = elevate_quintic_f64_to_nine([
        quintic[0].y,
        quintic[1].y,
        quintic[2].y,
        quintic[3].y,
        quintic[4].y,
        quintic[5].y,
    ]);
    let zs = elevate_quintic_f64_to_nine([
        quintic[0].z,
        quintic[1].z,
        quintic[2].z,
        quintic[3].z,
        quintic[4].z,
        quintic[5].z,
    ]);
    core::array::from_fn(|i| Point3::new(xs[i], ys[i], zs[i]))
}

fn elevate_quintic_f64_to_nine(quintic: [f64; 6]) -> [f64; 10] {
    let poly = BernsteinPoly::new(quintic.to_vec()).degree_elevate(4);
    let coeffs = poly.coeffs();
    core::array::from_fn(|i| coeffs[i])
}

/// Phase I — degree-elevate a 6×6 [`Point3`] grid (bi-5 tensor-
/// product) to a 10×10 grid (bi-9) by applying the univariate
/// elevation to every row, then every column of the intermediate
/// result.
///
/// Used by the K-P 2013 biquintic polar construction after the per-
/// sector macro-patches are assembled, to bring them onto the
/// rectangular-bi-9 patch type shared with the Phase H output.
pub fn elevate_bi5_grid_to_bi9(grid: [[Point3; 6]; 6]) -> [[Point3; 10]; 10] {
    // Pass 1: elevate each row (varying j, fixed i) from 6 CPs to 10.
    let mut intermediate: [[Point3; 10]; 6] = [[Point3::new(0.0, 0.0, 0.0); 10]; 6];
    for i in 0..6 {
        intermediate[i] = elevate_quintic_point_to_nine(grid[i]);
    }
    // Pass 2: for each j ∈ 0..10, elevate the column (varying i)
    // from 6 to 10 CPs.
    let mut out: [[Point3; 10]; 10] = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for j in 0..10 {
        let col: [Point3; 6] = core::array::from_fn(|i| intermediate[i][j]);
        let col9 = elevate_quintic_point_to_nine(col);
        for i in 0..10 {
            out[i][j] = col9[i];
        }
    }
    out
}

/// Phase I — degree-elevate a 7×7 [`Point3`] grid (bi-6 tensor-
/// product) to a 10×10 grid (bi-9). Row-then-column pass, same
/// pattern as [`elevate_bi5_grid_to_bi9`].
pub fn elevate_bi6_grid_to_bi9(grid: [[Point3; 7]; 7]) -> [[Point3; 10]; 10] {
    let mut intermediate: [[Point3; 10]; 7] = [[Point3::new(0.0, 0.0, 0.0); 10]; 7];
    for i in 0..7 {
        intermediate[i] = elevate_sextic_point_to_nine(grid[i]);
    }
    let mut out: [[Point3; 10]; 10] = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for j in 0..10 {
        let col: [Point3; 7] = core::array::from_fn(|i| intermediate[i][j]);
        let col9 = elevate_sextic_point_to_nine(col);
        for i in 0..10 {
            out[i][j] = col9[i];
        }
    }
    out
}

#[cfg(test)]
#[path = "kp_degree_elevate_tests.rs"]
mod tests;
