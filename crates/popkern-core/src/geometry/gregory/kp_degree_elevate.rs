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
/// formulas, which consume the guide's cross-derivative /
/// cross-second-derivative CPs at degree 9.
pub fn elevate_cubic_vec_to_nine(cubic: [Vec3; 4]) -> [Vec3; 10] {
    let poly = BernsteinPoly::new(cubic.to_vec()).degree_elevate(6);
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

#[cfg(test)]
#[path = "kp_degree_elevate_tests.rs"]
mod tests;
