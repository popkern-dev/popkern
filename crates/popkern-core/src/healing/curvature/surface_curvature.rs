//! Surface curvature computation from the shape operator.
//!
//! All functions accept any surface implementing [`SurfaceEval`] and
//! return `None` at degenerate points (poles, collapsed parameter lines)
//! where the shape operator is undefined.

use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::fundamental_forms::shape_operator;

/// Aggregate surface curvature at a single point.
///
/// Contains principal curvatures, Gaussian and mean curvature,
/// and the minimum radius of curvature.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceCurvature {
    /// First principal curvature (larger in absolute value).
    pub kappa_1: f64,
    /// Second principal curvature.
    pub kappa_2: f64,
    /// Gaussian curvature: K = kappa_1 * kappa_2.
    pub gaussian: f64,
    /// Mean curvature: H = (kappa_1 + kappa_2) / 2.
    pub mean: f64,
    /// Minimum radius of curvature: 1 / max(|kappa_1|, |kappa_2|).
    /// Infinite (f64::INFINITY) when both principal curvatures are zero.
    pub min_radius: f64,
}

/// Compute all surface curvature quantities at a point.
///
/// Returns `None` at degenerate points where the shape operator
/// cannot be computed (singular normal, collapsed parameter lines).
pub fn surface_curvature_at(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<SurfaceCurvature> {
    let so = shape_operator(surface, u, v)?;
    let k1 = so.kappa_1;
    let k2 = so.kappa_2;
    let max_abs = k1.abs().max(k2.abs());
    let min_r = if max_abs > 1e-15 {
        1.0 / max_abs
    } else {
        f64::INFINITY
    };
    Some(SurfaceCurvature {
        kappa_1: k1,
        kappa_2: k2,
        gaussian: k1 * k2,
        mean: (k1 + k2) / 2.0,
        min_radius: min_r,
    })
}

/// Principal curvatures at a surface point.
///
/// Returns `(kappa_1, kappa_2)` where `kappa_1 >= kappa_2` in absolute
/// value. Returns `None` at degenerate points.
pub fn principal_curvatures(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<(f64, f64)> {
    let so = shape_operator(surface, u, v)?;
    Some((so.kappa_1, so.kappa_2))
}

/// Gaussian curvature at a surface point.
///
/// K = kappa_1 * kappa_2 = det(II) / det(I).
///
/// - Positive on elliptic points (sphere-like).
/// - Zero on parabolic points (cylinder-like).
/// - Negative on hyperbolic points (saddle-like).
///
/// Returns `None` at degenerate points.
pub fn gaussian_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    Some(so.kappa_1 * so.kappa_2)
}

/// Mean curvature at a surface point.
///
/// H = (kappa_1 + kappa_2) / 2.
///
/// Zero on minimal surfaces (soap films). For a sphere of radius r,
/// H = 1/r (up to sign convention).
///
/// Returns `None` at degenerate points.
pub fn mean_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    Some((so.kappa_1 + so.kappa_2) / 2.0)
}

/// Minimum radius of curvature at a surface point.
///
/// Defined as `1 / max(|kappa_1|, |kappa_2|)`. Returns
/// `f64::INFINITY` when both principal curvatures are zero (flat surface).
///
/// Returns `None` at degenerate points.
pub fn min_radius_of_curvature(surface: &dyn SurfaceEval, u: f64, v: f64) -> Option<f64> {
    let so = shape_operator(surface, u, v)?;
    let max_abs = so.kappa_1.abs().max(so.kappa_2.abs());
    if max_abs > 1e-15 {
        Some(1.0 / max_abs)
    } else {
        Some(f64::INFINITY)
    }
}

#[cfg(test)]
#[path = "surface_curvature_tests.rs"]
mod tests;
