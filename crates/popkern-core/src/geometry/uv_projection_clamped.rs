//! Domain-clamped single-start Gauss-Newton projection (nested
//! module of `uv_projection` for the file-size budget). Used for
//! Gregory patches, whose global closest-point multistart can escape
//! the unit domain near degenerate pole corners (PK-087).

use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Single-start Gauss-Newton for any evaluable surface, iterates
/// clamped to the given parameter box.
pub(super) fn gauss_newton_single_clamped(
    surface: &SurfaceGeometry,
    point: &Point3,
    u0: f64,
    v0: f64,
    (u_lo, u_hi): (f64, f64),
    (v_lo, v_hi): (f64, f64),
) -> Option<(f64, f64)> {
    const MAX_ITER: usize = 20;
    const STEP_TOL: f64 = 1e-12;
    let mut u = u0.clamp(u_lo, u_hi);
    let mut v = v0.clamp(v_lo, v_hi);
    for _ in 0..MAX_ITER {
        let s = surface.point_at(u, v);
        let diff = *point - s;
        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);
        let a00 = su.dot(&su);
        let a01 = su.dot(&sv);
        let a11 = sv.dot(&sv);
        let det = a00 * a11 - a01 * a01;
        if det.abs() < 1e-30 {
            break;
        }
        let rhs_u = diff.dot(&su);
        let rhs_v = diff.dot(&sv);
        let du = (rhs_u * a11 - rhs_v * a01) / det;
        let dv = (rhs_v * a00 - rhs_u * a01) / det;
        u = (u + du).clamp(u_lo, u_hi);
        v = (v + dv).clamp(v_lo, v_hi);
        if du.abs() < STEP_TOL && dv.abs() < STEP_TOL {
            break;
        }
    }
    Some((u, v))
}

/// Single-start Gauss-Newton from a specific (u, v) hint.
pub(super) fn nurbs_gauss_newton_single(
    surface: &SurfaceGeometry,
    nurbs: &NurbsSurface,
    point: &Point3,
    u0: f64,
    v0: f64,
) -> Option<(f64, f64)> {
    const MAX_ITER: usize = 20;
    const STEP_TOL: f64 = 1e-12;
    const MAX_DIST: f64 = 1.0;

    let ((u_lo, u_hi), (v_lo, v_hi)) = nurbs.parameter_range();
    let mut u = u0.clamp(u_lo, u_hi);
    let mut v = v0.clamp(v_lo, v_hi);

    for _ in 0..MAX_ITER {
        let s = surface.point_at(u, v);
        let diff = *point - s;
        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);
        let rhs_u = diff.dot(&su);
        let rhs_v = diff.dot(&sv);
        let a00 = su.dot(&su);
        let a01 = su.dot(&sv);
        let a11 = sv.dot(&sv);
        let det = a00 * a11 - a01 * a01;
        if det.abs() < f64::EPSILON {
            break;
        }
        let du = (a11 * rhs_u - a01 * rhs_v) / det;
        let dv = (a00 * rhs_v - a01 * rhs_u) / det;
        u = (u + du).clamp(u_lo, u_hi);
        v = (v + dv).clamp(v_lo, v_hi);
        if du * du + dv * dv < STEP_TOL {
            break;
        }
    }

    let dist = point.distance_to(&surface.point_at(u, v));
    if dist <= MAX_DIST { Some((u, v)) } else { None }
}
