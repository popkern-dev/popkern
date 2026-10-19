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
    // A periodic direction has no edge to clamp at: evaluation wraps,
    // and the continuation consumer keeps the hint's branch. Clamped,
    // a closed rail whose hint sat 0.005 short of a periodic strip's
    // u = 0 stalled on that chart edge 0.21 off the rail (PK-891).
    let (u_periodic, v_periodic) = (nurbs.u_period().is_some(), nurbs.v_period().is_some());
    let clamp_u = |u: f64| if u_periodic { u } else { u.clamp(u_lo, u_hi) };
    let clamp_v = |v: f64| if v_periodic { v } else { v.clamp(v_lo, v_hi) };
    let mut u = clamp_u(u0);
    let mut v = clamp_v(v0);

    for _ in 0..MAX_ITER {
        // Fused analytic point+partials (PK-281): one basis
        // contraction per iteration.
        let (s, su, sv) = nurbs.point_and_first_partials(u, v);
        let diff = *point - s;
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
        u = clamp_u(u + du);
        v = clamp_v(v + dv);
        if du * du + dv * dv < STEP_TOL {
            break;
        }
    }

    let dist = point.distance_to(&surface.point_at(u, v));
    if dist <= MAX_DIST { Some((u, v)) } else { None }
}

/// Multi-start Gauss-Newton projection of a 3D point onto a NURBS surface.
///
/// Tries five starting points (parameter domain centre plus four corners) and
/// returns the `(u, v)` with the smallest final Euclidean distance, or `None`
/// if no starting point converges to within 1.0 units. Iterates CLAMP to the
/// chart domain — this projector's consumers (trim/clip/classification)
/// require in-chart parameters. Few-span charts only: multi-sheet charts
/// take the knot-seeded path above (PK-446/PK-449).
pub(super) fn nurbs_gauss_newton_project(
    surface: &SurfaceGeometry,
    nurbs: &NurbsSurface,
    point: &Point3,
) -> Option<(f64, f64)> {
    const MAX_ITER: usize = 20;
    const STEP_TOL: f64 = 1e-12;
    const MAX_DIST: f64 = 1.0;

    let ((u_lo, u_hi), (v_lo, v_hi)) = nurbs.parameter_range();
    let u_mid = 0.5 * (u_lo + u_hi);
    let v_mid = 0.5 * (v_lo + v_hi);

    let starts = [
        (u_mid, v_mid),
        (u_lo, v_lo),
        (u_lo, v_hi),
        (u_hi, v_lo),
        (u_hi, v_hi),
    ];

    let mut best_uv: Option<(f64, f64)> = None;
    let mut best_dist = f64::INFINITY;

    for (u0, v0) in starts {
        let mut u = u0;
        let mut v = v0;

        for _ in 0..MAX_ITER {
            // Fused analytic point+partials: one basis contraction per
            // iteration instead of the 5-evaluation point + FD pattern
            // (PK-281).
            let (s, su, sv) = nurbs.point_and_first_partials(u, v);
            let diff = *point - s;

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
        if dist < best_dist {
            best_dist = dist;
            best_uv = Some((u, v));
        }
    }

    if best_dist <= MAX_DIST { best_uv } else { None }
}

/// A Gregory patch's preimage within rounding (1e-12) of its unit
/// chart's edge, moved onto that edge. A dual-jet patch's boundary 2-jet
/// is its polynomial surrogate's only on the boundary itself, and the
/// certifier takes the surrogate only for pcurves lying on it exactly
/// (PK-541): a corner vertex projected to u = 1 − 1e-16, the pcurve
/// recompute fitted its line through that, and the rim seam then fell to
/// the rational form's 0/0 corner and stayed uncertified (PK-595).
pub(super) fn onto_unit_chart_edge((u, v): (f64, f64)) -> (f64, f64) {
    let snap = |x: f64| {
        if x.abs() <= 1e-12 {
            0.0
        } else if (x - 1.0).abs() <= 1e-12 {
            1.0
        } else {
            x
        }
    };
    (snap(u), snap(v))
}
