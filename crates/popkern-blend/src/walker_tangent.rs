//! Tangent stations (PK-764, `docs/architecture/fillet-vanishing-ends.md`):
//! where the spine point lies on both supports with their normals
//! parallel, the crease vanishes and the rolling ball touches both
//! supports at the spine point itself. The contact equations are
//! singular there (the Jacobian loses rank: a walk that warm-started
//! onto the point declared `NormalsParallel` at condition ~1e10, and a
//! finer station walk truncated at it), but the solution is known in
//! closed form: both contacts are the spine point, the centre lies `r`
//! along the common normal on the ball's side.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

use crate::contact::{BallSide, ContactSolution};
use crate::walker::WalkRequest;

/// `u` moved by whole periods to the branch nearest `near`.
fn near_branch(u: f64, near: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => u + p * ((near - u) / p).round(),
        _ => u,
    }
}

/// The exact ball at a tangent station, or `None` when the station is
/// regular (the supports cross there, or the spine point is off one of
/// them, or the two ball sides disagree).
pub(crate) fn tangent_station(
    request: &WalkRequest<'_>,
    spine_point: Point3,
    radius: f64,
    prev_uv0: (f64, f64),
    prev_uv1: (f64, f64),
    tolerance: &Tolerance,
) -> Option<ContactSolution> {
    let on = tolerance.linear().max(1e-9) * 10.0;
    let foot = |surface: &dyn popkern_core::geometry::surface::SurfaceEval,
                prev: (f64, f64),
                periodicity: (Option<f64>, Option<f64>)| {
        let (u, v) = surface.closest_point(&spine_point)?;
        let uv = (
            near_branch(u, prev.0, periodicity.0),
            near_branch(v, prev.1, periodicity.1),
        );
        (surface.point_at(uv.0, uv.1).distance_to(&spine_point) <= on).then_some(uv)
    };
    let uv0 = foot(request.surface_0, prev_uv0, request.surface_0_periodicity)?;
    let uv1 = foot(request.surface_1, prev_uv1, request.surface_1_periodicity)?;
    let side = |s: BallSide| if s == BallSide::Outside { 1.0 } else { -1.0 };
    let n0 = request.surface_0.normal_at(uv0.0, uv0.1).normalized()? * side(request.side_0);
    let n1 = request.surface_1.normal_at(uv1.0, uv1.1).normalized()? * side(request.side_1);
    // Parallel and on the same side: one ball touches both here.
    if n0.cross(&n1).length() > 1e-6 || n0.dot(&n1) <= 0.0 {
        return None;
    }
    Some(ContactSolution {
        uv0,
        uv1,
        center: spine_point + n0 * radius,
        contact0: spine_point,
        contact1: spine_point,
        iterations: 0,
        step_rejections: 0,
        residual: 0.0,
        condition: f64::INFINITY,
        jacobian: [[0.0; 4]; 4],
    })
}
