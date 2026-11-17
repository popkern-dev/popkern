//! Whether an overlap of a query with a face boundary edge ends where
//! it reports (PK-1190).
//!
//! A query running along a boundary edge for a stretch crosses the
//! face's trim only where it leaves that edge. The curve intersector
//! also reports an overlap "end" where the query's own parameter runs
//! out. A closed loop's closure point lying along the edge is one, and
//! so is a lift's end, or an end estimated a little short of either.
//! Taken as crossings, these split a window that runs across the
//! closure. A funnel collar's round end sat on the ring's top, its
//! bottom arc crossing both its cylinder's seam and the plane chart's
//! circle start. Its record kept 329 degrees of the circle instead of
//! the 180 between its two true corners, and the union declined
//! non-manifold. An end is an exit only when a step past it leaves the
//! edge: a closed query wraps round its closure, and a periodic chart
//! compares whole-period shifts.

use popkern_core::geometry::Vec3;
use popkern_core::topology::HalfEdgeHandle;

use super::clip_uv::UvCrossing;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;

/// The step past an overlap end, as a fraction of the query's range.
const STEP: f64 = 1e-3;
/// How far off the edge, against the step's length, a stepped point may
/// lie and still run along it.
const COINCIDENT: f64 = 1e-3;
/// How near a closed query's parameter end a point hit sits to be its
/// closure point, as a fraction of the query's range.
const AT_CLOSURE: f64 = 1e-6;

/// Whether the query `curve` leaves the boundary pcurve `boundary`
/// just past its parameter `t` (forward past the overlap's far end,
/// backward past its near end). An open query ending there leaves.
pub(super) fn leaves_boundary(
    (period_u, period_v): (Option<f64>, Option<f64>),
    boundary: &NurbsCurve,
    curve: &NurbsCurve,
    t: f64,
    forward: bool,
) -> bool {
    let (q0, q1) = curve.parameter_range();
    let span = q1 - q0;
    let fold = |d: Vec3| -> f64 {
        let reduce = |x: f64, p: Option<f64>| match p {
            Some(p) if p > 0.0 => x - p * (x / p).round(),
            _ => x,
        };
        Vec3::new(reduce(d.x, period_u), reduce(d.y, period_v), 0.0).length()
    };
    let closed = fold(curve.point_at(q1) - curve.point_at(q0)) <= span.abs() * 1e-9;
    let mut past = t + if forward { STEP * span } else { -STEP * span };
    if past > q1 || past < q0 {
        if !closed {
            return true;
        }
        past = if past > q1 { past - span } else { past + span };
    }
    let p = curve.point_at(past);
    let moved = fold(p - curve.point_at(t));
    // The edge's nearest point to each whole-period translate of the
    // step (a wrapped step lands a period from the edge's own chart).
    let (b0, b1) = boundary.parameter_range();
    let shifts = |p: Option<f64>| match p {
        Some(p) if p > 0.0 => vec![-p, 0.0, p],
        _ => vec![0.0],
    };
    let off = shifts(period_u)
        .into_iter()
        .flat_map(|du| shifts(period_v).into_iter().map(move |dv| (du, dv)))
        .map(|(du, dv)| {
            let q = popkern_core::geometry::Point3::new(p.x + du, p.y + dv, p.z);
            let near = popkern_core::geometry::curve::closest_point_in_window(boundary, &q, b0, b1);
            (boundary.point_at(near) - q).length()
        })
        .fold(f64::INFINITY, f64::min);
    // Still on the edge: coincident with it, off it by a thousandth of
    // the step (a tenth read a thread flank's shallow crossings of its
    // boundary helices as running along them).
    off > COINCIDENT * moved
}

/// A face surface's chart periods.
pub(super) fn chart_periods(surface: Option<&SurfaceGeometry>) -> (Option<f64>, Option<f64>) {
    surface.map_or((None, None), |s| {
        popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(s)
    })
}

/// The crossings an overlap of `curve` with the boundary pcurve of `he`
/// makes: each overlap end where the query leaves the edge.
pub(super) fn overlap_exits(
    periods: (Option<f64>, Option<f64>),
    (boundary, curve): (&NurbsCurve, &NurbsCurve),
    he: HalfEdgeHandle,
    t1_range: (f64, f64),
    t2_range: (f64, f64),
) -> Vec<UvCrossing> {
    // The overlap's near end (the lower query parameter) is left
    // stepping back, its far end stepping on: the intersector reports a
    // range either way round.
    let ascending = t2_range.0 <= t2_range.1;
    [
        (t1_range.0, t2_range.0, !ascending),
        (t1_range.1, t2_range.1, ascending),
    ]
    .into_iter()
    .filter(|&(tb, tq, forward)| {
        at_edge_end(boundary, tb) || leaves_boundary(periods, boundary, curve, tq, forward)
    })
    .map(|(tb, tq, _)| UvCrossing {
        t_curve: tq,
        boundary_he: he,
        t_boundary: tb,
    })
    .collect()
}

/// Whether a point hit at the query's parameter `t` crosses the
/// boundary pcurve: a closed loop's closure point reported where the
/// loop lies along the edge on both sides of it does not. Hits away
/// from the query's parameter ends are crossings as the intersector
/// found them.
pub(super) fn point_crosses(
    periods: (Option<f64>, Option<f64>),
    (boundary, curve): (&NurbsCurve, &NurbsCurve),
    t: f64,
) -> bool {
    let (q0, q1) = curve.parameter_range();
    let reach = (q1 - q0).abs() * AT_CLOSURE;
    if (t - q0).abs() > reach && (t - q1).abs() > reach {
        return true;
    }
    leaves_boundary(periods, boundary, curve, t, true)
        || leaves_boundary(periods, boundary, curve, t, false)
}

/// [`point_crosses`] on a face surface's chart.
pub(super) fn crosses_on(
    surface: Option<&SurfaceGeometry>,
    curves: (&NurbsCurve, &NurbsCurve),
    t: f64,
) -> bool {
    point_crosses(chart_periods(surface), curves, t)
}

/// Whether the boundary parameter `tb` is at either end of the boundary
/// pcurve: an overlap ending where its edge ends exits there in the
/// chart (a seam-split copy's rim ends at its seam, and a step past it
/// folded a period lands on the same circle's other end), so only an
/// end inside the edge is tested.
fn at_edge_end(boundary: &NurbsCurve, tb: f64) -> bool {
    let (b0, b1) = boundary.parameter_range();
    let reach = (b1 - b0).abs() * AT_CLOSURE;
    (tb - b0).abs() <= reach || (tb - b1).abs() <= reach
}
