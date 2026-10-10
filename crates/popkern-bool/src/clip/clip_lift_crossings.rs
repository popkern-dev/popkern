//! Lifted-loop crossings for periodic charts whose trim leaves the
//! loop's u-window: multi-turn charts (PK-673) and single-turn trims
//! straddling the window end (PK-704).
//!
//! A closed query that winds exactly once around a u-periodic
//! surface (a coaxial circle on a cylinder or cone) is a LOOP on the
//! surface: its preimage in the unrolled chart is the whole family of
//! period translates of one lift. A face whose trim unrolls several
//! periods (a helical thread band) lives in that unrolled cover, and
//! the stretch of the band the loop crosses can straddle the branch
//! window the query was normalized to — the band's stripe at the
//! loop's height runs across the branch boundary, so the plain pass
//! sees ONE of its two boundary crossings. Measured on the default
//! M3 nut's top mouth chamfer: the crest band's stripe at the cone
//! circle (v = 2.3069) spans u ∈ [−1.502, +0.0495] while the stripe
//! pick put the query on [−2π, 0]; the +0.0495 crossing of the lower
//! crest helix was off-window, one crossing read as "no hull ->
//! Outside", the crest band was never split by the chamfer and
//! classified whole as inside the tool, and the cellular body went
//! non-manifold along every crest edge.
//!
//! The honest crossing set of a loop is the union over ALL lifts:
//! translating the query by a whole period leaves its parameter
//! unchanged, so a crossing found on any lift is a crossing at the
//! same `t` of the loop. Containment downstream already folds lifts
//! (`point_inside_face_uv` probes u and u ± period), so the windows
//! these crossings bound classify consistently.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::intersection::{CurveIntersection, intersect_curves};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::clip_uv::UvCrossing;

/// Relative closure gate: the loop's end must sit one period from
/// its start in u and at its start's v within this fraction of the
/// period (chart units — the record's own uv fit slop is orders
/// below it, and an open record a full period long that does not
/// close is not a loop).
const LOOP_CLOSURE_FRACTION: f64 = 1e-6;

/// Whether `curve_uv` is a once-winding loop on a u-periodic chart of
/// period `period`: its ends one period apart in u, level in v.
pub(super) fn is_once_winding_loop(curve_uv: &NurbsCurve, period: f64) -> bool {
    let (t0, t1) = curve_uv.parameter_range();
    let a = curve_uv.point_at(t0);
    let b = curve_uv.point_at(t1);
    let gate = LOOP_CLOSURE_FRACTION * period;
    ((b.x - a.x).abs() - period).abs() < gate && (b.y - a.y).abs() < gate
}

/// Whether the boundary pcurves unroll more than one period in u
/// (a multi-turn chart).
/// The crossing a CLOSED query makes at its own closure point, when
/// that point lies on the face boundary and the crossing count is odd
/// (PK-775). The intersector misses a hit exactly at a curve's
/// parameter end. A torus band's cube-wall record starts and ends on
/// the band's v = 0 rim, so it saw its one interior crossing only; an
/// odd count gave no hull and the record dropped as Outside. A closed
/// curve crosses an even number of times, so the missing hit is the
/// closure. The query is closed, and its closure point on the
/// boundary, within 1e-5 of the boundary's uv diagonal: the record's uv
/// fit leaves its ends ~8e-7 apart on a 6.3-wide chart.
pub(super) fn closure_crossing(
    brep: &BRep,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    crossings: &[UvCrossing],
    period: Option<f64>,
) -> Option<UvCrossing> {
    if crossings.len().is_multiple_of(2) {
        return None;
    }
    let (t0, t1) = curve_uv.parameter_range();
    let (a, b) = (curve_uv.point_at(t0), curve_uv.point_at(t1));
    let du = b.x - a.x;
    let du = match period {
        Some(p) if p > 0.0 => du - p * (du / p).round(),
        _ => du,
    };
    let (u_lo, u_hi, v_lo, v_hi) = boundary
        .iter()
        .flat_map(|(_, pc)| pc.control_points.iter())
        .fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(a, b, c, d), p| (a.min(p.x), b.max(p.x), c.min(p.y), d.max(p.y)),
        );
    let diag = popkern_math::PortableMath::phypot(u_hi - u_lo, v_hi - v_lo);
    let band = (diag * 1e-5).max(1e-12);
    if du.abs() > band || (b.y - a.y).abs() > band {
        return None; // not closed in the chart
    }
    let mut best: Option<(f64, UvCrossing)> = None;
    for (he, pc) in boundary {
        let twin = brep.half_edges[*he].twin;
        if brep.half_edges.contains_key(twin)
            && brep.half_edges[twin].loop_ == brep.half_edges[*he].loop_
        {
            continue; // a seam use
        }
        for lift in [0.0, -1.0, 1.0] {
            let shift = period.map_or(0.0, |p| lift * p);
            let q = popkern_core::geometry::Point3::new(a.x + shift, a.y, 0.0);
            let Some(s) = pc.closest_point(&q) else {
                continue;
            };
            let d = pc.point_at(s).distance_to(&q);
            if d <= band && best.as_ref().is_none_or(|(bd, _)| d < *bd) {
                best = Some((
                    d,
                    UvCrossing {
                        t_curve: t0,
                        boundary_he: *he,
                        t_boundary: s,
                    },
                ));
            }
        }
    }
    best.map(|(_, c)| c)
}

#[cfg(test)]
fn boundary_is_multi_turn(boundary: &[(HalfEdgeHandle, NurbsCurve)], period: f64) -> bool {
    let (lo, hi) = boundary_u_range(boundary);
    hi - lo > period * (1.0 + 1e-9)
}

/// Whether the boundary pcurves reach outside the u-window the loop
/// `curve_uv` occupies — some of the trim lives only on another lift.
/// A multi-turn chart always does (its trim is wider than the
/// window); so does a single-turn face whose trim STRADDLES the
/// window's end (PK-704: a cylinder piece u ∈ [−0.262, 1.047] against
/// a plane circle on [0, 2π] — its u = −0.262 side met the query only
/// on the k = −1 lift, the plain pass saw one crossing and the record
/// dropped as Outside). The window gate is relative to the period.
pub(super) fn boundary_leaves_loop_window(
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    period: f64,
) -> bool {
    let (b_lo, b_hi) = boundary_u_range(boundary);
    let (q_lo, q_hi) = curve_uv
        .control_points
        .iter()
        .map(|p| p.x)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        });
    let gate = LOOP_CLOSURE_FRACTION * period;
    b_lo < q_lo - gate || b_hi > q_hi + gate
}

/// Whether a sub-period query leaves the window of a face that WRAPS a
/// whole period, so that some of its crossings exist only on another
/// lift (PK-775). The face's seam need not sit at u = 0: a rigid
/// rotation re-anchors a rotational chart and leaves the seam at the
/// shift angle. Two measured shapes are covered:
/// - a CLOSED record straddling the moved seam (a torus band's trim on
///   [−0.002, 6.281] against a cube-wall record on [5.83, 6.73], which
///   met the band's rim only on its k = −1 lift);
/// - a record lying ENTIRELY outside the window (a cylinder wall on
///   [−0.89, 5.39] against a cube-wall line at u = 5.5).
///
/// An OPEN record straddling a seam is the ±π dance's tuned class:
/// lifting it added phantom crossings at an equal-radius tee's
/// tangential singular points, and the union stranded record sides.
pub(super) fn query_leaves_wrapping_window(
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    period: f64,
) -> bool {
    let (b_lo, b_hi) = boundary_u_range(boundary);
    let (q_lo, q_hi) = curve_uv
        .control_points
        .iter()
        .map(|p| p.x)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        });
    let gate = LOOP_CLOSURE_FRACTION * period;
    let wraps = ((b_hi - b_lo) - period).abs() < gate;
    if !wraps || q_hi - q_lo >= period - gate {
        return false;
    }
    let (t0, t1) = curve_uv.parameter_range();
    let closed = curve_uv.point_at(t0).distance_to(&curve_uv.point_at(t1)) < 1e-5 * period;
    let straddles = q_lo < b_lo - gate || q_hi > b_hi + gate;
    let outside = q_lo > b_hi + gate || q_hi < b_lo - gate;
    (closed && straddles) || outside
}

fn boundary_u_range(boundary: &[(HalfEdgeHandle, NurbsCurve)]) -> (f64, f64) {
    boundary
        .iter()
        .flat_map(|(_, pc)| pc.control_points.iter().map(|p| p.x))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        })
}

/// Crossings of every whole-period lift of the loop `curve_uv`
/// (other than the lift itself) with the boundary pcurves, reported
/// at the loop's own parameter. Seam uses (an edge whose twin runs in
/// the same loop) are skipped as in the plain pass. Lifts are taken
/// over the boundary's whole unrolled u-range, so no turn of the
/// chart is missed; a lift that overlaps no boundary costs one
/// curve×curve call per pcurve and finds nothing.
pub(super) fn lifted_loop_crossings(
    brep: &BRep,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    period: f64,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    let (b_lo, b_hi) = boundary_u_range(boundary);
    let (q_lo, q_hi) = curve_uv
        .control_points
        .iter()
        .map(|p| p.x)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        });
    if !(b_lo.is_finite() && b_hi.is_finite() && q_lo.is_finite() && q_hi.is_finite()) {
        return Vec::new();
    }
    let k_min = ((b_lo - q_hi) / period).floor() as i64;
    let k_max = ((b_hi - q_lo) / period).ceil() as i64;
    let curve_range = curve_uv.parameter_range();
    let mut found = Vec::new();
    for k in k_min..=k_max {
        if k == 0 {
            continue;
        }
        let mut lift = curve_uv.clone();
        for p in &mut lift.control_points {
            p.x += k as f64 * period;
        }
        for (he, pc) in boundary {
            let range = pc.parameter_range();
            if (range.1 - range.0).abs() < tolerance.numerical_floor() {
                continue;
            }
            let twin = brep.half_edges[*he].twin;
            if brep.half_edges.contains_key(twin)
                && brep.half_edges[twin].loop_ == brep.half_edges[*he].loop_
            {
                continue;
            }
            let Ok(result) = intersect_curves(
                &CurveGeometry::Nurbs(pc.clone()),
                range,
                &CurveGeometry::Nurbs(lift.clone()),
                curve_range,
                tolerance,
            ) else {
                continue;
            };
            for hit in &result.intersections {
                match hit {
                    CurveIntersection::Point { t1, t2, .. } => found.push(UvCrossing {
                        t_curve: *t2,
                        boundary_he: *he,
                        t_boundary: *t1,
                    }),
                    CurveIntersection::Overlap { t1_range, t2_range } => {
                        found.push(UvCrossing {
                            t_curve: t2_range.0,
                            boundary_he: *he,
                            t_boundary: t1_range.0,
                        });
                        found.push(UvCrossing {
                            t_curve: t2_range.1,
                            boundary_he: *he,
                            t_boundary: t1_range.1,
                        });
                    }
                }
            }
        }
    }
    if popkern_core::trace::trace_enabled("clip") && !found.is_empty() {
        eprintln!(
            "    [uvclip-lift] {} crossing(s) on lifts k in [{k_min},{k_max}]",
            found.len()
        );
    }
    found
}

/// Add the loop's lifted crossings (see the module docs) to the plain
/// pass's `crossings` when `curve_uv` is a once-winding loop whose
/// face trim reaches outside the loop's u-window (a multi-turn chart,
/// or a single-turn trim straddling the window end); a no-op when the
/// trim lies within the window, where the plain pass already sees
/// every crossing and the tuned crossing sets stand. A lifted crossing at the boundary point
/// an existing crossing already reports (same half-edge, same
/// boundary parameter within the intersector's agreement) is the same
/// crossing seen from the next branch and is not added twice.
pub(super) fn extend_with_lifted_loop_crossings(
    crossings: &mut Vec<UvCrossing>,
    brep: &BRep,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    period: f64,
    tolerance: &Tolerance,
) {
    if !(is_once_winding_loop(curve_uv, period)
        && boundary_leaves_loop_window(boundary, curve_uv, period)
        || query_leaves_wrapping_window(boundary, curve_uv, period))
    {
        return;
    }
    for c in lifted_loop_crossings(brep, boundary, curve_uv, period, tolerance) {
        let duplicate = crossings.iter().any(|p| {
            p.boundary_he == c.boundary_he && {
                let span = boundary
                    .iter()
                    .find(|(h, _)| *h == c.boundary_he)
                    .map(|(_, pc)| {
                        let (a, b) = pc.parameter_range();
                        (b - a).abs()
                    })
                    .unwrap_or(1.0);
                (p.t_boundary - c.t_boundary).abs() < span * 1e-6
            }
        });
        if !duplicate {
            crossings.push(c);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{boundary_is_multi_turn, boundary_leaves_loop_window, is_once_winding_loop};
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    use popkern_core::topology::HalfEdgeHandle;
    use std::f64::consts::TAU;

    fn line(a: (f64, f64), b: (f64, f64)) -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(a.0, a.1, 0.0), Point3::new(b.0, b.1, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn loop_and_multi_turn_classification() {
        // A coaxial circle's chart image: one period in u, level in v.
        assert!(is_once_winding_loop(&line((-TAU, 2.3), (0.0, 2.3)), TAU));
        // Not a loop: short of a period, or climbing in v (a helix turn).
        assert!(!is_once_winding_loop(&line((0.0, 2.3), (6.0, 2.3)), TAU));
        assert!(!is_once_winding_loop(&line((0.0, 2.3), (TAU, 2.8)), TAU));
        // The crest band's unrolled trim spans ~4.4 periods.
        let band = vec![
            (HalfEdgeHandle::default(), line((-23.27, 0.57), (0.0, 2.43))),
            (HalfEdgeHandle::default(), line((0.0, 2.43), (4.21, 2.64))),
        ];
        assert!(boundary_is_multi_turn(&band, TAU));
        let single = vec![(HalfEdgeHandle::default(), line((0.0, 1.0), (6.0, 1.0)))];
        assert!(!boundary_is_multi_turn(&single, TAU));
        // A single-turn trim inside the loop's window: nothing on other
        // lifts. One straddling the window end (PK-704): its far side
        // meets the loop only on the k = -1 lift.
        let circle = line((0.0, 10.0), (TAU, 10.0));
        assert!(!boundary_leaves_loop_window(&single, &circle, TAU));
        let straddle = vec![(HalfEdgeHandle::default(), line((-0.262, 0.0), (1.047, 0.0)))];
        assert!(boundary_leaves_loop_window(&straddle, &circle, TAU));
        assert!(boundary_leaves_loop_window(&band, &circle, TAU));
    }
}
