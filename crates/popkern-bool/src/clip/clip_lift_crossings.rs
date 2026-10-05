//! Lifted-loop crossings for MULTI-TURN periodic charts (PK-673).
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
pub(super) fn boundary_is_multi_turn(
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    period: f64,
) -> bool {
    let (lo, hi) = boundary_u_range(boundary);
    hi - lo > period * (1.0 + 1e-9)
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
/// pass's `crossings` when `curve_uv` is a once-winding loop on a
/// multi-turn chart; a no-op otherwise, so single-period charts keep
/// their tuned crossing sets. A lifted crossing at the boundary point
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
    if !(is_once_winding_loop(curve_uv, period) && boundary_is_multi_turn(boundary, period)) {
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
    use super::{boundary_is_multi_turn, is_once_winding_loop};
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
    }
}
