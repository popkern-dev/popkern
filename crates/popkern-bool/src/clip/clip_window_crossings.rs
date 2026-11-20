//! The clip's crossings on the chart service's window (W2c, shadow).
//!
//! A record's crossings with a face's boundary on a u-periodic chart are
//! the crossings of every whole-period lift of the record with the
//! boundary pcurves as stored, each kept once at its place in space: a
//! lift by a period is the same curve. One rule for every record, closed
//! or open, wherever the face's seam sits, in place of the branch
//! normalization, the half-period shift and the lifted-loop conditions.
//! In shadow: `find_uv_crossings` logs where it disagrees with them.

use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::clip_uv::UvCrossing;

/// The window crossings of `query` with `boundary` on `face` (chart
/// period `period` along u), filtered as the plain pass's are.
pub(super) fn window_crossings(
    brep: &BRep,
    face: FaceHandle,
    (surface, boundary): (Option<&SurfaceGeometry>, &[(HalfEdgeHandle, NurbsCurve)]),
    query: &NurbsCurve,
    period: f64,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    // Each lift judged in its own frame (contact runs and grazes read the
    // record's sides against the boundary in uv: a lift's crossings read
    // against the unlifted record took a cone dimple's touches along its
    // edge for crossings), then merged.
    // The boundary on the branch its edges lie on in space (a pcurve stored
    // the long way round re-lifted from its edge).
    let read = surface.map(|s| super::clip_window_boundary::on_its_edges(brep, s, boundary));
    let boundary = read.as_deref().unwrap_or(boundary);
    let Some((k_min, k_max)) = super::clip_lift_crossings::lift_range(boundary, query, period)
    else {
        return Vec::new();
    };
    let mut crossings: Vec<UvCrossing> = Vec::new();
    for k in k_min..=k_max {
        let mut lift = query.clone();
        for p in &mut lift.control_points {
            p.x += k as f64 * period;
        }
        let mut found =
            super::clip_lift_crossings::lift_crossings(brep, boundary, &lift, period, tolerance);
        if found.is_empty() {
            continue;
        }
        super::clip_crossing_gap::keep_one_place(
            surface, boundary, &lift, &mut found, 0, tolerance,
        );
        found.sort_by(|a, b| a.t_curve.total_cmp(&b.t_curve));
        let found = super::clip_crossing_reports::collapse_crossing_reports(
            brep, face, &lift, found, tolerance,
        );
        let found = super::clip_crossing_reports::resolve_contact_runs_in_frame(
            brep,
            face,
            &lift,
            found,
            (tolerance, true),
        );
        crossings.extend(super::clip_crossing_graze::drop_grazes(
            (surface, boundary),
            &lift,
            found,
            tolerance,
        ));
    }
    // A crossing found on two lifts (at a lift's window edge) is one: the
    // same boundary edge at the same parameters. Reports at a closed
    // record's closure (at its start and its end, or on two edges meeting
    // there) are kept as the plain pass keeps them: the clip's arms read
    // them as its pair.
    let (t0, t1) = query.parameter_range();
    let eps = 1e-6 * (t1 - t0).abs();
    crossings.sort_by(|a, b| a.t_curve.total_cmp(&b.t_curve));
    // One place in space within the reach is one crossing (two lifts read
    // a crossing 1.15e-6 apart in the record's parameter, over a
    // parameter test of a millionth).
    let reach = brep.tolerance.linear_reach();
    let place = |c: &UvCrossing| {
        use popkern_core::geometry::surface::SurfaceEval as _;
        let q = query.point_at(c.t_curve);
        surface.map(|s| s.point_at(q.x, q.y))
    };
    crossings.dedup_by(|x, y| {
        x.boundary_he == y.boundary_he
            && match (place(x), place(y)) {
                (Some(a), Some(b)) => a.distance_to(&b) <= reach,
                _ => {
                    (x.t_curve - y.t_curve).abs() <= eps
                        && (x.t_boundary - y.t_boundary).abs() <= eps.max(1e-9)
                }
            }
    });
    merge_one_pass(brep, surface, query, &mut crossings);
    // A record running beside one boundary edge within the reach, its two
    // crossings of that edge either end of the stretch, touches it: two
    // curves within the tolerance meet once (tangency design §4.8; a
    // ball-ball circle grazing a rim crossed it 3.7e-3 apart, 7e-7 deep).
    drop_touch_excursions(brep, surface, query, &mut crossings);
    // A closed record on the boundary at its closure point (PK-775).
    let period_v = surface
        .and_then(|s| popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(s).1);
    if let Some(c) = super::clip_lift_crossings::closure_crossing(
        brep,
        boundary,
        query,
        &crossings,
        (Some(period), period_v),
    ) {
        crossings.insert(0, c);
    } else if let Some(pair) =
        closure_pair(brep, boundary, query, &crossings, (Some(period), period_v))
    {
        let [first, last] = pair;
        crossings.insert(0, first);
        crossings.push(last);
    }
    crossings
}

/// Whether `face` lies on a rotational analytic chart (cylinder, sphere,
/// torus, cone: the ones the half-period dance serves).
pub(super) fn rotational_chart(brep: &BRep, face: FaceHandle) -> bool {
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|s| {
            matches!(
                s.geometry,
                SurfaceGeometry::Cylinder { .. }
                    | SurfaceGeometry::Sphere { .. }
                    | SurfaceGeometry::Torus { .. }
                    | SurfaceGeometry::Cone { .. }
            )
        })
}

/// The closure pair of a closed record whose closure point lies on the
/// boundary and no crossing found anywhere (the window's per-lift filters
/// read a record meeting the boundary at its closure alone as a touch):
/// one crossing at each end of the record's range, on the boundary edge
/// the closure lies on. The clip's arms read this pair as the record kept
/// whole between its ends, as they read the plain pass's two reports at
/// the closure. `None` when any crossing was found, the record is open,
/// or its closure is off the boundary.
fn closure_pair(
    brep: &BRep,
    boundary: &[(HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    crossings: &[UvCrossing],
    periods: (Option<f64>, Option<f64>),
) -> Option<[UvCrossing; 2]> {
    // Only a record meeting the boundary at its closure alone: one that
    // crosses or rides the boundary elsewhere is windowed by those
    // crossings (a slot collar's rim riding its edge across the closure
    // took a spurious pair there).
    if !crossings.is_empty() {
        return None;
    }
    let (_, t1) = curve_uv.parameter_range();
    // The closure's crossing as the odd-count rule finds it (closed in the
    // chart, its point on a non-seam boundary edge within the band).
    let one = super::clip_lift_crossings::closure_crossing(
        brep,
        boundary,
        curve_uv,
        &[UvCrossing {
            t_curve: f64::NAN,
            boundary_he: HalfEdgeHandle::default(),
            t_boundary: f64::NAN,
        }],
        periods,
    )?;
    let end = UvCrossing {
        t_curve: t1,
        ..one.clone()
    };
    Some([one, end])
}

/// Remove the one pair of crossings of a record (its only ones), on one
/// boundary edge, whose
/// stretch of the record between stays within the linear reach of that
/// edge in space, plus the crossings' own fit deviation (nine stations
/// across it): a touch, not two crossings.
fn drop_touch_excursions(
    brep: &BRep,
    surface: Option<&SurfaceGeometry>,
    query: &NurbsCurve,
    crossings: &mut Vec<UvCrossing>,
) {
    use popkern_core::geometry::curve::closest_point_in_window;
    use popkern_core::geometry::surface::SurfaceEval as _;
    // Only a record whose crossings are that one pair: a cluster of
    // near-tangent crossings (a record passing through a rim's extreme)
    // holds real crossings among its pairs.
    let Some(surface) = surface.filter(|_| crossings.len() == 2) else {
        return;
    };
    let reach = brep.tolerance.linear_reach();
    let within = |a: &UvCrossing, b: &UvCrossing| -> bool {
        if a.boundary_he != b.boundary_he || b.t_curve <= a.t_curve {
            return false;
        }
        let Some(curve) = brep
            .half_edges
            .get(a.boundary_he)
            .and_then(|he| brep.edges.get(he.edge))
            .and_then(|e| e.curve.as_ref())
        else {
            return false;
        };
        let gap = |t: f64| {
            let q = query.point_at(t);
            let p = surface.point_at(q.x, q.y);
            let s = closest_point_in_window(&curve.geometry, &p, curve.t_start, curve.t_end);
            curve.geometry.point_at(s).distance_to(&p)
        };
        // The reach plus the crossings' own fit deviation (their distance
        // from the edge in space: a uv crossing a pcurve's fit puts off
        // the edge), as two records are read to meet (PK-1165).
        // Two places, farther apart than the reach: one crossing reported
        // twice at a point (a rim extreme's duplicate) is not an excursion.
        let (qa, qb) = (query.point_at(a.t_curve), query.point_at(b.t_curve));
        if surface
            .point_at(qa.x, qa.y)
            .distance_to(&surface.point_at(qb.x, qb.y))
            <= reach
        {
            return false;
        }
        let allowance = reach + gap(a.t_curve).max(gap(b.t_curve));
        (1..10).all(|k| gap(a.t_curve + (b.t_curve - a.t_curve) * k as f64 / 10.0) <= allowance)
    };
    // A stretch riding an edge from corner to corner (its ends crossings of
    // the neighbouring edges too) is the boundary-riding class, not a
    // touch: only a pair whose crossings no other edge shares is one.
    let (q0, q1) = query.parameter_range();
    let eps = 1e-6 * (q1 - q0).abs();
    let at_corner = |c: &UvCrossing, all: &[UvCrossing]| {
        all.iter()
            .any(|o| o.boundary_he != c.boundary_he && (o.t_curve - c.t_curve).abs() <= eps)
    };
    let mut i = 0;
    while i + 1 < crossings.len() {
        if within(&crossings[i], &crossings[i + 1])
            && !at_corner(&crossings[i], crossings)
            && !at_corner(&crossings[i + 1], crossings)
        {
            crossings.drain(i..i + 2);
        } else {
            i += 1;
        }
    }
}

/// Merge consecutive crossings of one boundary edge where the record stays
/// within the linear reach of that edge between them (three stations), as
/// one pass: kept at the crossing lying on the edge in space (two lifts
/// read one shallow crossing of a rim apart along it, one of them a
/// crossing of the pcurve's fit off the edge: a pass meets once).
fn merge_one_pass(
    brep: &BRep,
    surface: Option<&SurfaceGeometry>,
    query: &NurbsCurve,
    crossings: &mut Vec<UvCrossing>,
) {
    use popkern_core::geometry::curve::closest_point_in_window;
    use popkern_core::geometry::surface::SurfaceEval as _;
    let Some(surface) = surface else {
        return;
    };
    let reach = brep.tolerance.linear_reach();
    let at = |t: f64| {
        let q = query.point_at(t);
        surface.point_at(q.x, q.y)
    };
    let edge_of = |c: &UvCrossing| {
        brep.half_edges
            .get(c.boundary_he)
            .and_then(|he| brep.edges.get(he.edge))
            .and_then(|e| e.curve.as_ref())
    };
    let edge_gap = |c: &UvCrossing, t: f64| -> f64 {
        let Some(curve) = edge_of(c) else {
            return f64::INFINITY;
        };
        let p = at(t);
        let s = closest_point_in_window(&curve.geometry, &p, curve.t_start, curve.t_end);
        curve.geometry.point_at(s).distance_to(&p)
    };
    // The record's side of the edge at `t`: the sign of its offset from
    // the edge across the surface (normal × edge tangent).
    let side = |c: &UvCrossing, t: f64| -> f64 {
        let Some(curve) = edge_of(c) else {
            return 0.0;
        };
        let p = at(t);
        let s = closest_point_in_window(&curve.geometry, &p, curve.t_start, curve.t_end);
        let foot = curve.geometry.point_at(s);
        let Some((u, v)) = surface.closest_point(&foot) else {
            return 0.0;
        };
        let n = surface.normal_at(u, v);
        (p - foot)
            .dot(&n.cross(&curve.geometry.derivative_at(s, 1)))
            .signum()
    };
    let (t0, t1) = query.parameter_range();
    let short = 1e-3 * (t1 - t0).abs();
    let mut i = 0;
    while i + 1 < crossings.len() {
        let (a, b) = (&crossings[i], &crossings[i + 1]);
        // One pass changes sides once: a touch (in and back out within the
        // reach) is the excursion rule's, never merged here.
        let step = (b.t_curve - a.t_curve).abs().max(1e-6 * (t1 - t0).abs());
        // The two are one pass where one of them lies on the edge in space
        // (the other a crossing of the pcurve's fit, off the edge: 2e-5
        // along a shallow crossing, 3e-6 from the edge, 3.2e-7 from the
        // pcurve's image).
        let one = a.boundary_he == b.boundary_he
            && (b.t_curve - a.t_curve).abs() <= short
            && edge_gap(a, a.t_curve).min(edge_gap(b, b.t_curve)) <= reach
            && side(a, a.t_curve - step) * side(a, b.t_curve + step) < 0.0;
        if one {
            let drop = if edge_gap(a, a.t_curve) <= edge_gap(b, b.t_curve) {
                i + 1
            } else {
                i
            };
            crossings.remove(drop);
        } else {
            i += 1;
        }
    }
}
