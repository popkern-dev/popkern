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
    crossings.dedup_by(|x, y| {
        x.boundary_he == y.boundary_he
            && (x.t_curve - y.t_curve).abs() <= eps
            && (x.t_boundary - y.t_boundary).abs() <= eps.max(1e-9)
    });
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

/// Log a disagreement between the production crossings and the window's
/// (parameters compared within a millionth of the record's span).
pub(super) fn shadow_compare(
    face: FaceHandle,
    query: &NurbsCurve,
    (production, window): (&[UvCrossing], &[UvCrossing]),
    gap: &dyn Fn(&UvCrossing) -> f64,
) {
    let (t0, t1) = query.parameter_range();
    let eps = 1e-6 * (t1 - t0).abs().max(1e-12);
    // A crossing at the record's own end bounds nothing new (either
    // finder may report or omit it): compared off the ends.
    let ts = |c: &[UvCrossing]| {
        let mut v: Vec<f64> = c
            .iter()
            .map(|x| x.t_curve)
            .filter(|&t| t > t0 + eps && t < t1 - eps)
            .collect();
        v.sort_by(f64::total_cmp);
        v.dedup_by(|a, b| (*a - *b).abs() < eps);
        v
    };
    let (p, w) = (ts(production), ts(window));
    let same = p.len() == w.len() && p.iter().zip(&w).all(|(a, b)| (a - b).abs() < eps);
    if !same {
        let a = query.point_at(t0);
        let gaps = |c: &[UvCrossing]| -> Vec<String> {
            c.iter()
                .map(|x| format!("{:.6}@{:.1e}", x.t_curve, gap(x)))
                .collect()
        };
        crate::complex_shadow::append_log_line(&format!(
            "ClipWindowShadow face {face:?} query_start ({:.4},{:.4}) production {p:?} window {w:?} gaps production {:?} window {:?}",
            a.x,
            a.y,
            gaps(production),
            gaps(window)
        ));
    }
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

/// With the shadow log set, compare the window crossings with the
/// production ones for this record.
pub(super) fn shadow(
    brep: &BRep,
    face: FaceHandle,
    on: (Option<&SurfaceGeometry>, &[(HalfEdgeHandle, NurbsCurve)]),
    (query, period): (&NurbsCurve, f64),
    production: &[UvCrossing],
    tolerance: &Tolerance,
) {
    if std::env::var_os("POPKERN_COMPLEX_SHADOW_LOG").is_none() {
        return;
    }
    let window = window_crossings(brep, face, on, query, period, tolerance);
    // Each crossing's gap in space between the record's point and the
    // boundary's: a crossing far past the reach is a phantom.
    let gap = |c: &UvCrossing| -> f64 {
        use popkern_core::geometry::surface::SurfaceEval as _;
        let (Some(s), Some((_, pc))) = (on.0, on.1.iter().find(|(h, _)| *h == c.boundary_he))
        else {
            return f64::NAN;
        };
        let (qa, qb) = (query.point_at(c.t_curve), pc.point_at(c.t_boundary));
        s.point_at(qa.x, qa.y).distance_to(&s.point_at(qb.x, qb.y))
    };
    shadow_compare(face, query, (production, &window), &gap);
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

/// Remove each pair of consecutive crossings of one boundary edge whose
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
    let Some(surface) = surface else {
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
        let allowance = reach + gap(a.t_curve).max(gap(b.t_curve));
        (1..10).all(|k| gap(a.t_curve + (b.t_curve - a.t_curve) * k as f64 / 10.0) <= allowance)
    };
    let mut i = 0;
    while i + 1 < crossings.len() {
        if within(&crossings[i], &crossings[i + 1]) {
            crossings.drain(i..i + 2);
        } else {
            i += 1;
        }
    }
}
