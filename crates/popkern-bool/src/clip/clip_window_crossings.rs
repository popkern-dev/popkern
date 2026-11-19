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
use popkern_math::PortableMath as _;

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
    // A closed record's closure is one place: a crossing found at (or a
    // rounding past) its end parameter is the one at its start; and a
    // crossing found on two lifts (at a lift's window edge) is one.
    let (t0, t1) = query.parameter_range();
    let (a, b) = (query.point_at(t0), query.point_at(t1));
    let du = (b.x - a.x) - period * ((b.x - a.x) / period).round();
    let eps = 1e-6 * (t1 - t0).abs();
    if du.phypot(b.y - a.y) < 1e-6 * period {
        let at_start = crossings.iter().any(|c| (c.t_curve - t0).abs() <= eps);
        if at_start {
            crossings.retain(|c| c.t_curve < t1 - eps);
        } else {
            for c in &mut crossings {
                if c.t_curve >= t1 - eps {
                    c.t_curve = t0;
                }
            }
        }
    }
    crossings.sort_by(|a, b| a.t_curve.total_cmp(&b.t_curve));
    crossings.dedup_by(|x, y| (x.t_curve - y.t_curve).abs() <= eps);
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
