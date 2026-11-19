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
    let mut crossings = super::clip_lift_crossings::lifted_crossings(
        brep, boundary, query, period, tolerance, true,
    );
    // A closed record's closure is one place: a crossing found at (or a
    // rounding past) its end parameter is the one at its start.
    let (t0, t1) = query.parameter_range();
    // Closed modulo the period: a loop winding the chart once ends a
    // period from its start.
    let (a, b) = (query.point_at(t0), query.point_at(t1));
    let du = (b.x - a.x) - period * ((b.x - a.x) / period).round();
    if du.phypot(b.y - a.y) < 1e-6 * period {
        let eps = 1e-6 * (t1 - t0).abs();
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
    super::clip_crossing_gap::keep_one_place(
        surface,
        boundary,
        query,
        &mut crossings,
        0,
        tolerance,
    );
    crossings.sort_by(|a, b| a.t_curve.total_cmp(&b.t_curve));
    let crossings = super::clip_crossing_reports::collapse_crossing_reports(
        brep, face, query, crossings, tolerance,
    );
    let crossings =
        super::clip_crossing_reports::resolve_contact_runs(brep, face, query, crossings, tolerance);
    let mut crossings =
        super::clip_crossing_graze::drop_grazes((surface, boundary), query, crossings, tolerance);
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
    production: &[UvCrossing],
    window: &[UvCrossing],
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
        let b = query.point_at(t1);
        crate::complex_shadow::append_log_line(&format!(
            "ClipWindowShadow face {face:?} query_start ({:.4},{:.4}) end ({:.6},{:.6}) range ({t0},{t1}) production {p:?} window {w:?}",
            a.x, a.y, b.x, b.y
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
    shadow_compare(face, query, production, &window);
}
