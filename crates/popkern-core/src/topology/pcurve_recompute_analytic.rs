//! Exact analytic pcurve candidates for `recompute_face_pcurves`
//! (PK-210). Split out for the parent's frozen file-size ratchet.

use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;

/// Exact pcurve candidate for a curve lying in a Plane face —
/// delegates to [`crate::geometry::chart_pcurve::exact_plane_pcurve`]
/// (moved there so pcurve AUTHORS outside the recompute path can use
/// the same exact constructions; dependency direction: topology may
/// call geometry, never the reverse).
pub(super) fn exact_plane_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    crate::geometry::chart_pcurve::exact_plane_pcurve(surface, curve, t_lo, t_hi)
}

/// A placeholder candidate for an edge whose exact plane pcurve exists:
/// `install_pcurve_with_survival` replaces any candidate with that exact
/// pcurve, so the dense projection fit the recompute would otherwise
/// build is thrown away (PK-804: every planar rim circle was refit from
/// 64 samples on every boolean). The placeholder is the UV chord between
/// `ends`; it never survives the install.
pub(super) fn exact_plane_placeholder(
    surface: &SurfaceGeometry,
    curve: &crate::geometry::curve::CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
    ends: [(f64, f64); 2],
) -> Option<crate::geometry::NurbsCurve> {
    exact_plane_pcurve(surface, &curve.geometry, t_lo, t_hi)?;
    Some(crate::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: ends
            .iter()
            .map(|&(u, v)| crate::geometry::Point3::new(u, v, 0.0))
            .collect(),
        weights: vec![1.0, 1.0],
    })
}

/// Whether the edge's 3-D curve is straight (chord test at 8 interior
/// samples). The plane fast path in `recompute_face_pcurves` is only
/// valid for straight edges: a curved edge lying IN a plane (fillet
/// end-runoff arc) has a curved UV image, and a straight pcurve
/// misstates it by the full sagitta (PK-096: 2.4e-2 on runoff arcs —
/// STEP consumers flag the edge and refuse to stitch).
pub(super) fn edge_is_straight(edge_curve: &Option<crate::geometry::curve::CurveOnEdge>) -> bool {
    match edge_curve {
        None => true,
        Some(c) => {
            let eval: &dyn crate::geometry::curve::CurveEval = &c.geometry;
            let a = eval.point_at(c.t_start);
            let b = eval.point_at(c.t_end);
            let chord = b - a;
            let len2 = chord.dot(&chord);
            for i in 1..8 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                let p = eval.point_at(t);
                let ap = p - a;
                let f = if len2 < 1e-300 {
                    0.0
                } else {
                    ap.dot(&chord) / len2
                };
                let q = Point3::new(a.x + chord.x * f, a.y + chord.y * f, a.z + chord.z * f);
                if p.distance_to(&q) > 1.0e-9 {
                    return false;
                }
            }
            true
        }
    }
}

#[cfg(test)]
#[path = "pcurve_recompute_analytic_tests.rs"]
mod pcurve_recompute_analytic_tests;
