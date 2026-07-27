//! Exact analytic pcurve candidates for `recompute_face_pcurves`
//! (PK-210). Split out for the parent's frozen file-size ratchet.

use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;

/// Exact angle-parameterized elliptical pcurve for an Arc lying in a
/// Plane face — the boss-ring rail class.
///
/// `Arc::point_at` IS angle-linear, so an angle-parameterized
/// elliptical pcurve is rate-consistent with the curve; the generic
/// 23-sample cubic interpolation leaves ~`r·(Δθ)⁴/384` of sagitta
/// instead (measured 2.5e-4..9.2e-4 on the corpus boss-ring rails —
/// above the 1e-5 commercial trim gate, which read three watertight
/// parts as achieved-G0 with silent requested-G1 downgrades).
/// Projection uses the same 2×2 SVD decomposition as the plane×torus
/// analytic edge triple.
pub(super) fn exact_plane_arc_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    let (
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        },
    ) = (surface, curve)
    else {
        return None;
    };
    let (bu, bv) = crate::geometry::analytic_curves::circle_basis(*normal);
    let ul2 = u_axis.dot(u_axis);
    let vl2 = v_axis.dot(v_axis);
    if ul2 < 1e-30 || vl2 < 1e-30 {
        return None;
    }
    let cen: crate::geometry::Vec3 = *center - Point3::new(origin.x, origin.y, origin.z);
    let (a_uu, a_vu) = (radius * bu.dot(u_axis) / ul2, radius * bu.dot(v_axis) / vl2);
    let (a_uv, a_vv) = (radius * bv.dot(u_axis) / ul2, radius * bv.dot(v_axis) / vl2);
    let (u_radius, v_radius, rotation, phase) =
        crate::geometry::edge_triple::decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    Some(Pcurve::elliptical(
        (cen.dot(u_axis) / ul2, cen.dot(v_axis) / vl2),
        u_radius,
        v_radius,
        rotation,
        t_lo - phase,
        t_hi - phase,
        t_lo,
        t_hi,
    ))
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
