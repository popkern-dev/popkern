//! Exact analytic pcurve candidates for `recompute_face_pcurves`
//! (PK-210). Split out for the parent's frozen file-size ratchet.

use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;

/// Exact pcurve candidate for a curve lying in a Plane face. Dispatches
/// per curve form; returns `None` when no exact construction applies
/// (the caller falls back to the sampled fit, and the survival rule
/// measures whichever candidate is produced — exactness is verified,
/// not assumed).
pub(super) fn exact_plane_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    exact_plane_arc_pcurve(surface, curve, t_lo, t_hi)
        .or_else(|| exact_plane_nurbs_pcurve(surface, curve, t_lo, t_hi))
}

/// Exact angle-parameterized elliptical pcurve for an Arc or full
/// Circle lying in a Plane face — the boss-ring rail class.
///
/// `Arc::point_at` IS angle-linear, so an angle-parameterized
/// elliptical pcurve is rate-consistent with the curve; the generic
/// 23-sample cubic interpolation leaves ~`r·(Δθ)⁴/384` of sagitta
/// instead (measured 2.5e-4..9.2e-4 on the corpus boss-ring rails —
/// above the 1e-5 commercial trim gate, which read three watertight
/// parts as achieved-G0 with silent requested-G1 downgrades).
/// Projection uses the same 2×2 SVD decomposition as the plane×torus
/// analytic edge triple.
fn exact_plane_arc_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = surface
    else {
        return None;
    };
    let (center, normal, radius) = match curve {
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        }
        | CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => (center, normal, radius),
        _ => return None,
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

/// Exact NURBS pcurve for a NURBS curve lying in a Plane face — the
/// boolean hole-rim class (PK-287).
///
/// A plane's UV map is affine, and affine maps commute with rational
/// combinations: mapping every control point into UV while keeping
/// degree, knots, and weights yields a pcurve whose trace and
/// parameterization BOTH match the 3-D curve exactly. The generic
/// 23-sample projection walk + cubic interpolation leaves ~1.4e-4 of
/// interpolation sagitta on a r=0.3 rim circle instead — above the
/// stale-pcurve gate, so every repair pass re-flagged the rim forever
/// (`fix_stale_pcurves` could never measure its own recompute better).
/// Same orthogonal-axes projection convention as the arc path; the
/// caller's survival rule verifies the candidate by measurement, so a
/// non-orthogonal plane frame degrades to the sampled fit rather than
/// shipping a wrong exact-looking pcurve.
fn exact_plane_nurbs_pcurve(
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
        CurveGeometry::Nurbs(nurbs),
    ) = (surface, curve)
    else {
        return None;
    };
    let ul2 = u_axis.dot(u_axis);
    let vl2 = v_axis.dot(v_axis);
    if ul2 < 1e-30 || vl2 < 1e-30 {
        return None;
    }
    let mapped: Vec<Point3> = nurbs
        .control_points
        .iter()
        .map(|p| {
            let rel: crate::geometry::Vec3 = *p - Point3::new(origin.x, origin.y, origin.z);
            Point3::new(rel.dot(u_axis) / ul2, rel.dot(v_axis) / vl2, 0.0)
        })
        .collect();
    let full_domain = Pcurve {
        geometry: crate::geometry::pcurve::PcurveGeometry::Nurbs(
            crate::geometry::nurbs_curve::NurbsCurve {
                degree: nurbs.degree,
                knots: nurbs.knots.clone(),
                control_points: mapped,
                weights: nurbs.weights.clone(),
            },
        ),
        t_start: t_lo,
        t_end: t_hi,
    };
    // Normalize to the authored-pcurve structural contract: knot
    // domain == parameter range, forward along the half-edge. The
    // mapped NURBS spans the WHOLE curve domain and may run backward
    // (twin side); downstream consumers (UV walker, split matching)
    // are built against sliced-forward pcurves, and shipping the
    // full-circle domain broke the flush-union chain's face splits
    // (measured: flush_arm/flush_gusset SplitFailed). `to_nurbs`
    // slices and reverses by exact knot operations.
    Some(Pcurve::from_nurbs(full_domain.to_nurbs()))
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
