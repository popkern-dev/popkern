//! Exact, verified strip end-section curves for runout rewrites.
//!
//! A cylinder fillet strip terminating on a planar end face must end
//! in the exact cylinder ∩ plane section between its two tangency
//! feet: a circle arc when the end plane is perpendicular to the
//! strip axis, an ELLIPSE arc when it is oblique (PK-262 — the
//! oblique case was a recorded decline until pyramid3's base-edge
//! strips reached it: their far runouts kept an incoherent
//! complementary-side curve whose folded trim polygon manufactured a
//! SelfIntersection false positive).
//!
//! Both authors verify the candidate before returning it: endpoints
//! round-trip onto the snapped feet, and the mid-sample lies on the
//! cylinder AND on the end plane, all at 1e-9. Sweep selection is the
//! short (≤ π) side — the caller's runout detection already excludes
//! near-π spans, and the strip side of a sub-π foot pair is the short
//! way (the antipodal map is preserved by the ellipse
//! reparameterization, so the ≤ π rule carries over from the circle
//! case exactly).

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Round-trip acceptance for a verified section candidate.
const SECTION_ACCEPT: f64 = 1.0e-9;

/// Exact circle section at a perpendicular end plane: feet snapped
/// radially onto the intersection circle, both frame orientations
/// tried (a reflected frame silently traces the wrong sweep).
/// Returns the verified curve plus the snapped feet.
#[allow(clippy::too_many_arguments)]
pub(super) fn verified_circle_section(
    origin: &Point3,
    axis: &Vec3,
    radius: f64,
    c0: Point3,
    p1: Point3,
    p2: Point3,
    plane_origin: &Point3,
    n_hat: &Vec3,
) -> Option<(CurveOnEdge, Point3, Point3)> {
    let snap = |p: Point3| -> Option<Point3> {
        let r_vec = (p - c0) - *axis * (p - c0).dot(axis);
        Some(c0 + r_vec.normalized()? * radius)
    };
    let s1 = snap(p1)?;
    let s2 = snap(p2)?;
    let a = s1 - c0;
    let b = s2 - c0;
    for arc_normal in [*axis, *axis * -1.0] {
        let (f1, f2) = popkern_core::geometry::analytic_curves::circle_basis(arc_normal);
        let theta1 = a.dot(&f2).patan2(a.dot(&f1));
        let mut theta2 = b.dot(&f2).patan2(b.dot(&f1));
        while theta2 <= theta1 {
            theta2 += std::f64::consts::TAU;
        }
        if theta2 - theta1 > std::f64::consts::PI + SECTION_ACCEPT {
            continue;
        }
        let cand = CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center: c0,
                normal: arc_normal,
                radius,
                start_angle: theta1,
                end_angle: theta2,
            },
            t_start: theta1,
            t_end: theta2,
        };
        let ceval: &dyn CurveEval = &cand.geometry;
        let q0 = ceval.point_at(theta1);
        let q1 = ceval.point_at(theta2);
        let qm = ceval.point_at(0.5 * (theta1 + theta2));
        if q0.distance_to(&s1) < SECTION_ACCEPT
            && q1.distance_to(&s2) < SECTION_ACCEPT
            && (radial_dist(origin, axis, &qm) - radius).abs() < SECTION_ACCEPT
            && (qm - *plane_origin).dot(n_hat).abs() < SECTION_ACCEPT
        {
            return Some((cand, s1, s2));
        }
    }
    None
}

/// Exact ellipse section at an oblique end plane. The section of the
/// cylinder by a plane with unit normal `n̂` (with `axis·n̂ =
/// axis_cos`, `0 < |axis_cos| < 1`) is the ellipse centred on the
/// axis's piercing point `c0`, semi-major `r/|axis_cos|` along the
/// axis's in-plane projection, semi-minor `r` along the in-plane
/// perpendicular — every parameter point of that ellipse lies exactly
/// on the cylinder (the major elongation exactly cancels the axial
/// tilt). Feet are projected into the plane and normalized onto the
/// section in basis coefficients; both sweep orientations tried.
#[allow(clippy::too_many_arguments)]
pub(super) fn verified_ellipse_section(
    origin: &Point3,
    axis: &Vec3,
    radius: f64,
    c0: Point3,
    p1: Point3,
    p2: Point3,
    plane_origin: &Point3,
    n_hat: &Vec3,
    axis_cos: f64,
) -> Option<(CurveOnEdge, Point3, Point3)> {
    let major_dir = (*axis - *n_hat * axis_cos).normalized()?;
    let minor_dir = n_hat.cross(&major_dir);
    let u_basis = major_dir * (radius / axis_cos.abs());
    // Unit-circle coefficients of a foot against the section basis,
    // after projecting the foot into the end plane.
    let coefficients = |p: Point3| -> Option<(f64, f64)> {
        let in_plane = p - *n_hat * (p - c0).dot(n_hat);
        let d = in_plane - c0;
        let a_u = d.dot(&major_dir) * (axis_cos.abs() / radius);
        let a_v = d.dot(&minor_dir) / radius;
        let magnitude = (a_u * a_u + a_v * a_v).sqrt();
        (magnitude > 1.0e-12).then(|| (a_u / magnitude, a_v / magnitude))
    };
    let (cu1, cv1) = coefficients(p1)?;
    let (cu2, cv2) = coefficients(p2)?;
    let s1 = c0 + u_basis * cu1 + minor_dir * (radius * cv1);
    let s2 = c0 + u_basis * cu2 + minor_dir * (radius * cv2);
    for flip in [1.0, -1.0] {
        let v_basis = minor_dir * (radius * flip);
        let theta1 = (cv1 * flip).patan2(cu1);
        let mut theta2 = (cv2 * flip).patan2(cu2);
        while theta2 <= theta1 {
            theta2 += std::f64::consts::TAU;
        }
        if theta2 - theta1 > std::f64::consts::PI + SECTION_ACCEPT {
            continue;
        }
        let cand = CurveOnEdge {
            geometry: CurveGeometry::Ellipse {
                center: c0,
                u_basis,
                v_basis,
                start_angle: theta1,
                end_angle: theta2,
            },
            t_start: theta1,
            t_end: theta2,
        };
        let ceval: &dyn CurveEval = &cand.geometry;
        let q0 = ceval.point_at(theta1);
        let q1 = ceval.point_at(theta2);
        let qm = ceval.point_at(0.5 * (theta1 + theta2));
        if q0.distance_to(&s1) < SECTION_ACCEPT
            && q1.distance_to(&s2) < SECTION_ACCEPT
            && (radial_dist(origin, axis, &qm) - radius).abs() < SECTION_ACCEPT
            && (qm - *plane_origin).dot(n_hat).abs() < SECTION_ACCEPT
        {
            return Some((cand, s1, s2));
        }
    }
    None
}

fn radial_dist(origin: &Point3, axis: &Vec3, p: &Point3) -> f64 {
    let d = *p - *origin;
    (d - *axis * d.dot(axis)).length()
}

#[cfg(test)]
#[path = "strip_end_section_curves_tests.rs"]
mod tests;
