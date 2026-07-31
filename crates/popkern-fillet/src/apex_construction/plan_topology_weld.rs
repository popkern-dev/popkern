//! Weld-point geometry for the shared-corner plan classification
//! (PK-219): the contact-tangent-line closest approach and the
//! clamped-curve endpoint snap.

use popkern_core::geometry::{NurbsCurve, Point3};

/// Closest-approach midpoint of the two contact tangent lines, or
/// `None` when a tangent is missing, the lines are near-parallel, or
/// the weld point strays beyond the short-bridge band from either
/// termination (the contact lines meet at the support's corner
/// angle, so |X - T| ≈ gap / (2·sin(angle/2)) legitimately exceeds
/// the gap itself; the BAND bounds how far a re-termination may
/// move). The corner stays Bridged on decline — never a guessed
/// weld.
pub(super) fn contact_line_weld_point(
    t_b: Point3,
    dir_b: Option<popkern_core::geometry::Vec3>,
    t_a: Point3,
    dir_a: Option<popkern_core::geometry::Vec3>,
    gap: f64,
    band: f64,
) -> Option<Point3> {
    let d1 = dir_b?;
    let d2 = dir_a?;
    let w = t_b - t_a;
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d2.dot(&d2);
    let d = d1.dot(&w);
    let e = d2.dot(&w);
    let denom = a * c - b * b;
    if denom.abs() <= 1.0e-9 {
        return None;
    }
    let s = (b * e - c * d) / denom;
    let t = (a * e - b * d) / denom;
    let p1 = t_b + d1 * s;
    let p2 = t_a + d2 * t;
    let point = Point3::new(
        0.5 * (p1.x + p2.x),
        0.5 * (p1.y + p2.y),
        0.5 * (p1.z + p2.z),
    );
    let stray = band.max(4.0 * gap) + 1.0e-12;
    (point.distance_to(&t_b) <= stray && point.distance_to(&t_a) <= stray).then_some(point)
}

/// Snap a clamped NURBS curve's first/last control points to the
/// given endpoints so the curve agrees bitwise with the (possibly
/// welded) corner vertices. The displacement is bounded by the weld
/// gate (≤ 10·tol / 2 per endpoint) — below every downstream
/// endpoint-match gate; the arc interior is unchanged.
pub(super) fn snap_clamped_curve_endpoints(curve: &mut NurbsCurve, start: &Point3, end: &Point3) {
    if let Some(first) = curve.control_points.first_mut() {
        *first = *start;
    }
    if let Some(last) = curve.control_points.last_mut() {
        *last = *end;
    }
}
