//! Oblique plane-cone section: exact rational conic construction.
//!
//! Split from `plane_cone.rs` (structural 500-line cap). See
//! [`case_oblique`] for the shoulder-point construction notes.

use super::*;

// ── Oblique (ellipse) ───────────────────────────────────────────────────

/// Compute the 3D point on the cone-plane intersection at cone angle `u`.
///
/// The cone parameterization is:
///   P(u,v) = apex + axis*v + perp1*(v*tan_ha*cos(u)) + perp2*(v*tan_ha*sin(u))
///
/// The plane constraint is: (P - po) · pn = 0
///
/// Solving for v:  v(u) = c / (a + b*cos(u) + d*sin(u))
/// where:
///   c = (po - apex) · pn
///   a = ca · pn
///   b = perp1 · pn * tan_ha
///   d = perp2 · pn * tan_ha
#[allow(clippy::too_many_arguments)]
fn oblique_point(
    apex: Point3,
    perp1: Vec3,
    perp2: Vec3,
    ca: Vec3,
    tan_ha: f64,
    c: f64,
    a: f64,
    b: f64,
    d: f64,
    u: f64,
) -> Point3 {
    let denom = a + b * u.cos() + d * u.sin();
    let v = c / denom;
    let r = v * tan_ha;
    apex + ca * v + perp1 * (r * u.cos()) + perp2 * (r * u.sin())
}

/// General oblique case: intersection is an ellipse on the cone.
///
/// The curve keeps the cardinal cone-angle phase (quarters at u = 0,
/// π/2, π, 3π/2, winding about the cone axis — the same start/seam and
/// orientation as the perpendicular circle), and each quarter is an
/// EXACT rational conic arc: the middle control point is the
/// intersection of the endpoint tangent lines, and the weight comes
/// from the shoulder point where the true section crosses the median
/// from the chord midpoint to that intersection (w = s/(1−s); for a
/// circle quarter this reproduces 1/√2). The earlier `P + 1·T`
/// construction was not a valid rational conic and deviated from the
/// cone by more than a millimetre on shallow countersink
/// configurations — the PK-154 intersecting-wires defect.
///
/// Returns `None` for hyperbola/parabola (the denominator can reach
/// zero) and for degenerate arc geometry, signaling fallback to
/// marching.
#[allow(clippy::too_many_arguments)]
pub(super) fn case_oblique(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    tol: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (perp1, perp2) = orthonormal_basis(ca);
    let tan_ha = half_angle.tan();

    let c = (po - apex).dot(&pn);
    let a = ca.dot(&pn);
    let b = perp1.dot(&pn) * tan_ha;
    let d = perp2.dot(&pn) * tan_ha;

    // Ellipse condition: a² > b² + d² (denominator never zero).
    let disc = a * a - b * b - d * d;
    if disc < tol.angular() * tol.angular() {
        // Parabola or hyperbola: fall through to marching.
        return None;
    }

    // Check that v is positive for all u (intersection on the real cone half).
    // v(u) = c / (a + b*cos(u) + d*sin(u)).
    // The denominator ranges in [a - sqrt(b²+d²), a + sqrt(b²+d²)].
    let amp = (b * b + d * d).sqrt();
    let v_at_min = c / (a - amp);
    let v_at_max = c / (a + amp);
    if v_at_min < tol.linear() || v_at_max < tol.linear() {
        // The intersection is on the "wrong" side of the apex or
        // partially behind it.  Return empty.
        return Some(vec![]);
    }

    let pt = |u: f64| oblique_point(apex, perp1, perp2, ca, tan_ha, c, a, b, d, u);

    // Tangent at cone angle u (lies in the cutting plane):
    // dP/du = dv/du * (axis + perp1*tan_ha*cos(u) + perp2*tan_ha*sin(u))
    //       + v * (perp1*(-tan_ha*sin(u)) + perp2*(tan_ha*cos(u)))
    let tang = |u: f64| -> Vec3 {
        let denom = a + b * u.cos() + d * u.sin();
        let v = c / denom;
        let dv_du = c * (b * u.sin() - d * u.cos()) / (denom * denom);
        let from_v_change = (ca + perp1 * (tan_ha * u.cos()) + perp2 * (tan_ha * u.sin())) * dv_du;
        let from_u_change = (perp1 * (-tan_ha * u.sin()) + perp2 * (tan_ha * u.cos())) * v;
        from_v_change + from_u_change
    };

    // In-plane 2D frame for the tangent-line intersection.
    let e1 = pn.cross(&ca).normalized().unwrap_or_else(|| {
        // pn ∥ ca is the perpendicular case, unreachable here, but stay
        // total: any in-plane direction works.
        orthonormal_basis(pn).0
    });
    let e2 = pn.cross(&e1);
    let to_2d = |x: Point3| -> (f64, f64) {
        let r = x - po;
        (r.dot(&e1), r.dot(&e2))
    };

    let cos2 = half_angle.cos() * half_angle.cos();
    // One root of the double-cone quadratic along X(s) = m + s·dir,
    // s ∈ (0, 1): f > 0 inside the section, f < 0 outside, so the
    // bracketed sign change has exactly one root.
    let shoulder_fraction = |m: Point3, dir: Vec3| -> Option<f64> {
        let wm = m - apex;
        let alpha = wm.dot(&ca);
        let beta = dir.dot(&ca);
        let qa = beta * beta - cos2 * dir.dot(&dir);
        let qb = 2.0 * (alpha * beta - cos2 * wm.dot(&dir));
        let qc = alpha * alpha - cos2 * wm.dot(&wm);
        let discriminant = qb * qb - 4.0 * qa * qc;
        if discriminant < 0.0 {
            return None;
        }
        let sqrt_disc = discriminant.sqrt();
        let candidates = if qa.abs() < 1.0e-300 {
            [(-qc / qb), f64::NAN]
        } else {
            [
                (-qb + sqrt_disc) / (2.0 * qa),
                (-qb - sqrt_disc) / (2.0 * qa),
            ]
        };
        candidates
            .into_iter()
            .find(|s| s.is_finite() && *s > 0.0 && *s < 1.0)
    };

    use std::f64::consts::FRAC_PI_2;
    let mut control_points: Vec<Point3> = Vec::with_capacity(9);
    let mut weights: Vec<f64> = Vec::with_capacity(9);
    control_points.push(pt(0.0));
    weights.push(1.0);
    for k in 0..4u32 {
        let u0 = f64::from(k) * FRAC_PI_2;
        let u2 = u0 + FRAC_PI_2;
        let p0 = pt(u0);
        let p2 = if k == 3 { control_points[0] } else { pt(u2) };
        let t0 = tang(u0);
        let t2 = tang(u2);
        // Tangent-line intersection in the cutting plane.
        let (x0, y0) = to_2d(p0);
        let (x2, y2) = to_2d(p2);
        let (dx0, dy0) = (t0.dot(&e1), t0.dot(&e2));
        let (dx2, dy2) = (t2.dot(&e1), t2.dot(&e2));
        let det = dx0 * (-dy2) - (-dx2) * dy0;
        if det.abs() < 1.0e-14 {
            // Degenerate quarter (near-parabolic): let marching serve it.
            return None;
        }
        let rhs = (x2 - x0, y2 - y0);
        let s0 = (rhs.0 * (-dy2) - (-dx2) * rhs.1) / det;
        let p1 = p0 + t0 * s0;
        let m = Point3::new(
            0.5 * (p0.x + p2.x),
            0.5 * (p0.y + p2.y),
            0.5 * (p0.z + p2.z),
        );
        let s = shoulder_fraction(m, p1 - m)?;
        let w = s / (1.0 - s);
        control_points.push(p1);
        weights.push(w);
        control_points.push(p2);
        weights.push(1.0);
    }

    let curve_3d = NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points,
        weights,
    };

    // UV on plane: project all control points (affine map — the
    // rational conic maps exactly).
    let uv_plane = {
        let cps: Vec<Point3> = curve_3d
            .control_points
            .iter()
            .map(|p| plane_uv_pt(*p, po, pu, pv))
            .collect();
        NurbsCurve {
            degree: 2,
            knots: curve_3d.knots.clone(),
            control_points: cps,
            weights: curve_3d.weights.clone(),
        }
    };

    let uv_cone = cone_uv_curve(&curve_3d, apex, ca);

    let (uv1, uv2) = if swapped {
        (uv_cone, uv_plane)
    } else {
        (uv_plane, uv_cone)
    };

    Some(vec![SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1: uv1,
        curve_uv2: uv2,
        is_closed: true,
    })])
}
