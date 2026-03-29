//! Integration tests for conic constraint scenarios: point-on-conic,
//! tangent-line, tangent-arc, eccentricity, mixed sketches, and
//! serialization round-trips.

use popkern_constraint::{Constraint, Sketch, SolveStatus};

/// Tolerance for floating-point comparisons in solver output.
const TOL: f64 = 1e-6;

/// Assert that two floats are within [`TOL`] of each other.
fn approx(a: f64, b: f64, label: &str) {
    assert!(
        (a - b).abs() < TOL,
        "{}: expected {}, got {} (delta {})",
        label,
        b,
        a,
        (a - b).abs()
    );
}

// ── 1. Constrained ellipse: point moves onto ellipse ───────────────────────

#[test]
fn constrained_ellipse_point_on_ellipse() {
    let mut sketch = Sketch::new();

    // Ellipse at origin with a=3, b=2, no rotation.
    let ellipse = sketch.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ellipse));

    // Point initially off the ellipse — start close to the right-hand side
    // so the solver converges to a unique solution.
    let pt = sketch.add_point(3.5, 0.5);
    sketch.add_constraint(Constraint::PointOnEllipse(pt, ellipse));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "solve failed: {:?}",
        result.status
    );

    let (x, y) = sketch.point_position(pt).unwrap();
    let (_, _, a, b, _) = sketch.ellipse_params(ellipse).unwrap();

    // Verify the point lies on the ellipse: (x/a)^2 + (y/b)^2 = 1.
    let on_ellipse = (x / a).powi(2) + (y / b).powi(2);
    approx(on_ellipse, 1.0, "(x/a)^2 + (y/b)^2");
}

// ── 2. Parabolic mirror: point constrained to parabola ────────────────────

#[test]
fn parabolic_mirror_point_on_parabola() {
    let mut sketch = Sketch::new();

    // Parabola with vertex at origin, focal distance 1.0, no rotation.
    // Equation in local frame: v^2 = 4*f*u.
    let parabola = sketch.add_parabola(0.0, 0.0, 1.0, 0.0);
    sketch.add_constraint(Constraint::FixedParabola(parabola));

    // Point not on the parabola (at (2, 1), v^2 != 4fu). Start near (2, 2.5).
    let pt = sketch.add_point(2.0, 2.5);
    sketch.add_constraint(Constraint::PointOnParabola(pt, parabola));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    let (x, y) = sketch.point_position(pt).unwrap();
    let (_, _, f, _) = sketch.parabola_params(parabola).unwrap();

    // Verify v^2 = 4*f*u (parabola has zero rotation, so global == local).
    approx(y * y, 4.0 * f * x, "v^2 == 4*f*u");
}

// ── 3. Tangent line to ellipse ─────────────────────────────────────────────

#[test]
fn tangent_line_to_ellipse() {
    let mut sketch = Sketch::new();

    // Fixed ellipse at origin: a=3, b=2, no rotation.
    let ellipse = sketch.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ellipse));

    // Line starting roughly horizontal near the top of the ellipse (y~2).
    let p1 = sketch.add_point(-4.0, 2.1);
    let p2 = sketch.add_point(4.0, 1.9);
    let line = sketch.add_line(p1, p2);
    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::TangentLineEllipse(line, ellipse));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify the dual conic tangency condition (rotation=0, so global == local).
    // Line l*u + m*v + n = 0 is tangent to x^2/a^2 + y^2/b^2 = 1 iff
    // a^2*l^2 + b^2*m^2 = n^2.
    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    let (_, _, a, b, _) = sketch.ellipse_params(ellipse).unwrap();

    let dx = x2 - x1;
    let dy = y2 - y1;
    // Line: l*u + m*v + n = 0 where l=dy, m=-dx, n=dx*y1-dy*x1.
    let l = dy;
    let m = -dx;
    let n = dx * y1 - dy * x1;

    // Normalise by line length for numerical stability.
    let len = (dx * dx + dy * dy).sqrt();
    let (ln, mn, nn) = (l / len, m / len, n / len);

    let dual = a * a * ln * ln + b * b * mn * mn - nn * nn;
    approx(dual, 0.0, "dual conic tangency");
}

// ── 4. Arc tangent to ellipse ──────────────────────────────────────────────

#[test]
fn arc_tangent_to_ellipse() {
    let mut sketch = Sketch::new();

    // Fixed ellipse at origin: a=2, b=1, no rotation.
    let ellipse = sketch.add_ellipse(0.0, 0.0, 2.0, 1.0, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ellipse));

    // Arc of radius 0.5 initially near (2.5, 0) — close to the right tangent
    // point. TangentArcEllipse provides 2 equations; with cx and r free the
    // system is fully determined.
    let arc = sketch.add_arc(2.5, 0.0, 0.5, 0.0, std::f64::consts::TAU);
    sketch.add_constraint(Constraint::TangentArcEllipse(arc, ellipse));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify: minimum distance from arc centre to the ellipse equals r.
    let (acx, acy, r, _, _) = sketch.arc_params(arc).unwrap();
    let (ecx, ecy, a, b, _) = sketch.ellipse_params(ellipse).unwrap();

    let min_dist = (0..1000)
        .map(|i| {
            let t = i as f64 * std::f64::consts::TAU / 1000.0;
            let ex = ecx + a * t.cos();
            let ey = ecy + b * t.sin();
            ((acx - ex).powi(2) + (acy - ey).powi(2)).sqrt()
        })
        .fold(f64::MAX, f64::min);

    assert!(
        (min_dist - r).abs() < 1e-3,
        "expected min_dist ≈ r={}, got min_dist={}",
        r,
        min_dist
    );
}

// ── 5. Eccentricity constraint on ellipse ─────────────────────────────────

#[test]
fn ellipse_eccentricity_constraint() {
    let mut sketch = Sketch::new();

    // Ellipse at origin with a=2, b=1.
    let ellipse = sketch.add_ellipse(0.0, 0.0, 2.0, 1.0, 0.0);

    // EllipseEccentricity (1 eq) + FocalDistance (1 eq) together determine
    // a and b uniquely given e=0.5 and focal_dist=1.
    let e = 0.5_f64;
    sketch.add_constraint(Constraint::EllipseEccentricity(ellipse, e));
    sketch.add_constraint(Constraint::FocalDistance(ellipse, 1.0));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    let (_, _, a, b, _) = sketch.ellipse_params(ellipse).unwrap();

    // Verify a^2*(1 - e^2) = b^2.
    approx(b * b, a * a * (1.0 - e * e), "a^2*(1-e^2) == b^2");
}

// ── 6. Hyperbola eccentricity constraint ───────────────────────────────────

#[test]
fn hyperbola_eccentricity_constraint() {
    let mut sketch = Sketch::new();

    // Hyperbola at origin with a=2, b=1.
    let hyperbola = sketch.add_hyperbola(0.0, 0.0, 2.0, 1.0, 0.0);

    // Target eccentricity e=1.5 (valid: e > 1).
    // Constraint: b^2 = a^2*(e^2 - 1).
    let target_e = 1.5_f64;
    sketch.add_constraint(Constraint::HyperbolaEccentricity(hyperbola, target_e));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    let (_, _, a, b, _) = sketch.hyperbola_params(hyperbola).unwrap();

    // Verify b^2 = a^2*(e^2 - 1).
    approx(
        b * b,
        a * a * (target_e * target_e - 1.0),
        "b^2 == a^2*(e^2-1)",
    );
}

// ── 7. Mixed sketch: point on ellipse + tangent line ──────────────────────

#[test]
fn mixed_point_on_ellipse_and_tangent_line() {
    let mut sketch = Sketch::new();

    // Fixed ellipse at origin, a=3, b=2.
    let ellipse = sketch.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ellipse));

    // Point constrained to ellipse — start near the top.
    let pt = sketch.add_point(0.1, 2.1);
    sketch.add_constraint(Constraint::PointOnEllipse(pt, ellipse));

    // Line tangent to the ellipse; fix one endpoint to anchor direction.
    let lp1 = sketch.add_point(-5.0, 2.0);
    let lp2 = sketch.add_point(5.0, 2.0);
    let line = sketch.add_line(lp1, lp2);
    sketch.add_constraint(Constraint::Fixed(lp1));
    sketch.add_constraint(Constraint::TangentLineEllipse(line, ellipse));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify the point is on the ellipse.
    let (px, py) = sketch.point_position(pt).unwrap();
    let (_, _, a, b, _) = sketch.ellipse_params(ellipse).unwrap();
    approx((px / a).powi(2) + (py / b).powi(2), 1.0, "point on ellipse");

    // Verify the line satisfies the dual conic condition.
    let (x1, y1) = sketch.point_position(lp1).unwrap();
    let (x2, y2) = sketch.point_position(lp2).unwrap();
    let dx = x2 - x1;
    let dy = y2 - y1;
    let l = dy;
    let m = -dx;
    let n = dx * y1 - dy * x1;
    let len = (dx * dx + dy * dy).sqrt();
    let (ln, mn, nn) = (l / len, m / len, n / len);
    approx(
        a * a * ln * ln + b * b * mn * mn - nn * nn,
        0.0,
        "tangent line dual conic",
    );
}

// ── 8. Serialization round-trip with conic constraints ────────────────────

#[test]
fn serialize_roundtrip_conic_sketch() {
    let mut sketch = Sketch::new();

    let ellipse = sketch.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ellipse));
    let pt_e = sketch.add_point(2.5, 0.8);
    sketch.add_constraint(Constraint::PointOnEllipse(pt_e, ellipse));

    let parabola = sketch.add_parabola(5.0, 0.0, 1.0, 0.0);
    sketch.add_constraint(Constraint::FixedParabola(parabola));
    let pt_p = sketch.add_point(6.0, 2.1);
    sketch.add_constraint(Constraint::PointOnParabola(pt_p, parabola));

    let hyperbola = sketch.add_hyperbola(-6.0, 0.0, 2.0, 1.0, 0.0);
    sketch.add_constraint(Constraint::FixedHyperbola(hyperbola));
    let pt_h = sketch.add_point(-4.5, 0.1);
    sketch.add_constraint(Constraint::PointOnHyperbola(pt_h, hyperbola));

    sketch.solve().expect("solve original");

    // Serialize and deserialize.
    let json = serde_json::to_string(&sketch).expect("serialize");
    let mut loaded: Sketch = serde_json::from_str(&json).expect("deserialize");

    // Verify handles are still valid (element counts preserved via serde).
    assert!(
        loaded.ellipse_params(ellipse).is_some(),
        "ellipse handle invalid after round-trip"
    );
    assert!(
        loaded.parabola_params(parabola).is_some(),
        "parabola handle invalid after round-trip"
    );
    assert!(
        loaded.hyperbola_params(hyperbola).is_some(),
        "hyperbola handle invalid after round-trip"
    );
    assert!(loaded.point_position(pt_e).is_some(), "pt_e invalid");
    assert!(loaded.point_position(pt_p).is_some(), "pt_p invalid");
    assert!(loaded.point_position(pt_h).is_some(), "pt_h invalid");

    // Solve the deserialized sketch — should converge.
    let result = loaded.solve().expect("solve loaded");
    assert!(
        result.residual_norm < 1e-8,
        "loaded sketch residual too large: {}",
        result.residual_norm
    );
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "loaded sketch solve failed: {:?}",
        result.status
    );

    // Spot-check: the point on the ellipse is still on the ellipse.
    let (ex, ey) = loaded.point_position(pt_e).unwrap();
    let (_, _, a, b, _) = loaded.ellipse_params(ellipse).unwrap();
    approx(
        (ex / a).powi(2) + (ey / b).powi(2),
        1.0,
        "loaded: point on ellipse",
    );
}
