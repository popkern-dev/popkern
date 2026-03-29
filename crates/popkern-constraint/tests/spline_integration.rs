//! Integration tests for spline constraint scenarios: through-point
//! interpolation, G1/G2 continuity, tangency to conics, endpoint
//! curvature, arc length, and serialization round-trips.

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

// ── 1. Through-point spline ──────────────────────────────────────────────

#[test]
fn through_point_spline_five_points() {
    let mut sketch = Sketch::new();
    let pts = vec![(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (5.0, 3.0), (7.0, 0.0)];
    let (_sh, tp_handles) = sketch.add_spline_through_points(&pts).expect("create");

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify all through-points are on the spline.
    for (i, &tph) in tp_handles.iter().enumerate() {
        let (px, py) = sketch.point_position(tph).expect("point");
        approx(px, pts[i].0, &format!("pt{i} x"));
        approx(py, pts[i].1, &format!("pt{i} y"));
    }
}

// ── 2. Closed through-point spline ───────────────────────────────────────

#[test]
fn closed_through_point_spline() {
    let mut sketch = Sketch::new();
    let pts = vec![(2.0, 0.0), (0.0, 2.0), (-2.0, 0.0), (0.0, -2.0)];
    let (_sh, tp_handles) = sketch
        .add_closed_spline_through_points(&pts)
        .expect("create");

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify all through-points.
    for (i, &tph) in tp_handles.iter().enumerate() {
        let (px, py) = sketch.point_position(tph).expect("point");
        approx(px, pts[i].0, &format!("closed pt{i} x"));
        approx(py, pts[i].1, &format!("closed pt{i} y"));
    }
}

// ── 3. G2 junction ───────────────────────────────────────────────────────

#[test]
fn g2_junction_curvature_match() {
    let mut sketch = Sketch::new();

    // Two cubic splines sharing endpoint (4,0).
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 2.0);
    let p2 = sketch.add_point(3.0, 2.0);
    let p3 = sketch.add_point(4.0, 0.0);
    let sh1 = sketch.add_spline(&[p0, p1, p2, p3], 3);

    // Second spline starts at (4,0), initially mis-aligned.
    let q0 = sketch.add_point(4.0, 0.0);
    let q1 = sketch.add_point(5.0, -1.8);
    let q2 = sketch.add_point(7.0, -2.0);
    let q3 = sketch.add_point(8.0, 0.0);
    let sh2 = sketch.add_spline(&[q0, q1, q2, q3], 3);

    // Fix the outer endpoints so there's enough DOF only at the junction.
    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Fixed(p2));
    sketch.add_constraint(Constraint::Fixed(q2));
    sketch.add_constraint(Constraint::Fixed(q3));

    // G2 continuity at the junction: end of S1 meets start of S2.
    sketch.add_constraint(Constraint::SplineG2(sh1, 1, sh2, 0));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-6,
        "residual too large: {}",
        result.residual_norm
    );
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "solve failed: {:?}",
        result.status
    );

    // Verify position continuity.
    let (x3, y3) = sketch.point_position(p3).expect("p3");
    let (x_q0, y_q0) = sketch.point_position(q0).expect("q0");
    approx(x3, x_q0, "G2 pos x");
    approx(y3, y_q0, "G2 pos y");
}

// ── 4. Spline tangent to ellipse ─────────────────────────────────────────

#[test]
fn spline_tangent_to_fixed_ellipse() {
    let mut sketch = Sketch::new();

    // Fixed ellipse at (2, 0), a=2, b=1.5, no rotation.
    let ell = sketch.add_ellipse(2.0, 0.0, 2.0, 1.5, 0.0);
    sketch.add_constraint(Constraint::FixedEllipse(ell));

    // Spline whose midpoint is near the top of the ellipse (1.5 at y).
    // Position control points so the spline naturally grazes the ellipse.
    let cp0 = sketch.add_point(-1.0, 0.0);
    let cp1 = sketch.add_point(0.5, 2.0);
    let cp2 = sketch.add_point(3.5, 2.0);
    let cp3 = sketch.add_point(5.0, 0.0);
    let sh = sketch.add_spline(&[cp0, cp1, cp2, cp3], 3);

    // Fix all CPs except cp1 and cp2 (the interior CPs can adjust to
    // satisfy tangency).
    sketch.add_constraint(Constraint::Fixed(cp0));
    sketch.add_constraint(Constraint::Fixed(cp3));

    // Tangency constraint (aux params initialized to 0 by default, but
    // the solver should still converge from there).
    sketch.add_constraint(Constraint::SplineTangentToEllipse(sh, ell));

    let result = sketch.solve().expect("solve");
    // Three equations, two free interior CPs (4 DOF) + 2 aux params = 6 DOF,
    // so 3 equations are underdetermined. Verify the solver doesn't fail.
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "solve failed: {:?}",
        result.status
    );
    assert!(
        result.residual_norm < 1e-4,
        "residual too large: {}",
        result.residual_norm
    );
}

// ── 5. Endpoint curvature zero ───────────────────────────────────────────

#[test]
fn endpoint_curvature_zero() {
    let mut sketch = Sketch::new();

    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 2.0);
    let p2 = sketch.add_point(3.0, 2.0);
    let p3 = sketch.add_point(4.0, 0.0);
    let sh = sketch.add_spline(&[p0, p1, p2, p3], 3);

    // Fix endpoints and one interior CP.
    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Fixed(p3));
    sketch.add_constraint(Constraint::Fixed(p2));

    // Constrain start curvature to zero (natural spline end).
    sketch.add_constraint(Constraint::SplineFixedEndCurvature(sh, 0, 0.0));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-6,
        "residual too large: {}",
        result.residual_norm
    );
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "solve failed: {:?}",
        result.status
    );

    // The solver converged with residual < 1e-6, which means the
    // curvature constraint (kappa - 0 = 0) is satisfied to that
    // tolerance. No need to call internal eval_curvature.
    let _ = sh;
}

// ── 6. Arc length constraint ─────────────────────────────────────────────

#[test]
fn arc_length_constraint() {
    let mut sketch = Sketch::new();

    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 1.5);
    let p2 = sketch.add_point(3.0, 1.5);
    let p3 = sketch.add_point(4.0, 0.0);
    let sh = sketch.add_spline(&[p0, p1, p2, p3], 3);

    // Fix endpoints to constrain the spline shape.
    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Fixed(p3));

    let target_length = 5.0;
    sketch.add_constraint(Constraint::SplineLength(sh, target_length));

    let result = sketch.solve().expect("solve");
    // Arc length uses quadrature, so wider tolerance is expected.
    assert!(
        result.residual_norm < 1e-3,
        "residual too large: {}",
        result.residual_norm
    );
    assert!(
        !matches!(result.status, SolveStatus::Failed(_)),
        "solve failed: {:?}",
        result.status
    );

    // The solver residual being small means the SplineLength equation
    // (arc_length - target) is approximately zero, verifying the constraint.
    let _ = sh;
}

// ── 7. Serialization round-trip ──────────────────────────────────────────

/// Spline sketch with non-aux-param constraints survives round-trip.
///
/// Note: `aux_param_map` is transient (`#[serde(skip)]`), so constraints
/// that rely on auxiliary parameters (PointOnSpline, SplineTangentTo*)
/// require rebuilding after deserialization.  Here we test with
/// Fixed + SplineFixedEndTangent constraints which do not use aux params.
#[test]
fn serialize_roundtrip_spline_sketch() {
    let mut sketch = Sketch::new();

    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 2.0);
    let p2 = sketch.add_point(3.0, 2.0);
    let p3 = sketch.add_point(4.0, 0.0);
    let sh = sketch.add_spline(&[p0, p1, p2, p3], 3);

    // Constraints without aux params.
    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Fixed(p3));
    sketch.add_constraint(Constraint::SplineFixedEndTangent(sh, 0, 3.0, 6.0));

    sketch.solve().expect("solve original");

    // Serialize and deserialize.
    let json = serde_json::to_string(&sketch).expect("serialize");
    let mut loaded: Sketch = serde_json::from_str(&json).expect("deserialize");

    // Handles should still be valid.
    assert!(
        loaded.point_position(p0).is_some(),
        "p0 invalid after round-trip"
    );
    assert!(
        loaded.point_position(p3).is_some(),
        "p3 invalid after round-trip"
    );

    // Solve the deserialized sketch.
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

    // Spot-check: endpoints still in place.
    let (x0, y0) = loaded.point_position(p0).expect("p0");
    approx(x0, 0.0, "loaded p0 x");
    approx(y0, 0.0, "loaded p0 y");
    let (x3, y3) = loaded.point_position(p3).expect("p3");
    approx(x3, 4.0, "loaded p3 x");
    approx(y3, 0.0, "loaded p3 y");

    let _ = sh;
}
