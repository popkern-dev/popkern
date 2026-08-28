//! Integration tests: groups, rigid sets, configurations, and expressions
//! working together.

use popkern_constraint::{
    Sketch,
    constraint::Constraint,
    elements::ElementHandle,
    error::SolveStatus,
    expression::{BinOp, Expression},
};

const TOL: f64 = 1e-6;

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

// ── Test 1: Group + solve ─────────────────────────────────────────────────────

/// Two points with a Horizontal constraint in a "layout" group plus a Distance
/// constraint outside the group.  The Distance alone does not fix orientation;
/// without Horizontal the line can rotate, leaving 1 DOF.  With Horizontal it
/// is fully constrained.
#[test]
fn group_and_solve_rectangle() {
    let mut s = Sketch::new();

    let p0 = s.add_point(0.0, 0.0);
    let p1 = s.add_point(5.0, 0.0);

    // Fix the origin.
    s.add_constraint(Constraint::Fixed(p0));

    // Distance is always active — it constrains length but not direction.
    s.add_constraint(Constraint::Distance(p0, p1, 5.0));

    // Horizontal constraint lives in the "layout" group.
    let l = s.add_line(p0, p1);
    let grp = s.create_group("layout");
    let horiz_ch = s.add_constraint(Constraint::Horizontal(l));
    s.add_to_group(grp, horiz_ch).unwrap();

    // With group disabled the line can rotate — sketch is under-constrained.
    s.disable_group(grp);
    let result = s.solve().unwrap();
    assert!(
        matches!(result.status, SolveStatus::UnderConstrained(_)),
        "expected UnderConstrained with group disabled, got {:?}",
        result.status
    );

    // Re-enable the group — Horizontal fixes orientation, now fully constrained.
    s.enable_group(grp);
    let result2 = s.solve().unwrap();
    assert!(
        matches!(result2.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained with group enabled, got {:?}",
        result2.status
    );

    // Verify geometry: p1 should be at (5, 0).
    let (x1, y1) = s.point_position(p1).unwrap();
    assert!((x1 - 5.0).abs() < TOL, "p1.x={x1}, expected 5.0");
    assert!(y1.abs() < TOL, "p1.y={y1}, expected 0.0");
}

// ── Test 2: Rigid set + solve ─────────────────────────────────────────────────

/// 3-point right triangle made rigid. External distance drives the whole set.
/// All inter-point distances must be preserved.
#[test]
fn rigid_set_and_solve() {
    let mut s = Sketch::new();

    // Right triangle at origin.
    let a = s.add_point(0.0, 0.0);
    let b = s.add_point(3.0, 0.0);
    let c = s.add_point(0.0, 4.0);

    let d_ab = dist(s.point_position(a).unwrap(), s.point_position(b).unwrap());
    let d_ac = dist(s.point_position(a).unwrap(), s.point_position(c).unwrap());
    let d_bc = dist(s.point_position(b).unwrap(), s.point_position(c).unwrap());

    let _rsh = s
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    // External fixed point at (10, 0).
    let ext = s.add_point(10.0, 0.0);
    s.add_constraint(Constraint::Fixed(ext));

    // Distance from ext to set point a = 5.0 (places a at x=5 or x=15).
    s.add_constraint(Constraint::Distance(ext, a, 5.0));

    // Fix rotation by pinning a's y coordinate.
    s.add_constraint(Constraint::OrdinateY(a, 0.0));

    let result = s.solve().unwrap();
    assert!(
        result.residual_norm < 1e-7,
        "residual too large: {}",
        result.residual_norm
    );

    // Inter-point distances must be unchanged.
    let d_ab_new = dist(s.point_position(a).unwrap(), s.point_position(b).unwrap());
    let d_ac_new = dist(s.point_position(a).unwrap(), s.point_position(c).unwrap());
    let d_bc_new = dist(s.point_position(b).unwrap(), s.point_position(c).unwrap());

    assert!(
        (d_ab_new - d_ab).abs() < TOL,
        "AB dist changed: {d_ab_new} vs {d_ab}"
    );
    assert!(
        (d_ac_new - d_ac).abs() < TOL,
        "AC dist changed: {d_ac_new} vs {d_ac}"
    );
    assert!(
        (d_bc_new - d_bc).abs() < TOL,
        "BC dist changed: {d_bc_new} vs {d_bc}"
    );
}

// ── Test 3: Configuration switching ──────────────────────────────────────────

/// Two configurations "narrow" (dist=5) and "wide" (dist=10). Switching
/// between them restores geometry correctly.
#[test]
fn configuration_switching() {
    let mut s = Sketch::new();

    let p1 = s.add_point(0.0, 0.0);
    let p2 = s.add_point(5.0, 0.0);
    s.add_constraint(Constraint::Fixed(p1));
    let l = s.add_line(p1, p2);
    s.add_constraint(Constraint::Horizontal(l));
    let dist_c = s.add_constraint(Constraint::Distance(p1, p2, 5.0));

    s.solve().unwrap();

    // Save "narrow" config at dist=5.
    let cfg_narrow = s.create_config("narrow");
    s.save_config(cfg_narrow).unwrap();

    // Change to dist=10.
    s.set_constraint_value(dist_c, 10.0).unwrap();
    s.solve().unwrap();

    let (x2_wide, _) = s.point_position(p2).unwrap();
    assert!(
        (x2_wide - 10.0).abs() < TOL,
        "wide: p2.x={x2_wide}, expected 10.0"
    );

    // Save "wide" config at dist=10.
    let cfg_wide = s.create_config("wide");
    s.save_config(cfg_wide).unwrap();

    // Activate "narrow" — distance should revert to 5.0.
    s.activate_config(cfg_narrow).unwrap();
    let (x2_narrow, _) = s.point_position(p2).unwrap();
    assert!(
        (x2_narrow - 5.0).abs() < TOL,
        "narrow: p2.x={x2_narrow}, expected 5.0"
    );

    // Activate "wide" again — distance should be 10.0.
    s.activate_config(cfg_wide).unwrap();
    let (x2_w2, _) = s.point_position(p2).unwrap();
    assert!(
        (x2_w2 - 10.0).abs() < TOL,
        "wide (re-activate): p2.x={x2_w2}, expected 10.0"
    );
}

// ── Test 4: Expression-driven dimensions ─────────────────────────────────────

/// span=20, gap=Ref(span)/4=5. Changing span to 40 makes gap=10.
#[test]
fn expression_driven_dimensions() {
    let mut s = Sketch::new();

    // Two separate point pairs, both horizontal from origin.
    let origin = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(origin));

    let span_end = s.add_point(20.0, 0.0);
    let gap_end = s.add_point(5.0, 1.0);

    let l_span = s.add_line(origin, span_end);
    let l_gap = s.add_line(origin, gap_end);
    s.add_constraint(Constraint::Horizontal(l_span));
    s.add_constraint(Constraint::Horizontal(l_gap));

    let span_h = s.add_constraint(Constraint::Distance(origin, span_end, 20.0));
    let gap_h = s.add_constraint(Constraint::Distance(origin, gap_end, 1.0));

    // gap = Ref(span) / 4
    s.set_dimension_expression(
        gap_h,
        Expression::BinOp(
            Box::new(Expression::Ref(span_h)),
            BinOp::Div,
            Box::new(Expression::Literal(4.0)),
        ),
    )
    .unwrap();

    s.solve().unwrap();

    let (xg, _) = s.point_position(gap_end).unwrap();
    assert!(
        (xg - 5.0).abs() < TOL,
        "gap_end.x={xg}, expected 5.0 (span=20 / 4)"
    );

    // Change span to 40, re-solve — gap should become 10.
    s.set_constraint_value(span_h, 40.0).unwrap();
    s.solve().unwrap();

    let (xg2, _) = s.point_position(gap_end).unwrap();
    assert!(
        (xg2 - 10.0).abs() < TOL,
        "gap_end.x={xg2}, expected 10.0 (span=40 / 4)"
    );
}

// ── Test 5: Rigid set + group combo ──────────────────────────────────────────

/// Rigid set of 3 points with an external Distance constraint inside a
/// group. Disabling the group leaves the set with 3 DOF. Enabling it and
/// solving drives the set to the target distance.
#[test]
fn rigid_set_and_group_combo() {
    let mut s = Sketch::new();

    let a = s.add_point(0.0, 0.0);
    let b = s.add_point(1.0, 0.0);
    let c = s.add_point(0.0, 1.0);

    let _rsh = s
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    // External anchor point.
    let anchor = s.add_point(10.0, 0.0);
    s.add_constraint(Constraint::Fixed(anchor));

    // External distance inside a group.
    let grp = s.create_group("ext_dist");
    let dist_ch = s.add_constraint(Constraint::Distance(anchor, a, 3.0));
    // Also constrain y so rotation is determined.
    let ory_ch = s.add_constraint(Constraint::OrdinateY(a, 0.0));
    s.add_to_group(grp, dist_ch).unwrap();
    s.add_to_group(grp, ory_ch).unwrap();

    // Group disabled — rigid set should have 3 DOF (tx, ty, theta).
    s.disable_group(grp);
    let result = s.solve().unwrap();
    assert!(
        matches!(result.status, SolveStatus::UnderConstrained(_)),
        "expected UnderConstrained with group disabled, got {:?}",
        result.status
    );

    // Enable group — solver should drive the rigid body to satisfy the distance.
    s.enable_group(grp);
    let result2 = s.solve().unwrap();
    assert!(
        result2.residual_norm < 1e-7,
        "residual after enable: {}",
        result2.residual_norm
    );

    // Verify distance anchor -> a ≈ 3.0.
    let actual_dist = dist(
        s.point_position(anchor).unwrap(),
        s.point_position(a).unwrap(),
    );
    assert!(
        (actual_dist - 3.0).abs() < TOL,
        "dist(anchor, a)={actual_dist}, expected 3.0"
    );
}

// ── Test 6: Config + group interaction ───────────────────────────────────────

/// Saving one config with group enabled, another with group disabled.
/// Switching between them restores the group state.
#[test]
fn config_and_group_interaction() {
    let mut s = Sketch::new();

    let p0 = s.add_point(0.0, 0.0);
    let p1 = s.add_point(3.0, 0.0);
    s.add_constraint(Constraint::Fixed(p0));
    let l = s.add_line(p0, p1);
    s.add_constraint(Constraint::Horizontal(l));

    let grp = s.create_group("toggle");
    let dist_ch = s.add_constraint(Constraint::Distance(p0, p1, 8.0));
    s.add_to_group(grp, dist_ch).unwrap();

    // Save "enabled" config (group is still enabled at this point).
    s.solve().unwrap();
    let cfg_enabled = s.create_config("enabled");
    s.save_config(cfg_enabled).unwrap();

    // Disable group; save "disabled" config.
    s.disable_group(grp);
    s.solve().unwrap();
    let cfg_disabled = s.create_config("disabled");
    s.save_config(cfg_disabled).unwrap();

    // Activate enabled config — group should be enabled again.
    s.activate_config(cfg_enabled).unwrap();
    assert_eq!(
        s.is_group_enabled(grp),
        Some(true),
        "group should be enabled after activating cfg_enabled"
    );

    // Activate disabled config — group should be disabled.
    s.activate_config(cfg_disabled).unwrap();
    assert_eq!(
        s.is_group_enabled(grp),
        Some(false),
        "group should be disabled after activating cfg_disabled"
    );

    // Activate enabled again to verify toggle works both ways.
    s.activate_config(cfg_enabled).unwrap();
    assert_eq!(
        s.is_group_enabled(grp),
        Some(true),
        "group should be enabled again after second activation of cfg_enabled"
    );
    let (x1, _) = s.point_position(p1).unwrap();
    assert!(
        (x1 - 8.0).abs() < TOL,
        "p1.x={x1}, expected 8.0 after cfg_enabled re-activate"
    );
}

// ── Test 7: Expression + config ───────────────────────────────────────────────

/// A = 2 * Ref(B). Save a config. Change B manually. Save another config.
/// Switching back restores B.target, and A re-evaluates correctly on solve.
#[test]
fn expression_and_config() {
    let mut s = Sketch::new();

    let origin = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(origin));

    let pb = s.add_point(5.0, 0.0);
    let pa = s.add_point(10.0, 0.0);

    let lb = s.add_line(origin, pb);
    let la = s.add_line(origin, pa);
    s.add_constraint(Constraint::Horizontal(lb));
    s.add_constraint(Constraint::Horizontal(la));

    let b_h = s.add_constraint(Constraint::Distance(origin, pb, 5.0));
    let a_h = s.add_constraint(Constraint::Distance(origin, pa, 1.0));

    // A = 2 * Ref(B)
    s.set_dimension_expression(
        a_h,
        Expression::BinOp(
            Box::new(Expression::Literal(2.0)),
            BinOp::Mul,
            Box::new(Expression::Ref(b_h)),
        ),
    )
    .unwrap();

    s.solve().unwrap();

    // Verify initial state: A = 2*5 = 10.
    let (xa, _) = s.point_position(pa).unwrap();
    assert!((xa - 10.0).abs() < TOL, "initial pa.x={xa}, expected 10.0");

    // Save "cfg_b5" at B=5.
    let cfg_b5 = s.create_config("b5");
    s.save_config(cfg_b5).unwrap();

    // Change B to 8 and solve — A should become 16.
    s.set_constraint_value(b_h, 8.0).unwrap();
    s.solve().unwrap();
    let (xa2, _) = s.point_position(pa).unwrap();
    assert!(
        (xa2 - 16.0).abs() < TOL,
        "after B=8: pa.x={xa2}, expected 16.0"
    );

    // Save "cfg_b8".
    let cfg_b8 = s.create_config("b8");
    s.save_config(cfg_b8).unwrap();

    // Activate cfg_b5 — B.target should revert to 5, A re-evaluates to 10.
    s.activate_config(cfg_b5).unwrap();
    let (xa3, _) = s.point_position(pa).unwrap();
    assert!(
        (xa3 - 10.0).abs() < TOL,
        "after cfg_b5 re-activate: pa.x={xa3}, expected 10.0"
    );

    // Activate cfg_b8 — B.target should revert to 8, A re-evaluates to 16.
    s.activate_config(cfg_b8).unwrap();
    let (xa4, _) = s.point_position(pa).unwrap();
    assert!(
        (xa4 - 16.0).abs() < TOL,
        "after cfg_b8 re-activate: pa.x={xa4}, expected 16.0"
    );
}
