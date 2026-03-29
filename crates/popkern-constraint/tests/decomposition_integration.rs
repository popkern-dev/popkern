//! Integration tests for graph decomposition: verifying that independent
//! constraint groups are solved correctly when decomposed into sub-problems,
//! that bridging constraints merge components, that disabled groups cause
//! disconnection, and that component failure is isolated.

use popkern_constraint::elements::ElementHandle;
use popkern_constraint::{Constraint, Sketch, SolveStatus};

/// Floating-point comparison tolerance.
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

/// Euclidean distance between two 2D points.
fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

// ── Test 1: Independent groups solve correctly ────────────────────────────────

/// Two fully-constrained rectangles with no shared constraints.
/// Each has Fixed + Horizontal + Vertical + Distance constraints.
/// Both should solve to their correct geometry after a single solve call.
#[test]
fn independent_groups_solve_correctly() {
    let mut sketch = Sketch::new();

    // ── Rectangle A at origin, 6 × 4 ─────────────────────────────────────
    // Initial positions are exact: Fixed captures current coords at call time.
    let a0 = sketch.add_point(0.0, 0.0);
    let a1 = sketch.add_point(5.9, 0.0);
    let a2 = sketch.add_point(6.1, 3.9);
    let a3 = sketch.add_point(0.0, 4.0);

    let a_bottom = sketch.add_line(a0, a1);
    let a_right = sketch.add_line(a1, a2);
    let a_top = sketch.add_line(a2, a3);
    let a_left = sketch.add_line(a3, a0);

    sketch.add_constraint(Constraint::Fixed(a0));
    sketch.add_constraint(Constraint::Horizontal(a_bottom));
    sketch.add_constraint(Constraint::Horizontal(a_top));
    sketch.add_constraint(Constraint::Vertical(a_right));
    sketch.add_constraint(Constraint::Vertical(a_left));
    sketch.add_constraint(Constraint::Distance(a0, a1, 6.0));
    sketch.add_constraint(Constraint::Distance(a0, a3, 4.0));

    // ── Rectangle B far away, 3 × 2 ──────────────────────────────────────
    let b0 = sketch.add_point(50.0, 50.0);
    let b1 = sketch.add_point(52.9, 50.0);
    let b2 = sketch.add_point(53.1, 51.9);
    let b3 = sketch.add_point(50.0, 52.0);

    let b_bottom = sketch.add_line(b0, b1);
    let b_right = sketch.add_line(b1, b2);
    let b_top = sketch.add_line(b2, b3);
    let b_left = sketch.add_line(b3, b0);

    sketch.add_constraint(Constraint::Fixed(b0));
    sketch.add_constraint(Constraint::Horizontal(b_bottom));
    sketch.add_constraint(Constraint::Horizontal(b_top));
    sketch.add_constraint(Constraint::Vertical(b_right));
    sketch.add_constraint(Constraint::Vertical(b_left));
    sketch.add_constraint(Constraint::Distance(b0, b1, 3.0));
    sketch.add_constraint(Constraint::Distance(b0, b3, 2.0));

    let result = sketch.solve().expect("solve");
    assert!(
        matches!(result.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained, got {:?}",
        result.status
    );

    // Verify rectangle A corners.
    let (ax0, ay0) = sketch.point_position(a0).unwrap();
    let (ax1, ay1) = sketch.point_position(a1).unwrap();
    let (ax2, ay2) = sketch.point_position(a2).unwrap();
    let (ax3, ay3) = sketch.point_position(a3).unwrap();

    approx(ax0, 0.0, "A p0.x");
    approx(ay0, 0.0, "A p0.y");
    approx(ax1, 6.0, "A p1.x");
    approx(ay1, 0.0, "A p1.y");
    approx(ax2, 6.0, "A p2.x");
    approx(ay2, 4.0, "A p2.y");
    approx(ax3, 0.0, "A p3.x");
    approx(ay3, 4.0, "A p3.y");

    // Verify rectangle B corners.
    let (bx0, by0) = sketch.point_position(b0).unwrap();
    let (bx1, by1) = sketch.point_position(b1).unwrap();
    let (bx2, by2) = sketch.point_position(b2).unwrap();
    let (bx3, by3) = sketch.point_position(b3).unwrap();

    approx(bx0, 50.0, "B p0.x");
    approx(by0, 50.0, "B p0.y");
    approx(bx1, 53.0, "B p1.x");
    approx(by1, 50.0, "B p1.y");
    approx(bx2, 53.0, "B p2.x");
    approx(by2, 52.0, "B p2.y");
    approx(bx3, 50.0, "B p3.x");
    approx(by3, 52.0, "B p3.y");
}

// ── Test 2: Bridging constraint merges components ─────────────────────────────

/// Two independent fixed points.  Without a bridging constraint they form
/// two separate components.  Adding a Distance constraint between them
/// merges the components; the distance is satisfied after solving.
///
/// To keep the scenario solvable (not over-constrained), one point is fixed
/// and one is free — the free point moves to satisfy the distance.
#[test]
fn bridging_constraint_merges_components() {
    // ── Before bridge: two independent fixed points ───────────────────────
    let mut sketch = Sketch::new();
    let p_fixed = sketch.add_point(0.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(p_fixed));

    let p_free = sketch.add_point(10.0, 0.0);
    sketch.add_constraint(Constraint::OrdinateY(p_free, 0.0)); // keep on x-axis

    // Before bridge, p_free is under-constrained in X.
    let pre = sketch.solve().expect("pre-bridge solve");
    assert!(
        matches!(pre.status, SolveStatus::UnderConstrained(_)),
        "expected UnderConstrained before bridge, got {:?}",
        pre.status
    );

    // ── After bridge: add distance constraint linking both points ─────────
    sketch.add_constraint(Constraint::Distance(p_fixed, p_free, 7.0));

    let result = sketch.solve().expect("bridged solve");
    assert!(
        matches!(result.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained after bridge, got {:?}",
        result.status
    );

    // The bridge must be satisfied: |p_free - p_fixed| == 7.0.
    let d = dist(
        sketch.point_position(p_fixed).unwrap(),
        sketch.point_position(p_free).unwrap(),
    );
    assert!(
        (d - 7.0).abs() < TOL,
        "bridge distance: expected 7.0, got {}",
        d
    );
}

// ── Test 3: Disabled group creates disconnection ──────────────────────────────

/// Two groups of points connected by a bridging constraint placed inside a
/// constraint group.  Disabling the group removes the bridge, causing the two
/// groups to decompose into independent sub-problems.
/// One group is fully constrained; the other is under-constrained after the
/// bridge is removed.
#[test]
fn disabled_group_creates_disconnection() {
    let mut sketch = Sketch::new();

    // Group A: fully constrained fixed point.
    let pa = sketch.add_point(0.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(pa));

    // Group B: a free point with no constraints of its own.
    let pb = sketch.add_point(5.0, 5.0);

    // Bridge: Distance between pa and pb placed in a constraint group.
    let bridge_group = sketch.create_group("bridge");
    let bridge_ch = sketch.add_constraint(Constraint::Distance(pa, pb, 5.0));
    sketch.add_to_group(bridge_group, bridge_ch).unwrap();

    // With bridge enabled the two points form one component.  pb is partially
    // constrained by the distance but still has 1 rotational DOF.
    let result_enabled = sketch.solve().expect("solve with bridge enabled");
    assert!(
        matches!(
            result_enabled.status,
            SolveStatus::FullyConstrained | SolveStatus::UnderConstrained(_)
        ),
        "unexpected status with bridge enabled: {:?}",
        result_enabled.status
    );

    // Disable the bridge group.
    sketch.disable_group(bridge_group);

    // Now pa is a fixed point (its own component, fully constrained) and pb
    // is an unconstrained free point.  The sketch should report UnderConstrained.
    let result_disabled = sketch.solve().expect("solve with bridge disabled");
    assert!(
        matches!(result_disabled.status, SolveStatus::UnderConstrained(_)),
        "expected UnderConstrained with bridge disabled, got {:?}",
        result_disabled.status
    );

    // pa should still be at (0, 0) — its component is unaffected.
    let (xa, ya) = sketch.point_position(pa).unwrap();
    approx(xa, 0.0, "pa.x after bridge disabled");
    approx(ya, 0.0, "pa.y after bridge disabled");
}

// ── Test 4: Regression — single-component identical ──────────────────────────

/// The standard fully-constrained 5×3 rectangle from canonical_sketches.
/// Solving via the decomposed path (single component) must produce the same
/// corner positions and FullyConstrained status.
#[test]
fn single_component_regression() {
    let mut sketch = Sketch::new();

    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(4.0, 0.1);
    let p2 = sketch.add_point(4.1, 2.9);
    let p3 = sketch.add_point(-0.1, 3.1);

    let bottom = sketch.add_line(p0, p1);
    let right = sketch.add_line(p1, p2);
    let top = sketch.add_line(p2, p3);
    let left = sketch.add_line(p3, p0);

    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Horizontal(bottom));
    sketch.add_constraint(Constraint::Horizontal(top));
    sketch.add_constraint(Constraint::Vertical(right));
    sketch.add_constraint(Constraint::Vertical(left));
    sketch.add_constraint(Constraint::Distance(p0, p1, 5.0));
    sketch.add_constraint(Constraint::Distance(p0, p3, 3.0));

    let result = sketch.solve().expect("solve");
    assert!(
        matches!(result.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained, got {:?}",
        result.status
    );

    let (x0, y0) = sketch.point_position(p0).unwrap();
    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    let (x3, y3) = sketch.point_position(p3).unwrap();

    approx(x0, 0.0, "p0.x");
    approx(y0, 0.0, "p0.y");
    approx(x1, 5.0, "p1.x");
    approx(y1, 0.0, "p1.y");
    approx(x2, 5.0, "p2.x");
    approx(y2, 3.0, "p2.y");
    approx(x3, 0.0, "p3.x");
    approx(y3, 3.0, "p3.y");

    // Iteration count should be reasonable (not inflated by decomposition).
    assert!(
        result.iterations <= 20,
        "iteration count suspiciously high: {}",
        result.iterations
    );
}

// ── Test 5: Component failure isolation ───────────────────────────────────────

/// Two independent groups.  Group B is well-constrained (Fixed point).
/// Group A has a conflicting pair of constraints (Horizontal + Vertical on
/// the same non-zero-length line) that makes it over-constrained.
/// After solving, group B's point should still be at its fixed position.
#[test]
fn component_failure_isolation() {
    let mut sketch = Sketch::new();

    // Group A: over-constrained — horizontal + vertical + distance on one line.
    let a0 = sketch.add_point(0.0, 0.0);
    let a1 = sketch.add_point(3.0, 0.0);
    let la = sketch.add_line(a0, a1);
    sketch.add_constraint(Constraint::Fixed(a0));
    sketch.add_constraint(Constraint::Horizontal(la));
    sketch.add_constraint(Constraint::Vertical(la));
    sketch.add_constraint(Constraint::Distance(a0, a1, 3.0));

    // Group B: well-constrained — single fixed point far from group A.
    let b0 = sketch.add_point(100.0, 200.0);
    sketch.add_constraint(Constraint::Fixed(b0));

    let result = sketch.solve().expect("solve should not error");

    // Group A is over/failed; group B is still correctly constrained.
    // The overall result may be OverConstrained or Failed (A drives the outcome).
    match &result.status {
        SolveStatus::OverConstrained(_) | SolveStatus::Failed(_) => {}
        SolveStatus::FullyConstrained => {}
        SolveStatus::UnderConstrained(_) => {
            panic!(
                "unexpected UnderConstrained when A is over-constrained: {:?}",
                result.status
            );
        }
    }

    // Regardless of A's failure, group B's point must remain at (100, 200).
    let (bx, by) = sketch.point_position(b0).unwrap();
    approx(bx, 100.0, "b0.x after A failure");
    approx(by, 200.0, "b0.y after A failure");
}

// ── Test 6: Rigid set spans one component ─────────────────────────────────────

/// A rigid set of 3 points forms one component (3 DOF: tx, ty, theta).
/// An isolated unconstrained point is a separate component (2 DOF).
/// After solving with no constraints, the sketch is under-constrained with
/// at least 2 distinct free elements.
#[test]
fn rigid_set_spans_one_component() {
    let mut sketch = Sketch::new();

    // Rigid set: right triangle.
    let ra = sketch.add_point(0.0, 0.0);
    let rb = sketch.add_point(3.0, 0.0);
    let rc = sketch.add_point(0.0, 4.0);

    sketch
        .create_rigid_set(&[
            ElementHandle::Point(ra),
            ElementHandle::Point(rb),
            ElementHandle::Point(rc),
        ])
        .expect("create rigid set");

    // Isolated free point — no constraints, separate component.
    let _iso = sketch.add_point(50.0, 50.0);

    // Solve with no driving constraints — everything is free.
    let result = sketch.solve().expect("solve");

    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            // The rigid set and the isolated point are separate free elements.
            // Must have at least 2 distinct free elements reported.
            assert!(
                freedoms.len() >= 2,
                "expected at least 2 free elements, got {}",
                freedoms.len()
            );
            // Total reported DOF should be at least 4 (rigid set: tx+ty at
            // minimum; isolated point: 2).
            let total_dof: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert!(
                total_dof >= 4,
                "expected at least 4 DOF (rigid + isolated), got {}",
                total_dof
            );
        }
        other => {
            // No constraints at all → should be UnderConstrained.
            panic!("expected UnderConstrained, got {:?}", other);
        }
    }

    // Rigid set: add external constraints to fully constrain the rigid body
    // (tx, ty, theta = 3 DOF), then verify it solves.
    let mut sketch2 = Sketch::new();
    let sa = sketch2.add_point(0.0, 0.0);
    let sb = sketch2.add_point(3.0, 0.0);
    let sc = sketch2.add_point(0.0, 4.0);

    sketch2
        .create_rigid_set(&[
            ElementHandle::Point(sa),
            ElementHandle::Point(sb),
            ElementHandle::Point(sc),
        ])
        .expect("create rigid set");

    // Fix tx and ty by fixing sa, fix theta by pinning sb.y.
    sketch2.add_constraint(Constraint::Fixed(sa));
    sketch2.add_constraint(Constraint::OrdinateY(sb, 0.0));

    let result2 = sketch2.solve().expect("solve rigid set");
    assert!(
        matches!(result2.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained for rigid set with 3 external constraints, got {:?}",
        result2.status
    );

    // Verify inter-point distances are preserved.
    let d_ab = dist(
        sketch2.point_position(sa).unwrap(),
        sketch2.point_position(sb).unwrap(),
    );
    let d_ac = dist(
        sketch2.point_position(sa).unwrap(),
        sketch2.point_position(sc).unwrap(),
    );
    assert!(
        (d_ab - 3.0).abs() < TOL,
        "rigid AB dist: expected 3.0, got {}",
        d_ab
    );
    assert!(
        (d_ac - 4.0).abs() < TOL,
        "rigid AC dist: expected 4.0, got {}",
        d_ac
    );
}
