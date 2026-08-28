//! Adversarial tests for cycle detection in expression-driven dimensions.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::expression::{BinOp, Expression, UnaryFunc};

/// Direct self-referential expression: A = Ref(A).
#[test]
fn self_referential_cycle_detected() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p1, 1.0));

    // A = Ref(A) — direct self-reference.
    let result = sketch.set_dimension_expression(dist_a, Expression::Ref(dist_a));
    assert!(result.is_err(), "expected circular dependency error");
}

/// Longer cycle: A -> B -> C -> D -> A.
#[test]
fn long_chain_cycle_detected() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);
    let pc = sketch.add_point(3.0, 0.0);
    let pd = sketch.add_point(4.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 2.0));
    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 3.0));
    let dist_d = sketch.add_constraint(Constraint::Distance(p0, pd, 4.0));

    // Build chain: A -> B -> C -> D
    sketch
        .set_dimension_expression(dist_a, Expression::Ref(dist_b))
        .unwrap();
    sketch
        .set_dimension_expression(dist_b, Expression::Ref(dist_c))
        .unwrap();
    sketch
        .set_dimension_expression(dist_c, Expression::Ref(dist_d))
        .unwrap();

    // Now close the cycle: D -> A. Should fail.
    let result = sketch.set_dimension_expression(dist_d, Expression::Ref(dist_a));
    assert!(
        result.is_err(),
        "expected circular dependency error for A->B->C->D->A"
    );
}

/// Multiple independent subgraphs, one containing a cycle.
/// A -> B (no cycle), C -> D -> C (cycle). Setting D -> C should fail
/// without affecting the A -> B chain.
#[test]
fn independent_component_with_cycle() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);
    let pc = sketch.add_point(3.0, 0.0);
    let pd = sketch.add_point(4.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 2.0));
    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 3.0));
    let dist_d = sketch.add_constraint(Constraint::Distance(p0, pd, 4.0));

    // Independent chain: A = Ref(B) — no cycle.
    sketch
        .set_dimension_expression(dist_a, Expression::Ref(dist_b))
        .unwrap();

    // Start second chain: C = Ref(D) — no cycle yet.
    sketch
        .set_dimension_expression(dist_c, Expression::Ref(dist_d))
        .unwrap();

    // Close cycle in second component: D = Ref(C) — should fail.
    let result = sketch.set_dimension_expression(dist_d, Expression::Ref(dist_c));
    assert!(
        result.is_err(),
        "expected circular dependency error for C->D->C"
    );

    // The A -> B chain should still be intact and solvable.
    assert!(
        sketch.expressions.contains_key(&dist_a),
        "A -> B expression should still exist"
    );
}

/// Diamond dependency (A->B, A->C, B->D, C->D) is NOT a cycle.
/// This tests the `has_cycle` function with a DAG that has reconvergent
/// paths, ensuring the visited-set short-circuit doesn't false-positive.
#[test]
fn diamond_dependency_is_not_a_cycle() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);
    let pc = sketch.add_point(3.0, 0.0);
    let pd = sketch.add_point(6.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let la = sketch.add_line(p0, pa);
    let lb = sketch.add_line(p0, pb);
    let lc = sketch.add_line(p0, pc);
    let ld = sketch.add_line(p0, pd);
    sketch.add_constraint(Constraint::Horizontal(la));
    sketch.add_constraint(Constraint::Horizontal(lb));
    sketch.add_constraint(Constraint::Horizontal(lc));
    sketch.add_constraint(Constraint::Horizontal(ld));

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 1.0));
    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 1.0));
    let dist_d = sketch.add_constraint(Constraint::Distance(p0, pd, 1.0));

    // B = Ref(D), C = Ref(D), A = Ref(B) + Ref(C)
    // Graph: D <- B <- A, D <- C <- A
    // This is a diamond: A -> {B, C} -> D. No cycle.
    sketch
        .set_dimension_expression(dist_b, Expression::Ref(dist_d))
        .unwrap();
    sketch
        .set_dimension_expression(dist_c, Expression::Ref(dist_d))
        .unwrap();
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_b)),
                BinOp::Add,
                Box::new(Expression::Ref(dist_c)),
            ),
        )
        .expect("diamond dependency should NOT be detected as a cycle");

    // Solve should work: D=1, B=1, C=1, A=1+1=2.
    sketch.solve().unwrap();
    let (xa, _) = sketch.point_position(pa).unwrap();
    assert!((xa - 2.0).abs() < 1e-6, "expected pa.x=2, got {xa}");
}

/// Cycle hidden inside a BinOp: A = Ref(B) + Ref(A).
#[test]
fn cycle_hidden_in_binop() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 2.0));

    // A = Ref(B) + Ref(A) — self-reference through a BinOp.
    let result = sketch.set_dimension_expression(
        dist_a,
        Expression::BinOp(
            Box::new(Expression::Ref(dist_b)),
            BinOp::Add,
            Box::new(Expression::Ref(dist_a)),
        ),
    );
    assert!(
        result.is_err(),
        "expected circular dependency error for self-ref in BinOp"
    );
}

/// Cycle hidden inside a Func: A = Sqrt(Ref(B)), B = Neg(Ref(A)).
#[test]
fn cycle_hidden_in_func() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 2.0));

    // A = Sqrt(Ref(B))
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::Func(UnaryFunc::Sqrt, Box::new(Expression::Ref(dist_b))),
        )
        .unwrap();

    // B = Neg(Ref(A)) — closes cycle through Func nodes.
    let result = sketch.set_dimension_expression(
        dist_b,
        Expression::Func(UnaryFunc::Neg, Box::new(Expression::Ref(dist_a))),
    );
    assert!(
        result.is_err(),
        "expected circular dependency error for cycle through Func"
    );
}

/// Long linear chain with no cycle: A -> B -> C -> D -> E.
/// All should succeed.
#[test]
fn long_chain_no_cycle() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(2.0, 0.0);
    let pc = sketch.add_point(3.0, 0.0);
    let pd = sketch.add_point(4.0, 0.0);
    let pe = sketch.add_point(5.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let la = sketch.add_line(p0, pa);
    let lb = sketch.add_line(p0, pb);
    let lc = sketch.add_line(p0, pc);
    let ld = sketch.add_line(p0, pd);
    let le = sketch.add_line(p0, pe);
    sketch.add_constraint(Constraint::Horizontal(la));
    sketch.add_constraint(Constraint::Horizontal(lb));
    sketch.add_constraint(Constraint::Horizontal(lc));
    sketch.add_constraint(Constraint::Horizontal(ld));
    sketch.add_constraint(Constraint::Horizontal(le));

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 1.0));
    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 1.0));
    let dist_d = sketch.add_constraint(Constraint::Distance(p0, pd, 1.0));
    let dist_e = sketch.add_constraint(Constraint::Distance(p0, pe, 5.0));

    // Chain: A = Ref(B), B = Ref(C), C = Ref(D), D = Ref(E).
    sketch
        .set_dimension_expression(dist_a, Expression::Ref(dist_b))
        .unwrap();
    sketch
        .set_dimension_expression(dist_b, Expression::Ref(dist_c))
        .unwrap();
    sketch
        .set_dimension_expression(dist_c, Expression::Ref(dist_d))
        .unwrap();
    sketch
        .set_dimension_expression(dist_d, Expression::Ref(dist_e))
        .expect("linear chain should NOT be a cycle");

    // E=5, D=5, C=5, B=5, A=5.
    sketch.solve().unwrap();
    let (xa, _) = sketch.point_position(pa).unwrap();
    assert!((xa - 5.0).abs() < 1e-6, "expected pa.x=5, got {xa}");
}

/// Sqrt of negative value: A = Sqrt(Literal(-1.0)).
#[test]
fn sqrt_negative_fails() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p1, 1.0));

    sketch
        .set_dimension_expression(
            dist_a,
            Expression::Func(UnaryFunc::Sqrt, Box::new(Expression::Literal(-1.0))),
        )
        .unwrap();

    let result = sketch.solve();
    assert!(
        result.is_err(),
        "expected error from sqrt of negative number"
    );
}
