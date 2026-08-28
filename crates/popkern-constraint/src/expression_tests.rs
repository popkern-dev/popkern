//! Tests for expression-driven dimensions.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::expression::{BinOp, Expression, UnaryFunc};

/// Distance A has expression `2 * Ref(B)`. B.target=5 -> A.target=10
/// before solving, and geometry reflects it.
#[test]
fn basic_expression_mul_ref() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(5.0, 0.0);
    let p3 = sketch.add_point(10.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p1));
    let dist_b = sketch.add_constraint(Constraint::Distance(p1, p2, 5.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p1, p3, 1.0)); // placeholder
    let l1 = sketch.add_line(p1, p2);
    sketch.add_constraint(Constraint::Horizontal(l1));
    let l2 = sketch.add_line(p1, p3);
    sketch.add_constraint(Constraint::Horizontal(l2));

    // A = 2 * Ref(B)
    let expr = Expression::BinOp(
        Box::new(Expression::Literal(2.0)),
        BinOp::Mul,
        Box::new(Expression::Ref(dist_b)),
    );
    sketch
        .set_dimension_expression(dist_a, expr)
        .expect("should succeed");

    sketch.solve().expect("should solve");

    let (x3, _) = sketch.point_position(p3).unwrap();
    assert!((x3 - 10.0).abs() < 1e-6, "expected p3.x=10, got {x3}");
}

/// Chained: A = Ref(B) + 1, B = Ref(C) * 2. C.target=3 -> B=6 -> A=7.
#[test]
fn chained_expression() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(7.0, 0.0);
    let pb = sketch.add_point(6.0, 0.0);
    let pc = sketch.add_point(3.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let la = sketch.add_line(p0, pa);
    let lb = sketch.add_line(p0, pb);
    let lc = sketch.add_line(p0, pc);
    sketch.add_constraint(Constraint::Horizontal(la));
    sketch.add_constraint(Constraint::Horizontal(lb));
    sketch.add_constraint(Constraint::Horizontal(lc));

    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 3.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 1.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));

    // B = Ref(C) * 2
    sketch
        .set_dimension_expression(
            dist_b,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_c)),
                BinOp::Mul,
                Box::new(Expression::Literal(2.0)),
            ),
        )
        .unwrap();

    // A = Ref(B) + 1
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_b)),
                BinOp::Add,
                Box::new(Expression::Literal(1.0)),
            ),
        )
        .unwrap();

    sketch.solve().unwrap();

    let (xb, _) = sketch.point_position(pb).unwrap();
    let (xa, _) = sketch.point_position(pa).unwrap();
    assert!((xb - 6.0).abs() < 1e-6, "expected pb.x=6, got {xb}");
    assert!((xa - 7.0).abs() < 1e-6, "expected pa.x=7, got {xa}");
}

/// Circular reference: A = Ref(B), B = Ref(A) -> error on
/// set_dimension_expression.
#[test]
fn circular_reference_detected() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 0.0);
    let p2 = sketch.add_point(2.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p1, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, p2, 2.0));

    // A = Ref(B)
    sketch
        .set_dimension_expression(dist_a, Expression::Ref(dist_b))
        .unwrap();

    // B = Ref(A) -> should fail with circular dependency.
    let result = sketch.set_dimension_expression(dist_b, Expression::Ref(dist_a));
    assert!(result.is_err(), "expected circular dependency error");
}

/// Diamond dependency: A = Ref(C), B = Ref(C), D = Ref(A) + Ref(B).
/// No cycle — should work.
#[test]
fn diamond_dependency_no_cycle() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(1.0, 0.0);
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

    let dist_c = sketch.add_constraint(Constraint::Distance(p0, pc, 3.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 1.0));
    let dist_d = sketch.add_constraint(Constraint::Distance(p0, pd, 1.0));

    // A = Ref(C)
    sketch
        .set_dimension_expression(dist_a, Expression::Ref(dist_c))
        .unwrap();
    // B = Ref(C)
    sketch
        .set_dimension_expression(dist_b, Expression::Ref(dist_c))
        .unwrap();
    // D = Ref(A) + Ref(B)
    sketch
        .set_dimension_expression(
            dist_d,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_a)),
                BinOp::Add,
                Box::new(Expression::Ref(dist_b)),
            ),
        )
        .unwrap();

    sketch.solve().unwrap();

    let (xd, _) = sketch.point_position(pd).unwrap();
    // D = A + B = 3 + 3 = 6
    assert!((xd - 6.0).abs() < 1e-6, "expected pd.x=6, got {xd}");
}

/// Unary function: A = Sqrt(Ref(B)). B.target=25 -> A.target=5.
#[test]
fn unary_sqrt_expression() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(5.0, 0.0);
    let pb = sketch.add_point(25.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let la = sketch.add_line(p0, pa);
    let lb = sketch.add_line(p0, pb);
    sketch.add_constraint(Constraint::Horizontal(la));
    sketch.add_constraint(Constraint::Horizontal(lb));

    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 25.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));

    // A = Sqrt(Ref(B))
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::Func(UnaryFunc::Sqrt, Box::new(Expression::Ref(dist_b))),
        )
        .unwrap();

    sketch.solve().unwrap();

    let (xa, _) = sketch.point_position(pa).unwrap();
    assert!((xa - 5.0).abs() < 1e-6, "expected pa.x=5, got {xa}");
}

/// Division by zero: A = Ref(B) / Ref(C). C.target=0 -> error.
#[test]
fn division_by_zero_error() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let pa = sketch.add_point(1.0, 0.0);
    let pb = sketch.add_point(5.0, 0.0);

    let dist_a = sketch.add_constraint(Constraint::Distance(p0, pa, 1.0));
    let dist_b = sketch.add_constraint(Constraint::Distance(p0, pb, 5.0));

    // Create a dummy zero-value constraint. Use OrdinateX with value 0.
    let pz = sketch.add_point(0.0, 0.0);
    let dist_c = sketch.add_constraint(Constraint::OrdinateX(pz, 0.0));

    // A = Ref(B) / Ref(C) where C=0
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_b)),
                BinOp::Div,
                Box::new(Expression::Ref(dist_c)),
            ),
        )
        .unwrap();

    let result = sketch.solve();
    assert!(result.is_err(), "expected error from division by zero");
}

/// Ref on suppressed constraint: suppress B, set A = 2*Ref(B). Ref still
/// resolves (reads stored target, not affected by suppression).
#[test]
fn ref_on_suppressed_constraint() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(5.0, 0.0);
    let p2 = sketch.add_point(10.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let l1 = sketch.add_line(p0, p1);
    let l2 = sketch.add_line(p0, p2);
    sketch.add_constraint(Constraint::Horizontal(l1));
    sketch.add_constraint(Constraint::Horizontal(l2));

    let dist_b = sketch.add_constraint(Constraint::Distance(p0, p1, 5.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p2, 1.0));

    // Suppress dist_b by putting it in a disabled group.
    let grp = sketch.create_group("suppressed");
    sketch.add_to_group(grp, dist_b).unwrap();
    sketch.disable_group(grp);

    // A = 2 * Ref(B) — B is suppressed but the value 5.0 is still stored.
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Literal(2.0)),
                BinOp::Mul,
                Box::new(Expression::Ref(dist_b)),
            ),
        )
        .unwrap();

    sketch.solve().unwrap();

    let (x2, _) = sketch.point_position(p2).unwrap();
    assert!((x2 - 10.0).abs() < 1e-6, "expected p2.x=10, got {x2}");
}

/// Clear expression: set an expression, clear it, verify the constant
/// value persists.
#[test]
fn clear_expression_reverts_to_constant() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(3.0, 0.0);
    let p2 = sketch.add_point(1.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let l1 = sketch.add_line(p0, p1);
    let l2 = sketch.add_line(p0, p2);
    sketch.add_constraint(Constraint::Horizontal(l1));
    sketch.add_constraint(Constraint::Horizontal(l2));

    let dist_b = sketch.add_constraint(Constraint::Distance(p0, p1, 3.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p2, 1.0));

    // Set expression: A = Ref(B) * 2 = 6
    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Ref(dist_b)),
                BinOp::Mul,
                Box::new(Expression::Literal(2.0)),
            ),
        )
        .unwrap();

    sketch.solve().unwrap();
    let (xa, _) = sketch.point_position(p2).unwrap();
    assert!((xa - 6.0).abs() < 1e-6, "expected 6, got {xa}");

    // Clear expression. The value was written as 6.0 by evaluate_all.
    sketch.clear_dimension_expression(dist_a);

    // Now set the underlying value to 4.0 manually.
    sketch.set_constraint_value(dist_a, 4.0).unwrap();
    sketch.solve().unwrap();

    let (xa2, _) = sketch.point_position(p2).unwrap();
    assert!(
        (xa2 - 4.0).abs() < 1e-6,
        "after clear+set to 4, expected 4, got {xa2}"
    );
}

/// Serialization: Expression with Ref round-trips through serde_json.
#[test]
fn expression_serde_roundtrip() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(5.0, 0.0);
    let p2 = sketch.add_point(10.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(p0));
    let l1 = sketch.add_line(p0, p1);
    let l2 = sketch.add_line(p0, p2);
    sketch.add_constraint(Constraint::Horizontal(l1));
    sketch.add_constraint(Constraint::Horizontal(l2));

    let dist_b = sketch.add_constraint(Constraint::Distance(p0, p1, 5.0));
    let dist_a = sketch.add_constraint(Constraint::Distance(p0, p2, 1.0));

    sketch
        .set_dimension_expression(
            dist_a,
            Expression::BinOp(
                Box::new(Expression::Literal(2.0)),
                BinOp::Mul,
                Box::new(Expression::Ref(dist_b)),
            ),
        )
        .unwrap();

    // Serialize and deserialize the whole sketch.
    let json = serde_json::to_string(&sketch).expect("serialize");
    let mut sketch2: Sketch = serde_json::from_str(&json).expect("deserialize");

    // The expression should survive the round-trip.
    assert!(
        sketch2.expressions.contains_key(&dist_a),
        "expression should be present after deserialization"
    );

    // And solving the deserialized sketch should produce the same result.
    sketch2.solve().unwrap();
    let (x2, _) = sketch2.point_position(p2).unwrap();
    assert!((x2 - 10.0).abs() < 1e-6, "expected 10, got {x2}");
}

/// Setting expression on non-dimensional constraint should fail.
#[test]
fn expression_on_non_dimensional_fails() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(1.0, 0.0);
    let p2 = sketch.add_point(2.0, 0.0);

    let dist_b = sketch.add_constraint(Constraint::Distance(p0, p2, 2.0));
    let coincident = sketch.add_constraint(Constraint::Coincident(p0, p1));

    let result = sketch.set_dimension_expression(coincident, Expression::Ref(dist_b));
    assert!(result.is_err(), "expected NotDimensional error");
}

/// evaluate_expression preview without storing.
#[test]
fn evaluate_expression_preview() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(7.0, 0.0);

    let dist = sketch.add_constraint(Constraint::Distance(p0, p1, 7.0));

    let expr = Expression::BinOp(
        Box::new(Expression::Ref(dist)),
        BinOp::Add,
        Box::new(Expression::Literal(3.0)),
    );

    let result = sketch.evaluate_expression(&expr).unwrap();
    assert!((result - 10.0).abs() < 1e-10, "expected 10, got {result}");

    // No expression stored.
    assert!(
        sketch.expressions.is_empty(),
        "evaluate_expression should not store expression"
    );
}
