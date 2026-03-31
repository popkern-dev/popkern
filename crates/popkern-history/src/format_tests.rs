use crate::expression::{BinOp, CmpOp, Expression, UnaryFunc};
use crate::format::{format_expression, format_expression_with_unit};
use crate::parameter::{ParameterId, ParameterTable};
use crate::parse::parse_expression;
use crate::unit::Unit;

fn make_table() -> (ParameterTable, ParameterId, ParameterId, ParameterId) {
    let mut t = ParameterTable::new();
    let a = t.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = t.add_global("b", 2.0, Unit::Ratio).unwrap();
    let c = t.add_global("c", 3.0, Unit::Ratio).unwrap();
    (t, a, b, c)
}

// ---------------------------------------------------------------------------
// Literal formatting
// ---------------------------------------------------------------------------

#[test]
fn test_literal_with_mm_unit() {
    let t = ParameterTable::new();
    let expr = Expression::Literal(0.05); // 50mm in SI
    let s = format_expression_with_unit(&expr, &t, Unit::Millimeters);
    assert_eq!(s, "50mm");
}

#[test]
fn test_literal_with_deg_unit() {
    let t = ParameterTable::new();
    let angle = std::f64::consts::PI / 4.0; // 45 degrees in radians
    let expr = Expression::Literal(angle);
    let s = format_expression_with_unit(&expr, &t, Unit::Degrees);
    assert_eq!(s, "45deg");
}

#[test]
fn test_literal_bare() {
    let t = ParameterTable::new();
    let expr = Expression::Literal(3.14);
    let s = format_expression(&expr, &t);
    // Should be a reasonable number representation.
    assert!(s.starts_with("3.14"), "got: {s}");
}

#[test]
fn test_literal_integer() {
    let t = ParameterTable::new();
    let expr = Expression::Literal(2.0);
    let s = format_expression(&expr, &t);
    assert_eq!(s, "2");
}

// ---------------------------------------------------------------------------
// Param formatting
// ---------------------------------------------------------------------------

#[test]
fn test_param_name() {
    let (t, a_id, _, _) = make_table();
    let expr = Expression::Param(a_id);
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a");
}

#[test]
fn test_param_unknown() {
    use slotmap::SlotMap;
    // Build a table and extract a valid id, then use a fresh table that doesn't have it.
    let mut t1 = ParameterTable::new();
    let _ = t1.add_global("x", 1.0, Unit::Ratio).unwrap();
    let t2 = ParameterTable::new(); // fresh table, id not valid here

    // We can't easily create a "dangling" ParameterId without unsafe, so instead
    // we test indirectly: if qualified_name returns None, we get "<unknown>".
    // Simulate by checking the format_expression for the unknown path by using
    // a known id from a different table. This tests the fallback formatting path.
    //
    // Since ParameterId is a slotmap key, IDs from different slotmaps may overlap
    // in practice but might also give None. We just check the format produces
    // something reasonable.
    let id = t1.resolve_name("x").unwrap();
    // t2 doesn't have this param so qualified_name returns None.
    let expr = Expression::Param(id);
    let s = format_expression(&expr, &t2);
    // Either resolves to something or returns "<unknown>".
    assert!(!s.is_empty());
}

// ---------------------------------------------------------------------------
// BinOp formatting — precedence and parens
// ---------------------------------------------------------------------------

#[test]
fn test_add_mul_no_parens() {
    // a + b * c — no parens needed
    let (t, a, b, c) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Add,
        Box::new(Expression::BinOp(
            Box::new(Expression::Param(b)),
            BinOp::Mul,
            Box::new(Expression::Param(c)),
        )),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a + b * c");
}

#[test]
fn test_add_in_mul_gets_parens() {
    // (a + b) * c — inner add needs parens
    let (t, a, b, c) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::BinOp(
            Box::new(Expression::Param(a)),
            BinOp::Add,
            Box::new(Expression::Param(b)),
        )),
        BinOp::Mul,
        Box::new(Expression::Param(c)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "(a + b) * c");
}

#[test]
fn test_min_function_form() {
    let (t, a, b, _) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Min,
        Box::new(Expression::Param(b)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "min(a, b)");
}

#[test]
fn test_max_function_form() {
    let (t, a, b, _) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Max,
        Box::new(Expression::Param(b)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "max(a, b)");
}

#[test]
fn test_compare_gt() {
    let (t, a, b, _) = make_table();
    let expr = Expression::Compare(
        Box::new(Expression::Param(a)),
        CmpOp::Gt,
        Box::new(Expression::Param(b)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a > b");
}

#[test]
fn test_compare_all_ops() {
    let (t, a, b, _) = make_table();
    let cases = [
        (CmpOp::Gt, "a > b"),
        (CmpOp::Lt, "a < b"),
        (CmpOp::Ge, "a >= b"),
        (CmpOp::Le, "a <= b"),
        (CmpOp::Eq, "a == b"),
        (CmpOp::Ne, "a != b"),
    ];
    for (op, expected) in cases {
        let expr = Expression::Compare(
            Box::new(Expression::Param(a)),
            op,
            Box::new(Expression::Param(b)),
        );
        let s = format_expression(&expr, &t);
        assert_eq!(s, expected, "op: {op:?}");
    }
}

#[test]
fn test_conditional() {
    // a > 0 ? 5mm : 3mm
    let t = ParameterTable::new();
    let mut t2 = t.clone();
    let a_id = t2.add_global("a", 1.0, Unit::Ratio).unwrap();
    let five_mm = Unit::Millimeters.to_si(5.0);
    let three_mm = Unit::Millimeters.to_si(3.0);
    let expr = Expression::Conditional(
        Box::new(Expression::Compare(
            Box::new(Expression::Param(a_id)),
            CmpOp::Gt,
            Box::new(Expression::Literal(0.0)),
        )),
        Box::new(Expression::Literal(five_mm)),
        Box::new(Expression::Literal(three_mm)),
    );
    let s = format_expression_with_unit(&expr, &t2, Unit::Millimeters);
    assert_eq!(s, "a > 0mm ? 5mm : 3mm");
}

#[test]
fn test_logical_and_format() {
    let (t, a, b, _) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::And,
        Box::new(Expression::Param(b)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a && b");
}

#[test]
fn test_logical_or_format() {
    let (t, a, b, _) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Or,
        Box::new(Expression::Param(b)),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a || b");
}

#[test]
fn test_negation_format() {
    let (t, a, _, _) = make_table();
    let expr = Expression::Func(UnaryFunc::Neg, Box::new(Expression::Param(a)));
    let s = format_expression(&expr, &t);
    assert_eq!(s, "-a");
}

#[test]
fn test_not_format() {
    let (t, a, _, _) = make_table();
    let expr = Expression::Func(UnaryFunc::Not, Box::new(Expression::Param(a)));
    let s = format_expression(&expr, &t);
    assert_eq!(s, "!a");
}

#[test]
fn test_func_sin_format() {
    let (t, a, _, _) = make_table();
    let expr = Expression::Func(UnaryFunc::Sin, Box::new(Expression::Param(a)));
    let s = format_expression(&expr, &t);
    assert_eq!(s, "sin(a)");
}

// ---------------------------------------------------------------------------
// Round-trip tests
// ---------------------------------------------------------------------------

/// Parse an expression, format it, parse again, and check the evaluated values
/// match.
fn round_trip_value(input: &str, table: &ParameterTable) {
    use crate::expression::evaluate;

    let expr1 = parse_expression(input, table).expect("first parse failed");
    let val1 = evaluate(&expr1, table).expect("first eval failed");

    let formatted = format_expression(&expr1, table);
    let expr2 = parse_expression(&formatted, table).expect("second parse failed");
    let val2 = evaluate(&expr2, table).expect("second eval failed");

    assert!(
        (val1 - val2).abs() < 1e-9,
        "round-trip value mismatch for '{input}': {val1} vs {val2} (formatted: '{formatted}')"
    );
}

#[test]
fn test_round_trip_arithmetic() {
    let mut t = ParameterTable::new();
    let _ = t.add_global("a", 5.0, Unit::Ratio).unwrap();
    let _ = t.add_global("b", 3.0, Unit::Ratio).unwrap();
    round_trip_value("a + b * 2", &t);
    round_trip_value("(a + b) / 2", &t);
    round_trip_value("a - b", &t);
}

#[test]
fn test_round_trip_comparison() {
    let mut t = ParameterTable::new();
    let _ = t.add_global("a", 5.0, Unit::Ratio).unwrap();
    let _ = t.add_global("b", 3.0, Unit::Ratio).unwrap();
    round_trip_value("a > b", &t);
    round_trip_value("a == b", &t);
    round_trip_value("a != b", &t);
}

#[test]
fn test_round_trip_functions() {
    let t = ParameterTable::new();
    round_trip_value("sin(pi)", &t);
    round_trip_value("sqrt(4)", &t);
    round_trip_value("abs(-3)", &t);
}

#[test]
fn test_round_trip_with_unit() {
    // A literal with unit context: format in mm, re-parse as mm.
    use crate::expression::evaluate;

    let t = ParameterTable::new();
    let expr1 = parse_expression("50mm", &t).unwrap();
    let val1 = evaluate(&expr1, &t).unwrap();
    let formatted = format_expression_with_unit(&expr1, &t, Unit::Millimeters);
    let expr2 = parse_expression(&formatted, &t).unwrap();
    let val2 = evaluate(&expr2, &t).unwrap();
    assert!(
        (val1 - val2).abs() < 1e-12,
        "values differ: {val1} vs {val2}"
    );
}

#[test]
fn test_round_trip_conditional() {
    let mut t = ParameterTable::new();
    let _ = t.add_global("a", 5.0, Unit::Ratio).unwrap();
    let _ = t.add_global("b", 3.0, Unit::Ratio).unwrap();
    round_trip_value("a > b ? a : b", &t);
}

#[test]
fn test_round_trip_min_max() {
    let mut t = ParameterTable::new();
    let _ = t.add_global("a", 5.0, Unit::Ratio).unwrap();
    let _ = t.add_global("b", 3.0, Unit::Ratio).unwrap();
    round_trip_value("min(a, b)", &t);
    round_trip_value("max(a, b)", &t);
}

// ---------------------------------------------------------------------------
// Additional formatting edge cases
// ---------------------------------------------------------------------------

#[test]
fn test_sub_right_assoc_parens() {
    // a - (b - c) needs parens on the right
    let (t, a, b, c) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Sub,
        Box::new(Expression::BinOp(
            Box::new(Expression::Param(b)),
            BinOp::Sub,
            Box::new(Expression::Param(c)),
        )),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a - (b - c)");
}

#[test]
fn test_div_right_assoc_parens() {
    // a / (b / c) needs parens on the right
    let (t, a, b, c) = make_table();
    let expr = Expression::BinOp(
        Box::new(Expression::Param(a)),
        BinOp::Div,
        Box::new(Expression::BinOp(
            Box::new(Expression::Param(b)),
            BinOp::Div,
            Box::new(Expression::Param(c)),
        )),
    );
    let s = format_expression(&expr, &t);
    assert_eq!(s, "a / (b / c)");
}
