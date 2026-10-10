use std::f64::consts::PI;

use crate::expression::{BinOp, CmpOp, Expression, UnaryFunc};
use crate::parameter::ParameterTable;
use crate::unit::Unit;

use super::parse_expression;

fn table_with_params() -> (
    ParameterTable,
    crate::parameter::ParameterId,
    crate::parameter::ParameterId,
    crate::parameter::ParameterId,
    crate::parameter::ParameterId,
    crate::parameter::ParameterId,
) {
    let mut t = ParameterTable::new();
    let width = t.add_global("width", 0.01, Unit::Millimeters).unwrap();
    let a = t.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = t.add_global("b", 2.0, Unit::Ratio).unwrap();
    let cond = t.add_global("cond", 1.0, Unit::Ratio).unwrap();
    // Feature-scoped param via add_feature_param (internal) — use add_global with dotted name
    // Since add_global stores the whole name as key, "extrude1.distance" is the qualified name.
    // We use the internal add_feature_param via the node module. Since we can't easily call it,
    // add a global named "extrude1.distance" via the name_index only. Instead, let's just check
    // resolution: add_global doesn't support dots, so we'll use a simpler test approach.
    let _ = cond;
    let _ = b;
    let _ = a;
    let _ = width;
    // Return fresh table with all params.
    let mut t2 = ParameterTable::new();
    let width2 = t2.add_global("width", 0.01, Unit::Millimeters).unwrap();
    let a2 = t2.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b2 = t2.add_global("b", 2.0, Unit::Ratio).unwrap();
    let cond2 = t2.add_global("cond", 1.0, Unit::Ratio).unwrap();
    let pi2 = t2.resolve_name("pi").unwrap();
    (t2, width2, a2, b2, cond2, pi2)
}

#[test]
fn test_literal_bare() {
    let t = ParameterTable::new();
    let expr = parse_expression("50", &t).unwrap();
    assert!(matches!(expr, Expression::Literal(v) if (v - 50.0).abs() < 1e-12));
}

#[test]
fn test_literal_mm() {
    let t = ParameterTable::new();
    let expr = parse_expression("50mm", &t).unwrap();
    let expected = Unit::Millimeters.to_si(50.0);
    assert!(matches!(expr, Expression::Quantity { value: v, .. } if (v - expected).abs() < 1e-12));
}

#[test]
fn test_literal_deg() {
    let t = ParameterTable::new();
    let expr = parse_expression("45deg", &t).unwrap();
    let expected = Unit::Degrees.to_si(45.0);
    assert!(matches!(expr, Expression::Quantity { value: v, .. } if (v - expected).abs() < 1e-12));
}

#[test]
fn test_literal_in() {
    let t = ParameterTable::new();
    let expr = parse_expression("3.5in", &t).unwrap();
    let expected = Unit::Inches.to_si(3.5);
    assert!(matches!(expr, Expression::Quantity { value: v, .. } if (v - expected).abs() < 1e-12));
}

#[test]
fn test_literal_cm() {
    let t = ParameterTable::new();
    let expr = parse_expression("2cm", &t).unwrap();
    let expected = Unit::Centimeters.to_si(2.0);
    assert!(matches!(expr, Expression::Quantity { value: v, .. } if (v - expected).abs() < 1e-12));
}

#[test]
fn test_literal_rad() {
    let t = ParameterTable::new();
    let expr = parse_expression("1rad", &t).unwrap();
    assert!(
        matches!(expr, Expression::Quantity { value, unit: Unit::Radians } if (value - 1.0).abs() < 1e-12)
    );
}

#[test]
fn test_literal_m() {
    let t = ParameterTable::new();
    let expr = parse_expression("1m", &t).unwrap();
    assert!(
        matches!(expr, Expression::Quantity { value, unit: Unit::Meters } if (value - 1.0).abs() < 1e-12)
    );
}

#[test]
fn test_param_simple() {
    let (t, width_id, _, _, _, _) = table_with_params();
    let expr = parse_expression("width", &t).unwrap();
    assert!(matches!(expr, Expression::Param(id) if id == width_id));
}

#[test]
fn test_param_builtin_pi() {
    let t = ParameterTable::new();
    let pi_id = t.resolve_name("pi").unwrap();
    let expr = parse_expression("pi", &t).unwrap();
    assert!(matches!(expr, Expression::Param(id) if id == pi_id));
}

#[test]
fn test_dotted_param() {
    // We can test dotted name resolution by adding a param via add_global (which stores the
    // full string as key). add_global supports "extrude1.distance" as name.
    let mut t = ParameterTable::new();
    let dist_id = t
        .add_global("extrude1.distance", 0.01, Unit::Millimeters)
        .unwrap();
    let expr = parse_expression("extrude1.distance", &t).unwrap();
    assert!(matches!(expr, Expression::Param(id) if id == dist_id));
}

#[test]
fn test_arithmetic_add_mul() {
    // "width * 2 + 3mm" => BinOp(BinOp(Param(width), Mul, Literal(2)), Add, Literal(0.003))
    let (t, width_id, _, _, _, _) = table_with_params();
    let expr = parse_expression("width * 2 + 3mm", &t).unwrap();
    let expected_mm = Unit::Millimeters.to_si(3.0);
    match expr {
        Expression::BinOp(left, BinOp::Add, right) => {
            match *left {
                Expression::BinOp(ll, BinOp::Mul, rr) => {
                    assert!(matches!(*ll, Expression::Param(id) if id == width_id));
                    assert!(matches!(*rr, Expression::Literal(v) if (v - 2.0).abs() < 1e-12));
                }
                _ => panic!("expected inner BinOp(Mul)"),
            }
            assert!(
                matches!(*right, Expression::Quantity { value: v, .. } if (v - expected_mm).abs() < 1e-12)
            );
        }
        _ => panic!("expected BinOp(Add)"),
    }
}

#[test]
fn test_sin_pi_over_4() {
    // "sin(pi / 4)" => Func(Sin, BinOp(Param(pi), Div, Literal(4)))
    let t = ParameterTable::new();
    let pi_id = t.resolve_name("pi").unwrap();
    let expr = parse_expression("sin(pi / 4)", &t).unwrap();
    match expr {
        Expression::Func(UnaryFunc::Sin, arg) => match *arg {
            Expression::BinOp(left, BinOp::Div, right) => {
                assert!(matches!(*left, Expression::Param(id) if id == pi_id));
                assert!(matches!(*right, Expression::Literal(v) if (v - 4.0).abs() < 1e-12));
            }
            _ => panic!("expected BinOp(Div)"),
        },
        _ => panic!("expected Func(Sin)"),
    }
}

#[test]
fn test_min() {
    let (t, _, a_id, b_id, _, _) = table_with_params();
    let expr = parse_expression("min(a, b)", &t).unwrap();
    match expr {
        Expression::BinOp(left, BinOp::Min, right) => {
            assert!(matches!(*left, Expression::Param(id) if id == a_id));
            assert!(matches!(*right, Expression::Param(id) if id == b_id));
        }
        _ => panic!("expected BinOp(Min)"),
    }
}

#[test]
fn test_max() {
    let (t, _, a_id, b_id, _, _) = table_with_params();
    let expr = parse_expression("max(a, b)", &t).unwrap();
    match expr {
        Expression::BinOp(left, BinOp::Max, right) => {
            assert!(matches!(*left, Expression::Param(id) if id == a_id));
            assert!(matches!(*right, Expression::Param(id) if id == b_id));
        }
        _ => panic!("expected BinOp(Max)"),
    }
}

#[test]
fn test_ternary() {
    // "a > 10mm ? 5mm : 3mm" => Conditional(Compare, Literal, Literal)
    let (t, _, a_id, _, _, _) = table_with_params();
    let expr = parse_expression("a > 10mm ? 5mm : 3mm", &t).unwrap();
    match expr {
        Expression::Conditional(cond, then_e, else_e) => {
            match *cond {
                Expression::Compare(left, CmpOp::Gt, right) => {
                    assert!(matches!(*left, Expression::Param(id) if id == a_id));
                    let expected = Unit::Millimeters.to_si(10.0);
                    assert!(
                        matches!(*right, Expression::Quantity { value: v, .. } if (v - expected).abs() < 1e-12)
                    );
                }
                _ => panic!("expected Compare"),
            }
            let expected_5 = Unit::Millimeters.to_si(5.0);
            let expected_3 = Unit::Millimeters.to_si(3.0);
            assert!(
                matches!(*then_e, Expression::Quantity { value: v, .. } if (v - expected_5).abs() < 1e-12)
            );
            assert!(
                matches!(*else_e, Expression::Quantity { value: v, .. } if (v - expected_3).abs() < 1e-12)
            );
        }
        _ => panic!("expected Conditional"),
    }
}

#[test]
fn test_logical_and() {
    // "a > 0 && b < 10" => BinOp(Compare(a,Gt,0), And, Compare(b,Lt,10))
    let (t, _, a_id, b_id, _, _) = table_with_params();
    let expr = parse_expression("a > 0 && b < 10", &t).unwrap();
    match expr {
        Expression::BinOp(left, BinOp::And, right) => {
            assert!(
                matches!(*left, Expression::Compare(ref la, CmpOp::Gt, _) if matches!(**la, Expression::Param(id) if id == a_id))
            );
            assert!(
                matches!(*right, Expression::Compare(ref lb, CmpOp::Lt, _) if matches!(**lb, Expression::Param(id) if id == b_id))
            );
        }
        _ => panic!("expected BinOp(And)"),
    }
}

#[test]
fn test_paren_grouping() {
    // "(a + b) / 2" => BinOp(BinOp(a, Add, b), Div, Literal(2))
    let (t, _, a_id, b_id, _, _) = table_with_params();
    let expr = parse_expression("(a + b) / 2", &t).unwrap();
    match expr {
        Expression::BinOp(left, BinOp::Div, right) => {
            match *left {
                Expression::BinOp(ll, BinOp::Add, rr) => {
                    assert!(matches!(*ll, Expression::Param(id) if id == a_id));
                    assert!(matches!(*rr, Expression::Param(id) if id == b_id));
                }
                _ => panic!("expected BinOp(Add)"),
            }
            assert!(matches!(*right, Expression::Literal(v) if (v - 2.0).abs() < 1e-12));
        }
        _ => panic!("expected BinOp(Div)"),
    }
}

#[test]
fn test_unary_neg() {
    let (t, _, a_id, _, _, _) = table_with_params();
    let expr = parse_expression("-a", &t).unwrap();
    match expr {
        Expression::Func(UnaryFunc::Neg, arg) => {
            assert!(matches!(*arg, Expression::Param(id) if id == a_id));
        }
        _ => panic!("expected Func(Neg)"),
    }
}

#[test]
fn test_unary_not() {
    let (t, _, _, _, cond_id, _) = table_with_params();
    let expr = parse_expression("!cond", &t).unwrap();
    match expr {
        Expression::Func(UnaryFunc::Not, arg) => {
            assert!(matches!(*arg, Expression::Param(id) if id == cond_id));
        }
        _ => panic!("expected Func(Not)"),
    }
}

#[test]
fn test_unknown_param_error() {
    let t = ParameterTable::new();
    let result = parse_expression("unknown_param", &t);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(err.message.contains("unknown_param"));
}

#[test]
fn test_empty_input_error() {
    let t = ParameterTable::new();
    let result = parse_expression("", &t);
    assert!(result.is_err());
}

#[test]
fn test_incomplete_expression_error() {
    let t = ParameterTable::new();
    let result = parse_expression("3 +", &t);
    assert!(result.is_err());
    let err = result.unwrap_err();
    // Column should be at end of input.
    assert!(err.column > 0);
}

#[test]
fn test_double_plus_error() {
    let t = ParameterTable::new();
    let result = parse_expression("3 + + 3", &t);
    assert!(result.is_err());
}

#[test]
fn test_min_single_arg_error() {
    let (t, _, a_id, _, _, _) = table_with_params();
    let _ = a_id;
    let result = parse_expression("min(a)", &t);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(err.message.contains("2 argument"));
}

#[test]
fn test_pi_value() {
    let t = ParameterTable::new();
    let expr = parse_expression("pi", &t).unwrap();
    if let Expression::Param(id) = expr {
        let val = t.get(id).unwrap().value;
        assert!((val - PI).abs() < 1e-12);
    } else {
        panic!("expected Param");
    }
}

#[test]
fn test_whitespace_handling() {
    let t = ParameterTable::new();
    let expr = parse_expression("  42  ", &t).unwrap();
    assert!(matches!(expr, Expression::Literal(v) if (v - 42.0).abs() < 1e-12));
}

#[test]
fn test_logical_or() {
    let (t, _, a_id, b_id, _, _) = table_with_params();
    let expr = parse_expression("a > 0 || b < 10", &t).unwrap();
    match expr {
        Expression::BinOp(_, BinOp::Or, _) => {}
        _ => panic!("expected BinOp(Or)"),
    }
    let _ = (a_id, b_id);
}

#[test]
fn test_comparison_operators() {
    let (t, _, a_id, _, _, _) = table_with_params();
    let ops = [
        ("a >= 0", CmpOp::Ge),
        ("a <= 0", CmpOp::Le),
        ("a == 0", CmpOp::Eq),
        ("a != 0", CmpOp::Ne),
    ];
    for (input, expected_op) in ops {
        let expr = parse_expression(input, &t).unwrap();
        match expr {
            Expression::Compare(left, op, _) => {
                assert_eq!(op, expected_op, "input: {input}");
                assert!(matches!(*left, Expression::Param(id) if id == a_id));
            }
            _ => panic!("expected Compare for: {input}"),
        }
    }
}

#[test]
fn test_nested_unary_neg() {
    // --a should parse as Neg(Neg(a))
    let (t, _, a_id, _, _, _) = table_with_params();
    let expr = parse_expression("--a", &t).unwrap();
    match expr {
        Expression::Func(UnaryFunc::Neg, inner) => match *inner {
            Expression::Func(UnaryFunc::Neg, innermost) => {
                assert!(matches!(*innermost, Expression::Param(id) if id == a_id));
            }
            _ => panic!("expected inner Neg"),
        },
        _ => panic!("expected Func(Neg)"),
    }
}

#[test]
fn test_trailing_token_error() {
    let t = ParameterTable::new();
    // "3 3" has a trailing token after the first number
    let result = parse_expression("3 3", &t);
    assert!(result.is_err());
}
