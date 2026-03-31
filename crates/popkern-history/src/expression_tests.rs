use super::*;
use crate::measure_query::MeasureContext;
use crate::parameter::ParameterTable;
use crate::unit::Unit;

// ── Helpers ──────────────────────────────────────────────────────────────────

fn lit(v: f64) -> Expression {
    Expression::Literal(v)
}

fn param(id: ParameterId) -> Expression {
    Expression::Param(id)
}

fn binop(l: Expression, op: BinOp, r: Expression) -> Expression {
    Expression::BinOp(Box::new(l), op, Box::new(r))
}

fn func(f: UnaryFunc, arg: Expression) -> Expression {
    Expression::Func(f, Box::new(arg))
}

fn cmp(l: Expression, op: CmpOp, r: Expression) -> Expression {
    Expression::Compare(Box::new(l), op, Box::new(r))
}

fn cond(c: Expression, t: Expression, e: Expression) -> Expression {
    Expression::Conditional(Box::new(c), Box::new(t), Box::new(e))
}

/// Shorthand: evaluate with an empty MeasureContext.
fn eval(expr: &Expression, table: &ParameterTable) -> Result<f64, String> {
    evaluate(expr, table, &MeasureContext::empty())
}

// ── Evaluation: Literal ──────────────────────────────────────────────────────

#[test]
fn eval_literal() {
    let table = ParameterTable::new();
    assert_eq!(eval(&lit(42.0), &table).unwrap(), 42.0);
}

// ── Evaluation: Param reference ──────────────────────────────────────────────

#[test]
fn eval_param_reference() {
    let mut table = ParameterTable::new();
    let id = table.add_global("x", 7.5, Unit::Meters).unwrap();
    assert_eq!(eval(&param(id), &table).unwrap(), 7.5);
}

#[test]
fn eval_param_not_found() {
    let table = ParameterTable::new();
    let bad_id = ParameterId::default();
    assert!(eval(&param(bad_id), &table).is_err());
}

// ── Evaluation: BinOp ────────────────────────────────────────────────────────

#[test]
fn eval_add() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(3.0), BinOp::Add, lit(4.0)), &table).unwrap(),
        7.0
    );
}

#[test]
fn eval_sub() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(10.0), BinOp::Sub, lit(3.0)), &table).unwrap(),
        7.0
    );
}

#[test]
fn eval_mul() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(3.0), BinOp::Mul, lit(4.0)), &table).unwrap(),
        12.0
    );
}

#[test]
fn eval_div() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(10.0), BinOp::Div, lit(4.0)), &table).unwrap(),
        2.5
    );
}

#[test]
fn eval_div_by_zero() {
    let table = ParameterTable::new();
    let err = eval(&binop(lit(1.0), BinOp::Div, lit(0.0)), &table).unwrap_err();
    assert!(err.contains("division by zero"), "got: {err}");
}

#[test]
fn eval_and_both_true() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(1.0), BinOp::And, lit(2.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_and_one_false() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(1.0), BinOp::And, lit(0.0)), &table).unwrap(),
        0.0
    );
}

#[test]
fn eval_or_both_false() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(0.0), BinOp::Or, lit(0.0)), &table).unwrap(),
        0.0
    );
}

#[test]
fn eval_or_one_true() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(0.0), BinOp::Or, lit(5.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_min() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(3.0), BinOp::Min, lit(7.0)), &table).unwrap(),
        3.0
    );
}

#[test]
fn eval_max() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&binop(lit(3.0), BinOp::Max, lit(7.0)), &table).unwrap(),
        7.0
    );
}

// ── Evaluation: UnaryFunc ────────────────────────────────────────────────────

#[test]
fn eval_sin() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Sin, lit(0.0)), &table).unwrap();
    assert!((v - 0.0).abs() < 1e-12);
}

#[test]
fn eval_cos() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Cos, lit(0.0)), &table).unwrap();
    assert!((v - 1.0).abs() < 1e-12);
}

#[test]
fn eval_tan() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Tan, lit(0.0)), &table).unwrap();
    assert!((v - 0.0).abs() < 1e-12);
}

#[test]
fn eval_asin() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Asin, lit(1.0)), &table).unwrap();
    assert!((v - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
}

#[test]
fn eval_acos() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Acos, lit(1.0)), &table).unwrap();
    assert!((v - 0.0).abs() < 1e-12);
}

#[test]
fn eval_atan() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Atan, lit(1.0)), &table).unwrap();
    assert!((v - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
}

#[test]
fn eval_sqrt() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Sqrt, lit(9.0)), &table).unwrap();
    assert!((v - 3.0).abs() < 1e-12);
}

#[test]
fn eval_sqrt_negative() {
    let table = ParameterTable::new();
    let err = eval(&func(UnaryFunc::Sqrt, lit(-1.0)), &table).unwrap_err();
    assert!(err.contains("sqrt of negative"), "got: {err}");
}

#[test]
fn eval_abs() {
    let table = ParameterTable::new();
    assert_eq!(eval(&func(UnaryFunc::Abs, lit(-5.0)), &table).unwrap(), 5.0);
}

#[test]
fn eval_neg() {
    let table = ParameterTable::new();
    assert_eq!(eval(&func(UnaryFunc::Neg, lit(3.0)), &table).unwrap(), -3.0);
}

#[test]
fn eval_not_zero() {
    let table = ParameterTable::new();
    assert_eq!(eval(&func(UnaryFunc::Not, lit(0.0)), &table).unwrap(), 1.0);
}

#[test]
fn eval_not_nonzero() {
    let table = ParameterTable::new();
    assert_eq!(eval(&func(UnaryFunc::Not, lit(5.0)), &table).unwrap(), 0.0);
}

#[test]
fn eval_floor() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&func(UnaryFunc::Floor, lit(3.7)), &table).unwrap(),
        3.0
    );
}

#[test]
fn eval_ceil() {
    let table = ParameterTable::new();
    assert_eq!(eval(&func(UnaryFunc::Ceil, lit(3.2)), &table).unwrap(), 4.0);
}

#[test]
fn eval_round() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&func(UnaryFunc::Round, lit(3.5)), &table).unwrap(),
        4.0
    );
}

#[test]
fn eval_ln() {
    let table = ParameterTable::new();
    let v = eval(&func(UnaryFunc::Ln, lit(std::f64::consts::E)), &table).unwrap();
    assert!((v - 1.0).abs() < 1e-12);
}

#[test]
fn eval_ln_zero() {
    let table = ParameterTable::new();
    let err = eval(&func(UnaryFunc::Ln, lit(0.0)), &table).unwrap_err();
    assert!(err.contains("ln of non-positive"), "got: {err}");
}

#[test]
fn eval_ln_negative() {
    let table = ParameterTable::new();
    let err = eval(&func(UnaryFunc::Ln, lit(-1.0)), &table).unwrap_err();
    assert!(err.contains("ln of non-positive"), "got: {err}");
}

// ── Evaluation: CmpOp ───────────────────────────────────────────────────────

#[test]
fn eval_cmp_gt_true() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(5.0), CmpOp::Gt, lit(3.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_cmp_gt_false() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(3.0), CmpOp::Gt, lit(5.0)), &table).unwrap(),
        0.0
    );
}

#[test]
fn eval_cmp_lt() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(3.0), CmpOp::Lt, lit(5.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_cmp_ge() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(5.0), CmpOp::Ge, lit(5.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_cmp_le() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(5.0), CmpOp::Le, lit(5.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_cmp_eq() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(3.0), CmpOp::Eq, lit(3.0)), &table).unwrap(),
        1.0
    );
}

#[test]
fn eval_cmp_ne() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cmp(lit(3.0), CmpOp::Ne, lit(5.0)), &table).unwrap(),
        1.0
    );
}

// ── Evaluation: Conditional ──────────────────────────────────────────────────

#[test]
fn eval_conditional_true_branch() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cond(lit(1.0), lit(10.0), lit(20.0)), &table).unwrap(),
        10.0
    );
}

#[test]
fn eval_conditional_false_branch() {
    let table = ParameterTable::new();
    assert_eq!(
        eval(&cond(lit(0.0), lit(10.0), lit(20.0)), &table).unwrap(),
        20.0
    );
}

// ── Evaluation: Nested expression ────────────────────────────────────────────

#[test]
fn eval_nested_three_levels() {
    let table = ParameterTable::new();
    // (2 + 3) * sqrt(16) = 5 * 4 = 20
    let inner_add = binop(lit(2.0), BinOp::Add, lit(3.0));
    let sqrt_16 = func(UnaryFunc::Sqrt, lit(16.0));
    let expr = binop(inner_add, BinOp::Mul, sqrt_16);
    assert_eq!(eval(&expr, &table).unwrap(), 20.0);
}

// ── Evaluation: NaN/Inf detection ────────────────────────────────────────────

#[test]
fn eval_overflow_to_inf_detected() {
    let table = ParameterTable::new();
    // 1e308 + 1e308 overflows to +inf.
    let expr = binop(lit(1e308), BinOp::Add, lit(1e308));
    let err = eval(&expr, &table).unwrap_err();
    assert!(err.contains("NaN or infinite"), "got: {err}");
}
