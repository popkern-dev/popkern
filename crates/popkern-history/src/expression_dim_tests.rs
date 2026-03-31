use super::*;
use crate::error::HistoryError;
use crate::parameter::ParameterTable;
use crate::unit::{DimensionKind, Unit};

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

// ── Dimensional type checking: Literal ───────────────────────────────────────

#[test]
fn dim_literal_is_dimensionless() {
    let table = ParameterTable::new();
    assert_eq!(
        dimension_kind(&lit(42.0), &table).unwrap(),
        DimensionKind::Dimensionless
    );
}

// ── Dimensional type checking: Param ─────────────────────────────────────────

#[test]
fn dim_param_length() {
    let mut table = ParameterTable::new();
    let id = table.add_global("d", 1.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&param(id), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_param_angle() {
    let mut table = ParameterTable::new();
    let id = table.add_global("a", 1.0, Unit::Radians).unwrap();
    assert_eq!(
        dimension_kind(&param(id), &table).unwrap(),
        DimensionKind::Angle
    );
}

// ── Dimensional type checking: Add/Sub ───────────────────────────────────────

#[test]
fn dim_add_length_length() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 2.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(param(a), BinOp::Add, param(b)), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_add_length_angle_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 1.0, Unit::Radians).unwrap();
    assert!(matches!(
        dimension_kind(&binop(param(a), BinOp::Add, param(b)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

// ── Dimensional type checking: Mul/Div ───────────────────────────────────────

#[test]
fn dim_mul_length_dimensionless() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(param(a), BinOp::Mul, lit(2.0)), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_mul_dimensionless_length() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(lit(2.0), BinOp::Mul, param(a)), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_mul_length_length_gives_area() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 2.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(param(a), BinOp::Mul, param(b)), &table).unwrap(),
        DimensionKind::Area
    );
}

#[test]
fn dim_div_length_dimensionless() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(param(a), BinOp::Div, lit(2.0)), &table).unwrap(),
        DimensionKind::Length
    );
}

// ── Dimensional type checking: Trig functions ────────────────────────────────

#[test]
fn dim_sin_angle_to_dimensionless() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Radians).unwrap();
    assert_eq!(
        dimension_kind(&func(UnaryFunc::Sin, param(a)), &table).unwrap(),
        DimensionKind::Dimensionless
    );
}

#[test]
fn dim_sin_length_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert!(matches!(
        dimension_kind(&func(UnaryFunc::Sin, param(a)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

#[test]
fn dim_asin_dimensionless_to_angle() {
    let table = ParameterTable::new();
    assert_eq!(
        dimension_kind(&func(UnaryFunc::Asin, lit(0.5)), &table).unwrap(),
        DimensionKind::Angle
    );
}

#[test]
fn dim_asin_length_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert!(matches!(
        dimension_kind(&func(UnaryFunc::Asin, param(a)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

// ── Dimensional type checking: Sqrt ──────────────────────────────────────────

#[test]
fn dim_sqrt_dimensionless() {
    let table = ParameterTable::new();
    assert_eq!(
        dimension_kind(&func(UnaryFunc::Sqrt, lit(9.0)), &table).unwrap(),
        DimensionKind::Dimensionless
    );
}

#[test]
fn dim_sqrt_length_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert!(matches!(
        dimension_kind(&func(UnaryFunc::Sqrt, param(a)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

// ── Dimensional type checking: Abs/Neg preserve kind ─────────────────────────

#[test]
fn dim_abs_preserves_length() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", -1.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&func(UnaryFunc::Abs, param(a)), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_neg_preserves_angle() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Radians).unwrap();
    assert_eq!(
        dimension_kind(&func(UnaryFunc::Neg, param(a)), &table).unwrap(),
        DimensionKind::Angle
    );
}

// ── Dimensional type checking: Compare ───────────────────────────────────────

#[test]
fn dim_compare_length_length_is_dimensionless() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 2.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&cmp(param(a), CmpOp::Gt, param(b)), &table).unwrap(),
        DimensionKind::Dimensionless
    );
}

// ── Dimensional type checking: Conditional ───────────────────────────────────

#[test]
fn dim_conditional_both_length() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 2.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&cond(lit(1.0), param(a), param(b)), &table).unwrap(),
        DimensionKind::Length
    );
}

#[test]
fn dim_conditional_mixed_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 1.0, Unit::Radians).unwrap();
    assert!(matches!(
        dimension_kind(&cond(lit(1.0), param(a), param(b)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

// ── Dimensional type checking: Min/Max ───────────────────────────────────────

#[test]
fn dim_min_length_length() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let b = table.add_global("b", 2.0, Unit::Meters).unwrap();
    assert_eq!(
        dimension_kind(&binop(param(a), BinOp::Min, param(b)), &table).unwrap(),
        DimensionKind::Length
    );
}

// ── Dimensional type checking: And/Or ────────────────────────────────────────

#[test]
fn dim_and_dimensionless() {
    let table = ParameterTable::new();
    assert_eq!(
        dimension_kind(&binop(lit(1.0), BinOp::And, lit(0.0)), &table).unwrap(),
        DimensionKind::Dimensionless
    );
}

#[test]
fn dim_and_with_length_error() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    assert!(matches!(
        dimension_kind(&binop(param(a), BinOp::And, lit(1.0)), &table),
        Err(HistoryError::DimensionMismatch { .. })
    ));
}

// ── Topological sort ─────────────────────────────────────────────────────────

#[test]
fn topo_sort_chain() {
    let mut table = ParameterTable::new();
    let c = table.add_global("c", 10.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 0.0, Unit::Ratio).unwrap();
    let a = table.add_global("a", 0.0, Unit::Ratio).unwrap();
    table
        .set_expression(b, binop(param(c), BinOp::Add, lit(1.0)))
        .unwrap();
    table
        .set_expression(a, binop(param(b), BinOp::Mul, lit(2.0)))
        .unwrap();

    let order = topological_sort(&table).unwrap();
    let b_pos = order.iter().position(|&id| id == b).unwrap();
    let a_pos = order.iter().position(|&id| id == a).unwrap();
    assert!(b_pos < a_pos, "B must be evaluated before A");
}

#[test]
fn topo_sort_no_deps() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    table.set_expression(a, lit(5.0)).unwrap();
    table.set_expression(b, lit(6.0)).unwrap();

    let order = topological_sort(&table).unwrap();
    assert_eq!(order.len(), 2);
    assert!(order.contains(&a));
    assert!(order.contains(&b));
}

#[test]
fn topo_sort_cycle_detected() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    table.set_expression(a, param(b)).unwrap();
    table.set_expression(b, param(a)).unwrap();

    assert!(
        matches!(
            topological_sort(&table),
            Err(HistoryError::ExpressionCycle(_))
        ),
        "expected ExpressionCycle"
    );
}

// ── evaluate_all ─────────────────────────────────────────────────────────────

#[test]
fn evaluate_all_propagates_values() {
    let mut table = ParameterTable::new();
    let c = table.add_global("c", 10.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 0.0, Unit::Ratio).unwrap();
    let a = table.add_global("a", 0.0, Unit::Ratio).unwrap();
    table
        .set_expression(b, binop(param(c), BinOp::Add, lit(1.0)))
        .unwrap();
    table
        .set_expression(a, binop(param(b), BinOp::Mul, lit(2.0)))
        .unwrap();

    let warnings = evaluate_all(&mut table);
    assert!(warnings.is_empty(), "unexpected warnings: {warnings:?}");
    assert!((table.get(b).unwrap().value - 11.0).abs() < 1e-12);
    assert!((table.get(a).unwrap().value - 22.0).abs() < 1e-12);
}

#[test]
fn evaluate_all_math_error_keeps_last_value() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 42.0, Unit::Ratio).unwrap();
    table
        .set_expression(a, func(UnaryFunc::Sqrt, lit(-1.0)))
        .unwrap();

    let warnings = evaluate_all(&mut table);
    assert_eq!(warnings.len(), 1);
    assert!(matches!(
        &warnings[0].1,
        EvaluationWarning::MathError { .. }
    ));
    assert!((table.get(a).unwrap().value - 42.0).abs() < 1e-12);
}

#[test]
fn evaluate_all_bounds_clamped() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 0.0, Unit::Ratio).unwrap();
    table.set_bounds(a, 0.0, 10.0);
    table.set_expression(a, lit(20.0)).unwrap();

    let warnings = evaluate_all(&mut table);
    assert_eq!(warnings.len(), 1);
    assert!(matches!(
        &warnings[0].1,
        EvaluationWarning::BoundsClamped {
            original,
            clamped,
        } if (*original - 20.0).abs() < 1e-12 && (*clamped - 10.0).abs() < 1e-12
    ));
    assert!((table.get(a).unwrap().value - 10.0).abs() < 1e-12);
}

// ── expression_param_refs ────────────────────────────────────────────────────

#[test]
fn refs_literal_empty() {
    assert!(expression_param_refs(&lit(5.0)).is_empty());
}

#[test]
fn refs_param_single() {
    let id = ParameterId::default();
    assert_eq!(expression_param_refs(&param(id)), vec![id]);
}

#[test]
fn refs_nested_binop() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    let refs = expression_param_refs(&binop(param(a), BinOp::Add, param(b)));
    assert!(refs.contains(&a));
    assert!(refs.contains(&b));
    assert_eq!(refs.len(), 2);
}

#[test]
fn refs_conditional_all_branches() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    let c = table.add_global("c", 3.0, Unit::Ratio).unwrap();
    let refs = expression_param_refs(&cond(param(a), param(b), param(c)));
    assert_eq!(refs.len(), 3);
    assert!(refs.contains(&a));
    assert!(refs.contains(&b));
    assert!(refs.contains(&c));
}

// ── set_expression / clear_expression / dependents ───────────────────────────

#[test]
fn set_expression_on_readonly_errors() {
    let mut table = ParameterTable::new();
    let pi_id = table.resolve_name("pi").unwrap();
    assert!(matches!(
        table.set_expression(pi_id, lit(3.0)),
        Err(HistoryError::ReadOnlyParameter)
    ));
}

#[test]
fn clear_expression_removes_it() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    table.set_expression(a, lit(5.0)).unwrap();
    assert!(table.get(a).unwrap().expression.is_some());
    table.clear_expression(a);
    assert!(table.get(a).unwrap().expression.is_none());
}

#[test]
fn dependents_finds_references() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    table.set_expression(b, param(a)).unwrap();
    let deps = table.dependents(a);
    assert_eq!(deps, vec![b]);
}

#[test]
fn dependents_empty_when_no_refs() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    assert!(table.dependents(a).is_empty());
}

// ── remove with dependents ───────────────────────────────────────────────────

#[test]
fn remove_parameter_in_use_errors() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    table.set_expression(b, param(a)).unwrap();
    assert!(
        matches!(table.remove(a), Err(HistoryError::ParameterInUse(_))),
        "expected ParameterInUse"
    );
}

#[test]
fn remove_cascade_clears_dependent_expressions() {
    let mut table = ParameterTable::new();
    let a = table.add_global("a", 1.0, Unit::Ratio).unwrap();
    let b = table.add_global("b", 2.0, Unit::Ratio).unwrap();
    table.set_expression(b, param(a)).unwrap();
    table.remove_cascade(a);
    assert!(table.get(a).is_none());
    assert!(table.get(b).unwrap().expression.is_none());
}
