//! Expression evaluation for parametric parameters.
//!
//! An [`Expression`] encodes a formula that can drive a parameter's value from
//! other parameters and mathematical functions. Expressions are evaluated
//! against a [`ParameterTable`] and support cycle detection, dimensional type
//! checking, and topological ordering.

use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::error::HistoryError;
use crate::parameter::{EvaluationWarning, ParameterId, ParameterTable};
use crate::unit::DimensionKind;

/// A formula that computes a parameter's value from literals, other parameters,
/// arithmetic, functions, comparisons, and conditionals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Expression {
    /// A literal constant value (in SI units).
    Literal(f64),
    /// Reference to another parameter by stable ID.
    Param(ParameterId),
    /// Binary operation: left op right.
    BinOp(Box<Expression>, BinOp, Box<Expression>),
    /// Unary function application.
    Func(UnaryFunc, Box<Expression>),
    /// Comparison: left cmp right -> 1.0 (true) or 0.0 (false).
    Compare(Box<Expression>, CmpOp, Box<Expression>),
    /// Conditional: if condition != 0 then true_branch else false_branch.
    Conditional(Box<Expression>, Box<Expression>, Box<Expression>),
    /// A geometric measurement query evaluated against the current BRep.
    Measure(crate::measure_query::MeasureQuery),
}

/// Binary arithmetic / logical operators.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BinOp {
    /// Addition.
    Add,
    /// Subtraction.
    Sub,
    /// Multiplication.
    Mul,
    /// Division (errors on zero divisor).
    Div,
    /// Logical AND (nonzero = true).
    And,
    /// Logical OR (nonzero = true).
    Or,
    /// Minimum of two values.
    Min,
    /// Maximum of two values.
    Max,
}

/// Comparison operators producing 0.0 or 1.0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CmpOp {
    /// Greater than.
    Gt,
    /// Less than.
    Lt,
    /// Greater than or equal.
    Ge,
    /// Less than or equal.
    Le,
    /// Equal (within 1e-12 tolerance).
    Eq,
    /// Not equal (outside 1e-12 tolerance).
    Ne,
}

/// Unary math functions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum UnaryFunc {
    /// Sine (input in radians).
    Sin,
    /// Cosine (input in radians).
    Cos,
    /// Tangent (input in radians).
    Tan,
    /// Arcsine (output in radians).
    Asin,
    /// Arccosine (output in radians).
    Acos,
    /// Arctangent (output in radians).
    Atan,
    /// Square root (errors on negative input).
    Sqrt,
    /// Absolute value.
    Abs,
    /// Negation.
    Neg,
    /// Logical NOT (0 -> 1, nonzero -> 0).
    Not,
    /// Floor (round toward negative infinity).
    Floor,
    /// Ceiling (round toward positive infinity).
    Ceil,
    /// Round to nearest integer.
    Round,
    /// Natural logarithm (errors on non-positive input).
    Ln,
}

/// Evaluate an expression against the parameter table.
///
/// Pass [`MeasureContext::empty()`](crate::measure_query::MeasureContext::empty)
/// when no live model is available.
pub fn evaluate(
    expr: &Expression,
    table: &ParameterTable,
    ctx: &crate::measure_query::MeasureContext<'_>,
) -> Result<f64, String> {
    match expr {
        Expression::Literal(v) => Ok(*v),
        Expression::Param(id) => table
            .get(*id)
            .map(|p| p.value)
            .ok_or_else(|| format!("parameter {id:?} not found")),
        Expression::BinOp(left, op, right) => {
            let l = evaluate(left, table, ctx)?;
            let r = evaluate(right, table, ctx)?;
            let result = match op {
                BinOp::Add => l + r,
                BinOp::Sub => l - r,
                BinOp::Mul => l * r,
                BinOp::Div => {
                    if r == 0.0 {
                        return Err("division by zero".into());
                    }
                    l / r
                }
                BinOp::And => {
                    if l != 0.0 && r != 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                BinOp::Or => {
                    if l != 0.0 || r != 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                BinOp::Min => l.min(r),
                BinOp::Max => l.max(r),
            };
            check_finite(result)
        }
        Expression::Func(func, arg) => {
            let a = evaluate(arg, table, ctx)?;
            let result = match func {
                UnaryFunc::Sin => a.sin(),
                UnaryFunc::Cos => a.cos(),
                UnaryFunc::Tan => a.tan(),
                UnaryFunc::Asin => a.asin(),
                UnaryFunc::Acos => a.acos(),
                UnaryFunc::Atan => a.atan(),
                UnaryFunc::Sqrt => {
                    if a < 0.0 {
                        return Err("sqrt of negative number".into());
                    }
                    a.sqrt()
                }
                UnaryFunc::Abs => a.abs(),
                UnaryFunc::Neg => -a,
                UnaryFunc::Not => {
                    if a == 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                }
                UnaryFunc::Floor => a.floor(),
                UnaryFunc::Ceil => a.ceil(),
                UnaryFunc::Round => a.round(),
                UnaryFunc::Ln => {
                    if a <= 0.0 {
                        return Err("ln of non-positive number".into());
                    }
                    a.ln()
                }
            };
            check_finite(result)
        }
        Expression::Compare(left, op, right) => {
            let l = evaluate(left, table, ctx)?;
            let r = evaluate(right, table, ctx)?;
            let result = match op {
                CmpOp::Gt => l > r,
                CmpOp::Lt => l < r,
                CmpOp::Ge => l >= r,
                CmpOp::Le => l <= r,
                CmpOp::Eq => (l - r).abs() < 1e-12,
                CmpOp::Ne => (l - r).abs() >= 1e-12,
            };
            Ok(if result { 1.0 } else { 0.0 })
        }
        Expression::Conditional(cond, then_expr, else_expr) => {
            let c = evaluate(cond, table, ctx)?;
            if c != 0.0 {
                evaluate(then_expr, table, ctx)
            } else {
                evaluate(else_expr, table, ctx)
            }
        }
        Expression::Measure(query) => crate::measure_resolve::resolve_and_evaluate(query, ctx)
            .map_err(|e| format!("measurement failed: {e}")),
    }
}

/// Check that a value is finite (not NaN or infinite).
fn check_finite(v: f64) -> Result<f64, String> {
    if v.is_nan() || v.is_infinite() {
        Err("result is NaN or infinite".into())
    } else {
        Ok(v)
    }
}

/// Determine the [`DimensionKind`] that an expression produces.
///
/// Enforces dimensional consistency: e.g. adding a length to an angle is an
/// error, and `sin()` requires an angle input.
pub fn dimension_kind(
    expr: &Expression,
    table: &ParameterTable,
) -> Result<DimensionKind, HistoryError> {
    match expr {
        Expression::Literal(_) => Ok(DimensionKind::Dimensionless),
        Expression::Param(id) => table
            .get(*id)
            .map(|p| p.unit.kind())
            .ok_or(HistoryError::ExpressionUnresolved(*id)),
        Expression::BinOp(left, op, right) => {
            let lk = dimension_kind(left, table)?;
            let rk = dimension_kind(right, table)?;
            match op {
                BinOp::Add | BinOp::Sub => {
                    if lk == rk {
                        Ok(lk)
                    } else {
                        Err(dim_mismatch(lk, rk))
                    }
                }
                BinOp::Mul | BinOp::Div => mul_div_kind(lk, rk),
                BinOp::Min | BinOp::Max => {
                    if lk == rk {
                        Ok(lk)
                    } else {
                        Err(dim_mismatch(lk, rk))
                    }
                }
                BinOp::And | BinOp::Or => {
                    if lk == DimensionKind::Dimensionless && rk == DimensionKind::Dimensionless {
                        Ok(DimensionKind::Dimensionless)
                    } else {
                        Err(dim_mismatch(lk, rk))
                    }
                }
            }
        }
        Expression::Func(func, arg) => {
            let ak = dimension_kind(arg, table)?;
            match func {
                UnaryFunc::Sin | UnaryFunc::Cos | UnaryFunc::Tan => {
                    if ak == DimensionKind::Angle {
                        Ok(DimensionKind::Dimensionless)
                    } else {
                        Err(dim_mismatch(DimensionKind::Angle, ak))
                    }
                }
                UnaryFunc::Asin | UnaryFunc::Acos | UnaryFunc::Atan => {
                    if ak == DimensionKind::Dimensionless {
                        Ok(DimensionKind::Angle)
                    } else {
                        Err(dim_mismatch(DimensionKind::Dimensionless, ak))
                    }
                }
                UnaryFunc::Sqrt | UnaryFunc::Ln | UnaryFunc::Not => {
                    if ak == DimensionKind::Dimensionless {
                        Ok(DimensionKind::Dimensionless)
                    } else {
                        Err(dim_mismatch(DimensionKind::Dimensionless, ak))
                    }
                }
                UnaryFunc::Abs
                | UnaryFunc::Neg
                | UnaryFunc::Floor
                | UnaryFunc::Ceil
                | UnaryFunc::Round => Ok(ak),
            }
        }
        Expression::Compare(left, _op, right) => {
            let lk = dimension_kind(left, table)?;
            let rk = dimension_kind(right, table)?;
            if lk == rk {
                Ok(DimensionKind::Dimensionless)
            } else {
                Err(dim_mismatch(lk, rk))
            }
        }
        Expression::Conditional(_cond, then_expr, else_expr) => {
            let tk = dimension_kind(then_expr, table)?;
            let ek = dimension_kind(else_expr, table)?;
            if tk == ek {
                Ok(tk)
            } else {
                Err(dim_mismatch(tk, ek))
            }
        }
        Expression::Measure(query) => Ok(query.dimension_kind()),
    }
}

/// Build a [`HistoryError::DimensionMismatch`] with a dummy parameter ID.
fn dim_mismatch(expected: DimensionKind, got: DimensionKind) -> HistoryError {
    // Use a default ParameterId — the caller knows which parameter the
    // expression belongs to.
    HistoryError::DimensionMismatch {
        param: ParameterId::default(),
        expected,
        got,
    }
}

/// Dimensional rules for Mul/Div.
fn mul_div_kind(lk: DimensionKind, rk: DimensionKind) -> Result<DimensionKind, HistoryError> {
    use DimensionKind::*;
    match (lk, rk) {
        (Dimensionless, Dimensionless) => Ok(Dimensionless),
        (Length, Dimensionless) | (Dimensionless, Length) => Ok(Length),
        (Angle, Dimensionless) | (Dimensionless, Angle) => Ok(Angle),
        (Area, Dimensionless) | (Dimensionless, Area) => Ok(Area),
        (Volume, Dimensionless) | (Dimensionless, Volume) => Ok(Volume),
        (InverseLength, Dimensionless) | (Dimensionless, InverseLength) => Ok(InverseLength),
        (MomentOfArea, Dimensionless) | (Dimensionless, MomentOfArea) => Ok(MomentOfArea),
        (Length, Length) => Ok(Area),
        (Area, Length) | (Length, Area) => Ok(Volume),
        _ => Err(dim_mismatch(lk, rk)),
    }
}

/// Collect all [`ParameterId`] references in an expression.
pub(crate) fn expression_param_refs(expr: &Expression) -> Vec<ParameterId> {
    match expr {
        Expression::Literal(_) | Expression::Measure(_) => vec![],
        Expression::Param(id) => vec![*id],
        Expression::BinOp(l, _, r) | Expression::Compare(l, _, r) => {
            let mut refs = expression_param_refs(l);
            refs.extend(expression_param_refs(r));
            refs
        }
        Expression::Func(_, arg) => expression_param_refs(arg),
        Expression::Conditional(c, t, f) => {
            let mut refs = expression_param_refs(c);
            refs.extend(expression_param_refs(t));
            refs.extend(expression_param_refs(f));
            refs
        }
    }
}

/// Topological sort of all parameters with expressions.
///
/// Returns evaluation order (dependencies before dependents). Detects cycles
/// and returns [`HistoryError::ExpressionCycle`] with the participating IDs.
pub fn topological_sort(table: &ParameterTable) -> Result<Vec<ParameterId>, HistoryError> {
    // Build adjacency list: param → [params it depends on].
    let mut deps: HashMap<ParameterId, Vec<ParameterId>> = HashMap::new();
    for (id, param) in table.all_params() {
        if let Some(ref expr) = param.expression {
            deps.insert(id, expression_param_refs(expr));
        }
    }

    let mut visited = HashSet::new();
    let mut in_progress = HashSet::new();
    let mut order = Vec::new();

    for &id in deps.keys() {
        if !visited.contains(&id) {
            topo_visit(id, &deps, &mut visited, &mut in_progress, &mut order)?;
        }
    }

    Ok(order)
}

/// Recursive DFS visitor for topological sort.
fn topo_visit(
    id: ParameterId,
    deps: &HashMap<ParameterId, Vec<ParameterId>>,
    visited: &mut HashSet<ParameterId>,
    in_progress: &mut HashSet<ParameterId>,
    order: &mut Vec<ParameterId>,
) -> Result<(), HistoryError> {
    if visited.contains(&id) {
        return Ok(());
    }
    if in_progress.contains(&id) {
        return Err(HistoryError::ExpressionCycle(vec![id]));
    }
    in_progress.insert(id);

    if let Some(dep_ids) = deps.get(&id) {
        for &dep in dep_ids {
            // Only visit deps that themselves have expressions (are in the
            // deps map). Plain-value params don't need ordering.
            if deps.contains_key(&dep) {
                topo_visit(dep, deps, visited, in_progress, order)?;
            }
        }
    }

    in_progress.remove(&id);
    visited.insert(id);
    order.push(id);
    Ok(())
}

/// Evaluate all parameter expressions in topological order.
///
/// The [`MeasureContext`](crate::measure_query::MeasureContext) provides BRep
/// access for `Measure` expressions. Pass
/// [`MeasureContext::empty()`](crate::measure_query::MeasureContext::empty) when
/// no live model is available.
///
/// On math error the last good value is kept and an [`EvaluationWarning`] is
/// recorded. On bounds violation the value is clamped and a warning recorded.
/// Returns all warnings produced during this pass.
pub fn evaluate_all(
    table: &mut ParameterTable,
    ctx: &crate::measure_query::MeasureContext<'_>,
) -> Vec<(ParameterId, EvaluationWarning)> {
    let order = match topological_sort(table) {
        Ok(order) => order,
        Err(_) => return vec![], // Cycle — caller handles via topological_sort
    };

    let mut warnings = Vec::new();

    for id in order {
        let expr = match table.get(id).and_then(|p| p.expression.clone()) {
            Some(e) => e,
            None => continue,
        };

        match evaluate(&expr, table, ctx) {
            Ok(mut value) => {
                // Check bounds.
                if let Some((min, max)) = table.get(id).and_then(|p| p.bounds)
                    && (value < min || value > max)
                {
                    let clamped = value.clamp(min, max);
                    warnings.push((
                        id,
                        EvaluationWarning::BoundsClamped {
                            original: value,
                            clamped,
                        },
                    ));
                    value = clamped;
                }
                // Write resolved value.
                if let Some(param) = table.get_mut(id) {
                    param.value = value;
                    param.eval_warning = None;
                }
            }
            Err(msg) => {
                // Keep last good value, record warning.
                let warning = EvaluationWarning::MathError { message: msg };
                warnings.push((id, warning.clone()));
                if let Some(param) = table.get_mut(id) {
                    param.eval_warning = Some(warning);
                }
            }
        }
    }

    warnings
}

#[cfg(test)]
#[path = "expression_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "expression_dim_tests.rs"]
mod dim_tests;
