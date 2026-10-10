//! Expression evaluation for parametric parameters.
//!
//! An [`Expression`] encodes a formula that can drive a parameter's value from
//! other parameters and mathematical functions. Expressions are evaluated
//! against a [`ParameterTable`] and support cycle detection, dimensional type
//! checking, and topological ordering.

use popkern_math::PortableMath;
use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::error::HistoryError;
use crate::parameter::{EvaluationWarning, ParameterId, ParameterTable};

/// A formula that computes a parameter's value from literals, other parameters,
/// arithmetic, functions, comparisons, and conditionals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Expression {
    /// A literal constant value (in SI units).
    Literal(f64),
    /// A unit-suffixed literal (`180mm`, `10deg`): its SI value and the
    /// unit it was written in, which carries its dimension and its
    /// spelling (PK-754).
    Quantity {
        /// The value in SI units.
        value: f64,
        /// The unit the literal was written in.
        unit: crate::unit::Unit,
    },
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

/// Evaluate an expression against the parameter table and optional model context.
pub fn evaluate(
    expr: &Expression,
    table: &ParameterTable,
    ctx: &crate::measure_query::MeasureContext<'_>,
) -> Result<f64, String> {
    match expr {
        Expression::Literal(v) | Expression::Quantity { value: v, .. } => Ok(*v),
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
                UnaryFunc::Sin => a.psin(),
                UnaryFunc::Cos => a.pcos(),
                UnaryFunc::Tan => a.ptan(),
                UnaryFunc::Asin => a.pasin(),
                UnaryFunc::Acos => a.pacos(),
                UnaryFunc::Atan => a.patan(),
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
                    a.pln()
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
                // Scalar equality — not geometric tolerance; no Tolerance param.
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

/// Collect all [`ParameterId`] references in an expression.
pub(crate) fn expression_param_refs(expr: &Expression) -> Vec<ParameterId> {
    match expr {
        Expression::Literal(_) | Expression::Quantity { .. } | Expression::Measure(_) => vec![],
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

/// Topological sort of parameters with expressions (dependencies before dependents).
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
/// Pass [`MeasureContext::empty()`](crate::measure_query::MeasureContext::empty)
/// when no live model is available. On error, last good value is kept with a warning.
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

        // Detect whether this expression contains any MEASURE nodes.
        let has_measure = !crate::measure_query::expression_measure_refs(&expr).is_empty();

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
                // Write resolved value and update last_known for MEASURE expressions.
                if let Some(param) = table.get_mut(id) {
                    param.value = value;
                    param.eval_warning = None;
                    if has_measure {
                        param.last_known = Some(value);
                    }
                }
            }
            Err(msg) => {
                // For MEASURE failures, fall back to last_known if available.
                if has_measure && let Some(last) = table.get(id).and_then(|p| p.last_known) {
                    let warning = EvaluationWarning::MeasureStale {
                        message: msg.clone(),
                    };
                    warnings.push((id, warning.clone()));
                    if let Some(param) = table.get_mut(id) {
                        param.value = last;
                        param.eval_warning = Some(warning);
                    }
                    continue;
                }
                // Normal error path — keep last good value, record MathError.
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

#[path = "expression_dimension.rs"]
mod dimension;
pub use dimension::dimension_kind;

#[cfg(test)]
#[path = "expression_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "expression_dim_tests.rs"]
mod dim_tests;
