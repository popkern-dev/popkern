//! Expression-driven dimensions.
//!
//! Constraint dimension values (the `f64` in `Distance(p1, p2, target)`,
//! `Angle(l1, l2, target)`, etc.) are normally constants. This module lets
//! them be computed from formulas that reference other constraints' target
//! values.
//!
//! Expressions are stored in a `HashMap<ConstraintHandle, Expression>` on the
//! [`Sketch`](crate::Sketch), **not** inside the `Constraint` enum. The
//! [`evaluate_all`] pre-pass writes the computed values into each constraint's
//! `f64` field before the solver runs, so the Newton loop sees plain constants.

use crate::Sketch;
use crate::config::{get_constraint_value, set_constraint_value_on};
use crate::constraint::ConstraintHandle;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// A mathematical expression for computing dimension values.
///
/// Expressions form a tree of operations whose leaves are either literal
/// constants or references to the target value of another constraint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Expression {
    /// Reference another constraint's target value.
    Ref(ConstraintHandle),
    /// Literal constant.
    Literal(f64),
    /// Binary operation on two sub-expressions.
    BinOp(Box<Expression>, BinOp, Box<Expression>),
    /// Unary function applied to a sub-expression.
    Func(UnaryFunc, Box<Expression>),
}

/// Binary arithmetic operations.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum BinOp {
    /// Addition.
    Add,
    /// Subtraction.
    Sub,
    /// Multiplication.
    Mul,
    /// Division.
    Div,
}

/// Unary mathematical functions.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum UnaryFunc {
    /// Sine (radians).
    Sin,
    /// Cosine (radians).
    Cos,
    /// Tangent (radians).
    Tan,
    /// Square root (non-negative input only).
    Sqrt,
    /// Absolute value.
    Abs,
    /// Negation.
    Neg,
}

/// Evaluate an expression, resolving `Ref` to constraint target values.
///
/// Works regardless of group/suppression state -- the target value exists
/// even when the constraint is not participating in the solver.
pub(crate) fn evaluate(expr: &Expression, sketch: &Sketch) -> Result<f64, String> {
    match expr {
        Expression::Literal(v) => Ok(*v),
        Expression::Ref(handle) => {
            let entry = sketch.constraints.get(*handle).ok_or_else(|| {
                "expression references a constraint that does not exist".to_owned()
            })?;
            get_constraint_value(&entry.constraint)
                .ok_or_else(|| "referenced constraint does not have a dimensional value".to_owned())
        }
        Expression::BinOp(lhs, op, rhs) => {
            let l = evaluate(lhs, sketch)?;
            let r = evaluate(rhs, sketch)?;
            match op {
                BinOp::Add => Ok(l + r),
                BinOp::Sub => Ok(l - r),
                BinOp::Mul => Ok(l * r),
                BinOp::Div => {
                    if r == 0.0 {
                        return Err("division by zero in expression".to_owned());
                    }
                    Ok(l / r)
                }
            }
        }
        Expression::Func(func, inner) => {
            let v = evaluate(inner, sketch)?;
            match func {
                UnaryFunc::Sin => Ok(v.sin()),
                UnaryFunc::Cos => Ok(v.cos()),
                UnaryFunc::Tan => Ok(v.tan()),
                UnaryFunc::Sqrt => {
                    if v < 0.0 {
                        return Err("square root of negative number in expression".to_owned());
                    }
                    Ok(v.sqrt())
                }
                UnaryFunc::Abs => Ok(v.abs()),
                UnaryFunc::Neg => Ok(-v),
            }
        }
    }
}

/// Check if setting an expression on `target` would create a circular
/// dependency.
///
/// Walks the expression tree. For each `Ref(h)`: if `h == target`, there
/// is a cycle. If `h` has its own expression in `expressions`, that
/// expression is checked recursively. A `HashSet` tracks visited handles
/// to avoid re-traversal on diamond dependencies.
pub(crate) fn has_cycle(
    target: ConstraintHandle,
    expr: &Expression,
    expressions: &HashMap<ConstraintHandle, Expression>,
) -> bool {
    let mut visited = HashSet::new();
    has_cycle_inner(target, expr, expressions, &mut visited)
}

/// Recursive cycle-detection helper.
fn has_cycle_inner(
    target: ConstraintHandle,
    expr: &Expression,
    expressions: &HashMap<ConstraintHandle, Expression>,
    visited: &mut HashSet<ConstraintHandle>,
) -> bool {
    match expr {
        Expression::Literal(_) => false,
        Expression::Ref(h) => {
            if *h == target {
                return true;
            }
            if !visited.insert(*h) {
                // Already checked this handle -- no cycle through this path.
                return false;
            }
            if let Some(dep_expr) = expressions.get(h) {
                has_cycle_inner(target, dep_expr, expressions, visited)
            } else {
                false
            }
        }
        Expression::BinOp(lhs, _, rhs) => {
            has_cycle_inner(target, lhs, expressions, visited)
                || has_cycle_inner(target, rhs, expressions, visited)
        }
        Expression::Func(_, inner) => has_cycle_inner(target, inner, expressions, visited),
    }
}

/// Evaluate all expressions and write computed values into constraints.
///
/// Topologically sorts by dependency, then evaluates in order. Each
/// computed value is written into the constraint variant's first `f64`
/// field (via [`set_constraint_value_on`]) so the solver sees constant
/// values.
///
/// Note: for constraints with multiple `f64` fields (e.g.
/// `RelativePosition(dx, dy)`), expressions drive the **first** field
/// only.
pub(crate) fn evaluate_all(sketch: &mut Sketch) -> Result<(), String> {
    if sketch.expressions.is_empty() {
        return Ok(());
    }

    // Collect handles that have expressions.
    let expr_handles: Vec<ConstraintHandle> = sketch.expressions.keys().copied().collect();

    // Topological sort via depth-first post-order traversal.
    let order = topological_sort(&expr_handles, &sketch.expressions)?;

    // Evaluate in topological order and write back.
    for handle in order {
        // Clone expression to avoid borrow conflict.
        let expr = match sketch.expressions.get(&handle) {
            Some(e) => e.clone(),
            None => continue,
        };
        let computed = evaluate(&expr, sketch)?;
        let entry = sketch
            .constraints
            .get_mut(handle)
            .ok_or_else(|| "expression target constraint no longer exists".to_owned())?;
        if !set_constraint_value_on(&mut entry.constraint, computed) {
            return Err("expression target constraint is not dimensional".to_owned());
        }
    }

    Ok(())
}

/// Topological sort of expression handles by dependency order.
///
/// Returns handles in an order where each handle's `Ref` dependencies
/// (that also have expressions) come earlier in the list.
fn topological_sort(
    handles: &[ConstraintHandle],
    expressions: &HashMap<ConstraintHandle, Expression>,
) -> Result<Vec<ConstraintHandle>, String> {
    let handle_set: HashSet<ConstraintHandle> = handles.iter().copied().collect();
    let mut visited = HashSet::new();
    let mut in_stack = HashSet::new();
    let mut order = Vec::with_capacity(handles.len());

    for &handle in handles {
        if !visited.contains(&handle) {
            topo_visit(
                handle,
                expressions,
                &handle_set,
                &mut visited,
                &mut in_stack,
                &mut order,
            )?;
        }
    }

    Ok(order)
}

/// DFS visitor for topological sort.
fn topo_visit(
    handle: ConstraintHandle,
    expressions: &HashMap<ConstraintHandle, Expression>,
    handle_set: &HashSet<ConstraintHandle>,
    visited: &mut HashSet<ConstraintHandle>,
    in_stack: &mut HashSet<ConstraintHandle>,
    order: &mut Vec<ConstraintHandle>,
) -> Result<(), String> {
    if visited.contains(&handle) {
        return Ok(());
    }
    if !in_stack.insert(handle) {
        return Err("circular dependency detected during expression evaluation".to_owned());
    }

    if let Some(expr) = expressions.get(&handle) {
        let deps = collect_refs(expr);
        for dep in deps {
            if handle_set.contains(&dep) {
                topo_visit(dep, expressions, handle_set, visited, in_stack, order)?;
            }
        }
    }

    in_stack.remove(&handle);
    visited.insert(handle);
    order.push(handle);
    Ok(())
}

/// Collect all `Ref` handles from an expression tree.
fn collect_refs(expr: &Expression) -> Vec<ConstraintHandle> {
    let mut refs = Vec::new();
    collect_refs_inner(expr, &mut refs);
    refs
}

/// Recursive helper to collect `Ref` handles.
fn collect_refs_inner(expr: &Expression, refs: &mut Vec<ConstraintHandle>) {
    match expr {
        Expression::Literal(_) => {}
        Expression::Ref(h) => refs.push(*h),
        Expression::BinOp(lhs, _, rhs) => {
            collect_refs_inner(lhs, refs);
            collect_refs_inner(rhs, refs);
        }
        Expression::Func(_, inner) => collect_refs_inner(inner, refs),
    }
}

#[cfg(test)]
#[path = "expression_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "expression_cycle_tests.rs"]
mod cycle_tests;
