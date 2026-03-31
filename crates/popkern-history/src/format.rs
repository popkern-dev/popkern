//! Expression pretty-printer.
//!
//! Converts an [`Expression`] AST back to a human-readable string.  Minimal
//! parentheses are added based on operator precedence.

use crate::expression::{BinOp, CmpOp, Expression, UnaryFunc};
use crate::measure_query::{BoundExtent, Component, MeasureQuery};
use crate::parameter::ParameterTable;
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Precedence levels (lower number = lower precedence)
// ---------------------------------------------------------------------------

/// Operator precedence levels used for parenthesisation decisions.
const PREC_TERNARY: u8 = 1;
const PREC_LOGIC: u8 = 2;
const PREC_COMPARE: u8 = 3;
const PREC_ADD: u8 = 4;
const PREC_MUL: u8 = 5;
const PREC_UNARY: u8 = 6;
/// Atoms (literals, params, function calls) never need outer parens.
const PREC_ATOM: u8 = 7;

/// Return the precedence level for an expression.
fn expr_prec(expr: &Expression) -> u8 {
    match expr {
        Expression::Literal(_) | Expression::Param(_) | Expression::Measure(_) => PREC_ATOM,
        Expression::BinOp(_, op, _) => match op {
            BinOp::And | BinOp::Or => PREC_LOGIC,
            BinOp::Add | BinOp::Sub => PREC_ADD,
            BinOp::Mul | BinOp::Div => PREC_MUL,
            // Min/Max are rendered as function calls → atom precedence.
            BinOp::Min | BinOp::Max => PREC_ATOM,
        },
        Expression::Func(_, _) => PREC_ATOM,
        Expression::Compare(_, _, _) => PREC_COMPARE,
        Expression::Conditional(_, _, _) => PREC_TERNARY,
    }
}

// ---------------------------------------------------------------------------
// Internal recursive formatter
// ---------------------------------------------------------------------------

/// Format `expr`, wrapping in parentheses if its precedence is lower than
/// `min_prec`.
fn fmt_with_prec(
    expr: &Expression,
    table: &ParameterTable,
    min_prec: u8,
    unit: Option<Unit>,
) -> String {
    let s = fmt_inner(expr, table, unit);
    if expr_prec(expr) < min_prec {
        format!("({s})")
    } else {
        s
    }
}

/// Core formatter without parenthesisation wrapping.
fn fmt_inner(expr: &Expression, table: &ParameterTable, unit: Option<Unit>) -> String {
    match expr {
        Expression::Literal(v) => match unit {
            Some(u) if u != Unit::Ratio => {
                let display_val = u.from_si(*v);
                let suffix = u.suffix();
                // Format without trailing zeros but keep enough precision.
                format_number_with_unit(display_val, suffix)
            }
            _ => format_number(*v),
        },
        Expression::Param(id) => table
            .qualified_name(*id)
            .unwrap_or_else(|| "<unknown>".to_string()),

        Expression::BinOp(left, BinOp::Min, right) => {
            let l = fmt_with_prec(left, table, PREC_ATOM, unit);
            let r = fmt_with_prec(right, table, PREC_ATOM, unit);
            format!("min({l}, {r})")
        }
        Expression::BinOp(left, BinOp::Max, right) => {
            let l = fmt_with_prec(left, table, PREC_ATOM, unit);
            let r = fmt_with_prec(right, table, PREC_ATOM, unit);
            format!("max({l}, {r})")
        }
        Expression::BinOp(left, BinOp::And, right) => {
            let l = fmt_with_prec(left, table, PREC_LOGIC, unit);
            let r = fmt_with_prec(right, table, PREC_LOGIC + 1, unit);
            format!("{l} && {r}")
        }
        Expression::BinOp(left, BinOp::Or, right) => {
            let l = fmt_with_prec(left, table, PREC_LOGIC, unit);
            let r = fmt_with_prec(right, table, PREC_LOGIC + 1, unit);
            format!("{l} || {r}")
        }
        Expression::BinOp(left, op, right) => {
            let (prec, sym) = binop_prec_sym(*op);
            let l = fmt_with_prec(left, table, prec, unit);
            let r = fmt_with_prec(right, table, prec + 1, unit);
            format!("{l} {sym} {r}")
        }

        Expression::Func(func, arg) => match func {
            UnaryFunc::Neg => {
                let inner = fmt_with_prec(arg, table, PREC_UNARY, unit);
                format!("-{inner}")
            }
            UnaryFunc::Not => {
                let inner = fmt_with_prec(arg, table, PREC_UNARY, unit);
                format!("!{inner}")
            }
            _ => {
                let name = func_name(*func);
                let inner = fmt_with_prec(arg, table, 0, unit);
                format!("{name}({inner})")
            }
        },

        Expression::Compare(left, op, right) => {
            let sym = cmpop_sym(*op);
            let l = fmt_with_prec(left, table, PREC_COMPARE, unit);
            let r = fmt_with_prec(right, table, PREC_COMPARE + 1, unit);
            format!("{l} {sym} {r}")
        }

        Expression::Conditional(cond, then_expr, else_expr) => {
            let c = fmt_with_prec(cond, table, PREC_TERNARY + 1, unit);
            let t = fmt_with_prec(then_expr, table, PREC_TERNARY, unit);
            let e = fmt_with_prec(else_expr, table, PREC_TERNARY, unit);
            format!("{c} ? {t} : {e}")
        }

        Expression::Measure(query) => format!("MEASURE({})", format_measure_query(query)),
    }
}

// ---------------------------------------------------------------------------
// MEASURE query formatting
// ---------------------------------------------------------------------------

/// Format a [`MeasureQuery`] to the parser-accepted dotted-path syntax.
///
/// Body-level queries (e.g. `BodyVolume`) produce exact round-trip strings
/// like `"body.volume"`. Entity-level queries produce a deterministic
/// placeholder like `"op1_gen0.area"` that can be parsed back if an
/// [`EntityResolver`](crate::parse::EntityResolver) is provided.
pub(crate) fn format_measure_query(query: &MeasureQuery) -> String {
    match query {
        // Body queries
        MeasureQuery::BodyVolume => "body.volume".to_string(),
        MeasureQuery::BodySurfaceArea => "body.surface_area".to_string(),
        MeasureQuery::BodyCentroid(c) => format!("body.centroid.{}", comp_str(c)),
        MeasureQuery::BodyBound(e) => format!("body.bound.{}", bound_str(e)),
        MeasureQuery::BodyMomentOfInertia(a) => format!("body.moment_of_inertia.{}", axis_str(a)),

        // Edge queries
        MeasureQuery::EdgeLength(n) => format!("{}.length", entity_placeholder(n)),
        MeasureQuery::EdgeTangent(n, c) => {
            format!("{}.tangent.{}", entity_placeholder(n), comp_str(c))
        }

        // Face queries
        MeasureQuery::FaceArea(n) => format!("{}.area", entity_placeholder(n)),
        MeasureQuery::FacePerimeter(n) => format!("{}.perimeter", entity_placeholder(n)),
        MeasureQuery::FaceNormal(n, c) => {
            format!("{}.normal.{}", entity_placeholder(n), comp_str(c))
        }
        MeasureQuery::FaceCentroid(n, c) => {
            format!("{}.centroid.{}", entity_placeholder(n), comp_str(c))
        }
        MeasureQuery::FaceRadius(n) => format!("{}.radius", entity_placeholder(n)),
        MeasureQuery::FaceMinCurvature(n) => format!("{}.min_curvature", entity_placeholder(n)),
        MeasureQuery::FaceMaxCurvature(n) => format!("{}.max_curvature", entity_placeholder(n)),

        // Vertex queries
        MeasureQuery::VertexPosition(n, c) => {
            format!("{}.position.{}", entity_placeholder(n), comp_str(c))
        }

        // Pair queries
        MeasureQuery::MinDistance(a, b) => format!(
            "{}, {}.min_distance",
            entity_placeholder(a),
            entity_placeholder(b)
        ),
        MeasureQuery::DihedralAngle(a, b) => format!(
            "{}, {}.dihedral_angle",
            entity_placeholder(a),
            entity_placeholder(b)
        ),
        MeasureQuery::EdgeAngle(a, b) => format!(
            "{}, {}.edge_angle",
            entity_placeholder(a),
            entity_placeholder(b)
        ),
    }
}

fn comp_str(c: &Component) -> &'static str {
    match c {
        Component::X => "x",
        Component::Y => "y",
        Component::Z => "z",
    }
}

fn bound_str(e: &BoundExtent) -> &'static str {
    match e {
        BoundExtent::MinX => "min_x",
        BoundExtent::MinY => "min_y",
        BoundExtent::MinZ => "min_z",
        BoundExtent::MaxX => "max_x",
        BoundExtent::MaxY => "max_y",
        BoundExtent::MaxZ => "max_z",
        BoundExtent::SizeX => "size_x",
        BoundExtent::SizeY => "size_y",
        BoundExtent::SizeZ => "size_z",
    }
}

fn axis_str(a: &popkern_measure::Axis) -> &'static str {
    match a {
        popkern_measure::Axis::X => "x",
        popkern_measure::Axis::Y => "y",
        popkern_measure::Axis::Z => "z",
    }
}

/// Format an [`EntityName`] as a deterministic placeholder string.
///
/// Uses the pattern `op{id}_{category}{index}`, e.g. `op1_gen0`. This
/// placeholder can be parsed back when a matching
/// [`EntityResolver`](crate::parse::EntityResolver) is provided that maps
/// this string to the original [`EntityName`].
fn entity_placeholder(name: &popkern_core::naming::EntityName) -> String {
    use popkern_core::naming::RoleCategory;
    let cat = match name.role.category {
        RoleCategory::Modified => "mod",
        RoleCategory::Generated => "gen",
        RoleCategory::Lateral => "lat",
        RoleCategory::Cap => "cap",
        RoleCategory::Intersection => "int",
        RoleCategory::Junction => "jnc",
    };
    format!("op{}_{}{}", name.operation.0, cat, name.role.index)
}

/// Return the (precedence, symbol) for arithmetic binary operators.
fn binop_prec_sym(op: BinOp) -> (u8, &'static str) {
    match op {
        BinOp::Add => (PREC_ADD, "+"),
        BinOp::Sub => (PREC_ADD, "-"),
        BinOp::Mul => (PREC_MUL, "*"),
        BinOp::Div => (PREC_MUL, "/"),
        BinOp::And | BinOp::Or | BinOp::Min | BinOp::Max => {
            // These are handled before calling this function.
            unreachable!("And/Or/Min/Max handled separately")
        }
    }
}

/// Return the string symbol for a comparison operator.
fn cmpop_sym(op: CmpOp) -> &'static str {
    match op {
        CmpOp::Gt => ">",
        CmpOp::Lt => "<",
        CmpOp::Ge => ">=",
        CmpOp::Le => "<=",
        CmpOp::Eq => "==",
        CmpOp::Ne => "!=",
    }
}

/// Return the lowercase function name for a unary function.
fn func_name(func: UnaryFunc) -> &'static str {
    match func {
        UnaryFunc::Sin => "sin",
        UnaryFunc::Cos => "cos",
        UnaryFunc::Tan => "tan",
        UnaryFunc::Asin => "asin",
        UnaryFunc::Acos => "acos",
        UnaryFunc::Atan => "atan",
        UnaryFunc::Sqrt => "sqrt",
        UnaryFunc::Abs => "abs",
        UnaryFunc::Floor => "floor",
        UnaryFunc::Ceil => "ceil",
        UnaryFunc::Round => "round",
        UnaryFunc::Ln => "ln",
        UnaryFunc::Neg | UnaryFunc::Not => unreachable!("Neg/Not handled separately"),
    }
}

// ---------------------------------------------------------------------------
// Number formatting helpers
// ---------------------------------------------------------------------------

/// Format a raw f64 value with up to 10 significant digits, stripping trailing
/// zeros. Integers are formatted without a decimal point.
fn format_number(v: f64) -> String {
    // Use a fixed representation with enough precision but strip trailing zeros.
    let s = format!("{v:.10}");
    // Strip trailing zeros after decimal.
    if s.contains('.') {
        let s = s.trim_end_matches('0');
        let s = s.trim_end_matches('.');
        s.to_string()
    } else {
        s
    }
}

/// Format a display value with its unit suffix.
fn format_number_with_unit(display_val: f64, suffix: &str) -> String {
    let num = format_number(display_val);
    if suffix.is_empty() {
        num
    } else {
        format!("{num}{suffix}")
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Format an expression to a human-readable string.
///
/// Literals are formatted as plain numbers (no unit suffix). Parameter
/// references are rendered using their qualified name from `table`.
///
/// For unit-annotated output, use [`format_expression_with_unit`].
pub fn format_expression(expr: &Expression, table: &ParameterTable) -> String {
    fmt_inner(expr, table, None)
}

/// Format an expression with a specific display unit context for literals.
///
/// Literal values (stored in SI) are converted to `display_unit` and the unit
/// suffix is appended. For example, `Literal(0.05)` with
/// `Unit::Millimeters` produces `"50mm"`.
///
/// Passing `Unit::Ratio` (or any unit whose `suffix()` is `""`) produces the
/// same output as [`format_expression`].
pub fn format_expression_with_unit(
    expr: &Expression,
    table: &ParameterTable,
    display_unit: Unit,
) -> String {
    fmt_inner(expr, table, Some(display_unit))
}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
