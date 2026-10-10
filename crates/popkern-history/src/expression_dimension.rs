//! Dimensional type checking of [`Expression`]s (split from
//! `expression.rs` under the 500-line cap).

use super::{BinOp, Expression, UnaryFunc};
use crate::error::HistoryError;
use crate::parameter::{ParameterId, ParameterTable};
use crate::unit::DimensionKind;

/// Determine the [`DimensionKind`] that an expression produces.
pub fn dimension_kind(
    expr: &Expression,
    table: &ParameterTable,
) -> Result<DimensionKind, HistoryError> {
    match expr {
        Expression::Literal(_) => Ok(DimensionKind::Dimensionless),
        Expression::Quantity { unit, .. } => Ok(unit.kind()),
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
                    } else if let Some(dimension) = bare_literal_against(left, lk, right, rk) {
                        Err(HistoryError::UnitlessLiteralInSum { dimension })
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

/// The dimension a BARE number is summed with, when one side of a sum is
/// an unsuffixed literal and the other carries a dimension: `L - 180`
/// means `L - 180mm`, and the error says so (PK-754).
fn bare_literal_against(
    left: &Expression,
    lk: DimensionKind,
    right: &Expression,
    rk: DimensionKind,
) -> Option<DimensionKind> {
    match (left, right) {
        (Expression::Literal(_), _) if rk != DimensionKind::Dimensionless => Some(rk),
        (_, Expression::Literal(_)) if lk != DimensionKind::Dimensionless => Some(lk),
        _ => None,
    }
}
