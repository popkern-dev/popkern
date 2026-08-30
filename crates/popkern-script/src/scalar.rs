//! Scalar bindings (spec §1.3): every numeric argument of a feature op
//! is a literal, a reference to a named parameter, or an expression.
//! All three record as a feature-scoped parameter so the value stays
//! editable after the script has run.

use popkern_history::FeatureTree;
use popkern_history::expression::Expression;
use popkern_history::parameter::ParameterId;
use popkern_history::unit::Unit;

use crate::error::ScriptError;

/// A numeric argument.
#[derive(Debug, Clone, PartialEq)]
pub enum Scalar {
    /// A literal in document units.
    Literal(f64),
    /// The value of a named global parameter (by qualified name).
    Param(String),
    /// An expression over parameters, in the surface's expression
    /// language.
    Expr(String),
}

impl From<f64> for Scalar {
    fn from(v: f64) -> Self {
        Scalar::Literal(v)
    }
}

impl From<&str> for Scalar {
    /// A string is a parameter name when it resolves to one, else an
    /// expression.
    fn from(s: &str) -> Self {
        Scalar::Expr(s.to_string())
    }
}

/// Reference a global parameter by name.
pub fn p(name: &str) -> Scalar {
    Scalar::Param(name.to_string())
}

/// An expression argument.
pub fn expr(text: &str) -> Scalar {
    Scalar::Expr(text.to_string())
}

/// Create the parameter a feature op will own for `role`.
///
/// The parameter is added as a global under a unique pending name;
/// `FeatureTree::add_feature` re-scopes it to the feature and names it
/// `{feature}.{role}` (the history layer's registration rule). A
/// `Param` binding never re-scopes the referenced global — it records
/// as an expression `= <global>` so the global stays shared.
pub(crate) fn bind_scalar(
    tree: &mut FeatureTree,
    counter: &mut u64,
    role: &str,
    scalar: &Scalar,
    unit: Unit,
) -> Result<ParameterId, ScriptError> {
    *counter += 1;
    let pending = format!("__pending_{counter}_{role}");
    match scalar {
        Scalar::Literal(v) => {
            if !v.is_finite() {
                return Err(ScriptError::Argument {
                    argument: role.to_string(),
                    reason: format!("value must be finite, got {v}"),
                });
            }
            Ok(tree.add_global_parameter(&pending, *v, unit)?)
        }
        Scalar::Param(name) => {
            let target =
                tree.parameters()
                    .resolve_name(name)
                    .ok_or_else(|| ScriptError::Argument {
                        argument: role.to_string(),
                        reason: format!("unknown parameter '{name}'"),
                    })?;
            let value = tree
                .parameters()
                .get(target)
                .map(|p| p.value)
                .unwrap_or(0.0);
            let id = tree.add_global_parameter(&pending, value, unit)?;
            tree.set_parameter_expression(id, Expression::Param(target))?;
            Ok(id)
        }
        Scalar::Expr(text) => {
            let expression = popkern_history::parse::parse_expression(text, tree.parameters())
                .map_err(|e| ScriptError::Argument {
                    argument: role.to_string(),
                    reason: format!("expression '{text}': {e}"),
                })?;
            let id = tree.add_global_parameter(&pending, 0.0, unit)?;
            tree.set_parameter_expression(id, expression)?;
            Ok(id)
        }
    }
}
