//! Configuration override resolution and config-aware parameter evaluation.
//!
//! These helpers implement the three-level override chain:
//!   1. Walk the config inheritance chain for a given parameter/feature.
//!   2. Build the effective expression map (literal vs expression overrides).
//!   3. Evaluate all parameters for a specific config by cloning the table,
//!      applying overrides, calling [`expression::evaluate_all`], and returning
//!      the resulting value map.

use std::collections::HashMap;

use slotmap::SlotMap;

use crate::config::{ConfigId, Configuration, ParamOverride};
use crate::expression::Expression;
use crate::node::FeatureId;
use crate::parameter::{ParameterId, ParameterTable};

/// Resolve a parameter override by walking the config inheritance chain.
///
/// Returns `None` if no override is found anywhere in the chain (use default).
pub(crate) fn resolve_override(
    param: ParameterId,
    config: Option<ConfigId>,
    configs: &SlotMap<ConfigId, Configuration>,
) -> Option<&ParamOverride> {
    let mut current = config;
    while let Some(cid) = current {
        if let Some(cfg) = configs.get(cid) {
            if let Some(ov) = cfg.overrides.get(&param) {
                return Some(ov);
            }
            current = cfg.parent;
        } else {
            break;
        }
    }
    None
}

/// Resolve a feature suppression override by walking the config inheritance chain.
///
/// Returns `None` if no suppression override is found (use the node's default).
pub(crate) fn resolve_suppression(
    feature: FeatureId,
    config: Option<ConfigId>,
    configs: &SlotMap<ConfigId, Configuration>,
) -> Option<bool> {
    let mut current = config;
    while let Some(cid) = current {
        if let Some(cfg) = configs.get(cid) {
            if let Some(&sup) = cfg.suppression.get(&feature) {
                return Some(sup);
            }
            current = cfg.parent;
        } else {
            break;
        }
    }
    None
}

/// Build the effective expression map for a config.
///
/// Returns `(literal_overrides, expression_overrides)`:
/// - `literal_overrides`: params pinned to a value (skip expression eval).
/// - `expression_overrides`: params with a different expression than default.
pub(crate) fn build_effective_expressions(
    config: Option<ConfigId>,
    configs: &SlotMap<ConfigId, Configuration>,
    table: &ParameterTable,
) -> (HashMap<ParameterId, f64>, HashMap<ParameterId, Expression>) {
    let mut literals = HashMap::new();
    let mut expressions = HashMap::new();

    for (pid, _param) in table.all_params() {
        if let Some(ov) = resolve_override(pid, config, configs) {
            match ov {
                ParamOverride::Literal(v) => {
                    literals.insert(pid, *v);
                }
                ParamOverride::Expression(expr) => {
                    expressions.insert(pid, expr.clone());
                }
            }
        }
    }

    (literals, expressions)
}

/// Evaluate all parameters for a specific config.
///
/// Clones the parameter table, applies config overrides, evaluates all
/// expressions, and returns the resulting value map.  This is an expensive
/// operation — use [`FeatureTree::effective_value`] for the common case of
/// reading a single already-evaluated value.
pub(crate) fn evaluate_config(
    config: Option<ConfigId>,
    configs: &SlotMap<ConfigId, Configuration>,
    table: &ParameterTable,
) -> HashMap<ParameterId, f64> {
    if config.is_none() {
        // No config — return default values after evaluation.
        let mut clone = table.clone();
        let _ = crate::expression::evaluate_all(
            &mut clone,
            &crate::measure_query::MeasureContext::empty(),
        );
        return clone.all_params().map(|(pid, p)| (pid, p.value)).collect();
    }

    let (literals, expr_overrides) = build_effective_expressions(config, configs, table);

    // Clone and apply overrides.
    let mut clone = table.clone();
    for (pid, v) in &literals {
        if let Some(param) = clone.get_mut(*pid) {
            param.value = *v;
            // Clear expression so evaluate_all skips this param (value is
            // already pinned by the literal override).
            param.expression = None;
        }
    }
    for (pid, expr) in &expr_overrides {
        if let Some(param) = clone.get_mut(*pid) {
            param.expression = Some(expr.clone());
        }
    }

    // Evaluate in topological order.
    let _ =
        crate::expression::evaluate_all(&mut clone, &crate::measure_query::MeasureContext::empty());

    clone.all_params().map(|(pid, p)| (pid, p.value)).collect()
}

/// Check if adding `new_parent` as the parent of `config` would create a cycle.
///
/// Returns `true` if setting `config.parent = new_parent` would introduce a
/// cycle in the inheritance chain.
#[allow(dead_code)] // Used by Task 5 (config management API)
pub(crate) fn would_create_cycle(
    config: ConfigId,
    new_parent: ConfigId,
    configs: &SlotMap<ConfigId, Configuration>,
) -> bool {
    let mut current = Some(new_parent);
    while let Some(cid) = current {
        if cid == config {
            return true;
        }
        current = configs.get(cid).and_then(|c| c.parent);
    }
    false
}

#[cfg(test)]
#[path = "config_eval_tests.rs"]
mod tests;
