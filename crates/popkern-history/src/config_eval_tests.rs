// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//
// Tests for config_eval: resolve_override, resolve_suppression,
// build_effective_expressions, evaluate_config, would_create_cycle.

use slotmap::SlotMap;

use crate::config::{ConfigId, Configuration, ParamOverride};
use crate::config_eval::{
    build_effective_expressions, evaluate_config, resolve_override, resolve_suppression,
    would_create_cycle,
};
use crate::expression::Expression;
use crate::node::FeatureId;
use crate::parameter::{ParameterId, ParameterTable};
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn empty_config(name: &str, parent: Option<ConfigId>) -> Configuration {
    Configuration {
        name: name.into(),
        description: String::new(),
        parent,
        overrides: Default::default(),
        suppression: Default::default(),
    }
}

fn make_table_with_param(name: &str, value: f64) -> (ParameterTable, ParameterId) {
    let mut table = ParameterTable::new();
    let pid = table
        .add_global(name, value, Unit::Meters)
        .expect("add param");
    (table, pid)
}

// ---------------------------------------------------------------------------
// resolve_override: direct override on config
// ---------------------------------------------------------------------------

#[test]
fn resolve_override_direct_hit() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let (table, pid) = make_table_with_param("width", 1.0);

    let mut cfg = empty_config("Big", None);
    cfg.overrides.insert(pid, ParamOverride::Literal(2.0));
    let cid = configs.insert(cfg);

    let result = resolve_override(pid, Some(cid), &configs);
    assert!(result.is_some(), "should find direct override");
    match result.unwrap() {
        ParamOverride::Literal(v) => assert!((v - 2.0).abs() < 1e-12),
        _ => panic!("expected Literal"),
    }

    // Suppress unused variable warning
    let _ = table;
}

#[test]
fn resolve_override_no_match_returns_none() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let (table, pid) = make_table_with_param("height", 1.0);

    let cfg = empty_config("Empty", None);
    let cid = configs.insert(cfg);

    let result = resolve_override(pid, Some(cid), &configs);
    assert!(result.is_none(), "should return None when no override");

    let _ = table;
}

#[test]
fn resolve_override_none_config_returns_none() {
    let configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let (table, pid) = make_table_with_param("depth", 1.0);

    let result = resolve_override(pid, None, &configs);
    assert!(result.is_none(), "no config → no override");

    let _ = table;
}

// ---------------------------------------------------------------------------
// resolve_override: inheritance chain
// ---------------------------------------------------------------------------

#[test]
fn resolve_override_inherits_from_parent() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let (table, pid) = make_table_with_param("radius", 0.5);

    // Parent has override, child does not.
    let mut parent_cfg = empty_config("Parent", None);
    parent_cfg
        .overrides
        .insert(pid, ParamOverride::Literal(0.75));
    let parent_id = configs.insert(parent_cfg);

    let child_cfg = empty_config("Child", Some(parent_id));
    let child_id = configs.insert(child_cfg);

    let result = resolve_override(pid, Some(child_id), &configs);
    assert!(result.is_some(), "should inherit from parent");
    match result.unwrap() {
        ParamOverride::Literal(v) => assert!((v - 0.75).abs() < 1e-12),
        _ => panic!("expected Literal"),
    }

    let _ = table;
}

#[test]
fn resolve_override_child_overrides_parent() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let (table, pid) = make_table_with_param("angle", 1.0);

    let mut parent_cfg = empty_config("Parent", None);
    parent_cfg
        .overrides
        .insert(pid, ParamOverride::Literal(10.0));
    let parent_id = configs.insert(parent_cfg);

    let mut child_cfg = empty_config("Child", Some(parent_id));
    child_cfg
        .overrides
        .insert(pid, ParamOverride::Literal(20.0));
    let child_id = configs.insert(child_cfg);

    let result = resolve_override(pid, Some(child_id), &configs);
    match result.unwrap() {
        ParamOverride::Literal(v) => assert!((v - 20.0).abs() < 1e-12, "child wins"),
        _ => panic!("expected Literal"),
    }

    let _ = table;
}

// ---------------------------------------------------------------------------
// resolve_suppression
// ---------------------------------------------------------------------------

#[test]
fn resolve_suppression_direct_hit() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut features: slotmap::SlotMap<FeatureId, ()> = slotmap::SlotMap::with_key();
    let fid = features.insert(());

    let mut cfg = empty_config("Cfg", None);
    cfg.suppression.insert(fid, true);
    let cid = configs.insert(cfg);

    assert_eq!(resolve_suppression(fid, Some(cid), &configs), Some(true));
}

#[test]
fn resolve_suppression_inherits_from_parent() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut features: slotmap::SlotMap<FeatureId, ()> = slotmap::SlotMap::with_key();
    let fid = features.insert(());

    let mut parent_cfg = empty_config("Parent", None);
    parent_cfg.suppression.insert(fid, false);
    let parent_id = configs.insert(parent_cfg);

    let child_cfg = empty_config("Child", Some(parent_id));
    let child_id = configs.insert(child_cfg);

    assert_eq!(
        resolve_suppression(fid, Some(child_id), &configs),
        Some(false),
        "should inherit suppression=false from parent"
    );
}

#[test]
fn resolve_suppression_no_override_returns_none() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut features: slotmap::SlotMap<FeatureId, ()> = slotmap::SlotMap::with_key();
    let fid = features.insert(());

    let cfg = empty_config("Empty", None);
    let cid = configs.insert(cfg);

    assert_eq!(resolve_suppression(fid, Some(cid), &configs), None);
}

// ---------------------------------------------------------------------------
// build_effective_expressions
// ---------------------------------------------------------------------------

#[test]
fn build_effective_expressions_separates_literal_and_expression() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();

    let pid_a = table.add_global("a", 1.0, Unit::Meters).unwrap();
    let pid_b = table.add_global("b", 2.0, Unit::Meters).unwrap();

    let mut cfg = empty_config("Cfg", None);
    cfg.overrides.insert(pid_a, ParamOverride::Literal(5.0));
    cfg.overrides
        .insert(pid_b, ParamOverride::Expression(Expression::Literal(7.0)));
    let cid = configs.insert(cfg);

    let (literals, expressions) = build_effective_expressions(Some(cid), &configs, &table);

    assert_eq!(literals.len(), 1, "one literal override");
    assert_eq!(expressions.len(), 1, "one expression override");
    assert!((literals[&pid_a] - 5.0).abs() < 1e-12);
    assert!(matches!(expressions[&pid_b], Expression::Literal(v) if (v - 7.0).abs() < 1e-12));
}

#[test]
fn build_effective_expressions_no_config_returns_empty() {
    let configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();
    let _ = table.add_global("x", 3.0, Unit::Ratio).unwrap();

    let (literals, expressions) = build_effective_expressions(None, &configs, &table);
    assert!(literals.is_empty());
    assert!(expressions.is_empty());
}

// ---------------------------------------------------------------------------
// evaluate_config
// ---------------------------------------------------------------------------

#[test]
fn evaluate_config_no_config_returns_defaults() {
    let configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();
    let pid = table.add_global("x", 3.0, Unit::Ratio).unwrap();

    let values = evaluate_config(None, &configs, &table);
    assert!((values[&pid] - 3.0).abs() < 1e-12);
}

#[test]
fn evaluate_config_applies_literal_override() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();
    let pid = table.add_global("width", 1.0, Unit::Meters).unwrap();

    let mut cfg = empty_config("Big", None);
    cfg.overrides.insert(pid, ParamOverride::Literal(10.0));
    let cid = configs.insert(cfg);

    let values = evaluate_config(Some(cid), &configs, &table);
    assert!(
        (values[&pid] - 10.0).abs() < 1e-12,
        "literal override should pin value to 10.0, got {}",
        values[&pid]
    );
}

#[test]
fn evaluate_config_applies_expression_override() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();
    let pid_base = table.add_global("base", 5.0, Unit::Ratio).unwrap();
    let pid_derived = table.add_global("derived", 0.0, Unit::Ratio).unwrap();

    // Default: derived has no expression.
    // Config override: derived = Literal(42.0)
    let mut cfg = empty_config("Override", None);
    cfg.overrides.insert(
        pid_derived,
        ParamOverride::Expression(Expression::Literal(42.0)),
    );
    let cid = configs.insert(cfg);

    let values = evaluate_config(Some(cid), &configs, &table);
    assert!(
        (values[&pid_derived] - 42.0).abs() < 1e-12,
        "expression override should evaluate to 42.0, got {}",
        values[&pid_derived]
    );
    // Base param unaffected.
    assert!((values[&pid_base] - 5.0).abs() < 1e-12);
}

#[test]
fn evaluate_config_does_not_mutate_original_table() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let mut table = ParameterTable::new();
    let pid = table.add_global("val", 1.0, Unit::Ratio).unwrap();

    let mut cfg = empty_config("Cfg", None);
    cfg.overrides.insert(pid, ParamOverride::Literal(99.0));
    let cid = configs.insert(cfg);

    let _values = evaluate_config(Some(cid), &configs, &table);

    // Original table must be unchanged.
    let original_val = table.get(pid).map(|p| p.value).unwrap_or(0.0);
    assert!(
        (original_val - 1.0).abs() < 1e-12,
        "original table should not be mutated, got {original_val}"
    );
}

// ---------------------------------------------------------------------------
// would_create_cycle
// ---------------------------------------------------------------------------

#[test]
fn would_create_cycle_detects_direct_cycle() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let a = configs.insert(empty_config("A", None));
    // Setting a.parent = a would be a cycle.
    assert!(would_create_cycle(a, a, &configs));
}

#[test]
fn would_create_cycle_detects_indirect_cycle() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let a = configs.insert(empty_config("A", None));
    let b = configs.insert(empty_config("B", Some(a)));
    // Setting a.parent = b: the chain b→a→b would cycle.
    assert!(would_create_cycle(a, b, &configs));
}

#[test]
fn would_create_cycle_allows_valid_parent() {
    let mut configs: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let a = configs.insert(empty_config("A", None));
    let b = configs.insert(empty_config("B", None));
    // Setting b.parent = a is fine — no cycle.
    assert!(!would_create_cycle(b, a, &configs));
}
