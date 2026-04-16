//! Tests for configuration types and related error variants.

use slotmap::SlotMap;

use crate::FeatureTree;
use crate::config::{ConfigId, Configuration, DeleteConfigResult, ParamOverride};
use crate::error::HistoryError;
use crate::expression::Expression;

// ---------------------------------------------------------------------------
// ConfigId — slotmap key behaviour
// ---------------------------------------------------------------------------

#[test]
fn config_id_slotmap_insert_retrieve() {
    let mut map: SlotMap<ConfigId, &str> = SlotMap::with_key();
    let id = map.insert("hello");
    assert_eq!(map[id], "hello");
}

#[test]
fn config_id_slotmap_remove() {
    let mut map: SlotMap<ConfigId, u32> = SlotMap::with_key();
    let id = map.insert(42);
    assert_eq!(map.remove(id), Some(42));
    assert!(!map.contains_key(id));
}

// ---------------------------------------------------------------------------
// Configuration — serde round-trip
// ---------------------------------------------------------------------------

fn sample_config() -> Configuration {
    Configuration {
        name: "Large".into(),
        description: "Upscaled variant".into(),
        parent: None,
        overrides: Default::default(),
        suppression: Default::default(),
    }
}

#[test]
fn configuration_serde_roundtrip() {
    let cfg = sample_config();
    let json = serde_json::to_string(&cfg).expect("serialize");
    let back: Configuration = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.name, "Large");
    assert_eq!(back.description, "Upscaled variant");
    assert!(back.parent.is_none());
    assert!(back.overrides.is_empty());
    assert!(back.suppression.is_empty());
}

#[test]
fn configuration_with_parent_serde_roundtrip() {
    let mut map: SlotMap<ConfigId, Configuration> = SlotMap::with_key();
    let parent_id = map.insert(sample_config());

    let child = Configuration {
        name: "Large-Steel".into(),
        description: "".into(),
        parent: Some(parent_id),
        overrides: Default::default(),
        suppression: Default::default(),
    };

    let json = serde_json::to_string(&child).expect("serialize");
    let back: Configuration = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.parent, Some(parent_id));
}

// ---------------------------------------------------------------------------
// ParamOverride — serde round-trips
// ---------------------------------------------------------------------------

#[test]
fn param_override_literal_roundtrip() {
    let ov = ParamOverride::Literal(3.5);
    let json = serde_json::to_string(&ov).expect("serialize");
    let back: ParamOverride = serde_json::from_str(&json).expect("deserialize");
    match back {
        ParamOverride::Literal(v) => assert!((v - 3.5).abs() < 1e-12),
        _ => panic!("expected Literal"),
    }
}

#[test]
fn param_override_expression_roundtrip() {
    let expr = Expression::Literal(2.0);
    let ov = ParamOverride::Expression(expr);
    let json = serde_json::to_string(&ov).expect("serialize");
    let back: ParamOverride = serde_json::from_str(&json).expect("deserialize");
    match back {
        ParamOverride::Expression(Expression::Literal(v)) => {
            assert!((v - 2.0).abs() < 1e-12)
        }
        _ => panic!("expected Expression(Literal)"),
    }
}

// ---------------------------------------------------------------------------
// DeleteConfigResult
// ---------------------------------------------------------------------------

#[test]
fn delete_config_result_debug() {
    let r = DeleteConfigResult { reparented: vec![] };
    let s = format!("{r:?}");
    assert!(s.contains("DeleteConfigResult"));
}

// ---------------------------------------------------------------------------
// New HistoryError variants — display
// ---------------------------------------------------------------------------

#[test]
fn error_config_not_found_display() {
    let msg = HistoryError::ConfigNotFound.to_string();
    assert_eq!(msg, "configuration not found");
}

#[test]
fn error_inheritance_cycle_display() {
    let msg = HistoryError::InheritanceCycle.to_string();
    assert_eq!(msg, "configuration inheritance cycle detected");
}

#[test]
fn error_duplicate_config_name_display() {
    let msg = HistoryError::DuplicateConfigName("Foo".into()).to_string();
    assert_eq!(msg, "duplicate configuration name: Foo");
}

// ---------------------------------------------------------------------------
// FeatureTree::new() — configs field is empty, active_config is None
// ---------------------------------------------------------------------------

#[test]
fn feature_tree_new_has_empty_configs() {
    let tree = FeatureTree::new();
    assert!(tree.configs.is_empty());
}

#[test]
fn feature_tree_new_active_config_is_none() {
    let tree = FeatureTree::new();
    assert!(tree.active_config.is_none());
}
