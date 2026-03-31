use super::{ParamScope, ParameterTable};
use crate::error::HistoryError;
use crate::unit::{DimensionKind, Unit};
use slotmap::new_key_type;

// ── Built-in constants ────────────────────────────────────────────────────────

#[test]
fn new_has_three_builtins() {
    let table = ParameterTable::new();
    // pi, e, tau
    assert_eq!(table.all_params().count(), 3);
}

#[test]
fn pi_is_readonly() {
    let table = ParameterTable::new();
    let id = table.resolve_name("pi").expect("pi must exist");
    let param = table.get(id).unwrap();
    assert!(param.read_only);
}

#[test]
fn pi_value_approx() {
    let table = ParameterTable::new();
    let id = table.resolve_name("pi").unwrap();
    let param = table.get(id).unwrap();
    assert!((param.value - std::f64::consts::PI).abs() < 1e-15);
}

#[test]
fn e_value_approx() {
    let table = ParameterTable::new();
    let id = table.resolve_name("e").unwrap();
    let param = table.get(id).unwrap();
    assert!((param.value - std::f64::consts::E).abs() < 1e-15);
}

#[test]
fn tau_value_approx() {
    let table = ParameterTable::new();
    let id = table.resolve_name("tau").unwrap();
    let param = table.get(id).unwrap();
    assert!((param.value - std::f64::consts::TAU).abs() < 1e-15);
}

// ── add_global + resolve_name ─────────────────────────────────────────────────

#[test]
fn add_global_and_resolve_round_trip() {
    let mut table = ParameterTable::new();
    let id = table.add_global("width", 0.05, Unit::Meters).unwrap();
    let resolved = table.resolve_name("width");
    assert_eq!(resolved, Some(id));
    let param = table.get(id).unwrap();
    assert!((param.value - 0.05).abs() < 1e-15);
    assert_eq!(param.unit, Unit::Meters);
    assert_eq!(param.scope, ParamScope::Global);
}

#[test]
fn add_global_duplicate_name_errors() {
    let mut table = ParameterTable::new();
    table.add_global("radius", 1.0, Unit::Meters).unwrap();
    let result = table.add_global("radius", 2.0, Unit::Meters);
    assert!(
        matches!(result, Err(HistoryError::DuplicateParameterName(ref s)) if s == "radius"),
        "expected DuplicateParameterName, got {result:?}"
    );
}

#[test]
fn add_global_duplicate_with_builtin_errors() {
    let mut table = ParameterTable::new();
    let result = table.add_global("pi", 3.0, Unit::Ratio);
    assert!(matches!(
        result,
        Err(HistoryError::DuplicateParameterName(_))
    ));
}

// ── add_feature_param ────────────────────────────────────────────────────────

new_key_type! { struct TestFeatureId; }

// We need a FeatureId. Use the real one from the crate.
use crate::node::FeatureId;
use slotmap::SlotMap;

/// Allocate `n` distinct FeatureIds from a single SlotMap.
fn make_feature_ids(n: usize) -> Vec<FeatureId> {
    let mut sm: SlotMap<FeatureId, ()> = SlotMap::with_key();
    (0..n).map(|_| sm.insert(())).collect()
}

fn make_feature_id() -> FeatureId {
    make_feature_ids(1).remove(0)
}

#[test]
fn add_feature_param_and_feature_params_list() {
    let mut table = ParameterTable::new();
    let fid = make_feature_id();
    let id = table
        .add_feature_param("extrude1.distance", fid, 0.01, Unit::Meters)
        .unwrap();
    let list = table.feature_params(fid);
    assert_eq!(list, vec![id]);
    assert_eq!(table.get(id).unwrap().scope, ParamScope::Feature(fid));
}

#[test]
fn add_feature_param_short_name_extracted() {
    let mut table = ParameterTable::new();
    let fid = make_feature_id();
    let id = table
        .add_feature_param("extrude1.distance", fid, 0.01, Unit::Meters)
        .unwrap();
    assert_eq!(table.get(id).unwrap().name, "distance");
}

#[test]
fn add_feature_param_duplicate_errors() {
    let mut table = ParameterTable::new();
    let fid = make_feature_id();
    table
        .add_feature_param("feat1.h", fid, 1.0, Unit::Meters)
        .unwrap();
    let result = table.add_feature_param("feat1.h", fid, 2.0, Unit::Meters);
    assert!(matches!(
        result,
        Err(HistoryError::DuplicateParameterName(_))
    ));
}

// ── set_value ─────────────────────────────────────────────────────────────────

#[test]
fn set_value_on_readonly_errors() {
    let mut table = ParameterTable::new();
    let id = table.resolve_name("pi").unwrap();
    let result = table.set_value(id, 3.0, Unit::Ratio);
    assert!(matches!(result, Err(HistoryError::ReadOnlyParameter)));
}

#[test]
fn set_value_wrong_dimension_errors() {
    let mut table = ParameterTable::new();
    // Add a length parameter.
    let id = table.add_global("d", 1.0, Unit::Meters).unwrap();
    // Try setting it to an angle unit.
    let result = table.set_value(id, 45.0, Unit::Degrees);
    assert!(
        matches!(
            result,
            Err(HistoryError::DimensionMismatch {
                expected: DimensionKind::Length,
                got: DimensionKind::Angle,
                ..
            })
        ),
        "expected DimensionMismatch, got {result:?}"
    );
}

#[test]
fn set_value_converts_to_si() {
    let mut table = ParameterTable::new();
    let id = table.add_global("w", 1000.0, Unit::Millimeters).unwrap();
    // Now set to 500 mm — stored as 0.5 m.
    table.set_value(id, 500.0, Unit::Millimeters).unwrap();
    let param = table.get(id).unwrap();
    assert!((param.value - 0.5).abs() < 1e-12);
}

#[test]
fn set_value_accepts_same_kind_different_unit() {
    let mut table = ParameterTable::new();
    let id = table.add_global("d", 1.0, Unit::Meters).unwrap();
    // Set using millimeters (same dimension kind = Length).
    table.set_value(id, 25.4, Unit::Millimeters).unwrap();
    let param = table.get(id).unwrap();
    assert!((param.value - 0.0254).abs() < 1e-12);
}

// ── remove ────────────────────────────────────────────────────────────────────

#[test]
fn remove_ok_when_no_dependents() {
    let mut table = ParameterTable::new();
    let id = table.add_global("temp", 42.0, Unit::Ratio).unwrap();
    assert!(table.remove(id).is_ok());
    assert!(table.get(id).is_none());
    assert!(table.resolve_name("temp").is_none());
}

#[test]
fn remove_builtin_errors() {
    let mut table = ParameterTable::new();
    let id = table.resolve_name("e").unwrap();
    assert!(matches!(
        table.remove(id),
        Err(HistoryError::ReadOnlyParameter)
    ));
}

// ── remove_feature_params ─────────────────────────────────────────────────────

#[test]
fn remove_feature_params_removes_all_for_feature() {
    let mut table = ParameterTable::new();
    let mut fids = make_feature_ids(2);
    let fid = fids.remove(0);
    let fid2 = fids.remove(0);
    table
        .add_feature_param("f1.a", fid, 1.0, Unit::Meters)
        .unwrap();
    table
        .add_feature_param("f1.b", fid, 2.0, Unit::Meters)
        .unwrap();
    table
        .add_feature_param("f2.a", fid2, 3.0, Unit::Meters)
        .unwrap();

    table.remove_feature_params(fid);

    assert_eq!(table.feature_params(fid).len(), 0);
    // fid2's param should survive.
    assert_eq!(table.feature_params(fid2).len(), 1);
    // Built-ins should survive.
    assert!(table.resolve_name("pi").is_some());
}

// ── rename ────────────────────────────────────────────────────────────────────

#[test]
fn rename_updates_name_index() {
    let mut table = ParameterTable::new();
    let id = table.add_global("old_name", 1.0, Unit::Meters).unwrap();
    table.rename(id, "new_name").unwrap();
    assert_eq!(table.resolve_name("new_name"), Some(id));
    assert!(table.resolve_name("old_name").is_none());
}

#[test]
fn rename_to_existing_name_errors() {
    let mut table = ParameterTable::new();
    table.add_global("foo", 1.0, Unit::Meters).unwrap();
    let id2 = table.add_global("bar", 2.0, Unit::Meters).unwrap();
    let result = table.rename(id2, "foo");
    assert!(matches!(
        result,
        Err(HistoryError::DuplicateParameterName(_))
    ));
}

// ── qualified_name ────────────────────────────────────────────────────────────

#[test]
fn qualified_name_global() {
    let mut table = ParameterTable::new();
    let id = table.add_global("depth", 0.1, Unit::Meters).unwrap();
    assert_eq!(table.qualified_name(id), Some("depth".to_string()));
}

#[test]
fn qualified_name_feature_scoped() {
    let mut table = ParameterTable::new();
    let fid = make_feature_id();
    let id = table
        .add_feature_param("extrude1.distance", fid, 0.05, Unit::Meters)
        .unwrap();
    assert_eq!(
        table.qualified_name(id),
        Some("extrude1.distance".to_string())
    );
}

#[test]
fn qualified_name_invalid_id_returns_none() {
    let mut table = ParameterTable::new();
    // Add and then remove a global param so we have a ParameterId that was
    // once valid but is no longer in the table.  `qualified_name` should
    // return None for it.
    let id = table.add_global("ephemeral", 1.0, Unit::Ratio).unwrap();
    table.remove(id).unwrap();
    assert!(table.qualified_name(id).is_none());
}

// ── bounds ────────────────────────────────────────────────────────────────────

#[test]
fn set_and_get_bounds() {
    let mut table = ParameterTable::new();
    let id = table.add_global("angle", 0.0, Unit::Radians).unwrap();
    table.set_bounds(id, 0.0, std::f64::consts::PI);
    let param = table.get(id).unwrap();
    assert_eq!(param.bounds, Some((0.0, std::f64::consts::PI)));
}

#[test]
fn clear_bounds() {
    let mut table = ParameterTable::new();
    let id = table.add_global("scale", 1.0, Unit::Ratio).unwrap();
    table.set_bounds(id, 0.0, 10.0);
    table.clear_bounds(id);
    assert!(table.get(id).unwrap().bounds.is_none());
}

// ── all_params iterator ───────────────────────────────────────────────────────

#[test]
fn all_params_includes_globals_feature_params_and_builtins() {
    let mut table = ParameterTable::new();
    let fid = make_feature_id();
    table.add_global("g1", 1.0, Unit::Meters).unwrap();
    table
        .add_feature_param("feat.p1", fid, 2.0, Unit::Meters)
        .unwrap();
    // 3 built-ins + 1 global + 1 feature param = 5
    assert_eq!(table.all_params().count(), 5);
}
