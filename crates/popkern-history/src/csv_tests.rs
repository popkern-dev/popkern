//! Tests for CSV import and export of configurations.

use crate::FeatureTree;
use crate::config::ParamOverride;
use crate::expression::Expression;
use crate::feature::{Feature, SketchPlaneDesc};
use crate::node::FeatureNode;
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Add a dummy feature node with the given label to the tree's feature map and
/// timeline, returning its FeatureId.  Used for suppression-column tests.
fn add_dummy_feature(tree: &mut FeatureTree, label: &str) -> crate::node::FeatureId {
    use popkern_core::naming::OperationId;

    let node = FeatureNode {
        feature: Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(String::new()),
        },
        operation_id: OperationId::default(),
        label: label.to_string(),
        depends_on: vec![],
        suppressed: false,
        body: crate::node::BodyId(0),
        suppress_if: None,
        instance_op_ids: std::collections::HashMap::new(),
        last_result: None,
        selections: Vec::new(),
        provenance: Vec::new(),
    };
    let fid = tree.features.insert(node);
    tree.timeline.push(fid);
    fid
}

// ---------------------------------------------------------------------------
// Export tests
// ---------------------------------------------------------------------------

#[test]
fn export_empty_tree_headers_and_default_row() {
    let tree = FeatureTree::new();
    let csv = tree.export_configs_csv();
    let rows: Vec<&str> = csv.lines().collect();
    // Header + Default row only.
    assert_eq!(rows.len(), 2, "expected header + Default row, got:\n{csv}");
    assert!(
        rows[0].starts_with("name,description,parent"),
        "header: {}",
        rows[0]
    );
    assert!(rows[1].starts_with("Default"), "default row: {}", rows[1]);
}

#[test]
fn export_with_two_configs_three_rows() {
    let mut tree = FeatureTree::new();
    tree.add_config("Small", "Residential").unwrap();
    tree.add_config("Large", "Commercial").unwrap();
    let csv = tree.export_configs_csv();
    let rows: Vec<&str> = csv.lines().collect();
    // Header + Default + Small + Large = 4 rows.
    assert_eq!(rows.len(), 4, "rows:\n{csv}");
}

#[test]
fn export_literal_override_formatted_with_unit() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("width", 0.05, Unit::Millimeters)
        .unwrap();
    let cfg = tree.add_config("Small", "").unwrap();
    tree.set_config_override(cfg, pid, ParamOverride::Literal(0.025))
        .unwrap();
    let csv = tree.export_configs_csv();
    // Should contain "25mm" in the Small row.
    assert!(csv.contains("25mm"), "expected 25mm in:\n{csv}");
}

#[test]
fn export_expression_override_has_equals_prefix() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("width", 0.05, Unit::Millimeters)
        .unwrap();
    let cfg = tree.add_config("Derived", "").unwrap();
    let expr = Expression::Param(pid);
    tree.set_config_override(cfg, pid, ParamOverride::Expression(expr))
        .unwrap();
    let csv = tree.export_configs_csv();
    assert!(csv.contains("=width"), "expected =width in:\n{csv}");
}

#[test]
fn export_suppression_true_and_false() {
    let mut tree = FeatureTree::new();
    let fid = add_dummy_feature(&mut tree, "Chamfer1");
    let cfg = tree.add_config("NoFillet", "").unwrap();
    tree.set_config_suppression(cfg, fid, true);
    let csv = tree.export_configs_csv();
    assert!(
        csv.contains("Chamfer1.suppressed"),
        "suppression header missing:\n{csv}"
    );
    assert!(csv.contains("true"), "true missing:\n{csv}");
}

#[test]
fn export_suppression_false_value() {
    let mut tree = FeatureTree::new();
    let fid = add_dummy_feature(&mut tree, "Fillet1");
    let cfg = tree.add_config("WithFillet", "").unwrap();
    tree.set_config_suppression(cfg, fid, false);
    let csv = tree.export_configs_csv();
    // Row for WithFillet should contain "false" in the suppression column.
    let rows: Vec<&str> = csv.lines().collect();
    let config_row = rows.iter().find(|r| r.contains("WithFillet")).unwrap();
    assert!(config_row.contains("false"), "config row: {config_row}");
}

#[test]
fn export_parent_column_filled() {
    let mut tree = FeatureTree::new();
    let parent_id = tree.add_config("Base", "").unwrap();
    let _child = tree.add_derived_config("Child", "", parent_id).unwrap();
    let csv = tree.export_configs_csv();
    let rows: Vec<&str> = csv.lines().collect();
    let child_row = rows.iter().find(|r| r.contains("Child")).unwrap();
    assert!(
        child_row.contains("Base"),
        "child row should contain Base: {child_row}"
    );
}

#[test]
fn export_default_row_shows_base_param_value() {
    let mut tree = FeatureTree::new();
    // 50 mm stored as 0.05 m (SI).
    tree.parameter_table
        .add_global("height", 0.05, Unit::Millimeters)
        .unwrap();
    let csv = tree.export_configs_csv();
    let rows: Vec<&str> = csv.lines().collect();
    let default_row = rows[1];
    assert!(
        default_row.contains("50mm"),
        "expected 50mm in default row: {default_row}"
    );
}

// ---------------------------------------------------------------------------
// Round-trip test
// ---------------------------------------------------------------------------

#[test]
fn round_trip_export_import() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("radius", 0.01, Unit::Millimeters)
        .unwrap();
    let cfg1 = tree.add_config("Small", "Residential").unwrap();
    tree.set_config_override(cfg1, pid, ParamOverride::Literal(0.005))
        .unwrap();
    let cfg2 = tree.add_config("Large", "Commercial").unwrap();
    tree.set_config_override(cfg2, pid, ParamOverride::Literal(0.02))
        .unwrap();

    let csv = tree.export_configs_csv();

    let mut tree2 = FeatureTree::new();
    tree2
        .parameter_table
        .add_global("radius", 0.01, Unit::Millimeters)
        .unwrap();
    let result = tree2.import_configs_csv(&csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "unexpected warnings: {:?}",
        result.warnings
    );
    assert_eq!(result.created.len(), 2);

    // Verify overrides were imported.
    let configs2 = tree2.configs();
    let small_id = configs2
        .iter()
        .find(|&&id| tree2.config(id).map(|c| c.name == "Small").unwrap_or(false))
        .copied()
        .unwrap();
    let overrides = tree2.config_overrides(small_id).unwrap();
    let pid2 = tree2.parameter_table.resolve_name("radius").unwrap();
    assert!(overrides.contains_key(&pid2));
}

// ---------------------------------------------------------------------------
// Import-specific tests
// ---------------------------------------------------------------------------

#[test]
fn import_creates_new_config() {
    let mut tree = FeatureTree::new();
    let csv = "name,description,parent\nDefault,,\nFoo,desc,\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert_eq!(result.created.len(), 1);
    let configs = tree.configs();
    let foo = configs
        .iter()
        .find(|&&id| tree.config(id).map(|c| c.name == "Foo").unwrap_or(false));
    assert!(foo.is_some());
}

#[test]
fn import_updates_existing_config() {
    let mut tree = FeatureTree::new();
    let existing = tree.add_config("Existing", "old desc").unwrap();
    let csv = "name,description,parent\nExisting,new desc,\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert_eq!(result.updated, vec![existing]);
    assert_eq!(result.created.len(), 0);
}

#[test]
fn import_default_row_updates_base_param() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("gap", 0.01, Unit::Millimeters)
        .unwrap();
    let csv = "name,description,parent,gap\nDefault,,,5mm\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let param = tree.parameter_table.get(pid).unwrap();
    // 5mm in SI = 0.005 m
    assert!(
        (param.value - 0.005).abs() < 1e-12,
        "expected 0.005, got {}",
        param.value
    );
}

#[test]
fn import_expression_cell_parsed() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("width", 0.05, Unit::Millimeters)
        .unwrap();
    let csv = "name,description,parent,width\nDerived,,,=width\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let config_id = result.created[0];
    let overrides = tree.config_overrides(config_id).unwrap();
    match overrides.get(&pid) {
        Some(ParamOverride::Expression(_)) => {}
        other => panic!("expected Expression override, got {:?}", other),
    }
}

#[test]
fn import_literal_with_unit_converts_to_si() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("depth", 0.1, Unit::Millimeters)
        .unwrap();
    let csv = "name,description,parent,depth\nBig,,,200mm\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let cfg_id = result.created[0];
    let overrides = tree.config_overrides(cfg_id).unwrap();
    match overrides.get(&pid) {
        Some(ParamOverride::Literal(v)) => {
            // 200mm = 0.2 m
            assert!((*v - 0.2).abs() < 1e-12, "expected 0.2, got {v}");
        }
        other => panic!("expected Literal override, got {:?}", other),
    }
}

#[test]
fn import_empty_cell_means_no_override() {
    let mut tree = FeatureTree::new();
    tree.parameter_table
        .add_global("thickness", 0.003, Unit::Millimeters)
        .unwrap();
    let csv = "name,description,parent,thickness\nSame,,,\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let cfg_id = result.created[0];
    let overrides = tree.config_overrides(cfg_id).unwrap();
    assert!(overrides.is_empty(), "expected no overrides");
}

#[test]
fn import_suppression_column_true_false() {
    let mut tree = FeatureTree::new();
    let fid = add_dummy_feature(&mut tree, "Fillet1");
    let csv = "name,description,parent,Fillet1.suppressed\nNoFillet,,,true\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let cfg_id = result.created[0];
    let suppressions = tree.config_suppressions(cfg_id).unwrap();
    assert_eq!(suppressions.get(&fid), Some(&true));
}

#[test]
fn import_unknown_column_produces_warning() {
    let mut tree = FeatureTree::new();
    let csv = "name,description,parent,nonexistent_param\nFoo,,,42\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        !result.warnings.is_empty(),
        "expected a warning for unknown column"
    );
}

#[test]
fn import_bad_literal_cell_produces_warning() {
    let mut tree = FeatureTree::new();
    tree.parameter_table
        .add_global("width", 0.01, Unit::Millimeters)
        .unwrap();
    let csv = "name,description,parent,width\nFoo,,,notanumber\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        !result.warnings.is_empty(),
        "expected a warning for bad cell"
    );
    // Config should still be created.
    assert_eq!(result.created.len(), 1);
}

#[test]
fn import_with_parent_sets_derived_config() {
    let mut tree = FeatureTree::new();
    let csv = "name,description,parent\nBase,,\nChild,,Base\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    assert_eq!(result.created.len(), 2);

    let base_id = result
        .created
        .iter()
        .find(|&&id| tree.config(id).map(|c| c.name == "Base").unwrap_or(false))
        .copied()
        .unwrap();
    let child_id = result
        .created
        .iter()
        .find(|&&id| tree.config(id).map(|c| c.name == "Child").unwrap_or(false))
        .copied()
        .unwrap();

    let child_cfg = tree.config(child_id).unwrap();
    assert_eq!(child_cfg.parent, Some(base_id));
}

#[test]
fn import_quoted_cell_with_comma_handled() {
    let mut tree = FeatureTree::new();
    let pid = tree
        .parameter_table
        .add_global("width", 0.05, Unit::Millimeters)
        .unwrap();
    // Cell "=min(width, 10mm)" contains a comma, so it must be quoted.
    let csv = "name,description,parent,width\nQ,,,\"=min(width, 10mm)\"\n";
    let result = tree.import_configs_csv(csv).unwrap();
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    let cfg_id = result.created[0];
    let overrides = tree.config_overrides(cfg_id).unwrap();
    assert!(overrides.contains_key(&pid), "expected override for width");
}

#[test]
fn parse_csv_row_unquoted() {
    // Directly exercise the CSV row parser via the public export interface.
    let mut tree = FeatureTree::new();
    tree.parameter_table
        .add_global("a", 1.0, Unit::Ratio)
        .unwrap();
    tree.parameter_table
        .add_global("b", 2.0, Unit::Ratio)
        .unwrap();
    let cfg = tree.add_config("X", "desc").unwrap();
    let pid_b = tree.parameter_table.resolve_name("b").unwrap();
    tree.set_config_override(cfg, pid_b, ParamOverride::Literal(5.0))
        .unwrap();
    let csv = tree.export_configs_csv();
    // Import back and confirm.
    let mut tree2 = FeatureTree::new();
    tree2
        .parameter_table
        .add_global("a", 1.0, Unit::Ratio)
        .unwrap();
    tree2
        .parameter_table
        .add_global("b", 2.0, Unit::Ratio)
        .unwrap();
    let result = tree2.import_configs_csv(&csv).unwrap();
    assert_eq!(result.created.len(), 1);
}
