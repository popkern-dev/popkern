// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//
// Integration tests for config-aware replay: literal overrides, suppression
// overrides, and default-behavior preservation when no config is active.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::config::ParamOverride;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::unit::Unit;
use crate::{Configuration, FeatureResult, FeatureTree};

// ---------------------------------------------------------------------------
// Helpers (duplicated from replay_integration_tests to keep files independent)
// ---------------------------------------------------------------------------

fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Test 1: Replay with a literal config override uses the overridden value
// ---------------------------------------------------------------------------

#[test]
fn replay_with_literal_override_uses_config_value() {
    let mut tree = FeatureTree::new();

    // Add a "dist" parameter with default 1.0 m.
    let dist_id = tree
        .parameters_mut()
        .add_global("dist_cfg", 1.0, Unit::Meters)
        .expect("add dist_cfg");

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch",
    );
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_id },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );

    // Create a config that overrides dist_cfg → 0.5 m.
    let mut cfg = Configuration {
        name: "Half".into(),
        description: String::new(),
        parent: None,
        overrides: Default::default(),
        suppression: Default::default(),
    };
    cfg.overrides.insert(dist_id, ParamOverride::Literal(0.5));
    let config_id = tree.configs.insert(cfg);
    tree.active_config = Some(config_id);

    let result = tree.replay(Tolerance::default()).expect("replay");

    // Both features should succeed.
    for outcome in result.feature_results.values() {
        assert!(
            matches!(outcome, FeatureResult::Success { .. }),
            "all features should succeed, got {outcome:?}"
        );
    }

    // The resolved value in the table must reflect the config override.
    let resolved = tree
        .parameters()
        .get(dist_id)
        .map(|p| p.value)
        .unwrap_or(0.0);
    assert!(
        (resolved - 0.5).abs() < 1e-10,
        "config override should set dist_cfg to 0.5, got {resolved}"
    );
}

// ---------------------------------------------------------------------------
// Test 2: Replay with config suppression override skips the feature
// ---------------------------------------------------------------------------

#[test]
fn replay_with_config_suppression_skips_feature() {
    let mut tree = FeatureTree::new();

    let dist_id = tree
        .parameters_mut()
        .add_global("dist_sup", 1.0, Unit::Meters)
        .expect("add dist_sup");

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch",
    );
    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_id },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );

    // Config suppresses the extrude feature.
    let mut cfg = Configuration {
        name: "NoExtrude".into(),
        description: String::new(),
        parent: None,
        overrides: Default::default(),
        suppression: Default::default(),
    };
    cfg.suppression.insert(extrude_id, true);
    let config_id = tree.configs.insert(cfg);
    tree.active_config = Some(config_id);

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Suppressed)
        ),
        "config suppression should suppress the extrude, got {:?}",
        result.feature_results.get(&extrude_id)
    );

    // BRep should be a flat sketch (no extrusion).
    let brep = tree.current_brep().expect("brep");
    assert!(
        brep.faces.len() < 6,
        "suppressed extrude → fewer than 6 faces, got {}",
        brep.faces.len()
    );
}

// ---------------------------------------------------------------------------
// Test 3: Replay without a config uses defaults (regression guard)
// ---------------------------------------------------------------------------

#[test]
fn replay_without_config_uses_defaults() {
    let mut tree = FeatureTree::new();

    let dist_id = tree
        .parameters_mut()
        .add_global("dist_def", 1.0, Unit::Meters)
        .expect("add dist_def");

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch",
    );
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_id },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );

    // No active config.
    assert!(tree.active_config.is_none());

    let result = tree.replay(Tolerance::default()).expect("replay");
    for outcome in result.feature_results.values() {
        assert!(
            matches!(outcome, FeatureResult::Success { .. }),
            "default replay should succeed, got {outcome:?}"
        );
    }
    let brep = tree.current_brep().expect("brep");
    assert_eq!(brep.faces.len(), 6, "default replay → 6-face box");
}
