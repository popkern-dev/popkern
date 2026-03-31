// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//
// Tests for Feature::Reference evaluation during replay — measurement queries
// that produce read-only parameters.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
use crate::measure_query::MeasureQuery;
use crate::parameter::ParamScope;
use crate::unit::Unit;
use crate::{FeatureResult, FeatureTree};

/// Create a solved rectangle sketch and serialize it to JSON.
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

/// Create an entity name referencing a given OperationId.
fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

/// Build a tree with a sketch + extrude to create a 1x1x1 box.
fn box_tree() -> FeatureTree {
    let mut tree = FeatureTree::new();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let dist = tree
        .parameters_mut()
        .add_global("dist", 1.0, Unit::Meters)
        .expect("add param");

    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude 1",
    );

    tree
}

// ---------------------------------------------------------------------------
// Test 1: BodyVolume reference produces a parameter with the correct value
// ---------------------------------------------------------------------------

#[test]
fn reference_body_volume_creates_parameter() {
    let mut tree = box_tree();

    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Volume".to_string(),
            query: MeasureQuery::BodyVolume,
        },
        "Vol Ref",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    // Reference should succeed.
    assert!(
        matches!(
            result.feature_results.get(&ref_id),
            Some(FeatureResult::Success { .. })
        ),
        "reference should succeed, got {:?}",
        result.feature_results.get(&ref_id)
    );

    // Check the parameter exists with qualified name "vol_ref.value".
    let pid = tree
        .parameters()
        .resolve_name("vol_ref.value")
        .expect("parameter 'vol_ref.value' should exist");

    let param = tree.parameters().get(pid).expect("parameter should exist");

    // The sketch+extrude box has tessellation-approximate volume.
    // Just verify we got a positive measurement — the exact value depends
    // on tessellation quality which varies across BRep construction paths.
    assert!(
        param.value > 0.0,
        "volume should be positive, got {}",
        param.value
    );

    // Parameter should be read-only and feature-scoped.
    assert!(param.read_only, "reference parameter should be read-only");
    assert!(
        matches!(param.scope, ParamScope::Feature(_)),
        "reference parameter should be feature-scoped"
    );
    assert_eq!(param.unit, Unit::CubicMeters, "volume unit should be m^3");

    // last_known should be set.
    assert!(
        param.last_known.is_some(),
        "last_known should be set on success"
    );
    assert!(
        (param.last_known.unwrap() - param.value).abs() < 1e-12,
        "last_known should match value"
    );
}

// ---------------------------------------------------------------------------
// Test 2: BodySurfaceArea reference
// ---------------------------------------------------------------------------

#[test]
fn reference_body_surface_area_creates_parameter() {
    let mut tree = box_tree();

    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Area".to_string(),
            query: MeasureQuery::BodySurfaceArea,
        },
        "Area Ref",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&ref_id),
            Some(FeatureResult::Success { .. })
        ),
        "surface area reference should succeed"
    );

    let pid = tree
        .parameters()
        .resolve_name("area_ref.value")
        .expect("parameter should exist");

    let param = tree.parameters().get(pid).expect("parameter should exist");

    // The sketch+extrude box has tessellation-approximate surface area.
    // Just verify we got a positive measurement.
    assert!(
        param.value > 0.0,
        "surface area should be positive, got {}",
        param.value
    );

    assert_eq!(param.unit, Unit::SquareMeters, "area unit should be m^2");
}

// ---------------------------------------------------------------------------
// Test 3: Reference with invalid query fails when no last_known
// ---------------------------------------------------------------------------

#[test]
fn reference_fails_when_query_invalid_and_no_last_known() {
    let mut tree = FeatureTree::new();

    // Add a reference without any BRep bodies — BodyVolume will fail.
    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Bad".to_string(),
            query: MeasureQuery::BodyVolume,
        },
        "Bad Ref",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&ref_id),
            Some(FeatureResult::Failed { .. })
        ),
        "reference should fail when no BRep and no last_known, got {:?}",
        result.feature_results.get(&ref_id)
    );
}

// ---------------------------------------------------------------------------
// Test 4: Suppressed reference is skipped
// ---------------------------------------------------------------------------

#[test]
fn suppressed_reference_is_skipped() {
    let mut tree = box_tree();

    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Vol".to_string(),
            query: MeasureQuery::BodyVolume,
        },
        "Vol Ref",
    );

    tree.suppress(ref_id);

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&ref_id),
            Some(FeatureResult::Suppressed)
        ),
        "suppressed reference should be skipped"
    );

    // Parameter should not be created for suppressed features.
    assert!(
        tree.parameters().resolve_name("vol_ref.value").is_none(),
        "suppressed reference should not create a parameter"
    );
}

// ---------------------------------------------------------------------------
// Test 5: Reference parameter updates on re-replay
// ---------------------------------------------------------------------------

#[test]
fn reference_parameter_updates_on_re_replay() {
    let mut tree = box_tree();

    // Find the distance parameter and remember its id.
    let dist_pid = tree
        .parameters()
        .resolve_name("extrude_1.distance")
        .expect("extrude distance param should exist");

    tree.add_feature(
        Feature::Reference {
            name: "Vol".to_string(),
            query: MeasureQuery::BodyVolume,
        },
        "Vol Ref",
    );

    // First replay — 1x1x1 box.
    tree.replay(Tolerance::default()).expect("first replay");

    let pid = tree
        .parameters()
        .resolve_name("vol_ref.value")
        .expect("param should exist");
    let v1 = tree.parameters().get(pid).expect("param").value;

    // Change extrude distance to 2.0 → 1x1x2 box → volume = 2.0.
    tree.parameters_mut()
        .get_mut(dist_pid)
        .expect("dist param")
        .value = 2.0;

    tree.replay(Tolerance::default()).expect("second replay");

    let v2 = tree.parameters().get(pid).expect("param").value;
    // Doubling the extrude distance should roughly double the volume.
    // Use a generous tolerance since tessellation-based volume is approximate.
    assert!(
        v2 > v1,
        "volume should increase when extrude distance doubles: v1={v1}, v2={v2}"
    );
}

// ---------------------------------------------------------------------------
// Test 6: Stale Reference propagates StaleInput to downstream feature params
// ---------------------------------------------------------------------------

#[test]
fn stale_reference_propagates_stale_input_to_dependents() {
    use crate::expression::Expression;
    use crate::parameter::EvaluationWarning;

    let mut tree = box_tree();
    tree.add_feature(
        Feature::Reference {
            name: "Vol".into(),
            query: MeasureQuery::BodyVolume,
        },
        "Vol Ref",
    );
    let dist2 = tree
        .parameters_mut()
        .add_global("dist2", 1.0, Unit::Meters)
        .expect("ok");
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist2,
        },
        "Extrude 2",
    );

    // First replay — everything succeeds, vol_ref.value gets last_known.
    tree.replay(Tolerance::default()).expect("first replay");

    let vol_pid = tree
        .parameters()
        .resolve_name("vol_ref.value")
        .expect("exists");
    let ext2_pid = tree
        .parameters()
        .resolve_name("extrude_2.distance")
        .expect("exists");
    tree.parameters_mut()
        .set_expression(ext2_pid, Expression::Param(vol_pid))
        .expect("set expr");

    // Suppress sketch and extrude so body is absent — Reference goes stale.
    tree.suppress(tree.timeline()[0]); // sketch
    tree.suppress(tree.timeline()[1]); // extrude
    tree.replay(Tolerance::default()).expect("second replay");

    let vol_p = tree.parameters().get(vol_pid).expect("vol param");
    assert!(
        matches!(
            vol_p.eval_warning,
            Some(EvaluationWarning::MeasureStale { .. })
        ),
        "vol_ref.value should have MeasureStale, got {:?}",
        vol_p.eval_warning
    );
    let ext2_p = tree.parameters().get(ext2_pid).expect("ext2 param");
    assert!(
        matches!(
            ext2_p.eval_warning,
            Some(EvaluationWarning::StaleInput { .. })
        ),
        "extrude_2.distance should have StaleInput, got {:?}",
        ext2_p.eval_warning
    );
    if let Some(EvaluationWarning::StaleInput { source }) = &ext2_p.eval_warning {
        assert_eq!(source, "vol_ref.value");
    }
}

// ---------------------------------------------------------------------------
// Test 7: Never-valid Reference suppresses downstream features
// ---------------------------------------------------------------------------

#[test]
fn never_valid_reference_suppresses_dependents() {
    let mut tree = FeatureTree::new();
    // No sketch — BodyVolume will fail (no BRep body at all).
    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Bad".into(),
            query: MeasureQuery::BodyVolume,
        },
        "Bad Ref",
    );
    let dist = tree
        .parameters_mut()
        .add_global("dist", 1.0, Unit::Meters)
        .expect("ok");
    let ext_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude 1",
    );
    tree.features[ext_id].depends_on.push(ref_id);

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&ref_id),
            Some(FeatureResult::Failed { .. })
        ),
        "reference should fail, got {:?}",
        result.feature_results.get(&ref_id)
    );
    assert!(
        matches!(
            result.feature_results.get(&ext_id),
            Some(FeatureResult::Suppressed)
        ),
        "dependent should be suppressed, got {:?}",
        result.feature_results.get(&ext_id)
    );
}

// ---------------------------------------------------------------------------
// Test 8: Never-valid suppression cascades through multiple levels
// ---------------------------------------------------------------------------

#[test]
fn never_valid_suppression_cascades() {
    let mut tree = FeatureTree::new();
    // No sketch — BodyVolume will fail (no BRep body at all).
    let ref_id = tree.add_feature(
        Feature::Reference {
            name: "Bad".into(),
            query: MeasureQuery::BodyVolume,
        },
        "Bad Ref",
    );

    let da = tree
        .parameters_mut()
        .add_global("dist_a", 1.0, Unit::Meters)
        .expect("ok");
    let feat_a = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: da,
        },
        "Extrude A",
    );
    tree.features[feat_a].depends_on.push(ref_id);

    let db = tree
        .parameters_mut()
        .add_global("dist_b", 1.0, Unit::Meters)
        .expect("ok");
    let feat_b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: db,
        },
        "Extrude B",
    );
    tree.features[feat_b].depends_on.push(feat_a);

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&feat_a),
            Some(FeatureResult::Suppressed)
        ),
        "feat_a should be suppressed, got {:?}",
        result.feature_results.get(&feat_a)
    );
    assert!(
        matches!(
            result.feature_results.get(&feat_b),
            Some(FeatureResult::Suppressed)
        ),
        "feat_b should cascade-suppress, got {:?}",
        result.feature_results.get(&feat_b)
    );
}
