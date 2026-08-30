use super::*;
use crate::FeatureTree;
use crate::feature::{Feature, SketchPlaneDesc};
use crate::unit::Unit;
use popkern_sketch::ProfileShape;

fn analytic_circle(radius: ScalarBinding) -> AnalyticSketchDesc {
    AnalyticProfile {
        shape: ProfileShape::Circle { radius },
        center: [1.0.into(), 0.0.into()],
    }
}

/// The constraint arm's wire form IS the bare JSON string v3 documents
/// carry in `sketch_data`, so existing documents load unchanged.
#[test]
fn constraint_arm_is_the_bare_string_wire_form() {
    let payload = crate::SketchPayload::Constraint("{\"points\":[]}".into());
    let json = serde_json::to_value(&payload).unwrap();
    assert_eq!(json, serde_json::Value::String("{\"points\":[]}".into()));
    let loaded: SketchPayload = serde_json::from_str("\"v3-json\"").unwrap();
    assert_eq!(loaded, crate::SketchPayload::Constraint("v3-json".into()));
    assert_eq!(loaded.constraint_json(), Some("v3-json"));
    assert!(loaded.analytic().is_none());
}

#[test]
fn analytic_arm_round_trips_with_literal_and_parameter_scalars() {
    let mut tree = FeatureTree::new();
    let r = tree
        .parameters_mut()
        .add_global("r", 0.25, Unit::Meters)
        .unwrap();
    let payload = SketchPayload::Analytic(analytic_circle(r.into()));
    let json = serde_json::to_string(&payload).unwrap();
    // A bare number is a literal; the parameter id is a key object.
    assert!(json.contains("\"center\":[1.0,0.0]"));
    assert!(json.contains("\"idx\""));
    let back: SketchPayload = serde_json::from_str(&json).unwrap();
    assert_eq!(back, payload);
    assert_eq!(back.param_ids(), vec![("profile", r)]);
    assert!(back.constraint_json().is_none());
    assert_eq!(
        back.analytic().unwrap().center[0],
        ScalarBinding::Literal(1.0)
    );

    let literal = SketchPayload::Analytic(analytic_circle(0.25.into()));
    assert!(literal.param_ids().is_empty());
}

#[test]
fn bindings_resolve_against_the_table_and_missing_parameters_decline() {
    let mut tree = FeatureTree::new();
    let r = tree
        .parameters_mut()
        .add_global("r", 0.25, Unit::Meters)
        .unwrap();
    let resolved = resolve_profile(&analytic_circle(r.into()), tree.parameters()).unwrap();
    assert_eq!(resolved.shape, ProfileShape::Circle { radius: 0.25 });
    assert_eq!(resolved.center, [1.0, 0.0]);

    let other = FeatureTree::new();
    let err = resolve_profile(&analytic_circle(r.into()), other.parameters()).unwrap_err();
    assert!(matches!(err, ReplayError::FeatureFailed { ref label, .. } if label == "sketch"));
}

/// The feature-level entry: an analytic sketch's parameters are
/// reported like any other feature parameter.
#[test]
fn feature_param_ids_include_analytic_bindings() {
    let mut tree = FeatureTree::new();
    let r = tree
        .parameters_mut()
        .add_global("r", 0.25, Unit::Meters)
        .unwrap();
    let feature = Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: analytic_circle(r.into()).into(),
    };
    let ids = crate::feature_params::feature_param_ids(&feature);
    assert_eq!(ids, vec![("profile", r)]);
}
