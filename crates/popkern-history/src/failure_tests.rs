//! Projection tests for the structured feature-failure payload.

use crate::failure::FeatureFailure;
use popkern_bool::BoolError;

#[test]
fn boolean_projection_carries_kind_and_merge_stage() {
    let e = BoolError::MergeFailed {
        stage: "stitch".into(),
        reason: "open fragment".into(),
    };
    let f = FeatureFailure::boolean("union", &e);
    match &f {
        FeatureFailure::BooleanFailed {
            operation,
            bool_kind,
            stage,
            detail,
        } => {
            assert_eq!(operation, "union");
            assert_eq!(bool_kind, "MergeFailed");
            assert_eq!(stage.as_deref(), Some("stitch"));
            assert_eq!(detail, &e.to_string());
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(f.kind(), "BooleanFailed");
    assert_eq!(
        f.to_string(),
        e.to_string(),
        "Display is the kernel's message"
    );
    let json = serde_json::to_value(&f).expect("serialise");
    assert_eq!(json["kind"], "BooleanFailed");
    assert_eq!(json["stage"], "stitch");
    let back: FeatureFailure = serde_json::from_value(json).expect("deserialise");
    assert_eq!(back, f);
}

#[test]
fn string_sites_project_to_other_and_keep_their_text() {
    let f: FeatureFailure = "parameter P not found".into();
    assert_eq!(f.kind(), "Other");
    assert!(f.contains("not found"));
    let e = crate::error::ReplayError::FeatureFailed {
        label: "x".into(),
        error: FeatureFailure::sweep("extrude", &"degenerate direction").into(),
    };
    let f = FeatureFailure::from(e);
    assert_eq!(f.kind(), "SweepFailed");
    let e = crate::error::ReplayError::HistoryError(crate::HistoryError::FeatureNotFound);
    assert_eq!(FeatureFailure::from(e).kind(), "Other");
}
