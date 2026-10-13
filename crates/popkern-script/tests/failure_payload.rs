//! PK-388: a declined feature surfaces the kernel's structured payload
//! (spec §10.2/§10.3) — not a flattened string — through
//! `ScriptError::FeatureFailed`, and the payload serialises with the
//! spec's field names (`kind`, `blend_kind`, `reason_kind`,
//! `remediation`), which is the same JSON the PyO3 probe presents.

use popkern_history::FilletRecipeDesc;
use popkern_script::prelude::*;
use popkern_script::{BlendDeclineKind, BlendFailureReason, FeatureFailure, ScriptError};

fn unit_cube(m: &mut Model) -> Body {
    m.profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("sketch")
        .extrude(m, End::Blind(1.0.into()))
        .expect("extrude")
}

fn oversized_fillet() -> ScriptError {
    let mut m = Model::new();
    let cube = unit_cube(&mut m);
    cube.fillet(&mut m, edges(), 0.9)
        .expect_err("r = 0.9 on every unit-cube edge cannot construct")
}

#[test]
fn declined_fillet_exposes_per_edge_failures_with_remediation() {
    let err = oversized_fillet();
    let failure = err.feature_failure().expect("FeatureFailed class");
    assert_eq!(failure.kind(), "BlendDeclined");
    let blend = failure.as_blend().expect("blend payload");
    // Whole-decline hint points at the per-edge list; each edge carries
    // its own ReduceRadius with a concrete value.
    assert_eq!(blend.remediation.action, "Precheck", "{blend:?}");
    let BlendDeclineKind::PartialFailure {
        failed,
        total,
        failures,
    } = &blend.decline
    else {
        panic!("expected PartialFailure, got {:?}", blend.decline);
    };
    assert_eq!(*total, 12);
    assert_eq!(*failed, failures.len());
    assert!(*failed >= 1 && *failed <= 12, "{failed}");
    for f in failures {
        let BlendFailureReason::RadiusTooLarge {
            requested,
            max_radius,
            ..
        } = &f.reason
        else {
            panic!("expected RadiusTooLarge, got {:?}", f.reason);
        };
        assert_eq!(*requested, 0.9);
        let max = max_radius.expect("engine estimates a fitting radius");
        assert!(max > 0.0 && max < 0.9, "max_radius {max}");
        assert_eq!(f.remediation.action, "ReduceRadius");
        let suggested = f.remediation.value.expect("concrete value");
        assert!(suggested < 0.9 && suggested > 0.0, "{suggested}");
        assert!(
            f.edge_name.is_some(),
            "extruded-cube edges carry stable names; record lost it: {f:?}"
        );
    }
    // Display is the kernel's message, redundant with the fields.
    assert!(err.to_string().contains("of 12 edges failed"), "{err}");
}

#[test]
fn failure_payload_serialises_with_the_spec_field_names() {
    let err = oversized_fillet();
    let failure = err.feature_failure().expect("FeatureFailed class");
    let json = serde_json::to_value(failure).expect("serialise");
    assert_eq!(json["kind"], "BlendDeclined");
    assert_eq!(json["blend_kind"], "PartialFailure");
    assert_eq!(json["total"], 12);
    assert!(json["message"].as_str().is_some_and(|m| !m.is_empty()));
    for key in ["action", "parameter", "value", "rationale"] {
        assert!(json["remediation"].get(key).is_some(), "remediation.{key}");
    }
    let first = &json["failures"][0];
    assert_eq!(first["reason_kind"], "RadiusTooLarge");
    for key in [
        "edge",
        "edge_name",
        "requested",
        "max_radius",
        "location",
        "remediation",
    ] {
        assert!(first.get(key).is_some(), "failures[0].{key}");
    }
    assert_eq!(first["remediation"]["action"], "ReduceRadius");
    // Floats cross bit-exactly through JSON (spec §10.4).
    let back: FeatureFailure = serde_json::from_value(json).expect("deserialise");
    assert_eq!(&back, failure);
}

#[test]
fn migration_edit_failures_carry_the_same_payload() {
    let mut m = Model::new();
    let cube = unit_cube(&mut m);
    cube.fillet_with(&mut m, edges(), 0.1, FilletRecipeDesc::default(), "rounds")
        .expect("r = 0.1 serves");
    let rounds = m.feature("rounds").expect("fillet feature handle");
    let err = m
        .modify(&rounds, Edit::new().radius(0.9))
        .expect_err("r = 0.9 declines on replay");
    let failure = err.feature_failure().expect("FeatureFailed class");
    assert_eq!(failure.kind(), "BlendDeclined");
    assert!(matches!(
        failure.as_blend().map(|b| &b.decline),
        Some(BlendDeclineKind::PartialFailure { .. })
    ));
    // The rejected edit is rolled back: the original fillet still serves.
    assert!(m.replay().expect("replay").failures().next().is_none());
}

#[test]
fn a_declined_feature_recorded_again_declines_again() {
    // PK-805: the recorder reads only the EXECUTED steps' results
    // (`replay_appended`). Re-recording a rolled-back failure verbatim
    // repeats the step's content but not its feature id, and cached
    // records resume only for their own id, so the second attempt
    // executes and declines again rather than replaying a stale serve.
    let mut m = Model::new();
    let cube = unit_cube(&mut m);
    for attempt in 0..2 {
        let err = cube
            .fillet(&mut m, edges(), 0.9)
            .expect_err("r = 0.9 on every unit-cube edge cannot construct");
        assert!(
            err.feature_failure().is_some(),
            "attempt {attempt}: {err:?}"
        );
    }
}
