//! Tests for [`super::DirectConstructionError`].

use super::*;

#[test]
fn fallback_classification_excludes_self_check_and_cancelled() {
    assert!(
        !DirectConstructionError::SelfCheckFailedBug {
            measured: 0.5,
            required: 0.99,
        }
        .should_fall_back_to_sampled()
    );
    assert!(!DirectConstructionError::Cancelled.should_fall_back_to_sampled());
}

#[test]
fn fallback_classification_covers_condition_errors() {
    assert!(
        DirectConstructionError::BoundaryTangentDegenerate {
            spoke_idx: 0,
            t: 0.5,
            magnitude: 1e-12,
            tolerance: 1e-10,
        }
        .should_fall_back_to_sampled()
    );
    assert!(
        DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: 1e-10,
            tolerance: 1e-8,
            edge_hint: "edge_42".into(),
        }
        .should_fall_back_to_sampled()
    );
}

#[test]
fn fallback_classification_excludes_radius_rejections() {
    // Radius errors should NOT fall back silently — they indicate
    // user inputs that can't be filleted at all.
    assert!(
        !DirectConstructionError::RadiusZero {
            edge_id: "edge_7".into(),
        }
        .should_fall_back_to_sampled()
    );
    assert!(
        !DirectConstructionError::RadiusTooSmall {
            edge_id: "edge_7".into(),
            radius: 1e-7,
            min_allowed: 1e-5,
            char_len: 1.0,
        }
        .should_fall_back_to_sampled()
    );
}

#[test]
fn tag_is_stable_machine_readable() {
    assert_eq!(DirectConstructionError::Cancelled.tag(), "cancelled");
    assert_eq!(
        DirectConstructionError::SelfCheckFailedBug {
            measured: 0.5,
            required: 0.99,
        }
        .tag(),
        "self_check_failed_bug"
    );
    assert_eq!(
        DirectConstructionError::TooFewSubPatches { got: 2 }.tag(),
        "too_few_sub_patches"
    );
}

#[test]
fn display_messages_name_the_offending_context() {
    // Each error's Display should include the specific corner /
    // spoke / tolerance — grep-able from a log.
    let err = DirectConstructionError::BoundaryTangentDegenerate {
        spoke_idx: 42,
        t: 0.25,
        magnitude: 1e-15,
        tolerance: 1e-10,
    };
    let s = format!("{err}");
    assert!(s.contains("42"), "should name the spoke: {s}");
    assert!(
        s.contains("0.25") || s.contains("0.25"),
        "should name the t: {s}"
    );
}
