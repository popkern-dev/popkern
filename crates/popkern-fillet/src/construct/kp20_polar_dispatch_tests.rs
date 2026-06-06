//! Tests for KP20 polar dispatch diagnostics.

use super::timing::Kp20DispatchStageTiming;

#[test]
fn dispatch_stage_timing_label_is_structured() {
    let timing = Kp20DispatchStageTiming {
        support_domain_ms: 1,
        support_transition_ms: 2,
        support_normal_lift_ms: 3,
        origin_inner_ms: 4,
        paired_inner_ms: 5,
        support_inner_ms: 6,
        support_reparameterization_ms: 7,
    };

    let detail = timing.label();

    for label in [
        "timing_total_ms=28",
        "timing_support_domain_ms=1",
        "timing_support_transition_ms=2",
        "timing_support_normal_lift_ms=3",
        "timing_origin_inner_ms=4",
        "timing_paired_inner_ms=5",
        "timing_support_inner_ms=6",
        "timing_support_reparameterization_ms=7",
    ] {
        assert!(detail.contains(label), "missing {label} in {detail}");
    }
}
