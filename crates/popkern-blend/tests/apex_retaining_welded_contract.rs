//! PK-223: apex-retaining route contract on welded shared-corner apexes.
//!
//! The welded (uniform repeated-apex) partition is Euler-balanced only
//! under the V-free contract: on the apex-retaining contract the
//! per-fillet far-side trim machinery converges every trim onto one
//! shared support face (measured on the pyramid3 trihedral base corner)
//! and the retired scaffold carries unrepayable Euler debt. The route
//! must decline that class SEMANTICALLY — no null-loop half-edges
//! minted, no Production-role debt records, carrier absorption never
//! gated off by a broken candidate — while every class that served
//! before keeps serving.
//!
//! Root fixes pinned here (PK-223):
//! - connector-closed residual cycles home their connector twins on the
//!   residual loop and keep connector edges retired
//!   (`materialize_support_boundary_residual_cycles`), instead of
//!   leaving a Production edge with a null-loop side that broke
//!   validation for every downstream repair;
//! - the bigon weld declines fragments whose far twin is the fragment's
//!   own homed connector record (a closed residual home, not a
//!   duplicate-span sliver);
//! - the partition classifies the welded apex-retaining decline instead
//!   of reporting a bare Euler postcondition crash.

#[allow(dead_code)]
#[path = "support/direct_apex_validate_helpers.rs"]
mod direct_apex_validate_helpers;
#[path = "kp_gate_helpers.rs"]
mod helpers;

use direct_apex_validate_helpers::{
    breakdown, build_under_direct_apex, direct_apex_fallback_reasons,
};
use popkern_core::topology::LoopHandle;

/// Check-9 failures on the served body.
fn check9_failures(brep: &popkern_core::topology::BRep) -> usize {
    brep.validate_all()
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::CurveSurfaceDeviation { .. }
            )
        })
        .count()
}

/// The first-pass rejections the serving retry recorded, as text.
fn retry_details(result: &popkern_blend::BlendResult) -> Vec<String> {
    result
        .repair_report
        .serving_self_intersection_retries
        .iter()
        .map(|retry| retry.detail.clone())
        .collect()
}

/// The welded trihedral base corner (pyramid3, uniform radii) declines
/// the apex-retaining partition with the recorded class. The
/// vertex-blend fallback it then builds on direct-apex anchors is off
/// its faces (21 edges, worst 1.06e-1), so the served-body gate
/// rejects that first pass. The exclusion retry serves the corner on
/// vertex-blend anchors (PK-533). The retry record keeps the welded
/// decline class, and the served body is check-9 clean, has no
/// null-loop half-edges, and is structurally sound.
#[test]
fn welded_n3_apex_retaining_declines_semantically() {
    let _g = popkern_blend::dispatch_overrides::V_FREE_REPEATED_APEX.set_guard(false);
    let result = build_under_direct_apex(true, &[0.10, 0.10, 0.10]);
    assert_eq!(
        check9_failures(&result.brep),
        0,
        "the served body must hold check 9"
    );

    let details = retry_details(&result);
    assert_eq!(
        details.len(),
        1,
        "welded N=3 apex-retaining must be retried exactly once, got {details:?}"
    );
    assert!(
        details[0].contains(
            "welded shared-corner repeated-apex partition declined on the apex-retaining contract"
        ) && details[0].contains("off their faces"),
        "the retry must record the welded apex-retaining class and the gate, got: {}",
        details[0]
    );
    assert!(
        direct_apex_fallback_reasons(&result).is_empty(),
        "the retried corner never re-enters direct-apex"
    );

    let null_loop_half_edges = result
        .brep
        .half_edges
        .iter()
        .filter(|(_, data)| {
            data.loop_ == LoopHandle::default() || !result.brep.loops.contains_key(data.loop_)
        })
        .count();
    assert_eq!(
        null_loop_half_edges, 0,
        "served body must not carry null-loop or dangling-loop half-edges"
    );

    let metrics = breakdown(&result.brep.validate_all());
    assert_eq!(
        metrics.structural, 0,
        "served body must be structurally sound after the welded apex-retaining decline"
    );
}

/// The non-welded and N=4 apex-retaining classes once "kept serving"
/// on that contract, but their served bodies were off their faces:
/// up to 26 edges, worst 9.2e-1, with inverted support faces on the
/// mixed pyramids. The served-body gate now rejects every such first
/// pass. The exclusion retry serves each apex check-9 clean through
/// the vertex-blend dispatch, and the retry is recorded (PK-533).
/// The default V-free contract serves all three through direct-apex
/// Coons without a retry.
#[test]
fn non_welded_and_n4_apex_retaining_classes_keep_serving() {
    let _g = popkern_blend::dispatch_overrides::V_FREE_REPEATED_APEX.set_guard(false);
    for (name, radii) in [
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05][..]),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07][..]),
        ("pyramid3_mixed", &[0.08, 0.10, 0.12][..]),
    ] {
        let result = build_under_direct_apex(true, radii);
        assert_eq!(
            check9_failures(&result.brep),
            0,
            "{name}: the served body must hold check 9"
        );
        let details = retry_details(&result);
        assert!(
            details.len() == 1 && details[0].contains("off their faces"),
            "{name}: the apex-retaining first pass must be retried under the gate, got {details:?}"
        );
        assert_eq!(
            breakdown(&result.brep.validate_all()).structural,
            0,
            "{name}: served body must be structurally sound"
        );
    }
}
