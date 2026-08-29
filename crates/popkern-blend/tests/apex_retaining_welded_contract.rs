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

/// The welded trihedral base corner (pyramid3, uniform radii) declines
/// the apex-retaining partition with the recorded class, serves the
/// vertex-blend fallback, and the served body carries no null-loop half-edges
/// and no structural validation errors.
#[test]
fn welded_n3_apex_retaining_declines_semantically() {
    let _g = popkern_blend::dispatch_overrides::V_FREE_REPEATED_APEX.set_guard(false);
    let result = build_under_direct_apex(true, &[0.10, 0.10, 0.10]);

    let reasons = direct_apex_fallback_reasons(&result);
    assert_eq!(
        reasons.len(),
        1,
        "welded N=3 apex-retaining must fall back exactly once, got {reasons:?}"
    );
    assert!(
        reasons[0].contains(
            "welded shared-corner repeated-apex partition declined on the apex-retaining contract"
        ),
        "decline must record the welded apex-retaining class, got: {}",
        reasons[0]
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

/// Every apex-retaining class that served before the welded-N=3 decline
/// landed keeps serving natively: welded N=4 reuses existing
/// target-cylinder paths (Euler delta zero) and the bridged/mixed
/// classes keep their balanced far-side complexes.
#[test]
fn non_welded_and_n4_apex_retaining_classes_keep_serving() {
    let _g = popkern_blend::dispatch_overrides::V_FREE_REPEATED_APEX.set_guard(false);
    for (name, radii) in [
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05][..]),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07][..]),
        ("pyramid3_mixed", &[0.08, 0.10, 0.12][..]),
    ] {
        let result = build_under_direct_apex(true, radii);
        let reasons = direct_apex_fallback_reasons(&result);
        assert!(
            reasons.is_empty(),
            "{name} must keep serving on the apex-retaining contract, fell back: {reasons:?}"
        );
    }
}
