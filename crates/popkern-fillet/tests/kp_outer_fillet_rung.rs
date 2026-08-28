//! End-to-end pin for the outer-fillet G1 rung (PK-275).
//!
//! The rung extends the bi-nonic phase-3 L-M objective with
//! cap-vs-fillet G1 samples along each cap sub-patch's outer
//! boundary. It is armed by `POPKERN_OUTER_FILLET_G1=1`
//! (`three_quad_kp_promote.rs`) and had ZERO coverage before PK-275
//! — no test or caller ever set the variable. It was also
//! unaffordable when it did engage: with the outer block active,
//! `BinonicEval::jacobian` refused and the L-M driver fell back to
//! full-DoF finite differences (one complete ring rebuild per DoF
//! per iteration). PK-275 gave the outer block an analytic Jacobian
//! (FD-agreement pinned in `kp_binonic_jacobian_fd_tests.rs`) and a
//! structured outcome fact (`kp_binonic.promotion.outer_fillet_rung`)
//! so arming it is observable instead of silently absorbed.
//!
//! Fixture: the PK-249 extreme-asymmetry cube corner — the ladder
//! class that actually routes through the bi-nonic promotion
//! (`kp_binonic_extreme_ladder.rs` proves the route via its
//! boundary-candidates facts). This file owns the env var for its
//! whole process (own integration binary), the standard pattern for
//! env-gated rungs.
//!
//! Pinned:
//! 1. The armed promotion records the rung outcome fact — the rung
//!    can never again be silently unreachable (blind-green guard).
//! 2. The armed serve stays `validate_all`-clean.
//! 3. Determinism fingerprint: two in-process runs produce
//!    bit-identical worst-case metrics.

use popkern_core::geometry::gregory::GregoryPatch;

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
use helpers::{cap_faces_from_result, collect_patches};
use kp_cube_fixture::fillet_cube_result;

fn rung_facts(result: &popkern_fillet::BlendResult) -> Vec<String> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(f)
                if f.stage == "kp_binonic.promotion.outer_fillet_rung" =>
            {
                Some(f.detail_code.clone().unwrap_or_default())
            }
            _ => None,
        })
        .collect()
}

fn worst_ring_g1(result: &popkern_fillet::BlendResult) -> f64 {
    let faces = cap_faces_from_result(result);
    let patches: Vec<GregoryPatch> = collect_patches(&result.brep, &faces);
    helpers::measure_ring(&patches, helpers::N_SAMPLES).worst_g1_coplanarity
}

#[test]
fn armed_outer_fillet_rung_is_observable_and_serves_clean() {
    // SAFETY: this test binary owns its process; no other test in
    // this file reads the environment concurrently.
    unsafe { std::env::set_var("POPKERN_OUTER_FILLET_G1", "1") };

    let radii = [0.01, 0.10, 0.20];
    let t0 = std::time::Instant::now();
    let r1 = fillet_cube_result(&radii, true);
    let first = t0.elapsed();

    let facts = rung_facts(&r1);
    println!("armed rung facts: {facts:?} (first-run {first:?})");
    assert!(
        !facts.is_empty(),
        "the armed promotion must record kp_binonic.promotion.outer_fillet_rung facts — \
         an empty list means the rung silently stopped being reachable"
    );

    let errors = r1.brep.validate_all();
    assert!(
        errors.is_empty(),
        "armed serve must be validate_all-clean, got {} issues: {:?}",
        errors.len(),
        errors.iter().take(3).collect::<Vec<_>>()
    );

    // Determinism fingerprint on the served ring's worst G1.
    let g1_1 = worst_ring_g1(&r1);
    let r2 = fillet_cube_result(&radii, true);
    let g1_2 = worst_ring_g1(&r2);
    println!("armed rung cube_extreme_20: worst g1 {g1_1:.3e}");
    assert_eq!(
        g1_1.to_bits(),
        g1_2.to_bits(),
        "armed rung serve must be deterministic ({g1_1:.17e} vs {g1_2:.17e})"
    );
}
