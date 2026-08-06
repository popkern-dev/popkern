//! PK-240 extreme radius-asymmetry ladder: the boundary-ring serve
//! with the fold-triggered balanced-spoke retry.
//!
//! The straight-spoke KP binonic fold-declines on these fixtures
//! (the guided construction's regular characteristic map vs the
//! angularly clustered ring — see
//! `docs/architecture/lopsided-ring-balanced-spoke-fan.md`); the
//! balanced retry then serves the binonic ring at the commercial G1
//! gate through 1:50, and records its per-strength decline reasons
//! at 1:100 (PK-247 owns the tail).
//!
//! Split from `kp_binonic_commercial_gates.rs` for the 500-line
//! structural cap.

use popkern_fillet::{ApexConstructionRoute, ExactBoundaryCandidateOutcome};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
use helpers::{
    N_SAMPLES, assert_all_binonic, cap_faces_from_result, collect_patches, measure_ring,
};
use kp_cube_fixture::fillet_cube_result;

// ── K-P/binonic dispatch stress fixtures ─────────────────────
//
// Wider-asymmetry fixtures are gated on whichever route dispatch
// actually selects: binonic continuity gates where the K-P ring
// ships (≤ 1:30), legacy-route ownership tests beyond that.

/// 1:20 — PK-240: the straight-spoke binonic still fold-declines
/// (recorded), and the fold-triggered balanced-spoke retry then
/// serves the binonic ring at the commercial G1 gate
/// (measured g1 7.59e-4 at strength 0.50, validate_all CLEAN).
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_20() {
    assert_extreme_cube_ring_serve(&[0.01, 0.10, 0.20], "cube_extreme_20");
}

/// 1:30 radius asymmetry — same balanced-retry serve class as 1:20.
#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_30() {
    assert_extreme_cube_ring_serve(&[0.0067, 0.10, 0.20], "cube_extreme_30");
}

/// 1:40+: previously LEGACY-owned with zero measurable patches,
/// then (PK-239) a measurable ring with a creased plain-Gregory cap
/// (g1 7.962e-2 / 1.952e-1 / 2.249e-1 at 1:40/50/100). The PK-240
/// balanced retry now serves the binonic ring at the commercial G1
/// gate across the whole tail (measured 2.66e-3 at 1:50,
/// strength 0.33).
#[test]
fn cube_extreme_40_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_serve(&[0.005, 0.10, 0.20], "cube_extreme_40");
}

#[test]
fn cube_extreme_50_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_serve(&[0.004, 0.10, 0.20], "cube_extreme_50");
}

/// 1:100 — the retry runs and DECLINES every strength (recorded
/// per-strength reasons); the straight ring serve stands at its
/// measured ceiling. The remaining tail class is PK-247.
#[test]
fn cube_extreme_100_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_retry_declined(&[0.002, 0.10, 0.20], "cube_extreme_100", 2.8e-1);
}

/// PK-240 ladder contract for extreme radius asymmetry: the ring
/// applied-fact is recorded; the STRAIGHT-spoke binonic fold-decline
/// is recorded (the chart-vs-geometry mismatch is real and stays
/// visible); the fold-triggered balanced-spoke retry then serves
/// the binonic ring — every committed patch `RectangularBiNonic`,
/// balanced-fan fact recorded, validate_all CLEAN, ring G0 at the
/// machine-zero gate and internal-spoke G1 at the commercial 5e-2
/// gate. G2 across the balanced ring is recorded (measured ~1.2e0)
/// but not yet gated — the balanced-fan G2 frontier is PK-247.
fn assert_extreme_cube_ring_serve(radii: &[f64], fixture: &str) {
    let result = fillet_cube_result(radii, true);
    let route = result
        .repair_report
        .apex_routes
        .first()
        .unwrap_or_else(|| panic!("{fixture} should record an apex construction route"));
    assert_eq!(
        route.route,
        ApexConstructionRoute::Legacy,
        "{fixture} rings inside the legacy route"
    );
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::NotApplicable,
        "{fixture} should not pretend direct-apex/exact-boundary was attempted"
    );
    let facts: Vec<&str> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(f) => {
                Some(f.stage.as_str())
            }
            _ => None,
        })
        .collect();
    assert!(
        facts.contains(&"boundary_ring.applied"),
        "{fixture} must serve through the boundary ring carve; facts: {facts:?}"
    );
    assert!(
        facts.contains(&"kp_binonic.promotion.self_intersection"),
        "{fixture} must record the straight-spoke binonic fold decline; facts: {facts:?}"
    );
    assert!(
        facts.contains(&"three_quad.balanced_spoke_fan"),
        "{fixture} must record the balanced-spoke retry provenance; facts: {facts:?}"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{fixture} ring serve must validate clean, got {errors:?}"
    );
    let faces = cap_faces_from_result(&result);
    let patches = collect_patches(&result.brep, &faces);
    assert_eq!(
        patches.len(),
        3,
        "{fixture} must ship a measurable 3-patch ring"
    );
    assert_all_binonic(&patches, fixture);
    let m = measure_ring(&patches, N_SAMPLES);
    eprintln!(
        "balanced-retry {fixture}: g0={:.3e} g1={:.3e} g2={:.3e}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.worst_g2_residual
    );
    assert!(
        m.worst_g0_gap < helpers::G0_GATE,
        "{fixture} ring G0 gap {:.3e} exceeds the machine-zero gate",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_coplanarity <= 5.0e-2,
        "{fixture} internal-spoke G1 {:.3e} exceeds the commercial gate 5e-2",
        m.worst_g1_coplanarity
    );
    assert!(
        m.worst_g2_residual <= 2.0,
        "{fixture} ring G2 {:.3e} exceeds the recorded balanced-fan ceiling 2.0 (PK-247)",
        m.worst_g2_residual
    );
}

/// Tail class: the balanced retry runs and declines every strength;
/// the straight ring serve stands, pinned at its measured ceiling
/// with the decline reasons recorded (PK-247 owns the tail).
fn assert_extreme_cube_ring_retry_declined(radii: &[f64], fixture: &str, g1_ceiling: f64) {
    let result = fillet_cube_result(radii, true);
    let facts: Vec<&str> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(f) => {
                Some(f.stage.as_str())
            }
            _ => None,
        })
        .collect();
    assert!(
        facts.contains(&"boundary_ring.applied"),
        "{fixture} must serve through the boundary ring carve; facts: {facts:?}"
    );
    assert!(
        facts.contains(&"kp_binonic.promotion.self_intersection"),
        "{fixture} must record the straight-spoke binonic fold decline; facts: {facts:?}"
    );
    assert!(
        facts.contains(&"three_quad.balanced_spoke_fan.declined"),
        "{fixture} must record the balanced-retry decline reasons; facts: {facts:?}"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{fixture} ring serve must validate clean, got {errors:?}"
    );
    let faces = cap_faces_from_result(&result);
    let patches = collect_patches(&result.brep, &faces);
    assert_eq!(
        patches.len(),
        3,
        "{fixture} must ship a measurable 3-patch ring"
    );
    let m = measure_ring(&patches, N_SAMPLES);
    eprintln!(
        "retry-declined {fixture}: g0={:.3e} g1={:.3e}",
        m.worst_g0_gap, m.worst_g1_coplanarity
    );
    assert!(
        m.worst_g0_gap < helpers::G0_GATE,
        "{fixture} ring G0 gap {:.3e} exceeds the machine-zero gate",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_coplanarity <= g1_ceiling,
        "{fixture} internal-spoke crease {:.3e} exceeds the recorded ceiling {g1_ceiling:.1e} \
         (PK-247 owns the balanced-retry tail)",
        m.worst_g1_coplanarity
    );
}
