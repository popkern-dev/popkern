//! PK-240/PK-249 extreme radius-asymmetry ladder: the boundary-ring
//! serve through the chart-orientation-repaired exact-injection
//! binonic.
//!
//! History of the class: the straight-spoke KP binonic fold-declined
//! on these fixtures (the guided construction's regular
//! characteristic map vs the angularly clustered ring — see
//! `docs/architecture/lopsided-ring-balanced-spoke-fan.md`), and
//! PK-240's balanced-spoke retry served 1:20–1:50 while 1:100
//! declined every strength. PK-249 closed the class at the root:
//! the exact-injection candidate's fold is the severe end of a
//! chart-orientation defect (its guided rows are authored with no
//! orientation constraint), and the transactional orientation
//! repair — localized TFI unfold + seam/spoke polish under the
//! promotion's commercial gates — now unfolds the STRAIGHT-fan
//! exact candidate on the whole ladder. Every fixture serves the
//! binonic ring with a machine-zero outer boundary, no fold
//! decline, no balanced retry, and chart-honest committed patches
//! (see `docs/architecture/chart-orientation-honest-binonic-rings.md`).
//!
//! Split from `kp_binonic_commercial_gates.rs` for the 500-line
//! structural cap.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_fillet::{ApexConstructionRoute, ExactBoundaryCandidateOutcome};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
use helpers::{
    N_SAMPLES, assert_all_binonic, cap_faces_from_result, collect_patches, measure_ring,
};
use kp_cube_fixture::fillet_cube_result;

#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_20() {
    assert_extreme_cube_ring_serve(&[0.01, 0.10, 0.20], "cube_extreme_20");
}

#[test]
fn binonic_passes_commercial_gates_on_cube_extreme_30() {
    assert_extreme_cube_ring_serve(&[0.0067, 0.10, 0.20], "cube_extreme_30");
}

#[test]
fn cube_extreme_40_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_serve(&[0.005, 0.10, 0.20], "cube_extreme_40");
}

#[test]
fn cube_extreme_50_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_serve(&[0.004, 0.10, 0.20], "cube_extreme_50");
}

/// 1:100 — historically declined at EVERY balanced-retry strength
/// (the tiny 0.0028 arc against 0.24 sides; the guide boundary sat
/// 1.323e-1 off the rim and the relock accepted zero improvement
/// steps). PK-249's orientation repair unfolds the exact candidate,
/// so the fixture serves the same contract as the rest of the
/// ladder.
#[test]
fn cube_extreme_100_serves_boundary_ring_instead_of_fallback_cap() {
    assert_extreme_cube_ring_serve(&[0.002, 0.10, 0.20], "cube_extreme_100");
}

/// Chart handedness-reversal bands over the committed cap patches —
/// a test-side copy of the serving metric (contiguous opposed runs
/// on the four 0.02/0.98 rails against the mid-face normal), so the
/// pin measures committed geometry independent of the kernel's
/// implementation.
fn chart_reversal_band_count(patches: &[GregoryPatch]) -> usize {
    let mut bands = 0usize;
    for patch in patches {
        let GregoryPatch::RectangularBiNonic(binonic) = patch else {
            continue;
        };
        let normal = |u: f64, v: f64| binonic.partial_u(u, v).cross(&binonic.partial_v(u, v));
        let reference = normal(0.5, 0.5);
        for rail in 0..4 {
            let mut previous_opposed = false;
            for k in 0..=24 {
                let t = 0.02 + 0.96 * (k as f64) / 24.0;
                let (u, v) = match rail {
                    0 => (0.02, t),
                    1 => (0.98, t),
                    2 => (t, 0.02),
                    _ => (t, 0.98),
                };
                let opposed = normal(u, v).dot(&reference) < 0.0;
                if opposed && !previous_opposed {
                    bands += 1;
                }
                previous_opposed = opposed;
            }
        }
    }
    bands
}

/// PK-249 ladder contract for extreme radius asymmetry: the ring
/// applied-fact is recorded; the promotion's boundary-candidates
/// fact shows the EXACT-injection candidate serving through the
/// orientation repair (`provenance=exact`, `repaired`, boundary
/// residual at machine zero); no fold decline and no balanced-fan
/// retry are recorded (the straight fan serves); validate_all is
/// CLEAN; the committed 3-patch ring is all-binonic, chart-honest
/// (zero reversal bands), and meets ring G0 at the machine-zero
/// gate with internal-spoke G1 and ring G2 at the commercial 5e-2
/// gates.
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
    let facts: Vec<(&str, &str)> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(f) => {
                Some((f.stage.as_str(), f.detail_code.as_deref().unwrap_or("")))
            }
            _ => None,
        })
        .collect();
    assert!(
        facts
            .iter()
            .any(|(stage, _)| *stage == "boundary_ring.applied"),
        "{fixture} must serve through the boundary ring carve; facts: {facts:?}"
    );
    let selection = facts
        .iter()
        .find(|(stage, _)| *stage == "kp_binonic.promotion.boundary_candidates")
        .unwrap_or_else(|| {
            panic!("{fixture} must record the boundary-candidates selection; facts: {facts:?}")
        });
    assert!(
        selection.1.contains("provenance=exact"),
        "{fixture} must serve the exact-injection candidate, got {}",
        selection.1
    );
    assert!(
        selection.1.contains("repaired") && selection.1.contains("accepted("),
        "{fixture} must serve through the chart-orientation repair, got {}",
        selection.1
    );
    assert!(
        !facts
            .iter()
            .any(|(stage, _)| *stage == "kp_binonic.promotion.self_intersection"),
        "{fixture} must not fold-decline anymore (PK-249); facts: {facts:?}"
    );
    assert!(
        !facts
            .iter()
            .any(|(stage, _)| stage.starts_with("three_quad.balanced_spoke_fan")),
        "{fixture} must serve the straight fan without the balanced retry; facts: {facts:?}"
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
    let bands = chart_reversal_band_count(&patches);
    assert_eq!(
        bands, 0,
        "{fixture} committed patches must be chart-honest (PK-249)"
    );
    let m = measure_ring(&patches, N_SAMPLES);
    eprintln!(
        "repaired-exact {fixture}: g0={:.3e} g1={:.3e} g2={:.3e} bands={bands}",
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
        m.worst_g2_residual <= 5.0e-2,
        "{fixture} ring G2 {:.3e} exceeds the commercial gate 5e-2",
        m.worst_g2_residual
    );
}
