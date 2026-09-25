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
//! exact candidate on the whole ladder (see
//! `docs/architecture/chart-orientation-honest-binonic-rings.md`).
//!
//! PK-539 re-measured the ladder with the orientation-aware spoke
//! normal sine (the promotion's coplanarity metric had passed 7°
//! creases, and the optimizer had reduced it by leaning cross tangents
//! onto the spokes). 1:20 and 1:30 then declined, because the
//! orientation repair's sine-driven polish converged into folds across
//! the spokes.
//!
//! PK-544 serves every rung:
//! - spoke folds are measured against the patches' reference
//!   orientation;
//! - the repair clears them;
//! - the polish minimizes the oriented chord, which is monotone through
//!   a fold.
//!
//! PK-547 made the served rings G1 over the whole spoke. The promotion
//! window `[0.05, 0.95]` had hidden a certified 0.33 crease at the rim
//! end, where the spoke's end leg left the rim tangent plane. The end
//! legs are now built in the rim and center tangent planes and held
//! there, and every ring measure samples the end zones. Each rung's cap
//! seams certify G1 over the whole edge.
//!
//! Split from `kp_binonic_commercial_gates.rs` for the 500-line
//! structural cap.

use popkern_blend::{ApexConstructionRoute, ExactBoundaryCandidateOutcome};
use popkern_core::geometry::gregory::GregoryPatch;

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "support/kp_cube_fixture.rs"]
mod kp_cube_fixture;
use helpers::{
    N_SAMPLES, assert_all_binonic, cap_faces_from_result, collect_patches, measure_ring,
};
use kp_cube_fixture::fillet_cube_result;

#[test]
fn cube_extreme_20_serves_a_commercial_g1_ring() {
    assert_extreme_cube_ring_serve(&[0.01, 0.10, 0.20], "cube_extreme_20");
}

#[test]
fn cube_extreme_30_serves_a_commercial_g1_ring() {
    assert_extreme_cube_ring_serve(&[0.0067, 0.10, 0.20], "cube_extreme_30");
}

#[test]
fn cube_extreme_40_serves_a_commercial_g1_ring() {
    assert_extreme_cube_ring_serve(&[0.005, 0.10, 0.20], "cube_extreme_40");
}

#[test]
fn cube_extreme_50_serves_a_commercial_g1_ring() {
    assert_extreme_cube_ring_serve(&[0.004, 0.10, 0.20], "cube_extreme_50");
}

/// 1:100 — historically declined at EVERY balanced-retry strength
/// (the tiny 0.0028 arc against 0.24 sides; the guide boundary sat
/// 1.323e-1 off the rim and the relock accepted zero improvement
/// steps). PK-249's orientation repair unfolds the exact candidate,
/// so the fixture serves the same contract as the rest of the
/// ladder.
#[test]
fn cube_extreme_100_serves_a_commercial_g1_ring() {
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

/// The promotion facts of `result` as `(stage, detail)`.
fn promotion_facts(result: &popkern_blend::BlendResult) -> Vec<(String, String)> {
    result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            popkern_blend::ApexConstructionDiagnostic::ExactBoundaryCap(f) => {
                Some((f.stage.clone(), f.detail_code.clone().unwrap_or_default()))
            }
            _ => None,
        })
        .collect()
}

/// Ladder contract (PK-539/544/547). The ring served for extreme radius
/// asymmetry is the exact-injection candidate through the orientation
/// repair; every cap seam certifies G1 over the whole edge; chart-honest (zero reversal bands), validates
/// clean, and meets ring G0 at the machine-zero gate and spoke G1 at
/// the commercial 5e-2 bar as the orientation-aware normal-angle sine
/// measures it (the former coplanarity metric read 7e-2..1.3e-1
/// creases as passing). G2 is not claimed and is pinned as measured
/// (1.83..1.96). The corner claims G1 from the committed ring's
/// measurement.
fn assert_extreme_cube_ring_serve(radii: &[f64], fixture: &str) {
    let result = fillet_cube_result(radii, true);
    let route = result
        .repair_report
        .apex_routes
        .first()
        .unwrap_or_else(|| panic!("{fixture} should record an apex construction route"));
    assert_eq!(
        route.route,
        ApexConstructionRoute::VertexBlend,
        "{fixture} rings inside the vertex-blend route"
    );
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::NotApplicable,
        "{fixture} should not pretend direct-apex/exact-boundary was attempted"
    );
    let facts = promotion_facts(&result);
    let served = facts.iter().any(|(stage, detail)| {
        stage == "kp_binonic.promotion.boundary_candidates"
            && detail.contains("provenance=exact")
            && detail.contains("repaired")
            && detail.contains("accepted(")
    });
    assert!(
        served,
        "{fixture} must serve the repaired exact-injection candidate; facts: {facts:?}"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "{fixture} must validate clean, got {errors:?}"
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
    helpers::assert_cap_seams_certify_g1(&result.brep, &faces, fixture);
    let m = measure_ring(&patches, N_SAMPLES);
    eprintln!(
        "served {fixture}: g0={:.3e} g1={:.3e} g2={:.3e}",
        m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual
    );
    assert!(
        m.worst_g0_gap < helpers::G0_GATE,
        "{fixture} ring G0 {:.3e}",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_normal <= helpers::G1_GATE,
        "{fixture} spoke G1 {:.3e} exceeds the commercial bar",
        m.worst_g1_normal
    );
    assert!(
        m.worst_g2_residual <= MEASURED_RING_G2_PIN,
        "{fixture} ring G2 {:.3e} regressed past its measured pin",
        m.worst_g2_residual
    );
    for patch in &result.corner_patches {
        assert_eq!(
            patch.achieved_continuity,
            popkern_blend::ContinuityLevel::G1,
            "{fixture} claims the measured G1"
        );
    }
}

/// Measured maximum ring G2 on the serving rungs (1.96), pinned: the
/// ring does not claim G2.
const MEASURED_RING_G2_PIN: f64 = 2.0;
