//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates_live, collect_patches,
    fillet_corner_full, incident_edges, measure_live_cap_seams, pick_vertex_with_valence,
};

/// Fillet the N-gonal pyramid apex and measure the committed cap
/// across its LIVE seam topology (the setback fan ships 4 bi-9
/// pieces per sector; the legacy `measure_ring` frame does not
/// describe that layout — PK-087).
fn run_n_gonal_pyramid_polar(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    fixture: &str,
) -> helpers::LiveSeamMetrics {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner_full(&pyr, &edges, radii, true, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_live_cap_seams(&brep, &faces, N_SAMPLES)
}

fn fillet_n_gonal_pyramid_result(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
) -> popkern_fillet::BlendResult {
    fillet_n_gonal_pyramid_result_with_continuity(
        n,
        base_radius,
        height,
        radii,
        popkern_fillet::ContinuityLevel::G1,
    )
}

fn fillet_n_gonal_pyramid_result_with_continuity(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    continuity: popkern_fillet::ContinuityLevel,
) -> popkern_fillet::BlendResult {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let specs = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    let opts = FilletOptions {
        continuity,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None).unwrap()
}

/// Blunt-apex positive control (apex half-angle 18°, typical CAD):
/// pentagon5 at h=3.0 is the REAL commercial-CAD case; the sharp-apex
/// fixtures below stress the polar/setback path specifically. A
/// gentle corner may legitimately ship a sphere-fit cap (single
/// non-Gregory face) instead of a fan, so the assertion is "valid
/// output, and IF a multi-patch cap shipped its live seams pass the
/// commercial gates" rather than a specific patch census.
#[test]
fn pentagon5_blunt_ships_valid_commercial_cap() {
    let result = fillet_n_gonal_pyramid_result(5, 1.0, 3.0, &[0.05; 5]);
    result
        .brep
        .validate()
        .expect("pentagon5_blunt must leave a valid BRep");
    let faces: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.n_sides == 5)
        .map(|p| p.face)
        .collect();
    let m = measure_live_cap_seams(&result.brep, &faces, N_SAMPLES);
    if m.seam_count > 0 {
        // G0/G1 commercial gates only: this control requests the
        // default G1 fillet class, so the cap's seam G2 is recorded
        // but not gated (measured at re-pin 2026-07-15: g0=1.6e-13,
        // g1=2.5e-2, g2=7.0e-2 — the G2 jump across seams of a
        // G1-request blend on the wide flat blunt fan).
        println!(
            "live-seam pentagon5_blunt: seams={} g0={:.3e} g1={:.3e} g2={:.3e}",
            m.seam_count, m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual
        );
        assert!(
            m.worst_g0_gap < helpers::G0_GATE,
            "[pentagon5_blunt] live-seam G0 gap {:.3e} exceeds {:.1e} gate",
            m.worst_g0_gap,
            helpers::G0_GATE
        );
        assert!(
            m.worst_g1_normal < helpers::G1_GATE,
            "[pentagon5_blunt] live-seam G1 misalignment {:.3e} exceeds {:.1e} gate",
            m.worst_g1_normal,
            helpers::G1_GATE
        );
    } else {
        println!(
            "pentagon5_blunt: single-face cap ({} apex patches), no live seams",
            faces.len()
        );
    }
}
/// Blunt-apex G2-request contract (PK-122): the same 18°-half-angle
/// pentagon requesting G2. The wide flat fan's seam curvature jump
/// (measured 7.0e-2 on the G1-request control) exceeds the 5e-2 G2
/// gate, so the setback surgery's G2 postcondition must DECLINE the
/// fan rather than report an achieved G2 the geometry does not have;
/// the corner then ships through a fallback path. Whatever shipped,
/// any patch claiming achieved G2 must measure within the gate on its
/// live seams, and the output stays valid.
#[test]
fn pentagon5_blunt_g2_request_is_honest() {
    let result = fillet_n_gonal_pyramid_result_with_continuity(
        5,
        1.0,
        3.0,
        &[0.05; 5],
        popkern_fillet::ContinuityLevel::G2,
    );
    result
        .brep
        .validate()
        .expect("pentagon5_blunt under a G2 request must leave a valid BRep");
    let g2_claimed: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.n_sides == 5 && p.achieved_continuity == popkern_fillet::ContinuityLevel::G2)
        .map(|p| p.face)
        .collect();
    if g2_claimed.is_empty() {
        println!(
            "pentagon5_blunt G2 request: no patch claims achieved G2 (honest decline; {} apex patches shipped)",
            result
                .corner_patches
                .iter()
                .filter(|p| p.n_sides == 5)
                .count()
        );
        return;
    }
    let m = measure_live_cap_seams(&result.brep, &g2_claimed, N_SAMPLES);
    println!(
        "live-seam pentagon5_blunt(G2 request): seams={} g0={:.3e} g1={:.3e} g2={:.3e}",
        m.seam_count, m.worst_g0_gap, m.worst_g1_normal, m.worst_g2_residual
    );
    if m.seam_count > 0 {
        assert!(
            m.worst_g2_residual < helpers::G2_GATE,
            "[pentagon5_blunt/G2] a patch claiming achieved G2 must measure within the \
             {:.1e} live-seam G2 gate; got {:.3e}",
            helpers::G2_GATE,
            m.worst_g2_residual,
        );
    }
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates_live("pentagon5_uniform", &m);
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_mixed() {
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05],
        "pentagon5_mixed",
    );
    assert_commercial_gates_live("pentagon5_mixed", &m);
}

#[test]
fn pentagon5_mixed_no_longer_records_polar_dispatch_decline() {
    let result = fillet_n_gonal_pyramid_result(5, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    result
        .brep
        .validate()
        .expect("pentagon5_mixed polar route must leave a valid BRep");
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .all(|diagnostic| !matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "rational_gregory.dispatch.declined"
            )),
        "pentagon5_mixed should reach N>=5 Gregory/polar dispatch without fallback diagnostics"
    );
}

#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates_live("hexagon6_uniform", &m);
}

/// Re-pinned 2026-07-15 (PK-087): the support-cell topology this
/// test originally asserted was replaced by the setback vertex fan
/// (4 bi-9 pieces per sector); the commercial claim is now the live
/// seam gate, same as the sibling fixtures. (An interim "skew-trim
/// containment" theory pinned g1 at 1.2e-1 — FALSIFIED: the 8.0e-2
/// twist was an artifact of a projection-walk outlier-drop repair
/// that was itself reverted; the honest reading is g1 = 2.7e-4.)
#[test]
fn hexagon6_mixed_ships_commercial_setback_fan() {
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
    assert_commercial_gates_live("hexagon6_mixed", &m);
}

#[test]
fn polar_passes_commercial_gates_on_pentagon5_variable() {
    // Phase I.22 variable-radius fixture. Wider 1:3 ratio spread
    // (0.02..0.06) than `pentagon5_mixed` (0.04..0.07, ratio ~1.75)
    // stresses the per-sector radius threading. Core-level tests
    // (`variable_radius_*` in `kp15_bi6_polar_bridge_tests`)
    // already verify the construction is homogeneous + finite on
    // asymmetric radii; this fixture verifies end-to-end through
    // the fillet pipeline + commercial G1/G2 gates.
    let m = run_n_gonal_pyramid_polar(
        5,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.04, 0.02],
        "pentagon5_variable",
    );
    assert_commercial_gates_live("pentagon5_variable", &m);
}

#[test]
fn polar_passes_commercial_gates_on_hexagon6_variable() {
    // Phase I.22 variable-radius — hexagonal apex with 1:3 spread.
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
        "hexagon6_variable",
    );
    assert_commercial_gates_live("hexagon6_variable", &m);
}

#[test]
fn polar_passes_commercial_gates_on_octagon8_mixed() {
    let m = run_n_gonal_pyramid_polar(
        8,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
        "octagon8_mixed",
    );
    assert_commercial_gates_live("octagon8_mixed", &m);
}

// ── Dispatch smoke test (flag flip must not crash the pipeline) ──
//
// Flipping `use_kp_polar_gregory = true` must not panic `fillet_edges`
// end-to-end. This proves the flag's plumbing is wired, NOT that the
// polar engine actually runs on the fixture — with the upstream
// N ≥ 5 gap (see module header) the polar promote function is never
// invoked, so this test is compatible with that gap. When the gap
// closes, this test will still pass and the commercial-gate fixtures
// above will start being flip-able.

#[test]
fn polar_flag_flip_does_not_crash_pipeline_at_n5() {
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).unwrap();
    let Some(vertex) = pick_vertex_with_valence(&pyr, 5) else {
        return;
    };
    let edges = incident_edges(&pyr, vertex);
    let (_brep, _faces) =
        fillet_corner_full(&pyr, &edges, &[0.04, 0.06, 0.05, 0.07, 0.05], true, true);
}
