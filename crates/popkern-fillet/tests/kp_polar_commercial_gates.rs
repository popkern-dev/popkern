//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.
//!
//! # Fixture matrix
//!
//! | Fixture             | Valence | Radii                                      | Status |
//! |---------------------|---------|--------------------------------------------|--------|
//! | pentagon5_uniform   | 5       | all 0.05                                       | `#[ignore]` |
//! | pentagon5_mixed     | 5       | [0.04, 0.06, 0.05, 0.07, 0.05]                 | polar commercial gates pass |
//! | pentagon5_variable  | 5       | [0.02, 0.04, 0.06, 0.04, 0.02] (1:3 ratio)     | `#[ignore]` I.22 |
//! | hexagon6_uniform    | 6       | all 0.05                                       | `#[ignore]` |
//! | hexagon6_mixed      | 6       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06]           | `#[ignore]` |
//! | hexagon6_variable   | 6       | [0.02, 0.04, 0.06, 0.06, 0.04, 0.02] (1:3)     | `#[ignore]` I.22 |
//! | octagon8_mixed      | 8       | [0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06] | `#[ignore]` |
//!
//! # Apex-angle regimes (important for fixture selection)
//!
//! The polar path is **not** the commercial-grade answer for every
//! N-valent corner — it's the answer specifically for corners where
//! the simpler paths (sphere-fit, polar-NURBS) fail. Empirical probe
//! (`kp_polar_pipeline_probe.rs`) shows three distinct regimes by
//! apex half-angle `θ = atan(base_radius / height)`:
//!
//! | θ     | Polar reachable? | Pipeline choice                  |
//! |-------|------------------|----------------------------------|
//! | ≤ 27° | No               | Sphere-fit / single-surface path |
//! | 27°–45° | Sometimes (depends on r) | Polar when polar-NURBS insufficient |
//! | > 45° | Yes (sharp apex) | Polar attempted; folds likely    |
//!
//! Blunt-apex fixtures (θ ≤ 27° — commercial-CAD typical) bypass
//! polar entirely and take sphere-fit. They're positive controls
//! for the sphere-fit path, not tests of polar. Tests below mark
//! this explicitly in their names (`_blunt` = expected sphere-fit;
//! `_sharp` or unqualified = polar-engaged).
//!
//! `pentagon5_mixed` is active after high-valence face-loop
//! normalization restored the N-edge topology prerequisite before
//! polar dispatch. Remaining ignored fixtures are either selected
//! non-polar paths or polar stress cases awaiting the same measured
//! promotion.
//!
//! The uniform-radii fixtures are deliberately included to catch
//! regressions where a sphere-fit shortcut or upstream dispatch
//! diverts the ring away from the polar rung; `assert_all_binonic`
//! fails loudly if a fixture stops producing measurable
//! `RectangularBiNonic` patches.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, collect_patches, fillet_corner_full,
    incident_edges, measure_ring, pick_vertex_with_valence,
};

fn run_n_gonal_pyramid_polar(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
    fixture: &str,
) -> helpers::RingMetrics {
    assert_eq!(radii.len(), n);
    let pyr = BRep::n_gonal_pyramid(n, base_radius, height, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).expect("n_gonal_pyramid apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner_full(&pyr, &edges, radii, true, true);
    let patches = collect_patches(&brep, &faces);
    assert_all_binonic(&patches, fixture);
    measure_ring(&patches, N_SAMPLES)
}

fn fillet_n_gonal_pyramid_result(
    n: usize,
    base_radius: f64,
    height: f64,
    radii: &[f64],
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

/// Blunt-apex positive control (apex half-angle 18°, typical CAD).
/// Expected to take sphere-fit path and bypass polar entirely.
/// The `assert_all_binonic` assertion in `run_n_gonal_pyramid_polar`
/// will fail at the moment because the sphere-fit path produces a
/// single non-Gregory face, not N Gregory sub-patches. This is the
/// intended commercial behaviour for gentle corners — the test is
/// `#[ignore]` until the probe is refactored to assert "valid
/// commercial output" (sphere-fit OR polar-bi9) rather than
/// specifically "N Gregory bi-nonic patches."
///
/// Documented here as a reminder: pentagon5 at h=3.0 (18° apex) is
/// the REAL commercial-CAD test case; the `_uniform` / `_mixed` /
/// `_variable` sharp-apex tests below stress the polar path
/// specifically.
#[ignore]
#[test]
fn pentagon5_blunt_takes_sphere_fit_path() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 3.0, &[0.05; 5], "pentagon5_blunt");
    assert_commercial_gates("pentagon5_blunt", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_pentagon5_uniform() {
    let m = run_n_gonal_pyramid_polar(5, 1.0, 0.8, &[0.05; 5], "pentagon5_uniform");
    assert_commercial_gates("pentagon5_uniform", &m);
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
    assert_commercial_gates("pentagon5_mixed", &m);
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

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_uniform() {
    let m = run_n_gonal_pyramid_polar(6, 1.0, 0.8, &[0.05; 6], "hexagon6_uniform");
    assert_commercial_gates("hexagon6_uniform", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_hexagon6_mixed() {
    let m = run_n_gonal_pyramid_polar(
        6,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        "hexagon6_mixed",
    );
    assert_commercial_gates("hexagon6_mixed", &m);
}

#[test]
fn hexagon6_mixed_records_structured_polar_promotion_decline() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    result
        .brep
        .validate()
        .expect("hexagon6_mixed fallback route must leave a valid BRep");
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp_polar_bi6.promotion.declined"
                        && cap.stage == "kp_polar_bi6.promotion.lm_rejected"
                        && cap.sector_name == "kp_polar_bi6_ring_n6"
                        && cap.detail_code.as_deref().is_some_and(|detail| {
                            detail.contains("L-M reason=full_deep_lm_no_fold_safe_iterate")
                        })
            )),
        "hexagon6_mixed should expose its polar promotion rejection as structured diagnostics"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp_polar_bi6.promotion.declined"
                        && cap.stage == "kp_polar_bi6.promotion.deep_fill_rejected"
                        && cap.sector_name == "kp_polar_bi6_ring_n6"
                        && cap.detail_code.as_deref() == Some("deep_fill_self_intersection")
            )),
        "hexagon6_mixed should classify boundary-owned deep-fill rejection structurally"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .all(|diagnostic| !matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.stage == "kp20_polar_macro.dispatch.piece_map_extrapolates"
            )),
        "bounded KP20 piece maps should no longer report triangle-domain extrapolation"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp20_polar_macro.dispatch.fact"
                        && cap.stage == "kp20_polar_macro.dispatch.piece_map_orientation"
                        && cap.sector_name == "kp20_polar_macro_ring_n6"
                        && cap.detail_code.as_deref().is_some_and(|detail| {
                            detail.contains("non_positive_orientation_maps=3")
                                && detail.contains("NearOriginRight:")
                                && detail.contains("NearMidpointRight:")
                                && detail.contains("NearMidpointLeft:")
                        })
            )),
        "hexagon6_mixed should expose KP20 piece-map orientation facts structurally"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp20_polar_macro.dispatch.fact"
                        && cap.stage == "kp20_polar_macro.dispatch.nonoverlap_coupled_lift"
                        && cap.sector_name == "kp20_polar_macro_ring_n6"
                        && cap.detail_code.as_deref().is_some_and(|detail| {
                            detail.contains("status=sampled")
                                && detail.contains("default_folds=22")
                                && detail.contains("folds=8")
                                && detail.contains("piece_folds=NearOriginLeft:2|NearOriginRight:2|NearMidpointRight:0|NearMidpointLeft:4")
                                && detail.contains("fold_zones=origin_corner:0|edge_shoulder:3|interior:5")
                                && detail.contains("fold_ownership=origin_pieces:4|support_pieces:4")
                                && detail.contains("dominant_fold_zone=interior")
                                && detail.contains("recommended_next_move=support_inner_interior_row")
                                && detail.contains("fold_site_count=8")
                                && detail.contains("fold_site_target_piece=NearMidpointLeft")
                                && detail.contains("fold_site_target_count=4")
                                && detail.contains("fold_site_worst=sector:3|piece:NearOriginRight|uv:(6.250e-2,6.250e-2)")
                                && detail.contains("fold_site_touching=NearOriginLeft:[0,3]|NearOriginRight:[1,3]|NearMidpointLeft:[1,2,4,5]")
                                && detail.contains("origin_inner_evaluated=12")
                                && detail.contains("origin_inner_accepted=1")
                                && detail.contains("origin_inner_rejected=11")
                                && detail.contains("origin_inner_start_folds=9")
                                && detail.contains("origin_inner_end_folds=8")
                                && detail.contains("origin_inner_mode=normal")
                                && detail.contains("origin_inner_scale=-2.000e-1")
                                && detail.contains("origin_inner_pre_evaluated=60")
                                && detail.contains("origin_inner_pre_accepted=0")
                                && detail.contains("origin_inner_pre_rejected=60")
                                && detail.contains("origin_inner_pre_mode=constrained_interior_solve_none|coupled_origin_support_none|origin_corner_none")
                                && detail.contains("paired_inner_evaluated=12")
                                && detail.contains("paired_inner_accepted=0")
                                && detail.contains("paired_inner_rejected=12")
                                && detail.contains("paired_inner_start_folds=8")
                                && detail.contains("paired_inner_end_folds=8")
                                && detail.contains("paired_inner_mode=none")
                                && detail.contains("paired_inner_scale=0.000e0")
                                && detail.contains("support_inner_evaluated=24")
                                && detail.contains("support_inner_accepted=0")
                                && detail.contains("support_inner_rejected=24")
                                && detail.contains("support_inner_start_folds=8")
                                && detail.contains("support_inner_end_folds=8")
                                && detail.contains("support_inner_mode=none")
                                && detail.contains("support_inner_scale=0.000e0")
                                && detail.contains("support_domain_evaluated=8")
                                && detail.contains("support_domain_accepted=1")
                                && detail.contains("support_domain_rejected=7")
                                && detail.contains("support_domain_start_folds=20")
                                && detail.contains("support_domain_end_folds=18")
                                && detail.contains("support_domain_target_piece=NearMidpointLeft")
                                && detail.contains("support_domain_decline_reason=none")
                                && detail.contains("support_domain_refine_evaluated=10")
                                && detail.contains("support_domain_refine_accepted=0")
                                && detail.contains("support_domain_refine_rejected=10")
                                && detail.contains("support_domain_refine_start_folds=18")
                                && detail.contains("support_domain_refine_end_folds=18")
                                && detail.contains("support_domain_refine_target_piece=NearMidpointLeft")
                                && detail.contains("support_domain_refine_decline_reason=support_domain_refine_no_fold_reduction")
                                && detail.contains("coupled_domain_evaluated=8")
                                && detail.contains("coupled_domain_accepted=0")
                                && detail.contains("coupled_domain_rejected=8")
                                && detail.contains("coupled_domain_start_folds=18")
                                && detail.contains("coupled_domain_end_folds=18")
                                && detail.contains("coupled_domain_origin_shift=0.000e0")
                                && detail.contains("coupled_domain_origin_scale=1.000e0")
                                && detail.contains("coupled_domain_decline_reason=coupled_domain_no_fold_reduction")
                                && detail.contains("guide_samples=90")
                                && detail.contains("guide_reversals=44")
                                && detail.contains("guide_min_normal_cos=-9.980e-1")
                                && detail.contains("guide_worst=sector:2|uv:(2.500e-1,2.500e-1)")
                                && detail.contains("guide_classification=guide_orientation_reversal")
                                && detail.contains("guide_repair_evaluated=4")
                                && detail.contains("guide_repair_accepted=0")
                                && detail.contains("guide_repair_rejected=4")
                                && detail.contains("guide_repair_start_reversals=44")
                                && detail.contains("guide_repair_end_reversals=44")
                                && detail.contains("guide_repair_blend=0.000e0")
                                && detail.contains("guide_repair_decline_reason=guide_repair_downstream_no_fold_reduction")
                                && detail.contains("split_domain_evaluated=10")
                                && detail.contains("split_domain_accepted=0")
                                && detail.contains("split_domain_rejected=10")
                                && detail.contains("split_domain_start_folds=18")
                                && detail.contains("split_domain_end_folds=18")
                                && detail.contains("split_domain_fraction=5.000e-1")
                                && detail.contains("split_domain_interface_pull=0.000e0")
                                && detail.contains("split_domain_decline_reason=split_domain_no_fold_reduction")
                                && detail.contains("preassembly_interior_evaluated=24")
                                && detail.contains("preassembly_interior_accepted=0")
                                && detail.contains("preassembly_interior_rejected=24")
                                && detail.contains("preassembly_interior_start_folds=18")
                                && detail.contains("preassembly_interior_end_folds=18")
                                && detail.contains("preassembly_interior_decline_reason=preassembly_warp_no_fold_reduction")
                                && detail.contains("support_transition_evaluated=8")
                                && detail.contains("support_transition_accepted=0")
                                && detail.contains("support_transition_rejected=8")
                                && detail.contains("support_transition_start_folds=18")
                                && detail.contains("support_transition_end_folds=18")
                                && detail.contains("support_transition_target_piece=NearMidpointLeft")
                                && detail.contains("support_transition_decline_reason=support_transition_")
                                && detail.contains("support_reparam_evaluated=8")
                                && detail.contains("support_reparam_fast_accepted=8")
                                && detail.contains("support_reparam_fast_rejected=0")
                                && detail.contains("support_reparam_lift_evaluated=1")
                                && detail.contains("support_reparam_accepted=0")
                                && detail.contains("support_reparam_rejected=8")
                                && detail.contains("support_reparam_start_folds=8")
                                && detail.contains("support_reparam_end_folds=8")
                                && detail.contains("support_reparam_target_piece=NearMidpointLeft")
                                && detail.contains("support_reparam_params=none")
                                && detail.contains("support_reparam_decline_reason=support_reparameterization_no_lift_improvement")
                                && detail.contains("tensor_border_residual=0.000e0")
                                && detail.contains("internal_c2_position=0.000e0")
                                && detail.contains("internal_c2_first=0.000e0")
                                && detail.contains("origin_g2_strip=0.000e0")
                                && detail.contains("optimizer_evaluated=")
                                && detail.contains("optimizer_accepted=5")
                                && detail.contains("accepted_steps=4")
                        })
            )),
        "hexagon6_mixed should expose the bounded nonoverlap-coupled lift candidate structurally"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp20_polar_macro.dispatch.declined"
                        && cap.stage == "kp20_polar_macro.dispatch.nonoverlap_coupled_lift_self_intersection"
                        && cap.sector_name == "kp20_polar_macro_ring_n6"
                        && cap.detail_code.as_deref().is_some_and(|detail| {
                            detail.contains("folds=8")
                                && detail.contains("piece_folds=NearOriginLeft:2|NearOriginRight:2|NearMidpointRight:0|NearMidpointLeft:4")
                                && detail.contains("origin_g2_strip=0.000e0")
                        })
            )),
        "hexagon6_mixed should decline the bounded lift candidate until it is fold-free"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|diagnostic| matches!(
                diagnostic,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                    if cap.code == "kp20_polar_macro.dispatch.declined"
                        && cap.stage == "kp20_polar_macro.dispatch.self_intersection"
                        && cap.sector_name == "kp20_polar_macro_ring_n6"
                        && cap.detail_code.as_deref().is_some_and(|detail| {
                            detail.contains("piece=") && detail.contains("jacobian_reversal=true")
                        })
            )),
        "hexagon6_mixed should expose the bounded KP20 candidate rejection structurally"
    );
}

#[ignore]
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
    assert_commercial_gates("pentagon5_variable", &m);
}

#[ignore]
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
    assert_commercial_gates("hexagon6_variable", &m);
}

#[ignore]
#[test]
fn polar_passes_commercial_gates_on_octagon8_mixed() {
    let m = run_n_gonal_pyramid_polar(
        8,
        1.0,
        0.8,
        &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
        "octagon8_mixed",
    );
    assert_commercial_gates("octagon8_mixed", &m);
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
