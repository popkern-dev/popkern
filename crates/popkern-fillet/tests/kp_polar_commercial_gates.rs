//! Phase I — commercial-grade continuity gates on the K-P 2015
//! bi-6 polar output (N ≥ 5-valent vertices).
//!
//! Sibling of `kp_binonic_commercial_gates.rs` (Phase H, N ∈ {3, 4});
//! this file exercises the polar rung behind
//! [`FilletOptions::use_kp_polar_gregory`]. Metric definitions and
//! gate thresholds (G0 ≤ 1e-10, G1 ≤ 5e-2, G2 ≤ 5e-2) come from
//! `kp_gate_helpers.rs`.
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

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, assert_all_binonic, assert_commercial_gates, collect_patches, fillet_corner_full,
    incident_edges, measure_ring, measure_shared_gregory_trims_with_sources,
    pick_vertex_with_valence,
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
fn hexagon6_mixed_accepts_support_cell_topology() {
    let result = fillet_n_gonal_pyramid_result(6, 1.0, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05, 0.06]);
    result
        .brep
        .validate()
        .expect("hexagon6_mixed support-cell route must leave a valid BRep");
    let support_piece_count =
        support_cell_piece_count(&result.corner_patches).expect("support-cell provenance");
    let support_cells_per_sector = support_piece_count / 6;
    assert_eq!(result.corner_patches.len(), support_piece_count);
    assert!(
        matches!(support_cells_per_sector, 4 | 8),
        "hexagon6_mixed support-cell route should accept a supported inventory, got piece_count={} cells_per_sector={}",
        support_piece_count,
        support_cells_per_sector
    );
    let support_cell_faces = result
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect::<Vec<_>>();
    let support_cell_sources = result
        .corner_patches
        .iter()
        .filter_map(|patch| {
            patch
                .sub_patch_provenance
                .map(|provenance| (patch.face, provenance.piece_index))
        })
        .collect::<Vec<_>>();
    let trim_metrics = measure_shared_gregory_trims_with_sources(
        &result.brep,
        &support_cell_faces,
        &support_cell_sources,
        N_SAMPLES,
    );
    assert!(
        trim_metrics.shared_trim_count >= 18,
        "hexagon6_mixed support-cell topology should expose internal shared Gregory trims, got {}",
        trim_metrics.shared_trim_count
    );
    assert!(
        trim_metrics.pcurve_out_of_unit_count == 0
            && trim_metrics.degenerate_pcurve_count == 0
            && trim_metrics.worst_abs_uv <= 1.0 + 1.0e-9,
        "hexagon6_mixed support-cell pcurves should stay in unit nondegenerate support domains; out_of_unit={} degenerate={} worst_abs_uv={:.3e} degenerate_source_pair={:?} degenerate_uv_edges={:?}",
        trim_metrics.pcurve_out_of_unit_count,
        trim_metrics.degenerate_pcurve_count,
        trim_metrics.worst_abs_uv,
        trim_metrics.degenerate_source_pair,
        trim_metrics.degenerate_uv_edges
    );
    assert!(
        trim_metrics.worst_g0_gap.is_finite()
            && trim_metrics.worst_normal_disagreement.is_finite()
            && trim_metrics.worst_source_pair.is_some(),
        "hexagon6_mixed support-cell shared-trim metrics should remain observable; g0={:.3e} normal={:.3e} worst_source_pair={:?}",
        trim_metrics.worst_g0_gap,
        trim_metrics.worst_normal_disagreement,
        trim_metrics.worst_source_pair
    );
    assert!(
        trim_metrics.worst_g0_gap <= 2.0e-3,
        "hexagon6_mixed support-cell shared-trim G0 should stay bounded after segment-aware reconciliation; g0={:.3e} normal={:.3e} best_boundary_gap={:.3e} max_row_trim_count={} worst_source_pair={:?} worst_uv_edges={:?}",
        trim_metrics.worst_g0_gap,
        trim_metrics.worst_normal_disagreement,
        trim_metrics.worst_pair_best_boundary_gap,
        trim_metrics.max_boundary_row_shared_trim_count,
        trim_metrics.worst_source_pair,
        trim_metrics.worst_uv_edges
    );
    assert!(
        trim_metrics.max_boundary_row_shared_trim_count >= 2,
        "hexagon6_mixed support-cell residual should expose multi-trim boundary-row ownership; max_row_trim_count={}",
        trim_metrics.max_boundary_row_shared_trim_count
    );
    assert!(result.corner_patches.iter().all(|patch| {
        patch.n_sides == 6
            && patch.sub_patch_provenance.is_some_and(|provenance| {
                provenance.sector_count == 6
                    && provenance.piece_count == support_piece_count
                    && matches!(
                        provenance.role,
                        popkern_fillet::CornerSubPatchRole::SupportCell
                    )
            })
    }));
    let diagnostics = &result.repair_report.apex_diagnostics;
    assert!(
        diagnostics.iter().all(|diagnostic| !matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "rational_gregory.dispatch.declined"
                    && cap.sector_name == "n6_corner_cap"
        )),
        "hexagon6_mixed support-cell route should not fall back through rational Gregory decline"
    );
    let support_cell_dispatch_details = support_cell_dispatch_details(diagnostics);
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "fillet.apply_blend.pipeline.fact"
                    && cap.stage == "fillet.apply_blend.pipeline.timing"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains(&format!("corner_patch_count={support_piece_count}"))
                            && detail.contains("support_seam_job_count=1")
                    })
        )),
        "hexagon6_mixed support-cell acceptance should remain observable through pipeline timing diagnostics"
    );
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp20_polar_macro.dispatch.fact"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail_contains_support_cell_inventory(
                            detail,
                            support_piece_count,
                            support_cells_per_sector,
                        )
                    })
        )),
        "hexagon6_mixed support-cell route should expose selected primitive inventory; piece_count={} cells_per_sector={} dispatch_details={:?}",
        support_piece_count,
        support_cells_per_sector,
        support_cell_dispatch_details
    );
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp20_polar_macro.dispatch.fact"
                    && cap.stage == "kp20_polar_macro.dispatch.support_cell_assembled_quality"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains("selected=promoted")
                            && detail.contains(&format!(
                                "selected_surface_count={support_piece_count}"
                            ))
                            && detail.contains(&format!(
                                "selected_cells_per_sector={support_cells_per_sector}"
                            ))
                            && detail.contains("validation_status=valid")
                    })
        )),
        "hexagon6_mixed support-cell route should emit post-apply assembled quality diagnostics"
    );
    assert!(
        support_cell_dispatch_details.iter().any(|detail| detail
            .contains("selected_candidate=pcurve_then_boundary")
            && (detail.contains("repair_status=accepted_non_regressing_final_repair")
                || detail.contains("repair_status=declined_to_preserve_final_quality"))),
        "hexagon6_mixed support-cell route should emit transactional shared-trim repair diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .any(|detail| detail.contains("candidate_status=accepted_selected")),
        "hexagon6_mixed support-cell selector should emit per-candidate score diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .any(|detail| detail.contains("candidate=pcurve_then_boundary_if_g0_safe")),
        "hexagon6_mixed support-cell selector should score the split boundary-safe candidate; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        support_cell_dispatch_details
            .iter()
            .any(|detail| detail.contains(
                "candidate=pcurve_then_boundary_if_g0_safe;candidate_status=accepted_not_selected"
            ) && detail.contains("operation_status=boundary_kept")),
        "hexagon6_mixed split candidate should explain its boundary substep outcome; support_cell_details={support_cell_dispatch_details:?}"
    );
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp20_polar_macro.dispatch.fact"
                    && cap.stage == "kp20_polar_macro.dispatch.support_cell_shared_trim_quality"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains("selected=promoted")
                            && detail.contains("measurement_scope=final_post_pipeline")
                            && detail
                                .contains("quality_status=current_budget_pass_commercial_g0_exceeded")
                            && detail.contains("current_g0_budget=2.000000000000000e-3")
                            && detail.contains("commercial_g0_target=1.000000000000000e-10")
                            && detail.contains(&format!(
                                "shared_trim_count={}",
                                trim_metrics.shared_trim_count
                            ))
                            && detail.contains("pcurve_out_of_unit_count=0")
                            && detail.contains("degenerate_pcurve_count=0")
                            && detail.contains(&format!(
                                "max_boundary_row_shared_trim_count={}",
                                trim_metrics.max_boundary_row_shared_trim_count
                            ))
                    })
        )),
        "hexagon6_mixed support-cell route should promote shared-trim quality into structured diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
}

fn support_cell_piece_count(patches: &[popkern_fillet::CornerPatchInfo]) -> Option<usize> {
    patches.iter().try_fold(None, |current, patch| {
        let provenance = patch.sub_patch_provenance?;
        if !matches!(
            provenance.role,
            popkern_fillet::CornerSubPatchRole::SupportCell
        ) {
            return None;
        }
        match current {
            Some(piece_count) if piece_count != provenance.piece_count => None,
            Some(piece_count) => Some(Some(piece_count)),
            None => Some(Some(provenance.piece_count)),
        }
    })?
}

fn detail_contains_support_cell_inventory(
    detail: &str,
    piece_count: usize,
    cells_per_sector: usize,
) -> bool {
    detail.contains(&format!("support_cell_primitive_patch_count={piece_count}"))
        && detail.contains("support_cell_primitive_selection_policy=")
        && match cells_per_sector {
            4 => support_cell_grid_detail_matches(detail, &[(1, 4), (4, 1), (2, 2)]),
            8 => support_cell_grid_detail_matches(detail, &[(2, 4), (4, 2)]),
            _ => false,
        }
}

fn support_cell_grid_detail_matches(detail: &str, grids: &[(usize, usize)]) -> bool {
    grids.iter().any(|&(cells_u, cells_v)| {
        detail.contains(&format!("support_cell_primitive_cells_u={cells_u}"))
            && detail.contains(&format!("support_cell_primitive_cells_v={cells_v}"))
    })
}

fn support_cell_dispatch_details(
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
) -> Vec<String> {
    diagnostics
        .iter()
        .filter_map(|diagnostic| {
            let popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic
            else {
                return None;
            };
            (cap.code == "kp20_polar_macro.dispatch.fact" && cap.stage.contains("support_cell"))
                .then(|| cap.detail_code.clone())
                .flatten()
        })
        .collect()
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
