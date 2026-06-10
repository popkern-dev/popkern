#![allow(dead_code)]

#[path = "support/kp_gate_shared_trims.rs"]
mod kp_gate_shared_trims;
#[path = "support/kp_support_cell_candidate_diagnostics.rs"]
mod kp_support_cell_candidate_diagnostics;
#[path = "kp_polar_support_cell_repair_facts.rs"]
mod support_cell_repair_facts;

const SUPPORT_CELL_SHARED_TRIM_SAMPLES: usize = 100;

pub(crate) fn assert_support_cell_quality_gate(
    fixture: &str,
    result: &popkern_fillet::BlendResult,
    n_sides: usize,
) {
    result
        .brep
        .validate()
        .unwrap_or_else(|_| panic!("{fixture} support-cell route must leave a valid BRep"));
    let support_piece_count =
        support_cell_piece_count(&result.corner_patches).expect("support-cell provenance");
    let support_cells_per_sector = support_piece_count / n_sides;
    assert_eq!(result.corner_patches.len(), support_piece_count);
    assert!(
        matches!(support_cells_per_sector, 4 | 8 | 16),
        "{fixture} support-cell route should accept a supported inventory, got piece_count={} cells_per_sector={}",
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
    let trim_metrics = kp_gate_shared_trims::measure_shared_gregory_trims_with_sources(
        &result.brep,
        &support_cell_faces,
        &support_cell_sources,
        SUPPORT_CELL_SHARED_TRIM_SAMPLES,
    );
    assert!(
        trim_metrics.shared_trim_count >= n_sides * 3,
        "{fixture} support-cell topology should expose internal shared Gregory trims, got {}",
        trim_metrics.shared_trim_count
    );
    assert!(
        trim_metrics.pcurve_out_of_unit_count == 0
            && trim_metrics.degenerate_pcurve_count == 0
            && trim_metrics.worst_abs_uv <= 1.0 + 1.0e-9,
        "{fixture} support-cell pcurves should stay in unit nondegenerate support domains; out_of_unit={} degenerate={} worst_abs_uv={:.3e} degenerate_source_pair={:?} degenerate_uv_edges={:?}",
        trim_metrics.pcurve_out_of_unit_count,
        trim_metrics.degenerate_pcurve_count,
        trim_metrics.worst_abs_uv,
        trim_metrics.degenerate_source_pair,
        trim_metrics.degenerate_uv_edges
    );
    assert!(
        trim_metrics.worst_g0_gap.is_finite()
            && trim_metrics.worst_normal_disagreement.is_finite()
            && trim_metrics.worst_g2_residual.is_finite()
            && trim_metrics.worst_source_pair.is_some(),
        "{fixture} support-cell shared-trim metrics should remain observable; g0={:.3e} normal={:.3e} g2={:.3e} worst_source_pair={:?}",
        trim_metrics.worst_g0_gap,
        trim_metrics.worst_normal_disagreement,
        trim_metrics.worst_g2_residual,
        trim_metrics.worst_source_pair
    );
    assert!(
        trim_metrics.g2_unmeasured_sample_count == 0,
        "{fixture} support-cell shared-trim curvature should be measurable at every sample (fail closed); unmeasured={} g2={:.3e}",
        trim_metrics.g2_unmeasured_sample_count,
        trim_metrics.worst_g2_residual
    );
    assert!(
        trim_metrics.worst_g0_gap <= 1.0e-10,
        "{fixture} support-cell shared trims should meet the commercial G0 target by construction (chain-curve network: every row an exact restriction of its chain's canonical curve); g0={:.3e} normal={:.3e} best_boundary_gap={:.3e} max_row_trim_count={} worst_source_pair={:?} worst_uv_edges={:?}",
        trim_metrics.worst_g0_gap,
        trim_metrics.worst_normal_disagreement,
        trim_metrics.worst_pair_best_boundary_gap,
        trim_metrics.max_boundary_row_shared_trim_count,
        trim_metrics.worst_source_pair,
        trim_metrics.worst_uv_edges
    );
    assert!(
        trim_metrics.max_boundary_row_shared_trim_count >= 2,
        "{fixture} support-cell residual should expose multi-trim boundary-row ownership; max_row_trim_count={}",
        trim_metrics.max_boundary_row_shared_trim_count
    );
    assert_support_cell_provenance(fixture, result, n_sides, support_piece_count);
    let diagnostics = &result.repair_report.apex_diagnostics;
    assert_no_rational_gregory_fallback(fixture, diagnostics, n_sides);
    let support_cell_dispatch_details = support_cell_dispatch_details(diagnostics);
    assert_pipeline_timing_reports_support_cells(
        fixture,
        diagnostics,
        support_piece_count,
        &support_cell_dispatch_details,
    );
    assert_selected_support_cell_inventory(
        fixture,
        diagnostics,
        support_piece_count,
        support_cells_per_sector,
        &support_cell_dispatch_details,
    );
    assert_assembled_quality_reports_valid_selection(
        fixture,
        diagnostics,
        support_piece_count,
        support_cells_per_sector,
    );
    kp_support_cell_candidate_diagnostics::assert_candidate_diagnostics(
        fixture,
        &support_cell_dispatch_details,
    );
    support_cell_repair_facts::assert_selected_repair_facts_are_consistent(
        fixture,
        &support_cell_dispatch_details,
    );
    assert_trim_ownership_plan_diagnostic(
        fixture,
        diagnostics,
        &trim_metrics,
        &support_cell_dispatch_details,
    );
    assert_final_shared_trim_quality_diagnostic(
        fixture,
        diagnostics,
        &trim_metrics,
        &support_cell_dispatch_details,
    );
}

fn assert_support_cell_provenance(
    fixture: &str,
    result: &popkern_fillet::BlendResult,
    n_sides: usize,
    support_piece_count: usize,
) {
    assert!(
        result.corner_patches.iter().all(|patch| {
            patch.n_sides == n_sides
                && patch.sub_patch_provenance.is_some_and(|provenance| {
                    provenance.sector_count == n_sides
                        && provenance.piece_count == support_piece_count
                        && matches!(
                            provenance.role,
                            popkern_fillet::CornerSubPatchRole::SupportCell
                        )
                })
        }),
        "{fixture} support-cell patches should preserve support-cell provenance"
    );
}

pub(crate) fn support_cell_piece_count(
    patches: &[popkern_fillet::CornerPatchInfo],
) -> Option<usize> {
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

pub(crate) fn detail_contains_support_cell_inventory(
    detail: &str,
    piece_count: usize,
    cells_per_sector: usize,
) -> bool {
    detail.contains(&format!("support_cell_primitive_patch_count={piece_count}"))
        && detail.contains("support_cell_primitive_selection_policy=")
        && match cells_per_sector {
            4 => support_cell_grid_detail_matches(detail, &[(1, 4), (4, 1), (2, 2)]),
            8 => support_cell_grid_detail_matches(detail, &[(2, 4), (4, 2)]),
            16 => support_cell_grid_detail_matches(detail, &[(4, 4)]),
            _ => false,
        }
}

pub(crate) fn support_cell_dispatch_details(
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

fn assert_no_rational_gregory_fallback(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    n_sides: usize,
) {
    let sector_name = format!("n{n_sides}_corner_cap");
    assert!(
        diagnostics.iter().all(|diagnostic| !matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "rational_gregory.dispatch.declined"
                    && cap.sector_name == sector_name
        )),
        "{fixture} support-cell route should not fall back through rational Gregory decline"
    );
}

fn assert_pipeline_timing_reports_support_cells(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    support_piece_count: usize,
    support_cell_dispatch_details: &[String],
) {
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
        "{fixture} support-cell acceptance should remain observable through pipeline timing diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
}

fn assert_selected_support_cell_inventory(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    support_piece_count: usize,
    support_cells_per_sector: usize,
    support_cell_dispatch_details: &[String],
) {
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
        "{fixture} support-cell route should expose selected primitive inventory; piece_count={} cells_per_sector={} dispatch_details={:?}",
        support_piece_count,
        support_cells_per_sector,
        support_cell_dispatch_details
    );
}

fn assert_assembled_quality_reports_valid_selection(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    support_piece_count: usize,
    support_cells_per_sector: usize,
) {
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
        "{fixture} support-cell route should emit post-apply assembled quality diagnostics"
    );
}

fn assert_trim_ownership_plan_diagnostic(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    trim_metrics: &kp_gate_shared_trims::SharedTrimMetrics,
    support_cell_dispatch_details: &[String],
) {
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp20_polar_macro.dispatch.fact"
                    && cap.stage == "kp20_polar_macro.dispatch.support_cell_trim_ownership_plan"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains("selected=promoted")
                            && detail.contains("measurement_scope=final_post_pipeline")
                            && detail.contains("ownership_status=piecewise_segment_ownership_required")
                            && detail.contains(&format!(
                                "shared_trim_count={}",
                                trim_metrics.shared_trim_count
                            ))
                            && detail.contains("full_row_owned_trim_count=")
                            && detail.contains("segment_owned_trim_count=")
                            && detail.contains(&format!(
                                "max_boundary_row_shared_trim_count={}",
                                trim_metrics.max_boundary_row_shared_trim_count
                            ))
                            && detail.contains("piecewise_boundary_committed=already_segmented")
                    })
        )),
        "{fixture} support-cell route should expose trim ownership planning diagnostics before piecewise boundary mutation; support_cell_details={support_cell_dispatch_details:?}"
    );
}

fn assert_final_shared_trim_quality_diagnostic(
    fixture: &str,
    diagnostics: &[popkern_fillet::ApexConstructionDiagnostic],
    trim_metrics: &kp_gate_shared_trims::SharedTrimMetrics,
    support_cell_dispatch_details: &[String],
) {
    assert!(
        diagnostics.iter().any(|diagnostic| matches!(
            diagnostic,
            popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(cap)
                if cap.code == "kp20_polar_macro.dispatch.fact"
                    && cap.stage == "kp20_polar_macro.dispatch.support_cell_shared_trim_quality"
                    && cap.detail_code.as_deref().is_some_and(|detail| {
                        detail.contains("selected=promoted")
                            && detail.contains("measurement_scope=final_post_pipeline")
                            && detail.contains("piecewise_boundary_committed=already_segmented")
                            && detail.contains("quality_status=commercial_g1_exceeded")
                            && detail.contains("current_g0_budget=5.000000000000000e-4")
                            && detail.contains("commercial_g0_target=1.000000000000000e-10")
                            && detail.contains("commercial_g1_target=1.250000000000000e-3")
                            && detail.contains("commercial_g2_target=5.000000000000000e-2")
                            && detail.contains("worst_g2_residual=")
                            && detail.contains("g2_unmeasured_sample_count=0")
                            && detail.contains(&format!(
                                "shared_trim_count={}",
                                trim_metrics.shared_trim_count
                            ))
                            && detail.contains("partial_shared_trim_count=")
                            && detail.contains("pcurve_out_of_unit_count=0")
                            && detail.contains("degenerate_pcurve_count=0")
                            && detail.contains(&format!(
                                "max_boundary_row_shared_trim_count={}",
                                trim_metrics.max_boundary_row_shared_trim_count
                            ))
                            && detail.contains("worst_trim=present")
                            && detail.contains("worst_trim_span=")
                            && detail.contains("worst_face_a=")
                            && detail.contains("worst_row_a=")
                            && detail.contains("worst_face_b=")
                            && detail.contains("worst_row_b=")
                            && detail.contains("worst_sample_fraction=")
                    })
        )),
        "{fixture} support-cell route should promote shared-trim quality into structured diagnostics; support_cell_details={support_cell_dispatch_details:?}"
    );
}

fn support_cell_grid_detail_matches(detail: &str, grids: &[(usize, usize)]) -> bool {
    grids.iter().any(|&(cells_u, cells_v)| {
        detail.contains(&format!("support_cell_primitive_cells_u={cells_u}"))
            && detail.contains(&format!("support_cell_primitive_cells_v={cells_v}"))
    })
}
