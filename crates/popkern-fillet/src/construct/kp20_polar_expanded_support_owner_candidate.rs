//! Support-owner candidate selection for expanded KP20 support topology.

use popkern_core::geometry::gregory::Kp20PieceId;

use super::super::kp20_polar_expanded_support_binding::expanded_support_commit_shape;
use super::super::kp20_polar_expanded_support_diagnostics::{
    Kp20ExpandedSupportFoldLocalization, localize_expanded_support_folds,
};
use super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::super::kp20_polar_support_surface_split::support_owner_piece;
use super::support_cell_primitive::{
    Kp20SupportCellPrimitiveReport, measure_support_cell_primitive,
};
use super::{
    Kp20ExpandedSupportMode, Kp20ExpandedSupportTopologyReport, Kp20ExpandedSupportTopologyTrial,
    emitted_candidate_patchwork, expanded_support_topology_trial, fragmented_macro_fold_count,
    retained_fold_count, support_lane_reauthor_label, support_sites_by_sector, worst_support_site,
};

pub(super) fn measure_support_owner_candidates(
    baseline: &Kp20SupportNormalLiftReport,
) -> Kp20ExpandedSupportTopologyReport {
    let mut reports = Vec::new();
    for target_piece in support_owner_candidate_order() {
        reports.push(measure_support_owner_candidate(baseline, target_piece));
    }
    let best = reports
        .iter()
        .cloned()
        .reduce(|current, candidate| {
            if report_is_better(&candidate, &current) {
                candidate
            } else {
                current
            }
        })
        .unwrap_or_else(|| {
            initial_report(
                baseline.summary.fold_count,
                Kp20PieceId::NearMidpointLeft,
                None,
            )
        });
    if let Some(combined) = combined_owner_support_cell_report(&reports, best.clone()) {
        return combined;
    }
    reports
        .into_iter()
        .find(preserves_local_support_cell_floor)
        .unwrap_or(best)
}

fn combined_owner_support_cell_report(
    reports: &[Kp20ExpandedSupportTopologyReport],
    mut base: Kp20ExpandedSupportTopologyReport,
) -> Option<Kp20ExpandedSupportTopologyReport> {
    let left = reports
        .iter()
        .find(|report| report.target_piece == Kp20PieceId::NearMidpointLeft)?;
    let right = reports
        .iter()
        .find(|report| report.target_piece == Kp20PieceId::NearMidpointRight)?;
    let left_cells = &left.support_cell_primitive;
    let right_cells = &right.support_cell_primitive;
    if !support_cell_primitive_is_combinable(left_cells)
        || !support_cell_primitive_is_combinable(right_cells)
    {
        return None;
    }
    let left_cells_per_target = left_cells.cells_u * left_cells.cells_v;
    let right_cells_per_target = right_cells.cells_u * right_cells.cells_v;
    if left_cells_per_target == 0 || left_cells_per_target != right_cells_per_target {
        return None;
    }
    if support_cell_targets_overlap(left_cells, right_cells) {
        return None;
    }

    let total_patch_count = left_cells.emitted_patches.len() + right_cells.emitted_patches.len();
    let mut emitted_patches = Vec::with_capacity(total_patch_count);
    let mut emitted_provenance = Vec::with_capacity(total_patch_count);
    append_owner_support_cells(
        &mut emitted_patches,
        &mut emitted_provenance,
        left_cells,
        total_patch_count,
    );
    append_owner_support_cells(
        &mut emitted_patches,
        &mut emitted_provenance,
        right_cells,
        total_patch_count,
    );
    if emitted_patches.len() != total_patch_count || emitted_provenance.len() != total_patch_count {
        return None;
    }

    base.evaluated = reports.iter().map(|report| report.evaluated).sum();
    base.rejected = reports.iter().map(|report| report.rejected).sum();
    base.target_sector = None;
    base.support_cell_primitive = Kp20SupportCellPrimitiveReport {
        evaluated: left_cells.evaluated + right_cells.evaluated,
        target_count: left_cells.target_count + right_cells.target_count,
        patch_count: total_patch_count,
        folded_patch_count: left_cells.folded_patch_count + right_cells.folded_patch_count,
        cells_u: left_cells.cells_u,
        cells_v: left_cells.cells_v,
        topology_path_g0: left_cells
            .topology_path_g0
            .max(right_cells.topology_path_g0),
        selection_policy: "multi_owner_support_cell_combined",
        construction_mode: combined_construction_mode(left_cells, right_cells),
        target_piece: Kp20PieceId::NearMidpointLeft,
        emitted_patches,
        emitted_provenance,
        committed_floor_patches: Vec::new(),
        decline_reason: "support_cell_primitive_multi_owner_fold_free",
    };
    base.decline_reason = "expanded_support_topology_multi_owner_support_cell";
    Some(base)
}

fn support_cell_primitive_is_combinable(cells: &Kp20SupportCellPrimitiveReport) -> bool {
    cells.folded_patch_count == 0
        && !cells.emitted_patches.is_empty()
        && !cells.emitted_provenance.is_empty()
        && cells.emitted_patches.len() == cells.emitted_provenance.len()
        && cells.topology_path_g0.is_finite()
}

fn support_cell_targets_overlap(
    left: &Kp20SupportCellPrimitiveReport,
    right: &Kp20SupportCellPrimitiveReport,
) -> bool {
    left.emitted_provenance.iter().any(|left_item| {
        right
            .emitted_provenance
            .iter()
            .any(|right_item| left_item.sector_index == right_item.sector_index)
    })
}

fn append_owner_support_cells(
    emitted_patches: &mut Vec<popkern_core::geometry::gregory::GregoryPatch>,
    emitted_provenance: &mut Vec<crate::types::CornerSubPatchProvenance>,
    cells: &Kp20SupportCellPrimitiveReport,
    total_patch_count: usize,
) {
    for (patch, provenance) in cells
        .emitted_patches
        .iter()
        .cloned()
        .zip(cells.emitted_provenance.iter().copied())
    {
        let mut rewritten = provenance;
        rewritten.piece_index = emitted_provenance.len();
        rewritten.piece_count = total_patch_count;
        emitted_patches.push(patch);
        emitted_provenance.push(rewritten);
    }
}

fn combined_construction_mode(
    left: &Kp20SupportCellPrimitiveReport,
    right: &Kp20SupportCellPrimitiveReport,
) -> &'static str {
    if left.construction_mode == right.construction_mode {
        left.construction_mode
    } else {
        "multi_owner_mixed_construction_modes"
    }
}

fn support_owner_candidate_order() -> [Kp20PieceId; 2] {
    [
        Kp20PieceId::NearMidpointLeft,
        Kp20PieceId::NearMidpointRight,
    ]
}

fn preserves_local_support_cell_floor(report: &Kp20ExpandedSupportTopologyReport) -> bool {
    let cells = &report.support_cell_primitive;
    !report.emitted_patches.is_empty()
        && !cells.emitted_patches.is_empty()
        && cells.folded_patch_count == 0
        && cells.topology_path_g0.is_finite()
}

fn measure_support_owner_candidate(
    baseline: &Kp20SupportNormalLiftReport,
    target_piece: Kp20PieceId,
) -> Kp20ExpandedSupportTopologyReport {
    let target = worst_support_site(&baseline.summary, target_piece);
    let mut report = initial_report(
        baseline.summary.fold_count,
        target_piece,
        target.map(|site| site.sector),
    );
    let Some(target) = target else {
        report.decline_reason = "expanded_support_topology_no_target_site";
        return report;
    };
    let targets = support_sites_by_sector(&baseline.summary, target_piece);
    if targets
        .iter()
        .any(|site| site.sector >= baseline.macro_patches.len())
    {
        report.decline_reason = "expanded_support_topology_target_out_of_range";
        return report;
    }

    let support_piece = support_owner_piece(target.piece);
    report.support_cell_primitive = measure_support_cell_primitive(
        &baseline.macro_patches,
        &baseline.summary.fold_sites,
        support_piece,
    );
    let retained_macro_folds = retained_fold_count(&baseline.summary, &targets, support_piece);
    report.retained_macro_folds = retained_macro_folds;
    let mut best: Option<Kp20ExpandedSupportTopologyTrial> = None;
    'search: {
        for split_parameter in [0.38, 0.50, 0.62] {
            for fragment_rebuild in [false, true] {
                report.evaluated += 1;
                let macro_fragment_folds = fragmented_macro_fold_count(
                    &baseline.macro_patches,
                    &targets,
                    split_parameter,
                    fragment_rebuild,
                );
                let trial = expanded_support_topology_trial(
                    &baseline.macro_patches,
                    &targets,
                    support_piece,
                    split_parameter,
                    fragment_rebuild,
                    None,
                    retained_macro_folds,
                    macro_fragment_folds,
                );
                if best.as_ref().is_none_or(|best| trial.better_than(best)) {
                    best = Some(trial);
                }
            }
        }
        for split_parameter in [0.42, 0.50, 0.58] {
            for fragment_rebuild in [false, true] {
                let macro_fragment_folds = fragmented_macro_fold_count(
                    &baseline.macro_patches,
                    &targets,
                    split_parameter,
                    fragment_rebuild,
                );
                for mode in Kp20ExpandedSupportMode::all() {
                    for scale in [0.65] {
                        report.evaluated += 1;
                        let trial = expanded_support_topology_trial(
                            &baseline.macro_patches,
                            &targets,
                            support_piece,
                            split_parameter,
                            fragment_rebuild,
                            Some((mode, scale)),
                            retained_macro_folds,
                            macro_fragment_folds,
                        );
                        if best.as_ref().is_none_or(|best| trial.better_than(best)) {
                            let improves = trial.measured_folds < baseline.summary.fold_count;
                            best = Some(trial);
                            if improves {
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
    }

    let Some(best) = best else {
        return report;
    };
    report.rejected = report.evaluated;
    let emitted = emitted_candidate_patchwork(
        &baseline.macro_patches,
        &targets,
        support_piece,
        best.split_parameter,
        best.fragment_rebuild,
        best.reauthor,
    );
    report.emitted_patches = emitted.patches;
    report.emitted_provenance = emitted.provenance;
    report.fold_localization = best.fold_localization;
    let emitted_fold_localization =
        localize_expanded_support_folds(&report.emitted_patches, &report.emitted_provenance);
    report.expanded_patch_folds = emitted_fold_localization.total_folds();
    report.measured_folds = report.expanded_patch_folds;
    report.fold_localization = emitted_fold_localization;
    report.emitted_patch_count = report.emitted_patches.len();
    report.measured_patch_count = targets.len() * 2;
    let shape = expanded_support_commit_shape(
        baseline.macro_patches.len(),
        targets.len(),
        report.emitted_patch_count,
    );
    report.topology_face_count = shape.face_count;
    report.topology_pentagon_count = shape.pentagon_count;
    report.topology_alignment = shape.alignment;
    report.split_parameter = best.split_parameter;
    report.fragment_rebuild = best.fragment_rebuild;
    report.support_lane_reauthor = support_lane_reauthor_label(best.reauthor);
    report.support_lane_reauthor_scale = best.reauthor.map(|(_, scale)| scale);
    report.decline_reason = if best.measured_folds < baseline.summary.fold_count {
        "expanded_support_topology_not_emittable"
    } else {
        "expanded_support_topology_no_fold_reduction"
    };
    report
}

fn initial_report(
    start_folds: usize,
    target_piece: Kp20PieceId,
    target_sector: Option<usize>,
) -> Kp20ExpandedSupportTopologyReport {
    Kp20ExpandedSupportTopologyReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds,
        measured_folds: start_folds,
        retained_macro_folds: start_folds,
        expanded_patch_folds: 0,
        target_piece,
        target_sector,
        emitted_patch_count: 0,
        measured_patch_count: 0,
        topology_face_count: 0,
        topology_pentagon_count: 0,
        topology_alignment: "expanded_support_topology_no_shape",
        split_parameter: 0.0,
        fragment_rebuild: false,
        support_lane_reauthor: "raw_support_split",
        support_lane_reauthor_scale: None,
        support_cell_primitive: Kp20SupportCellPrimitiveReport::default(),
        emitted_patches: Vec::new(),
        emitted_provenance: Vec::new(),
        fold_localization: Kp20ExpandedSupportFoldLocalization::default(),
        decline_reason: "expanded_support_topology_no_candidate",
    }
}

fn report_is_better(
    candidate: &Kp20ExpandedSupportTopologyReport,
    current: &Kp20ExpandedSupportTopologyReport,
) -> bool {
    report_score(candidate) < report_score(current)
}

fn report_score(report: &Kp20ExpandedSupportTopologyReport) -> (bool, usize, usize, usize, usize) {
    (
        report.emitted_patches.is_empty(),
        report.retained_macro_folds + report.expanded_patch_folds,
        report.fold_localization.support_lane_folds(),
        report.emitted_patch_count,
        owner_rank(report.target_piece),
    )
}

fn owner_rank(piece: Kp20PieceId) -> usize {
    match piece {
        Kp20PieceId::NearMidpointLeft => 0,
        Kp20PieceId::NearMidpointRight => 1,
        _ => 2,
    }
}
