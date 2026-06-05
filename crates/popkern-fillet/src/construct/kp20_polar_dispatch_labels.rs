//! Structured label helpers for KP20 polar dispatch diagnostics.

use popkern_core::geometry::gregory::{Kp20Bi6ValidationReport, Kp20PieceId};

use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::origin_inner_row::{
    Kp20OriginInnerRowReport, Kp20PairedInnerRowReport, Kp20SupportInnerRowReport,
};

pub(super) fn piece_counts_label(counts: [usize; 4]) -> String {
    let pieces = Kp20PieceId::all();
    format!(
        "{:?}:{}|{:?}:{}|{:?}:{}|{:?}:{}",
        pieces[0], counts[0], pieces[1], counts[1], pieces[2], counts[2], pieces[3], counts[3]
    )
}

pub(super) fn fold_ownership_label(summary: &Kp20SamplerFoldSummary) -> String {
    let zones = summary.fold_zone_counts();
    let origin_folds = summary.piece_fold_counts[Kp20PieceId::NearOriginLeft.index()]
        + summary.piece_fold_counts[Kp20PieceId::NearOriginRight.index()];
    let support_folds = summary.piece_fold_counts[Kp20PieceId::NearMidpointRight.index()]
        + summary.piece_fold_counts[Kp20PieceId::NearMidpointLeft.index()];
    format!(
        "fold_zones=origin_corner:{}|edge_shoulder:{}|interior:{};fold_ownership=origin_pieces:{}|support_pieces:{};dominant_fold_zone={};recommended_next_move={}",
        zones.origin_corner,
        zones.edge_shoulder,
        zones.interior,
        origin_folds,
        support_folds,
        dominant_fold_zone_label(zones.origin_corner, zones.edge_shoulder, zones.interior),
        recommended_next_move_label(origin_folds, support_folds, zones.interior),
    )
}

pub(super) fn validation_label(validation: Kp20Bi6ValidationReport) -> String {
    format!(
        "tensor_border_residual={:.3e};internal_c2_position={:.3e};internal_c2_first={:.3e};origin_g2_strip={:.3e}",
        validation.max_tensor_border_residual,
        validation.max_internal_c2_position_residual,
        validation.max_internal_c2_first_residual,
        validation.max_origin_g2_strip_residual,
    )
}

pub(super) fn nonoverlap_coupled_lift_detail(
    default: &Kp20SamplerFoldSummary,
    lifted: &Kp20SamplerFoldSummary,
    report: &Kp20SupportNormalLiftReport,
    origin_inner: &Kp20OriginInnerRowReport,
    paired_inner: &Kp20PairedInnerRowReport,
    support_inner: &Kp20SupportInnerRowReport,
) -> String {
    format!(
        "status=sampled;default_folds={};folds={};piece_folds={};{};{};origin_inner_evaluated={};origin_inner_accepted={};origin_inner_rejected={};origin_inner_start_folds={};origin_inner_end_folds={};origin_inner_mode={};origin_inner_scale={:.3e};paired_inner_evaluated={};paired_inner_accepted={};paired_inner_rejected={};paired_inner_start_folds={};paired_inner_end_folds={};paired_inner_mode={};paired_inner_scale={:.3e};support_inner_evaluated={};support_inner_accepted={};support_inner_rejected={};support_inner_start_folds={};support_inner_end_folds={};support_inner_mode={};support_inner_scale={:.3e};optimizer_active_fold_sites={};optimizer_row_variables={};optimizer_evaluated={};optimizer_accepted={};optimizer_rejected={};accepted_steps={}",
        default.fold_count,
        lifted.fold_count,
        piece_counts_label(lifted.piece_fold_counts),
        fold_ownership_label(lifted),
        validation_label(report.validation),
        origin_inner.evaluated,
        origin_inner.accepted,
        origin_inner.rejected,
        origin_inner.start_folds,
        origin_inner.end_folds,
        origin_inner.mode,
        origin_inner.scale,
        paired_inner.evaluated,
        paired_inner.accepted,
        paired_inner.rejected,
        paired_inner.start_folds,
        paired_inner.end_folds,
        paired_inner.mode,
        paired_inner.scale,
        support_inner.evaluated,
        support_inner.accepted,
        support_inner.rejected,
        support_inner.start_folds,
        support_inner.end_folds,
        support_inner.mode,
        support_inner.scale,
        report.optimizer_stats.active_fold_sites,
        report.optimizer_stats.row_variable_count,
        report.optimizer_stats.evaluated_candidates,
        report.optimizer_stats.accepted_candidates,
        report.optimizer_stats.rejected_candidates,
        report.accepted_steps.len(),
    )
}

fn dominant_fold_zone_label(
    origin_corner: usize,
    edge_shoulder: usize,
    interior: usize,
) -> &'static str {
    if interior >= origin_corner && interior >= edge_shoulder {
        "interior"
    } else if edge_shoulder >= origin_corner {
        "edge_shoulder"
    } else {
        "origin_corner"
    }
}

fn recommended_next_move_label(
    origin_folds: usize,
    support_folds: usize,
    interior_folds: usize,
) -> &'static str {
    if interior_folds > 0 && origin_folds > support_folds {
        "origin_inner_interior_row"
    } else if interior_folds > 0 {
        "support_inner_interior_row"
    } else {
        "edge_shoulder_row"
    }
}
