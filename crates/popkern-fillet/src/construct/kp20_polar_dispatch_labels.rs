//! Structured label helpers for KP20 polar dispatch diagnostics.

use popkern_core::geometry::gregory::{Kp20Bi6ValidationReport, Kp20PieceId};

use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::super::kp20_polar_support_reparameterization::{
    Kp20CoupledDomainReport, Kp20SplitDomainReport, Kp20SupportDomainRefineReport,
    Kp20SupportDomainReport, Kp20SupportReparameterizationReport, Kp20SupportTransitionReport,
};
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
        "fold_zones=origin_corner:{}|edge_shoulder:{}|interior:{};fold_ownership=origin_pieces:{}|support_pieces:{};dominant_fold_zone={};recommended_next_move={};{}",
        zones.origin_corner,
        zones.edge_shoulder,
        zones.interior,
        origin_folds,
        support_folds,
        dominant_fold_zone_label(zones.origin_corner, zones.edge_shoulder, zones.interior),
        recommended_next_move_label(origin_folds, support_folds, zones.interior),
        fold_site_localization_label(summary),
    )
}

fn fold_site_localization_label(summary: &Kp20SamplerFoldSummary) -> String {
    let target = Kp20PieceId::NearMidpointLeft;
    let target_sites = summary
        .fold_sites
        .iter()
        .filter(|site| site.piece == target)
        .count();
    let worst = summary
        .worst_fold_site()
        .map(|site| {
            format!(
                "sector:{}|piece:{:?}|uv:({:.3e},{:.3e})|cos:{:.3e}|ratio:{:.3e}",
                site.sector,
                site.piece,
                site.u,
                site.v,
                site.min_jacobian_normal_cos,
                site.min_jacobian_ratio
            )
        })
        .unwrap_or_else(|| "none".to_string());
    format!(
        "fold_site_count={};fold_site_target_piece={:?};fold_site_target_count={};fold_site_worst={};fold_site_touching={}",
        summary.fold_sites.len(),
        target,
        target_sites,
        worst,
        fold_site_touching_label(summary),
    )
}

fn fold_site_touching_label(summary: &Kp20SamplerFoldSummary) -> String {
    let mut touched = Vec::new();
    for piece in Kp20PieceId::all() {
        let sectors = summary
            .fold_sites
            .iter()
            .filter(|site| site.piece == piece)
            .map(|site| site.sector.to_string())
            .collect::<Vec<_>>();
        if !sectors.is_empty() {
            touched.push(format!("{piece:?}:[{}]", sectors.join(",")));
        }
    }
    if touched.is_empty() {
        "none".to_string()
    } else {
        touched.join("|")
    }
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
    support_reparameterization: &Kp20SupportReparameterizationReport,
    support_transition: &Kp20SupportTransitionReport,
    support_domain: &Kp20SupportDomainReport,
    support_domain_refine: &Kp20SupportDomainRefineReport,
    coupled_domain: &Kp20CoupledDomainReport,
    split_domain: &Kp20SplitDomainReport,
) -> String {
    format!(
        "status=sampled;default_folds={};folds={};piece_folds={};{};{};origin_inner_evaluated={};origin_inner_accepted={};origin_inner_rejected={};origin_inner_start_folds={};origin_inner_end_folds={};origin_inner_mode={};origin_inner_scale={:.3e};origin_inner_pre_evaluated={};origin_inner_pre_accepted={};origin_inner_pre_rejected={};origin_inner_pre_mode={};paired_inner_evaluated={};paired_inner_accepted={};paired_inner_rejected={};paired_inner_start_folds={};paired_inner_end_folds={};paired_inner_mode={};paired_inner_scale={:.3e};support_inner_evaluated={};support_inner_accepted={};support_inner_rejected={};support_inner_start_folds={};support_inner_end_folds={};support_inner_mode={};support_inner_scale={:.3e};support_domain_evaluated={};support_domain_accepted={};support_domain_rejected={};support_domain_start_folds={};support_domain_end_folds={};support_domain_target_piece={:?};support_domain_origin=({:.3e},{:.3e});support_domain_span=({:.3e},{:.3e});support_domain_skew={:.3e};support_domain_decline_reason={};support_domain_refine_evaluated={};support_domain_refine_accepted={};support_domain_refine_rejected={};support_domain_refine_start_folds={};support_domain_refine_end_folds={};support_domain_refine_target_piece={:?};support_domain_refine_origin=({:.3e},{:.3e});support_domain_refine_span=({:.3e},{:.3e});support_domain_refine_skew={:.3e};support_domain_refine_decline_reason={};coupled_domain_evaluated={};coupled_domain_accepted={};coupled_domain_rejected={};coupled_domain_start_folds={};coupled_domain_end_folds={};coupled_domain_support_origin=({:.3e},{:.3e});coupled_domain_support_span=({:.3e},{:.3e});coupled_domain_support_skew={:.3e};coupled_domain_origin_shift={:.3e};coupled_domain_origin_scale={:.3e};coupled_domain_decline_reason={};split_domain_evaluated={};split_domain_accepted={};split_domain_rejected={};split_domain_start_folds={};split_domain_end_folds={};split_domain_support_origin=({:.3e},{:.3e});split_domain_support_span=({:.3e},{:.3e});split_domain_support_skew={:.3e};split_domain_fraction={:.3e};split_domain_interface_pull={:.3e};split_domain_decline_reason={};support_transition_evaluated={};support_transition_accepted={};support_transition_rejected={};support_transition_start_folds={};support_transition_end_folds={};support_transition_target_piece={:?};support_transition_alpha={:.3e};support_transition_normal_scale={:.3e};support_transition_skew_scale={:.3e};support_transition_decline_reason={};support_reparam_evaluated={};support_reparam_fast_accepted={};support_reparam_fast_rejected={};support_reparam_lift_evaluated={};support_reparam_accepted={};support_reparam_rejected={};support_reparam_start_folds={};support_reparam_end_folds={};support_reparam_target_piece={:?};support_reparam_params={};support_reparam_decline_reason={};optimizer_active_fold_sites={};optimizer_row_variables={};optimizer_evaluated={};optimizer_accepted={};optimizer_rejected={};accepted_steps={}",
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
        origin_inner.pre_evaluated,
        origin_inner.pre_accepted,
        origin_inner.pre_rejected,
        origin_inner.pre_mode,
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
        support_domain.evaluated,
        support_domain.accepted,
        support_domain.rejected,
        support_domain.start_folds,
        support_domain.end_folds,
        support_domain.target_piece,
        support_domain.origin_u,
        support_domain.origin_v,
        support_domain.span_u,
        support_domain.span_v,
        support_domain.skew,
        support_domain.decline_reason,
        support_domain_refine.evaluated,
        support_domain_refine.accepted,
        support_domain_refine.rejected,
        support_domain_refine.start_folds,
        support_domain_refine.end_folds,
        support_domain_refine.target_piece,
        support_domain_refine.origin_u,
        support_domain_refine.origin_v,
        support_domain_refine.span_u,
        support_domain_refine.span_v,
        support_domain_refine.skew,
        support_domain_refine.decline_reason,
        coupled_domain.evaluated,
        coupled_domain.accepted,
        coupled_domain.rejected,
        coupled_domain.start_folds,
        coupled_domain.end_folds,
        coupled_domain.support_origin_u,
        coupled_domain.support_origin_v,
        coupled_domain.support_span_u,
        coupled_domain.support_span_v,
        coupled_domain.support_skew,
        coupled_domain.origin_shift,
        coupled_domain.origin_scale,
        coupled_domain.decline_reason,
        split_domain.evaluated,
        split_domain.accepted,
        split_domain.rejected,
        split_domain.start_folds,
        split_domain.end_folds,
        split_domain.support_origin_u,
        split_domain.support_origin_v,
        split_domain.support_span_u,
        split_domain.support_span_v,
        split_domain.support_skew,
        split_domain.split_fraction,
        split_domain.interface_pull,
        split_domain.decline_reason,
        support_transition.evaluated,
        support_transition.accepted,
        support_transition.rejected,
        support_transition.start_folds,
        support_transition.end_folds,
        support_transition.target_piece,
        support_transition.alpha,
        support_transition.normal_scale,
        support_transition.skew_scale,
        support_transition.decline_reason,
        support_reparameterization.evaluated,
        support_reparameterization.fast_accepted,
        support_reparameterization.fast_rejected,
        support_reparameterization.lift_evaluated,
        support_reparameterization.accepted,
        support_reparameterization.rejected,
        support_reparameterization.start_folds,
        support_reparameterization.end_folds,
        support_reparameterization.target_piece,
        support_reparameterization.params,
        support_reparameterization.decline_reason,
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
