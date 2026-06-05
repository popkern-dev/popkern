//! Structured label helpers for KP20 polar dispatch diagnostics.

use popkern_core::geometry::gregory::{Kp20Bi6ValidationReport, Kp20PieceId};

use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;

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
    if interior_folds > 0 && origin_folds >= support_folds {
        "origin_inner_interior_row"
    } else if interior_folds > 0 {
        "support_inner_interior_row"
    } else {
        "edge_shoulder_row"
    }
}
