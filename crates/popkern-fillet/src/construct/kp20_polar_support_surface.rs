//! Support-surface-owned KP20 polar candidate scaffold.

use popkern_core::geometry::gregory::{Kp20Bi6ValidationReport, Kp20PieceId};

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_support_image::Kp20SupportImageCandidate;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceReport {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) stats: Kp20SupportSurfaceStats,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceStats {
    pub(super) active_fold_pieces: usize,
    pub(super) orientation_rows: usize,
    pub(super) evaluated_candidates: usize,
    pub(super) accepted_candidates: usize,
    pub(super) rejected_candidates: usize,
}

/// Build the initial support-surface-owned candidate from the selected,
/// exact-validated support-image macro cap.
///
/// The current slice establishes semantic ownership and diagnostics for the
/// future construction primitive. It does not mutate rows or replace accepted
/// behavior until a candidate improves fold validation while preserving the
/// exact tensor-border and pole invariants.
pub(super) fn support_surface_candidate_summary(
    support_image: &Kp20SupportImageCandidate,
) -> Kp20SupportSurfaceReport {
    let active_fold_pieces = active_fold_piece_count(&support_image.summary);
    Kp20SupportSurfaceReport {
        summary: support_image.summary.clone(),
        validation: support_image.validation,
        stats: Kp20SupportSurfaceStats {
            active_fold_pieces,
            orientation_rows: active_fold_pieces,
            evaluated_candidates: 1,
            accepted_candidates: 0,
            rejected_candidates: 1,
        },
        decline_reason: "support_surface_orientation_rows_pending",
    }
}

fn active_fold_piece_count(summary: &Kp20SamplerFoldSummary) -> usize {
    let mut active = [false; 4];
    for site in &summary.fold_sites {
        active[site.piece.index()] = true;
    }
    Kp20PieceId::all()
        .into_iter()
        .filter(|piece| active[piece.index()])
        .count()
}
