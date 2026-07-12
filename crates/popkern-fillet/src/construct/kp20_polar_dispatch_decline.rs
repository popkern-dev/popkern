//! Structured decline recording for the KP20 polar dispatch
//! (split out of `kp20_polar_dispatch.rs` for the 500-line budget).

use super::*;

pub(super) fn record_candidate_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    summary: &Kp20SamplerFoldSummary,
    validation: Kp20Bi6ValidationReport,
) {
    record_decline(
        diagnostics,
        reason_code,
        n_sides,
        Some(format!(
            "folds={};piece_folds={};{}",
            summary.fold_count,
            piece_counts_label(summary.piece_fold_counts),
            validation_label(validation),
        )),
        Some(summary.min_jacobian_normal_cos),
    );
}

pub(super) fn record_self_intersection(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    idx: usize,
    detail: &popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    n: usize,
) {
    let label = detail.diagnostic_label();
    let sector = idx / 4;
    let piece = Kp20PieceId::all()[idx % 4];
    record_decline(
        diagnostics,
        "self_intersection",
        n,
        Some(format!(
            "sector={sector};piece={piece:?};patch_index={idx};{label}"
        )),
        Some(detail.min_jacobian_normal_cos),
    );
}

pub(super) fn record_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
) {
    diagnostics.push(crate::construct::kp20_polar_diagnostics::dispatch_declined(
        reason_code,
        n_sides,
        detail,
        residual,
        Some(SELF_INT_TOL),
    ));
}

pub(super) fn record_lift_budget_decline(
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
    n: usize,
    default_summary: &Kp20SamplerFoldSummary,
    bounded_summary: Option<&Kp20SamplerFoldSummary>,
) {
    let bounded_detail = bounded_summary
        .map(|summary| {
            format!(
                ";bounded_folds={};bounded_piece_folds={}",
                summary.fold_count,
                piece_counts_label(summary.piece_fold_counts),
            )
        })
        .unwrap_or_default();
    record_decline(
        diagnostics,
        "nonoverlap_coupled_lift_budget_exceeded",
        n,
        Some(format!(
            "folds={};piece_folds={};lift_solver_budget=limited_lift_only_for_valence_ge_6{}",
            default_summary.fold_count,
            piece_counts_label(default_summary.piece_fold_counts),
            bounded_detail,
        )),
        Some(default_summary.min_jacobian_normal_cos),
    );
}
