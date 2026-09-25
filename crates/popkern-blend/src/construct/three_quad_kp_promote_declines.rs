//! Decline reports of the K-P bi-nonic promotion. Split from
//! `three_quad_kp_promote.rs` for the 500-line structural cap.

use crate::construct::promotion_ladder::KpBinonicPromotionReport;

/// Decline report without selection facts.
pub(super) fn declined(
    reason_code: &'static str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
) -> KpBinonicPromotionReport {
    declined_with_selection(
        reason_code,
        n_sides,
        detail,
        residual,
        tolerance,
        Vec::new(),
    )
}

/// Decline report carrying the boundary-candidates selection fact
/// (and any PK-275 outer-fillet rung facts) alongside the decline
/// itself, so a declined promotion still records which candidate
/// paths were tried and why each fell short (PK-249 observability).
pub(super) fn declined_with_selection(
    reason_code: &'static str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
    selection_fact: Vec<crate::types::ApexConstructionDiagnostic>,
) -> KpBinonicPromotionReport {
    let mut diagnostics = vec![
        crate::construct::three_quad_kp_diagnostics::promotion_declined(
            reason_code,
            n_sides,
            detail,
            residual,
            tolerance,
        ),
    ];
    diagnostics.extend(selection_fact);
    KpBinonicPromotionReport {
        relaxed_coplanarity: None,
        diagnostics,
    }
}
