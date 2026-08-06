//! Structured diagnostics for K-P rectangular promotion transactions.

pub(super) fn outer_boundary_residual_exceeded(
    promoted_residual: f64,
    unpromoted_residual: Option<f64>,
    gate: f64,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "kp_binonic.promotion".to_string(),
            code: "kp_binonic.promotion.outer_boundary_residual_exceeded".to_string(),
            detail_code: unpromoted_residual.map(|residual| {
                format!(
                    "unpromoted_boundary_residual={residual:.3e};promoted_boundary_residual={promoted_residual:.3e}"
                )
            }),
            sector_name: "kp_binonic_ring".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: Some(promoted_residual),
            outer_side: None,
            sample_parameter: None,
            tolerance: Some(gate),
            shared_fact_ids: Vec::new(),
        },
    )
}

pub(super) fn promotion_declined(
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: format!("kp_binonic.promotion.{reason_code}"),
            code: "kp_binonic.promotion.declined".to_string(),
            detail_code: detail,
            sector_name: format!("kp_binonic_ring_n{n_sides}"),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual,
            outer_side: None,
            sample_parameter: None,
            tolerance,
            shared_fact_ids: Vec::new(),
        },
    )
}

pub(super) fn support_cell_candidate_nonselection(
    reason_code: &str,
    n_sides: usize,
    patch_count: usize,
) -> crate::types::ApexConstructionDiagnostic {
    super::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_candidate_nonselection",
        n_sides,
        Some(format!(
            "evaluated=0;target_count=0;patch_count={patch_count};folded_patch_count=0;cells_u=0;cells_v=0;topology_path_g0=inf;selection_policy=support_cell_not_evaluated;construction_mode=none;target_piece=NearMidpointLeft;committed_floor_patch_count=0;decline_reason=support_cell_primitive_not_evaluated_{reason_code}"
        )),
        None,
        None,
    )
}

/// PK-240: the balanced-spoke fan engaged — the natural sector split
/// at the ring center was degenerate and the spoke tangents were
/// redistributed to a regular 2π/N spacing so the fan matches the
/// regular characteristic map the KP promotion rungs assume.
pub(super) fn balanced_spoke_fan_engaged(
    natural_max_span: f64,
    balanced_max_span: f64,
    strength: f64,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "three_quad.balanced_spoke_fan".to_string(),
            code: "three_quad.balanced_spoke_fan".to_string(),
            detail_code: Some(format!(
                "natural_max_span_deg={:.1};balanced_max_span_deg={:.1};strength={strength:.2}",
                natural_max_span.to_degrees(),
                balanced_max_span.to_degrees()
            )),
            sector_name: "n_quad_ring".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: Some(natural_max_span),
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}

/// PK-240: the fold-triggered balanced-spoke retry ran and no
/// strength served — per-strength decline reasons recorded so the
/// class stays attributable (e.g. 1:100 declines every strength).
pub(super) fn balanced_spoke_fan_declined(
    reasons: &str,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "three_quad.balanced_spoke_fan.declined".to_string(),
            code: "three_quad.balanced_spoke_fan.declined".to_string(),
            detail_code: Some(reasons.to_string()),
            sector_name: "n_quad_ring".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: None,
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}
