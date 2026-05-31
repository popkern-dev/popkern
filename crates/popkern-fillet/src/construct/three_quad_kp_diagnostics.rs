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
