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

/// PK-249: the chart-orientation repair engaged on a binonic
/// candidate whose accepted charts carried handedness-reversal
/// bands (or whose exact-injection candidate folded). `outcome` is
/// `repaired` or `declined(<reason>)`; the detail carries the
/// before/after band counts and the G1 delta so serves stay
/// attributable either way.
pub(super) fn chart_orientation_fact(
    outcome: &str,
    pre_bands: usize,
    post_bands: usize,
    pre_g1: f64,
    post_g1: f64,
    post_g2: Option<f64>,
) -> crate::types::ApexConstructionDiagnostic {
    let post_g2 = post_g2
        .map(|g2| format!(";post_g2={g2:.3e}"))
        .unwrap_or_default();
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "kp_binonic.chart_orientation".to_string(),
            code: "kp_binonic.chart_orientation".to_string(),
            detail_code: Some(format!(
                "outcome={outcome};pre_bands={pre_bands};post_bands={post_bands};pre_g1={pre_g1:.3e};post_g1={post_g1:.3e}{post_g2}"
            )),
            sector_name: "kp_binonic_ring".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: Some(post_g1),
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}

/// PK-249: which boundary-fidelity candidate the binonic promotion
/// selected (`guide` / `exact` / `relocked`) and why the stronger
/// paths fell short when they did — the exact-injection and relock
/// outcomes were previously silently filtered, leaving the
/// extreme-asymmetry decline class (1:100 boundary residual
/// 1.323e-1) unattributable.
pub(super) fn boundary_candidate_selection(
    provenance: &str,
    exact_outcome: &str,
    relock_outcome: &str,
    boundary_residual: Option<f64>,
    n_sides: usize,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "kp_binonic.promotion.boundary_candidates".to_string(),
            code: "kp_binonic.promotion.boundary_candidates".to_string(),
            detail_code: Some(format!(
                "provenance={provenance};exact={exact_outcome};relock={relock_outcome}"
            )),
            sector_name: format!("kp_binonic_ring_n{n_sides}"),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: boundary_residual,
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
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
    chart_reversal_bands: usize,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "three_quad.balanced_spoke_fan".to_string(),
            code: "three_quad.balanced_spoke_fan".to_string(),
            detail_code: Some(format!(
                "natural_max_span_deg={:.1};balanced_max_span_deg={:.1};strength={strength:.2};chart_reversal_bands={chart_reversal_bands}",
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

/// PK-249: a banded balanced-retry acceptance served after the
/// full strength walk found no chart-clean candidate — the walk's
/// per-strength outcomes ride the served report so the preference
/// stays attributable.
pub(super) fn balanced_spoke_fan_walk(walk: &str) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "three_quad.balanced_spoke_fan.walk".to_string(),
            code: "three_quad.balanced_spoke_fan.walk".to_string(),
            detail_code: Some(walk.to_string()),
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
