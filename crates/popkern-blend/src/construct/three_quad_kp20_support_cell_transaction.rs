use popkern_core::error::ValidationError;

use crate::types::ApexConstructionDiagnostic;

#[derive(Debug)]
pub(super) enum Kp20SupportCellApplyTransactionResult<T> {
    Promoted(T),
    Fallback(T),
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportCellValidationFailure {
    kind: &'static str,
    class: &'static str,
    detail: String,
    edge: Option<String>,
    gap: Option<f64>,
}

impl Kp20SupportCellValidationFailure {
    fn construction(error: String) -> Self {
        Self {
            kind: "ConstructionError",
            class: "construction_commit_failure",
            detail: error,
            edge: None,
            gap: None,
        }
    }

    #[cfg(test)]
    fn shell_gap_for_test(edge: &str, gap: f64) -> Self {
        Self {
            kind: "ShellGap",
            class: "shared_boundary_geometry_mismatch",
            detail: format!("ShellGap {{ edge: {edge}, gap: {gap} }}"),
            edge: Some(edge.to_string()),
            gap: Some(gap),
        }
    }

    fn detail(&self) -> &str {
        &self.detail
    }
}

impl From<ValidationError> for Kp20SupportCellValidationFailure {
    fn from(error: ValidationError) -> Self {
        match error {
            ValidationError::ShellGap { edge, gap } => Self {
                kind: "ShellGap",
                class: "shared_boundary_geometry_mismatch",
                detail: format!("ShellGap {{ edge: {edge:?}, gap: {gap} }}"),
                edge: Some(format!("{edge:?}")),
                gap: Some(gap),
            },
            other => Self {
                kind: validation_error_kind(&other),
                class: validation_error_class(&other),
                detail: format!("{other:?}"),
                edge: None,
                gap: None,
            },
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run_kp20_support_cell_apply_transaction<S: Clone, T>(
    state: &mut S,
    n_sides: usize,
    promoted_surface_count: usize,
    promoted_cells_per_sector: usize,
    fallback_surface_count: Option<usize>,
    apply_promoted: impl FnOnce(&mut S) -> Result<T, String>,
    validate: impl Fn(&S) -> Result<(), Kp20SupportCellValidationFailure>,
    apply_fallback: impl FnOnce(&mut S, ApexConstructionDiagnostic) -> Result<T, String>,
) -> Result<Kp20SupportCellApplyTransactionResult<T>, String> {
    let before_promoted = state.clone();
    match apply_promoted(state) {
        Ok(result) => match validate(state) {
            Ok(()) => Ok(Kp20SupportCellApplyTransactionResult::Promoted(result)),
            Err(error) => {
                *state = before_promoted;
                run_support_cell_fallback(
                    state,
                    n_sides,
                    promoted_surface_count,
                    promoted_cells_per_sector,
                    fallback_surface_count,
                    error,
                    apply_fallback,
                    validate,
                )
            }
        },
        Err(error) => {
            *state = before_promoted;
            run_support_cell_fallback(
                state,
                n_sides,
                promoted_surface_count,
                promoted_cells_per_sector,
                fallback_surface_count,
                Kp20SupportCellValidationFailure::construction(error),
                apply_fallback,
                validate,
            )
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_support_cell_fallback<S: Clone, T>(
    state: &mut S,
    n_sides: usize,
    promoted_surface_count: usize,
    promoted_cells_per_sector: usize,
    fallback_surface_count: Option<usize>,
    error: Kp20SupportCellValidationFailure,
    apply_fallback: impl FnOnce(&mut S, ApexConstructionDiagnostic) -> Result<T, String>,
    validate: impl Fn(&S) -> Result<(), Kp20SupportCellValidationFailure>,
) -> Result<Kp20SupportCellApplyTransactionResult<T>, String> {
    if promoted_cells_per_sector <= 4 {
        return Err(format!(
            "KP20 support-cell commit left invalid BRep: {}",
            error.detail()
        ));
    }
    let Some(fallback_surface_count) = fallback_surface_count else {
        return Err(format!(
            "KP20 higher-resolution support-cell commit left invalid BRep and no four-cell floor was available: {}",
            error.detail()
        ));
    };
    let before_fallback = state.clone();
    let diagnostic = support_cell_transaction_decline(
        n_sides,
        promoted_surface_count,
        fallback_surface_count,
        error,
    );
    let result = apply_fallback(state, diagnostic)?;
    if validate(state).is_ok() {
        return Ok(Kp20SupportCellApplyTransactionResult::Fallback(result));
    }
    *state = before_fallback;
    Err("KP20 support-cell four-cell floor fallback left invalid BRep".to_string())
}

fn support_cell_transaction_decline(
    n_sides: usize,
    promoted_surface_count: usize,
    fallback_surface_count: usize,
    failure: Kp20SupportCellValidationFailure,
) -> ApexConstructionDiagnostic {
    let failure_gap = failure
        .gap
        .map(|gap| format!("{gap:.15e}"))
        .unwrap_or_else(|| "none".to_string());
    let failure_edge = failure.edge.as_deref().unwrap_or("none");
    crate::construct::kp20_polar_diagnostics::dispatch_declined(
        "support_cell_transaction_promoted_apply",
        n_sides,
        Some(format!(
            "promoted_surface_count={promoted_surface_count};fallback_surface_count={fallback_surface_count};failure_kind={};failure_class={};failure_edge={failure_edge};failure_gap={failure_gap};error={}",
            failure.kind,
            failure.class,
            failure.detail()
        )),
        None,
        None,
    )
}

fn validation_error_class(error: &ValidationError) -> &'static str {
    match error {
        ValidationError::ShellGap { .. } => "shared_boundary_geometry_mismatch",
        ValidationError::CurveSurfaceDeviation { .. } => "surface_curve_compatibility_mismatch",
        ValidationError::OpenLoop { .. } => "topology_loop_closure_mismatch",
        ValidationError::SelfIntersection { .. } | ValidationError::TwistedSurface { .. } => {
            "surface_shape_invalid"
        }
        ValidationError::LoopPeriodWinding { .. }
        | ValidationError::PcurveDeviation { .. }
        | ValidationError::CurveWindowOutsideSpan { .. } => "trim_domain_invalid",
        ValidationError::EulerPoincare { .. }
        | ValidationError::DanglingReference { .. }
        | ValidationError::BrokenTwin { .. }
        | ValidationError::BrokenCycle { .. }
        | ValidationError::InconsistentVertex { .. }
        | ValidationError::NonManifoldEdge { .. }
        | ValidationError::LowValence { .. }
        | ValidationError::InconsistentOrientation { .. }
        | ValidationError::WindingDishonest { .. }
        | ValidationError::InnerLoopWinding { .. }
        | ValidationError::CanonicalBackward { .. }
        | ValidationError::MissingPcurves { .. } => "topology_invariant_mismatch",
        ValidationError::Precision { .. } => "precision_budget_mismatch",
    }
}

fn validation_error_kind(error: &ValidationError) -> &'static str {
    match error {
        ValidationError::EulerPoincare { .. } => "EulerPoincare",
        ValidationError::DanglingReference { .. } => "DanglingReference",
        ValidationError::BrokenTwin { .. } => "BrokenTwin",
        ValidationError::BrokenCycle { .. } => "BrokenCycle",
        ValidationError::InconsistentVertex { .. } => "InconsistentVertex",
        ValidationError::NonManifoldEdge { .. } => "NonManifoldEdge",
        ValidationError::LowValence { .. } => "LowValence",
        ValidationError::InconsistentOrientation { .. } => "InconsistentOrientation",
        ValidationError::OpenLoop { .. } => "OpenLoop",
        ValidationError::CurveSurfaceDeviation { .. } => "CurveSurfaceDeviation",
        ValidationError::ShellGap { .. } => "ShellGap",
        ValidationError::SelfIntersection { .. } => "SelfIntersection",
        ValidationError::WindingDishonest { .. } => "WindingDishonest",
        ValidationError::InnerLoopWinding { .. } => "InnerLoopWinding",
        ValidationError::MissingPcurves { .. } => "MissingPcurves",
        ValidationError::TwistedSurface { .. } => "TwistedSurface",
        ValidationError::LoopPeriodWinding { .. } => "LoopPeriodWinding",
        ValidationError::PcurveDeviation { .. } => "PcurveDeviation",
        ValidationError::CanonicalBackward { .. } => "CanonicalBackward",
        ValidationError::CurveWindowOutsideSpan { .. } => "CurveWindowOutsideSpan",
        ValidationError::Precision { .. } => "Precision",
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn support_cell_assembled_quality_fact(
    n_sides: usize,
    selected: &'static str,
    selected_surface_count: usize,
    selected_cells_per_sector: usize,
    promoted_surface_count: usize,
    fallback_surface_count: Option<usize>,
    validation_status: &'static str,
    validation_error: Option<&str>,
) -> ApexConstructionDiagnostic {
    let fallback_surface_count = fallback_surface_count
        .map(|count| count.to_string())
        .unwrap_or_else(|| "none".to_string());
    let validation_error = validation_error.unwrap_or("none");
    crate::construct::kp20_polar_diagnostics::dispatch_fact(
        "support_cell_assembled_quality",
        n_sides,
        Some(format!(
            "selected={selected};selected_surface_count={selected_surface_count};selected_cells_per_sector={selected_cells_per_sector};promoted_surface_count={promoted_surface_count};fallback_surface_count={fallback_surface_count};validation_status={validation_status};validation_error={validation_error}"
        )),
        None,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_eight_cell_promoted_apply_restores_and_falls_back() {
        let mut state = 0;
        let result = run_kp20_support_cell_apply_transaction(
            &mut state,
            6,
            48,
            8,
            Some(24),
            |state| {
                *state = 1;
                Ok("promoted")
            },
            |state| {
                if *state == 1 {
                    Err(Kp20SupportCellValidationFailure::shell_gap_for_test(
                        "EdgeHandle(207v1)",
                        6.7e-3,
                    ))
                } else {
                    Ok(())
                }
            },
            |state, diagnostic| {
                *state = 2;
                let ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
                    panic!("expected exact-boundary diagnostic");
                };
                assert_eq!(cap.code, "kp20_polar_macro.dispatch.declined");
                assert_eq!(
                    cap.stage,
                    "kp20_polar_macro.dispatch.support_cell_transaction_promoted_apply"
                );
                let detail = cap.detail_code.expect("transaction detail");
                assert!(detail.contains("promoted_surface_count=48"));
                assert!(detail.contains("fallback_surface_count=24"));
                assert!(detail.contains("failure_kind=ShellGap"));
                assert!(detail.contains("failure_class=shared_boundary_geometry_mismatch"));
                assert!(detail.contains("failure_edge=EdgeHandle(207v1)"));
                assert!(detail.contains("failure_gap=6.700000000000000e-3"));
                Ok("fallback")
            },
        )
        .expect("transaction should fall back");

        assert!(matches!(
            result,
            Kp20SupportCellApplyTransactionResult::Fallback("fallback")
        ));
        assert_eq!(state, 2);
    }

    #[test]
    fn assembled_quality_fact_records_selected_inventory() {
        let ApexConstructionDiagnostic::ExactBoundaryCap(cap) = support_cell_assembled_quality_fact(
            6,
            "fallback_four_cell_floor",
            24,
            4,
            48,
            Some(24),
            "valid_after_promoted_decline",
            None,
        ) else {
            panic!("expected exact-boundary diagnostic");
        };

        assert_eq!(cap.code, "kp20_polar_macro.dispatch.fact");
        assert_eq!(
            cap.stage,
            "kp20_polar_macro.dispatch.support_cell_assembled_quality"
        );
        let detail = cap.detail_code.expect("quality detail");
        assert!(detail.contains("selected=fallback_four_cell_floor"));
        assert!(detail.contains("selected_surface_count=24"));
        assert!(detail.contains("selected_cells_per_sector=4"));
        assert!(detail.contains("promoted_surface_count=48"));
        assert!(detail.contains("fallback_surface_count=24"));
        assert!(detail.contains("validation_status=valid_after_promoted_decline"));
    }

    #[test]
    fn invalid_eight_cell_promoted_apply_without_floor_restores_and_declines() {
        let mut state = 0;
        let error = run_kp20_support_cell_apply_transaction(
            &mut state,
            6,
            48,
            8,
            None,
            |state| {
                *state = 1;
                Ok("promoted")
            },
            |_| {
                Err(Kp20SupportCellValidationFailure::shell_gap_for_test(
                    "EdgeHandle(207v1)",
                    6.7e-3,
                ))
            },
            |_, _| Ok("fallback"),
        )
        .expect_err("missing floor should decline");

        assert!(error.contains("no four-cell floor was available"));
        assert_eq!(state, 0);
    }

    #[test]
    fn invalid_sixteen_cell_promoted_apply_restores_and_falls_back() {
        let mut state = 0;
        let result = run_kp20_support_cell_apply_transaction(
            &mut state,
            6,
            96,
            16,
            Some(24),
            |state| {
                *state = 1;
                Ok("promoted")
            },
            |state| {
                if *state == 1 {
                    Err(Kp20SupportCellValidationFailure::shell_gap_for_test(
                        "EdgeHandle(311v1)",
                        2.5e-3,
                    ))
                } else {
                    Ok(())
                }
            },
            |state, diagnostic| {
                *state = 2;
                let ApexConstructionDiagnostic::ExactBoundaryCap(cap) = diagnostic else {
                    panic!("expected exact-boundary diagnostic");
                };
                let detail = cap.detail_code.expect("transaction detail");
                assert!(detail.contains("promoted_surface_count=96"));
                assert!(detail.contains("fallback_surface_count=24"));
                assert!(detail.contains("failure_edge=EdgeHandle(311v1)"));
                Ok("fallback")
            },
        )
        .expect("transaction should fall back");

        assert!(matches!(
            result,
            Kp20SupportCellApplyTransactionResult::Fallback("fallback")
        ));
        assert_eq!(state, 2);
    }
}
