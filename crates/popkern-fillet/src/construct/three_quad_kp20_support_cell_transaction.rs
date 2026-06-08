use crate::types::ApexConstructionDiagnostic;

#[derive(Debug)]
pub(super) enum Kp20SupportCellApplyTransactionResult<T> {
    Promoted(T),
    Fallback(T),
}

pub(super) fn run_kp20_support_cell_apply_transaction<S: Clone, T>(
    state: &mut S,
    n_sides: usize,
    promoted_surface_count: usize,
    promoted_cells_per_sector: usize,
    fallback_surface_count: Option<usize>,
    apply_promoted: impl FnOnce(&mut S) -> Result<T, String>,
    validate: impl Fn(&S) -> Result<(), String>,
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
                error,
                apply_fallback,
                validate,
            )
        }
    }
}

fn run_support_cell_fallback<S: Clone, T>(
    state: &mut S,
    n_sides: usize,
    promoted_surface_count: usize,
    promoted_cells_per_sector: usize,
    fallback_surface_count: Option<usize>,
    error: String,
    apply_fallback: impl FnOnce(&mut S, ApexConstructionDiagnostic) -> Result<T, String>,
    validate: impl Fn(&S) -> Result<(), String>,
) -> Result<Kp20SupportCellApplyTransactionResult<T>, String> {
    if promoted_cells_per_sector != 8 {
        return Err(format!(
            "KP20 support-cell commit left invalid BRep: {error}"
        ));
    }
    let Some(fallback_surface_count) = fallback_surface_count else {
        return Err(format!(
            "KP20 eight-cell support-cell commit left invalid BRep and no four-cell floor was available: {error}"
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
    error: String,
) -> ApexConstructionDiagnostic {
    crate::construct::kp20_polar_diagnostics::dispatch_declined(
        "support_cell_transaction_promoted_apply",
        n_sides,
        Some(format!(
            "promoted_surface_count={promoted_surface_count};fallback_surface_count={fallback_surface_count};error={error}"
        )),
        None,
        None,
    )
}

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
                    Err("ShellGap(6.7e-3)".to_string())
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
                assert!(detail.contains("ShellGap"));
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
            |_| Err("ShellGap(6.7e-3)".to_string()),
            |_, _| Ok("fallback"),
        )
        .expect_err("missing floor should decline");

        assert!(error.contains("no four-cell floor was available"));
        assert_eq!(state, 0);
    }
}
