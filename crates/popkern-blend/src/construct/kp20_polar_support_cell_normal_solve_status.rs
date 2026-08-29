//! Diagnostic formatting for support-cell normal-solve reports.

use super::worst_g1_normal_solve;

pub(super) fn report_was_kept(report: &worst_g1_normal_solve::WorstG1NormalSolveReport) -> bool {
    report.status == "kept"
}

pub(super) fn worst_g1_normal_solve_status(
    report: &worst_g1_normal_solve::WorstG1NormalSolveReport,
) -> String {
    format!(
        "worst_g1_normal_solve_status={};worst_g1_normal_decline_reason={};worst_g1_normal_control_authority={};worst_g1_normal_sample_count={};worst_g1_normal_variable_count={};worst_g1_normal_trial_count={};worst_g1_normal_validation_reject_count={};worst_g1_normal_quality_reject_count={};worst_g1_normal_limited_control_count={};max_worst_g1_normal_control_delta={:.6e};worst_g1_normal_damping={};best_worst_g1_normal_trial_damping={};best_worst_g1_normal_trial_g1={};best_worst_g1_normal_trial_local={};worst_g1_before={:.15e};worst_g1_after={:.15e};local_worst_g1_before={:.15e};local_worst_g1_after={:.15e}",
        report.status,
        report.decline_reason,
        report.control_authority,
        report.sample_count,
        report.variable_count,
        report.trial_count,
        report.validation_reject_count,
        report.quality_reject_count,
        report.limited_control_count,
        report.max_control_delta,
        report
            .damping
            .map(|damping| format!("{damping:.15e}"))
            .unwrap_or_else(|| "none".to_string()),
        optional_f64_status(report.best_trial_damping),
        optional_f64_status(report.best_trial_worst_g1),
        optional_f64_status(report.best_trial_local_disagreement),
        report.before_worst_g1,
        report.after_worst_g1,
        report.local_disagreement_before,
        report.local_disagreement_after
    )
}

fn optional_f64_status(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.15e}"))
        .unwrap_or_else(|| "none".to_string())
}

pub(super) fn declined_authority_summary(
    prefix: &str,
    report: &worst_g1_normal_solve::WorstG1NormalSolveReport,
) -> String {
    format!(
        "{prefix}_status={};{prefix}_decline_reason={};{prefix}_best_trial_g1={};{prefix}_best_trial_local={};{prefix}_variable_count={};{prefix}_limited_control_count={}",
        report.status,
        report.decline_reason,
        optional_f64_status(report.best_trial_worst_g1),
        optional_f64_status(report.best_trial_local_disagreement),
        report.variable_count,
        report.limited_control_count
    )
}
