//! Whole-candidate validation comparison for exact-boundary commit rollback.

use popkern_core::error::ValidationError;
use popkern_core::topology::FaceHandle;

use super::{ExactBoundaryCommitDeclineReason, ExactBoundaryCommitDiagnostic};

pub(crate) fn post_commit_validation_diagnostic(
    baseline_errors: &[ValidationError],
    post_errors: &[ValidationError],
    tolerance: f64,
) -> Option<ExactBoundaryCommitDiagnostic> {
    if post_errors.is_empty() {
        return None;
    }
    if baseline_errors.is_empty() {
        let err = most_informative_validation_error(post_errors)?;
        return Some(post_commit_validation_failed(
            validation_error_face(err),
            validation_error_residual(err),
        ));
    }
    let baseline_residual = validation_error_max_residual(baseline_errors).unwrap_or(0.0);
    let post_residual = validation_error_max_residual(post_errors).unwrap_or(0.0);
    if post_errors.len() > baseline_errors.len()
        || post_residual > baseline_residual + tolerance.max(1.0e-9)
    {
        let err = most_informative_validation_error(post_errors)?;
        return Some(post_commit_validation_failed(
            validation_error_face(err),
            validation_error_residual(err).or(Some(post_residual)),
        ));
    }
    None
}

fn most_informative_validation_error(errors: &[ValidationError]) -> Option<&ValidationError> {
    errors.iter().max_by(|a, b| {
        validation_error_residual(a)
            .unwrap_or(0.0)
            .total_cmp(&validation_error_residual(b).unwrap_or(0.0))
    })
}

fn validation_error_max_residual(errors: &[ValidationError]) -> Option<f64> {
    errors
        .iter()
        .filter_map(validation_error_residual)
        .max_by(|a, b| a.total_cmp(b))
}

fn post_commit_validation_failed(
    face: Option<FaceHandle>,
    residual: Option<f64>,
) -> ExactBoundaryCommitDiagnostic {
    ExactBoundaryCommitDiagnostic {
        reason: ExactBoundaryCommitDeclineReason::PostCommitValidationFailed,
        code: ExactBoundaryCommitDeclineReason::PostCommitValidationFailed.code(),
        sector_name: None,
        face,
        edge: None,
        cap_half_edge: None,
        support_half_edge: None,
        shared_trim_failure: None,
        residual,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
    }
}

fn validation_error_residual(err: &ValidationError) -> Option<f64> {
    match err {
        ValidationError::OpenLoop { gap, .. } | ValidationError::ShellGap { gap, .. } => Some(*gap),
        ValidationError::CurveSurfaceDeviation { max_deviation, .. } => Some(*max_deviation),
        _ => None,
    }
}

fn validation_error_face(err: &ValidationError) -> Option<FaceHandle> {
    match err {
        ValidationError::InconsistentOrientation { face }
        | ValidationError::TwistedSurface { face }
        | ValidationError::LoopPeriodWinding { face, .. } => Some(*face),
        ValidationError::SelfIntersection { face_a, .. } => Some(*face_a),
        _ => None,
    }
}
