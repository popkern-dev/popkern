//! Structured diagnostics for targeted KP20 support-side reparameterization.

use popkern_core::geometry::gregory::Kp20PieceId;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportReparameterizationReport {
    pub(super) evaluated: usize,
    pub(super) fast_accepted: usize,
    pub(super) fast_rejected: usize,
    pub(super) lift_evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) params: String,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportTransitionReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) alpha: f64,
    pub(super) normal_scale: f64,
    pub(super) skew_scale: f64,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportDomainReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) origin_u: f64,
    pub(super) origin_v: f64,
    pub(super) span_u: f64,
    pub(super) span_v: f64,
    pub(super) skew: f64,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportDomainRefineReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) origin_u: f64,
    pub(super) origin_v: f64,
    pub(super) span_u: f64,
    pub(super) span_v: f64,
    pub(super) skew: f64,
    pub(super) decline_reason: &'static str,
}
