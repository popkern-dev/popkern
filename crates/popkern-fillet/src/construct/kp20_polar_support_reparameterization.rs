//! Structured diagnostics for targeted KP20 support-side reparameterization.

use popkern_core::geometry::gregory::Kp20PieceId;

#[derive(Debug, Clone)]
pub(super) struct Kp20GuideQualityReport {
    pub(super) samples: usize,
    pub(super) reversals: usize,
    pub(super) min_normal_cos: f64,
    pub(super) worst_sector: usize,
    pub(super) worst_u: f64,
    pub(super) worst_v: f64,
    pub(super) classification: &'static str,
}

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

#[derive(Debug, Clone)]
pub(super) struct Kp20CoupledDomainReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) support_origin_u: f64,
    pub(super) support_origin_v: f64,
    pub(super) support_span_u: f64,
    pub(super) support_span_v: f64,
    pub(super) support_skew: f64,
    pub(super) origin_shift: f64,
    pub(super) origin_scale: f64,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SplitDomainReport {
    pub(super) guide_samples: usize,
    pub(super) guide_reversals: usize,
    pub(super) guide_min_normal_cos: f64,
    pub(super) guide_worst_sector: usize,
    pub(super) guide_worst_u: f64,
    pub(super) guide_worst_v: f64,
    pub(super) guide_classification: &'static str,
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) end_folds: usize,
    pub(super) support_origin_u: f64,
    pub(super) support_origin_v: f64,
    pub(super) support_span_u: f64,
    pub(super) support_span_v: f64,
    pub(super) support_skew: f64,
    pub(super) split_fraction: f64,
    pub(super) interface_pull: f64,
    pub(super) decline_reason: &'static str,
    pub(super) assembly_evaluated: usize,
    pub(super) assembly_accepted: usize,
    pub(super) assembly_rejected: usize,
    pub(super) assembly_start_folds: usize,
    pub(super) assembly_end_folds: usize,
    pub(super) assembly_normal_scale: f64,
    pub(super) assembly_radial_scale: f64,
    pub(super) assembly_decline_reason: &'static str,
}
