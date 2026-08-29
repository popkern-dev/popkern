//! Structured diagnostics for polar K-P bi-6 promotion transactions.

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

pub(super) fn promotion_declined(
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
    tolerance: Option<f64>,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: format!("kp_polar_bi6.promotion.{reason_code}"),
        code: "kp_polar_bi6.promotion.declined".to_string(),
        detail_code: detail,
        sector_name: format!("kp_polar_bi6_ring_n{n_sides}"),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: None,
        outer_side: None,
        sample_parameter: None,
        tolerance,
        shared_fact_ids: Vec::new(),
    })
}

/// Emit a diagnostic trace when the polar K-P promotion rejects.
/// Gated on the `kp-polar` trace topic so production builds stay
/// silent; set `=1` to see which stage rejected and why.
pub(super) fn trace_polar(stage: &str, n: usize, detail: Option<String>) {
    if popkern_core::trace::trace_enabled("kp-polar") {
        match detail {
            Some(d) => eprintln!("[kp_polar_bi6 N={n}] fallback @ {stage}: {d}"),
            None => eprintln!("[kp_polar_bi6 N={n}] fallback @ {stage}"),
        }
    }
}
