//! Structured diagnostics for KP20 polar macro-cap dispatch.

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

pub(super) fn dispatch_declined(
    reason_code: &str,
    n_sides: usize,
    detail: Option<String>,
    residual: Option<f64>,
    tolerance: Option<f64>,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: format!("kp20_polar_macro.dispatch.{reason_code}"),
        code: "kp20_polar_macro.dispatch.declined".to_string(),
        detail_code: detail,
        sector_name: format!("kp20_polar_macro_ring_n{n_sides}"),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual,
        outer_side: None,
        sample_parameter: None,
        tolerance,
        shared_fact_ids: Vec::new(),
    })
}
