//! Informational note constructors for the exact-boundary facade
//! (observability channel records — see PK-201).

use super::super::cross_spoke_g1::CrossSpokeOutcome;

/// Informational record of the committed cap's own internal-seam
/// normal misalignment (PK-197 criterion: measured on every commit).
pub(super) fn mid_seam_note(
    n_faces: usize,
    worst: f64,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "exact_boundary_cap.mid_seam_normal".to_string(),
            code: "exact_boundary_cap.mid_seam_normal".to_string(),
            detail_code: Some(format!("faces={n_faces};worst_mid_seam_normal={worst:.3e}")),
            sector_name: format!("n{n_faces}_corner_cap"),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: Some(worst),
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}

/// Informational record of the cross-spoke correction.
pub(super) fn cross_spoke_note(
    n_sectors: usize,
    outcome: &CrossSpokeOutcome,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "exact_boundary_cap.cross_spoke_g1".to_string(),
            code: "exact_boundary_cap.cross_spoke_g1".to_string(),
            detail_code: Some(format!(
                "spokes={};corrected_control_points={};dimension_mismatched_neighbours={}",
                outcome.spokes, outcome.corrected, outcome.dimension_mismatched
            )),
            sector_name: format!("n{n_sectors}_corner_cap"),
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
