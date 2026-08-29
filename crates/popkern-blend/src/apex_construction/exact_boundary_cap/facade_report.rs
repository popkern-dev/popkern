//! Exact-boundary facade report construction.

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::types::{
    ApexConstructionDiagnostic, CapContinuityReport, CapSupportContinuityDiagnostic,
    ExactBoundaryCandidateOutcome, ExactBoundaryCandidatePerformanceReport,
    ExactBoundaryCapDiagnostic, ExactBoundaryCapReport, ExactBoundaryReportDiagnostic,
    ExactBoundaryReportProvenance, ExactBoundaryReportSector, ExactBoundaryRouteDecisionStep,
    ExactBoundaryToleranceEvidence,
};

use super::super::sector::ExactBoundaryCapSectorInput;
use super::ExactBoundaryCapCandidateStatus;

/// Decline exit for a candidate whose boundary-curve preparation
/// failed before any construction: restores the snapshot and builds
/// the declined result with the same report/diagnostic ordering
/// rules as the main path (informational after the decision,
/// fact-ids attached, performance finalized).
#[allow(clippy::too_many_arguments)]
pub(super) fn declined_curve_preparation_result(
    brep: &mut BRep,
    snapshot: BRep,
    mut diagnostics: Vec<ApexConstructionDiagnostic>,
    informational: Vec<ApexConstructionDiagnostic>,
    candidate_sectors: Vec<ExactBoundaryCapSectorInput>,
    mut performance: ExactBoundaryCandidatePerformanceReport,
    started: std::time::Instant,
) -> super::ExactBoundaryCapCandidateResult {
    use super::super::facts::attach_shared_repair_fact_ids;
    use super::super::performance::finalize_performance_report;

    *brep = snapshot;
    diagnostics.extend(informational);
    attach_shared_repair_fact_ids(brep, &mut diagnostics);
    finalize_performance_report(&mut performance, started, &diagnostics);
    let report = build_exact_boundary_report(
        brep,
        ExactBoundaryCapCandidateStatus::Declined,
        &candidate_sectors,
        &[],
        &performance,
        &diagnostics,
    );
    super::ExactBoundaryCapCandidateResult {
        status: ExactBoundaryCapCandidateStatus::Declined,
        sectors: candidate_sectors,
        committed_faces: Vec::new(),
        performance,
        diagnostics,
        report,
        pole: None,
    }
}

pub(crate) fn build_exact_boundary_report(
    brep: &BRep,
    status: ExactBoundaryCapCandidateStatus,
    sectors: &[ExactBoundaryCapSectorInput],
    committed_faces: &[FaceHandle],
    performance: &ExactBoundaryCandidatePerformanceReport,
    diagnostics: &[ApexConstructionDiagnostic],
) -> ExactBoundaryCapReport {
    let report_diagnostics: Vec<_> = diagnostics
        .iter()
        .filter_map(exact_boundary_report_diagnostic)
        .collect();
    let source_fact_ids = unique_source_fact_ids(&report_diagnostics);
    let candidate_status = match status {
        ExactBoundaryCapCandidateStatus::Accepted => ExactBoundaryCandidateOutcome::Accepted,
        ExactBoundaryCapCandidateStatus::Declined => ExactBoundaryCandidateOutcome::Declined,
    };
    ExactBoundaryCapReport {
        schema_version: 1,
        route_label: "direct_apex".to_string(),
        candidate_status,
        decision_chain: vec![match candidate_status {
            ExactBoundaryCandidateOutcome::Accepted => {
                ExactBoundaryRouteDecisionStep::accepted("exact_boundary.one_patch")
            }
            ExactBoundaryCandidateOutcome::Declined => ExactBoundaryRouteDecisionStep::declined(
                "exact_boundary.one_patch",
                performance.deterministic_skip_code.clone(),
            ),
            ExactBoundaryCandidateOutcome::NotApplicable => {
                ExactBoundaryRouteDecisionStep::not_attempted("exact_boundary.one_patch")
            }
        }],
        provenance: ExactBoundaryReportProvenance::for_status(candidate_status),
        sectors: sectors
            .iter()
            .map(|sector| exact_boundary_report_sector(brep, sector))
            .collect(),
        committed_faces: committed_faces.to_vec(),
        preview_available: !sectors.is_empty(),
        tolerance_evidence: ExactBoundaryToleranceEvidence {
            model_linear_tolerance: Some(brep.tolerance.linear()),
            max_reported_residual: report_diagnostics
                .iter()
                .filter_map(|diagnostic| diagnostic.residual)
                .reduce(f64::max),
            deterministic_skip_code: performance.deterministic_skip_code.clone(),
            ..ExactBoundaryToleranceEvidence::default()
        },
        continuity: report_continuity(diagnostics),
        performance: performance.clone(),
        diagnostics: report_diagnostics,
        source_fact_ids,
    }
}

fn exact_boundary_report_sector(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
) -> ExactBoundaryReportSector {
    ExactBoundaryReportSector {
        sector_index: sector.sector_index,
        sector_name: sector.sector_name.clone(),
        fillet_index: sector.fillet_index,
        cap_face: sector.cap_face,
        fillet_end_edge: sector.fillet_end_edge,
        face_edge: sector.face_edge,
        fillet_support_face: sector.fillet_support_face,
        bridge_support_face: sector.bridge_support_face,
        cap_half_edges: [
            Some(sector.fillet_seam.cap_half_edge),
            sector.bridge_seam.as_ref().map(|seam| seam.cap_half_edge),
        ],
        support_half_edges: [
            Some(sector.fillet_seam.support_half_edge),
            sector
                .bridge_seam
                .as_ref()
                .map(|seam| seam.support_half_edge),
        ],
        cap_pcurves_present: [
            half_edge_has_pcurve(brep, sector.fillet_seam.cap_half_edge),
            sector
                .bridge_seam
                .as_ref()
                .is_some_and(|seam| half_edge_has_pcurve(brep, seam.cap_half_edge)),
        ],
        support_pcurves_present: [
            half_edge_has_pcurve(brep, sector.fillet_seam.support_half_edge),
            sector
                .bridge_seam
                .as_ref()
                .is_some_and(|seam| half_edge_has_pcurve(brep, seam.support_half_edge)),
        ],
    }
}

fn half_edge_has_pcurve(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    brep.half_edges
        .get(half_edge)
        .is_some_and(|half_edge| half_edge.pcurve.is_some())
}

fn exact_boundary_report_diagnostic(
    diagnostic: &ApexConstructionDiagnostic,
) -> Option<ExactBoundaryReportDiagnostic> {
    match diagnostic {
        ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) => {
            Some(report_diagnostic_from_exact_boundary(diagnostic))
        }
        ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => {
            Some(report_diagnostic_from_cap_support(diagnostic))
        }
        _ => None,
    }
}

fn report_diagnostic_from_exact_boundary(
    diagnostic: &ExactBoundaryCapDiagnostic,
) -> ExactBoundaryReportDiagnostic {
    ExactBoundaryReportDiagnostic {
        schema_version: 1,
        stage: diagnostic.stage.clone(),
        code: diagnostic.code.clone(),
        detail_code: diagnostic.detail_code.clone(),
        sector_name: diagnostic.sector_name.clone(),
        edge: diagnostic.edge,
        support_half_edge: diagnostic.support_half_edge,
        support_face: diagnostic.support_face,
        residual: diagnostic.residual,
        outer_side: diagnostic.outer_side.clone(),
        sample_parameter: diagnostic.sample_parameter,
        tolerance: diagnostic.tolerance,
        shared_fact_ids: diagnostic.shared_fact_ids.iter().map(|id| id.0).collect(),
    }
}

fn report_diagnostic_from_cap_support(
    diagnostic: &CapSupportContinuityDiagnostic,
) -> ExactBoundaryReportDiagnostic {
    let g0_declined =
        diagnostic.continuity.decline_code.as_deref() == Some("cap_outer_g0_exceeds_gate");
    ExactBoundaryReportDiagnostic {
        schema_version: 1,
        stage: "exact_boundary_cap.cap_support_continuity".to_string(),
        code: diagnostic
            .continuity
            .decline_code
            .clone()
            .unwrap_or_else(|| "cap_outer_g1_exceeds_gate".to_string()),
        detail_code: Some(cap_support_detail_code(diagnostic)),
        sector_name: String::new(),
        edge: Some(diagnostic.edge),
        support_half_edge: None,
        support_face: Some(diagnostic.support_face),
        residual: if g0_declined {
            diagnostic.continuity.worst_g0_gap
        } else {
            diagnostic.continuity.worst_g1_normal_deficit
        },
        outer_side: None,
        sample_parameter: diagnostic.sample_parameter,
        tolerance: Some(if g0_declined {
            diagnostic.continuity.target.g0_position_tolerance
        } else {
            diagnostic.continuity.target.g1_normal_deficit_tolerance
        }),
        shared_fact_ids: Vec::new(),
    }
}

fn cap_support_detail_code(diagnostic: &CapSupportContinuityDiagnostic) -> String {
    format!(
        "{};cap_uv={:?};support_uv={:?};cap_curve_gap={:?};support_curve_gap={:?}",
        diagnostic.measurement_scope,
        diagnostic.cap_uv,
        diagnostic.support_uv,
        diagnostic.cap_curve_gap,
        diagnostic.support_curve_gap
    )
}

fn report_continuity(diagnostics: &[ApexConstructionDiagnostic]) -> CapContinuityReport {
    diagnostics
        .iter()
        .find_map(|diagnostic| match diagnostic {
            ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => {
                Some(diagnostic.continuity.clone())
            }
            _ => None,
        })
        .unwrap_or_default()
}

fn unique_source_fact_ids(diagnostics: &[ExactBoundaryReportDiagnostic]) -> Vec<usize> {
    let mut fact_ids = Vec::new();
    for fact_id in diagnostics
        .iter()
        .flat_map(|diagnostic| diagnostic.shared_fact_ids.iter().copied())
    {
        if !fact_ids.contains(&fact_id) {
            fact_ids.push(fact_id);
        }
    }
    fact_ids
}
