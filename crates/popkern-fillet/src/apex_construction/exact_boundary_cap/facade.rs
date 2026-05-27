//! Semantic exact-boundary cap candidate facade.

use std::time::Instant;

use popkern_core::topology::{BRep, FaceHandle};

use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::types::{
    ApexConstructionDiagnostic, ExactBoundaryCandidateOutcome,
    ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapDiagnostic, ExactBoundaryCapReport,
    ExactBoundaryReportDiagnostic, ExactBoundaryReportProvenance, ExactBoundaryReportSector,
    ExactBoundaryToleranceEvidence,
};

use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::commit::{commit_exact_boundary_sector_surfaces, post_commit_validation_diagnostic};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes};
use super::facts::attach_shared_repair_fact_ids;
use super::pcurves::{ExactBoundaryPcurveBudget, prepare_exact_boundary_pcurves};
use super::performance::{
    finalize_performance_report, record_pcurve_report, record_support_family_counts,
    record_validation_report, support_sample_count, surface_control_point_count,
};
use super::ribbon::{ExactBoundaryRibbonBudget, prepare_exact_boundary_ribbons};
use super::sector::{ExactBoundaryCapSectorInput, extract_exact_boundary_cap_sectors};
use super::solve::{ExactBoundarySolveBudget, solve_exact_boundary_surfaces};
use super::topology::materialize_exact_boundary_sector_topology;
use super::validation::{
    ExactBoundaryValidationBudget, ExactBoundaryValidationInputs, validate_exact_boundary_candidate,
};

/// Status of one exact-boundary cap candidate attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryCapCandidateStatus {
    /// Candidate passed validation and committed exact-boundary sector faces.
    Accepted,
    /// Candidate declined and restored the input BRep snapshot.
    Declined,
}

/// Reusable exact-boundary candidate result for construction, diagnostics, and UX.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryCapCandidateResult {
    /// Candidate status after validation/commit or rollback.
    pub status: ExactBoundaryCapCandidateStatus,
    /// Sector ownership records. On accept these reference committed sector faces.
    pub sectors: Vec<ExactBoundaryCapSectorInput>,
    /// Faces committed by the accepted candidate.
    pub committed_faces: Vec<FaceHandle>,
    /// Candidate-only performance counters and budget evidence.
    pub performance: ExactBoundaryCandidatePerformanceReport,
    /// Versioned report for UI and programmatic consumers.
    pub report: ExactBoundaryCapReport,
    /// Structured diagnostics emitted during decline.
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

impl ExactBoundaryCapCandidateResult {
    /// Number of exact curves prepared without approximation.
    pub(crate) fn prepared_curve_count(&self) -> usize {
        self.performance.prepared_curve_count
    }

    /// Number of generated inner spokes prepared without topology mutation.
    pub(crate) fn inner_spoke_count(&self) -> usize {
        self.performance.inner_spoke_count
    }

    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub(crate) fn ribbon_count(&self) -> usize {
        self.performance.support_ribbon_count
    }

    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub(crate) fn solved_surface_count(&self) -> usize {
        self.performance.solved_surface_count
    }

    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub(crate) fn pcurve_pair_count(&self) -> usize {
        self.performance.pcurve_pair_count
    }
}

/// Run the exact-boundary cap candidate transactionally.
///
/// The function mutates `brep` only when all exact-boundary postconditions pass.
/// Any diagnostic or incomplete pcurve inventory restores the input snapshot and
/// returns a declined candidate result for the caller's containment path.
pub(crate) fn run_exact_boundary_cap_candidate(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) -> Result<ExactBoundaryCapCandidateResult, ApexBuildError> {
    let started = Instant::now();
    let snapshot = brep.clone();
    let baseline_validation_errors = snapshot.validate_all();
    let mut candidate_sectors = extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors)
        .map_err(|diagnostic| diagnostic.into_build_error())?;
    let mut performance = ExactBoundaryCandidatePerformanceReport {
        sector_count: candidate_sectors.len(),
        ..Default::default()
    };
    let mut diagnostics = Vec::new();
    let boundary_curves = match prepare_exact_boundary_curves(
        brep,
        &candidate_sectors,
        ExactBoundaryCurveBudget::default(),
    ) {
        Ok(prepared) => prepared,
        Err(diagnostic) => {
            diagnostics.push(diagnostic.into_apex_diagnostic());
            *brep = snapshot;
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
            return Ok(ExactBoundaryCapCandidateResult {
                status: ExactBoundaryCapCandidateStatus::Declined,
                sectors: candidate_sectors,
                committed_faces: Vec::new(),
                performance,
                diagnostics,
                report,
            });
        }
    };
    performance.prepared_curve_count = boundary_curves.curves.len();

    let inner_spokes = match prepare_exact_boundary_inner_spokes(
        &candidate_sectors,
        &boundary_curves,
        ExactBoundaryInnerSpokeBudget::default(),
    ) {
        Ok(prepared) => Some(prepared),
        Err(diagnostic) => {
            diagnostics.push(diagnostic.into_apex_diagnostic());
            None
        }
    };
    performance.inner_spoke_count = inner_spokes
        .as_ref()
        .map(|prepared| prepared.sectors.len() * 2)
        .unwrap_or(0);

    let topology = match inner_spokes.as_ref() {
        Some(inner_spokes) => match materialize_exact_boundary_sector_topology(
            brep,
            cap_face,
            &candidate_sectors,
            inner_spokes,
        ) {
            Ok(topology) => {
                candidate_sectors = topology.sectors.clone();
                Some(topology)
            }
            Err(diagnostic) => {
                diagnostics.push(diagnostic.into_apex_diagnostic());
                None
            }
        },
        None => None,
    };

    let ribbons = match topology.as_ref() {
        Some(_) => match prepare_exact_boundary_ribbons(
            brep,
            &candidate_sectors,
            &boundary_curves,
            ExactBoundaryRibbonBudget::default(),
        ) {
            Ok(prepared) => Some(prepared),
            Err(diagnostic) => {
                diagnostics.push(diagnostic.into_apex_diagnostic());
                None
            }
        },
        None => None,
    };
    performance.support_ribbon_count = ribbons
        .as_ref()
        .map(|prepared| prepared.ribbons.len())
        .unwrap_or(0);
    performance.support_sample_count = ribbons.as_ref().map(support_sample_count).unwrap_or(0);
    if let Some(ribbons) = ribbons.as_ref() {
        record_support_family_counts(&mut performance, ribbons);
    }

    let solved_surfaces = match (inner_spokes.as_ref(), ribbons.as_ref()) {
        (Some(inner_spokes), Some(ribbons)) => {
            match solve_exact_boundary_surfaces(
                inner_spokes,
                ribbons,
                ExactBoundarySolveBudget::default(),
            ) {
                Ok(solved) => Some(solved),
                Err(diagnostic) => {
                    diagnostics.push(diagnostic.into_apex_diagnostic());
                    None
                }
            }
        }
        _ => None,
    };
    performance.solved_surface_count = solved_surfaces
        .as_ref()
        .map(|solved| solved.validated_surface_count())
        .unwrap_or(0);
    performance.surface_control_point_count = solved_surfaces
        .as_ref()
        .map(surface_control_point_count)
        .unwrap_or(0);

    let mut committed = false;
    let mut committed_faces = Vec::new();
    if let (Some(topology), Some(inner_spokes), Some(ribbons), Some(solved_surfaces)) = (
        topology.as_ref(),
        inner_spokes.as_ref(),
        ribbons.as_ref(),
        solved_surfaces.as_ref(),
    ) {
        match prepare_exact_boundary_pcurves(
            inner_spokes,
            ribbons,
            solved_surfaces,
            ExactBoundaryPcurveBudget::default(),
        ) {
            Ok(pcurves) => {
                record_pcurve_report(&mut performance, &pcurves);
                let _topology_evidence = (topology.pole_vertex, topology.spoke_edges.len());
                match validate_exact_boundary_candidate(ExactBoundaryValidationInputs {
                    brep,
                    topology,
                    sectors: &candidate_sectors,
                    spokes: inner_spokes,
                    ribbons,
                    solved: solved_surfaces,
                    pcurves: &pcurves,
                    budget: ExactBoundaryValidationBudget::default(),
                }) {
                    Ok(validation_report) => {
                        record_validation_report(&mut performance, validation_report);
                        committed_faces = candidate_sectors
                            .iter()
                            .map(|sector| sector.cap_face)
                            .collect();
                        if let Err(diagnostic) = commit_exact_boundary_sector_surfaces(
                            brep,
                            cap_face,
                            &candidate_sectors,
                            solved_surfaces,
                            &pcurves,
                        ) {
                            diagnostics.push(diagnostic.into_apex_diagnostic());
                        } else if let Some(diagnostic) = post_commit_validation_diagnostic(
                            &baseline_validation_errors,
                            &brep.validate_all(),
                            brep.tolerance.linear(),
                        ) {
                            diagnostics.push(diagnostic.into_apex_diagnostic());
                        } else {
                            committed = true;
                        }
                    }
                    Err(diagnostic) => diagnostics.push(diagnostic.into_apex_diagnostic()),
                }
            }
            Err(diagnostic) => diagnostics.push(diagnostic.into_apex_diagnostic()),
        }
    }

    let accepted = committed && diagnostics.is_empty() && performance.pcurve_pair_count > 0;
    if !accepted {
        *brep = snapshot;
        committed_faces.clear();
    }
    attach_shared_repair_fact_ids(brep, &mut diagnostics);
    finalize_performance_report(&mut performance, started, &diagnostics);
    let status = if accepted {
        ExactBoundaryCapCandidateStatus::Accepted
    } else {
        ExactBoundaryCapCandidateStatus::Declined
    };
    let report = build_exact_boundary_report(
        brep,
        status,
        &candidate_sectors,
        &committed_faces,
        &performance,
        &diagnostics,
    );

    Ok(ExactBoundaryCapCandidateResult {
        status,
        sectors: candidate_sectors,
        committed_faces,
        performance,
        report,
        diagnostics,
    })
}

pub(super) fn build_exact_boundary_report(
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
        continuity: Default::default(),
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
        cap_half_edges: sector.cap_half_edges,
        support_half_edges: sector.support_half_edges,
        cap_pcurves_present: [
            half_edge_has_pcurve(brep, sector.cap_half_edges[0]),
            half_edge_has_pcurve(brep, sector.cap_half_edges[1]),
        ],
        support_pcurves_present: [
            half_edge_has_pcurve(brep, sector.support_half_edges[0]),
            half_edge_has_pcurve(brep, sector.support_half_edges[1]),
        ],
    }
}

fn half_edge_has_pcurve(brep: &BRep, half_edge: popkern_core::topology::HalfEdgeHandle) -> bool {
    brep.half_edges
        .get(half_edge)
        .is_some_and(|half_edge| half_edge.pcurve.is_some())
}

fn exact_boundary_report_diagnostic(
    diagnostic: &ApexConstructionDiagnostic,
) -> Option<ExactBoundaryReportDiagnostic> {
    let ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) = diagnostic else {
        return None;
    };
    Some(report_diagnostic_from_exact_boundary(diagnostic))
}

fn report_diagnostic_from_exact_boundary(
    diagnostic: &ExactBoundaryCapDiagnostic,
) -> ExactBoundaryReportDiagnostic {
    ExactBoundaryReportDiagnostic {
        schema_version: 1,
        stage: diagnostic.stage.clone(),
        code: diagnostic.code.clone(),
        sector_name: diagnostic.sector_name.clone(),
        edge: diagnostic.edge,
        support_half_edge: diagnostic.support_half_edge,
        support_face: diagnostic.support_face,
        residual: diagnostic.residual,
        tolerance: None,
        shared_fact_ids: diagnostic.shared_fact_ids.iter().map(|id| id.0).collect(),
    }
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
