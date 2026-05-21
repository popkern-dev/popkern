//! Semantic exact-boundary cap candidate facade.

use std::time::Instant;

use popkern_core::topology::{BRep, FaceHandle};

use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCandidatePerformanceReport};

use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::commit::commit_exact_boundary_sector_surfaces;
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
use super::validation::{ExactBoundaryValidationBudget, validate_exact_boundary_candidate};

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
    let mut candidate_sectors = extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors)
        .map_err(|diagnostic| diagnostic.into_build_error())?;
    let mut performance = ExactBoundaryCandidatePerformanceReport {
        sector_count: candidate_sectors.len(),
        ..Default::default()
    };
    let boundary_curves = prepare_exact_boundary_curves(
        brep,
        &candidate_sectors,
        ExactBoundaryCurveBudget::default(),
    )
    .map_err(|diagnostic| diagnostic.into_build_error())?;
    performance.prepared_curve_count = boundary_curves.curves.len();
    let mut diagnostics = Vec::new();

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

    let ribbons = match prepare_exact_boundary_ribbons(
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
    if let (Some(inner_spokes), Some(ribbons), Some(solved_surfaces)) = (
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
                match materialize_exact_boundary_sector_topology(
                    brep,
                    cap_face,
                    &candidate_sectors,
                    inner_spokes,
                ) {
                    Ok(topology) => {
                        let _topology_evidence = (topology.pole_vertex, topology.spoke_edges.len());
                        match validate_exact_boundary_candidate(
                            brep,
                            &topology,
                            &topology.sectors,
                            inner_spokes,
                            solved_surfaces,
                            &pcurves,
                            ExactBoundaryValidationBudget::default(),
                        ) {
                            Ok(validation_report) => {
                                record_validation_report(&mut performance, validation_report);
                                candidate_sectors = topology.sectors;
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

    Ok(ExactBoundaryCapCandidateResult {
        status: if accepted {
            ExactBoundaryCapCandidateStatus::Accepted
        } else {
            ExactBoundaryCapCandidateStatus::Declined
        },
        sectors: candidate_sectors,
        committed_faces,
        performance,
        diagnostics,
    })
}
