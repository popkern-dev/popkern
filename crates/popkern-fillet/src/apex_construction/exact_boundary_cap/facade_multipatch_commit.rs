//! Live commit bridge for exact-boundary multipatch route attempts.

use popkern_core::error::ValidationError;
use popkern_core::topology::{BRep, FaceHandle};

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryRouteDecisionStep};

use super::super::commit::{
    commit_exact_boundary_sector_surfaces, post_commit_validation_diagnostic,
};
use super::super::multipatch::ExactBoundaryMultiPatchCandidate;
use super::super::multipatch::commit::commit_multipatch_candidate_transactionally;
use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::sector::ExactBoundaryCapSectorInput;
use super::super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};
use super::multipatch_route::ExactBoundaryMultipatchRouteAttempt;
use crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic;

pub(super) enum ExactBoundaryMultipatchCommitOutcome {
    Accepted {
        faces: Vec<FaceHandle>,
    },
    Declined {
        step: ExactBoundaryRouteDecisionStep,
        diagnostic: ApexConstructionDiagnostic,
    },
}

pub(super) fn commit_multipatch_route_transactionally(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    baseline_validation_errors: &[ValidationError],
    attempt: &ExactBoundaryMultipatchRouteAttempt,
    candidate: &ExactBoundaryMultiPatchCandidate,
) -> ExactBoundaryMultipatchCommitOutcome {
    let multipatch_commit =
        match commit_multipatch_candidate_transactionally(brep, &attempt.sector, candidate) {
            Ok(report) => report,
            Err(diagnostic) => {
                return declined(diagnostic.code.clone(), diagnostic.into_apex_diagnostic());
            }
        };
    if let Err(diagnostic) = commit_unaffected_one_patch_sectors(
        brep,
        cap_face,
        sectors,
        solved_surfaces,
        pcurves,
        attempt,
    ) {
        return declined(
            diagnostic.code.to_string(),
            diagnostic.into_apex_diagnostic(),
        );
    }
    if let Some(diagnostic) = post_commit_validation_diagnostic(
        baseline_validation_errors,
        &brep.validate_all(),
        brep.tolerance.linear(),
    ) {
        return declined(
            diagnostic.code.to_string(),
            diagnostic.into_apex_diagnostic(),
        );
    }
    let mut faces = sectors
        .iter()
        .map(|sector| sector.cap_face)
        .collect::<Vec<_>>();
    faces.push(multipatch_commit.inner_face);
    if let Some(diagnostic) = cap_support_continuity_diagnostic(brep, &faces, 5.0e-3, 5.0e-2) {
        return declined(apex_diagnostic_code(&diagnostic), diagnostic);
    }
    ExactBoundaryMultipatchCommitOutcome::Accepted { faces }
}

fn commit_unaffected_one_patch_sectors(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    attempt: &ExactBoundaryMultipatchRouteAttempt,
) -> Result<(), super::super::commit::ExactBoundaryCommitDiagnostic> {
    let mut remaining_sectors = sectors
        .iter()
        .filter(|sector| sector.sector_name != attempt.sector.sector_name)
        .cloned()
        .collect::<Vec<_>>();
    for (index, sector) in remaining_sectors.iter_mut().enumerate() {
        sector.sector_index = index;
    }
    if remaining_sectors.is_empty() {
        return Ok(());
    }
    let remaining_solved = ExactBoundarySolvedSurfaces {
        surfaces: filter_surfaces_for_sectors(&solved_surfaces.surfaces, &remaining_sectors),
    };
    let remaining_pcurves = PreparedExactBoundaryPcurves {
        pairs: filter_pcurves_for_sectors(&pcurves.pairs, &remaining_sectors),
    };
    commit_exact_boundary_sector_surfaces(
        brep,
        cap_face,
        &remaining_sectors,
        &remaining_solved,
        &remaining_pcurves,
    )
    .map(|_| ())
}

fn filter_surfaces_for_sectors(
    surfaces: &[ExactBoundarySectorSurface],
    sectors: &[ExactBoundaryCapSectorInput],
) -> Vec<ExactBoundarySectorSurface> {
    let mut filtered = surfaces
        .iter()
        .filter(|surface| {
            sectors
                .iter()
                .any(|sector| sector.sector_name == surface.sector_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    for surface in &mut filtered {
        if let Some(index) = sectors
            .iter()
            .position(|sector| sector.sector_name == surface.sector_name)
        {
            surface.sector_index = index;
        }
    }
    filtered
}

fn filter_pcurves_for_sectors(
    pairs: &[PreparedExactBoundaryPcurvePair],
    sectors: &[ExactBoundaryCapSectorInput],
) -> Vec<PreparedExactBoundaryPcurvePair> {
    let mut filtered = pairs
        .iter()
        .filter(|pair| {
            sectors
                .iter()
                .any(|sector| sector.sector_name == pair.sector_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    for pair in &mut filtered {
        if let Some(index) = sectors
            .iter()
            .position(|sector| sector.sector_name == pair.sector_name)
        {
            pair.sector_index = index;
        }
    }
    filtered
}

fn declined(
    code: String,
    diagnostic: ApexConstructionDiagnostic,
) -> ExactBoundaryMultipatchCommitOutcome {
    ExactBoundaryMultipatchCommitOutcome::Declined {
        step: ExactBoundaryRouteDecisionStep::declined("exact_boundary.multipatch", Some(code)),
        diagnostic,
    }
}

fn apex_diagnostic_code(diagnostic: &ApexConstructionDiagnostic) -> String {
    match diagnostic {
        ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) => diagnostic.code.clone(),
        ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => diagnostic
            .continuity
            .decline_code
            .clone()
            .unwrap_or_else(|| "cap_outer_g1_exceeds_gate".to_string()),
        _ => "apex_construction_declined".to_string(),
    }
}
