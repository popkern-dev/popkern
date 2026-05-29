//! Candidate-only multipatch route bridge for exact-boundary reports.

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryRouteDecisionStep};

use super::super::boundary::PreparedExactBoundaryCurves;
use super::super::corner::PreparedExactBoundaryInnerSpokes;
use super::super::multipatch::budget::ExactBoundaryMultiPatchBudget;
use super::super::multipatch::report::ExactBoundaryMultiPatchCandidateOutcome;
use super::super::multipatch::strip::ExactBoundaryStripBoundaryDirection;
use super::super::multipatch::{
    ExactBoundaryMultiPatchCandidate, ExactBoundaryMultiPatchCandidateRequest,
    prepare_multipatch_candidate,
};
use super::super::ribbon::PreparedExactBoundaryRibbons;
use super::super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::super::validation::{
    ExactBoundaryValidationDeclineReason, ExactBoundaryValidationDiagnostic,
};

pub(super) struct ExactBoundaryMultipatchRouteAttempt {
    pub(super) step: ExactBoundaryRouteDecisionStep,
    pub(super) sector: ExactBoundaryCapSectorInput,
    pub(super) candidate: Option<ExactBoundaryMultiPatchCandidate>,
}

pub(super) fn attempt_multipatch_after_one_patch_decline(
    diagnostic: &ExactBoundaryValidationDiagnostic,
    sectors: &[ExactBoundaryCapSectorInput],
    boundary_curves: &PreparedExactBoundaryCurves,
    ribbons: &PreparedExactBoundaryRibbons,
    inner_spokes: &PreparedExactBoundaryInnerSpokes,
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) -> Option<ExactBoundaryMultipatchRouteAttempt> {
    if diagnostic.reason != ExactBoundaryValidationDeclineReason::ConstrainedStripBandUnsatisfied {
        return None;
    }
    let side = diagnostic_side(diagnostic)?;
    let sector = sectors
        .iter()
        .find(|sector| diagnostic.sector_name.as_deref() == Some(sector.sector_name.as_str()))?;
    let boundary = boundary_curves.curves.iter().find(|boundary| {
        boundary.sector_name == sector.sector_name
            && boundary.side == side
            && Some(boundary.edge) == diagnostic.edge
    })?;
    let ribbon = ribbons.ribbons.iter().find(|ribbon| {
        ribbon.sector_name == sector.sector_name
            && ribbon.side == side
            && Some(ribbon.support_half_edge) == diagnostic.support_half_edge
            && Some(ribbon.support_face) == diagnostic.support_face
    })?;
    let selected_sides = [side];
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector,
        selected_sides: &selected_sides,
        boundary: &boundary.curve,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &ribbon.samples,
        pole: inner_spokes.pole,
        budget: ExactBoundaryMultiPatchBudget::default(),
    });
    match attempt.report.outcome {
        ExactBoundaryMultiPatchCandidateOutcome::CommitEligible => {
            Some(ExactBoundaryMultipatchRouteAttempt {
                step: ExactBoundaryRouteDecisionStep::commit_eligible("exact_boundary.multipatch"),
                sector: sector.clone(),
                candidate: attempt.candidate,
            })
        }
        ExactBoundaryMultiPatchCandidateOutcome::Declined => {
            let diagnostic_code = attempt
                .report
                .diagnostics
                .first()
                .map(|diagnostic| diagnostic.code.clone());
            diagnostics.extend(
                attempt
                    .report
                    .diagnostics
                    .into_iter()
                    .map(|diagnostic| diagnostic.into_apex_diagnostic()),
            );
            Some(ExactBoundaryMultipatchRouteAttempt {
                step: ExactBoundaryRouteDecisionStep::declined(
                    "exact_boundary.multipatch",
                    diagnostic_code,
                ),
                sector: sector.clone(),
                candidate: None,
            })
        }
        ExactBoundaryMultiPatchCandidateOutcome::Pending => None,
    }
}

fn diagnostic_side(
    diagnostic: &ExactBoundaryValidationDiagnostic,
) -> Option<ExactBoundaryOuterSide> {
    match diagnostic.outer_side.as_deref() {
        Some("FilletEnd") => Some(ExactBoundaryOuterSide::FilletEnd),
        Some("FaceEdge") => Some(ExactBoundaryOuterSide::FaceEdge),
        _ => None,
    }
}
