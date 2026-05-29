//! Candidate-only multipatch route bridge for exact-boundary reports.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryRouteDecisionStep};

use super::super::boundary::{PreparedExactBoundaryCurve, PreparedExactBoundaryCurves};
use super::super::commit::commit_exact_boundary_sector_surfaces;
use super::super::corner::PreparedExactBoundaryInnerSpokes;
use super::super::multipatch::budget::ExactBoundaryMultiPatchBudget;
use super::super::multipatch::report::ExactBoundaryMultiPatchCandidateOutcome;
use super::super::multipatch::strip::ExactBoundaryStripBoundaryDirection;
use super::super::multipatch::{
    ExactBoundaryMultiPatchCandidate, ExactBoundaryMultiPatchCandidateRequest,
    prepare_multipatch_candidate,
};
use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::ribbon::{ExactBoundaryRibbonSample, PreparedExactBoundaryRibbons};
use super::super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::super::solve::ExactBoundarySolvedSurfaces;
use super::super::validation::{
    ExactBoundaryValidationDeclineReason, ExactBoundaryValidationDiagnostic,
};

pub(super) struct ExactBoundaryMultipatchRouteAttempt {
    pub(super) step: ExactBoundaryRouteDecisionStep,
    pub(super) sector: ExactBoundaryCapSectorInput,
    pub(super) candidate: Option<ExactBoundaryMultiPatchCandidate>,
}

pub(super) fn collect_live_one_patch_cap_support_failures(
    brep: &BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
) -> Vec<ExactBoundaryValidationDiagnostic> {
    let mut probe = brep.clone();
    if commit_exact_boundary_sector_surfaces(
        &mut probe,
        cap_face,
        sectors,
        solved_surfaces,
        pcurves,
    )
    .is_err()
    {
        return Vec::new();
    }
    let Some(diagnostic) =
        cap_support_continuity_diagnostic(&probe, &sector_faces(sectors), 5.0e-3, 5.0e-2)
    else {
        return Vec::new();
    };
    let Some(pair) = pcurve_pair_for_cap_support_diagnostic(pcurves, &diagnostic) else {
        return Vec::new();
    };
    vec![ExactBoundaryValidationDiagnostic::constrained_strip_band(
        pair,
        diagnostic_continuity_residual(&diagnostic),
        0.5,
    )]
}

pub(super) fn attempt_multipatches_after_one_patch_decline(
    diagnostics_to_repair: &[ExactBoundaryValidationDiagnostic],
    sectors: &[ExactBoundaryCapSectorInput],
    boundary_curves: &PreparedExactBoundaryCurves,
    ribbons: &PreparedExactBoundaryRibbons,
    pcurves: &PreparedExactBoundaryPcurves,
    inner_spokes: &PreparedExactBoundaryInnerSpokes,
    diagnostics: &mut Vec<ApexConstructionDiagnostic>,
) -> Vec<ExactBoundaryMultipatchRouteAttempt> {
    let mut attempts = Vec::new();
    for diagnostic in diagnostics_to_repair {
        let Some(attempt) = attempt_multipatch_after_one_patch_decline(
            diagnostic,
            sectors,
            boundary_curves,
            ribbons,
            pcurves,
            inner_spokes,
            diagnostics,
        ) else {
            continue;
        };
        if attempts
            .iter()
            .any(|existing: &ExactBoundaryMultipatchRouteAttempt| {
                existing.sector.sector_name == attempt.sector.sector_name
            })
        {
            attempts.push(ExactBoundaryMultipatchRouteAttempt {
                step: ExactBoundaryRouteDecisionStep::declined(
                    "exact_boundary.multipatch",
                    Some("exact_boundary_cap.multipatch.multiple_sides_required".to_string()),
                ),
                sector: attempt.sector,
                candidate: None,
            });
        } else {
            attempts.push(attempt);
        }
    }
    attempts
}

pub(super) fn attempt_multipatch_after_one_patch_decline(
    diagnostic: &ExactBoundaryValidationDiagnostic,
    sectors: &[ExactBoundaryCapSectorInput],
    boundary_curves: &PreparedExactBoundaryCurves,
    ribbons: &PreparedExactBoundaryRibbons,
    pcurves: &PreparedExactBoundaryPcurves,
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
    let support_samples =
        dense_support_samples(boundary, pcurves, ribbon).unwrap_or_else(|| ribbon.samples.clone());
    let selected_sides = [side];
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector,
        selected_sides: &selected_sides,
        boundary: &boundary.curve,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples,
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

fn sector_faces(sectors: &[ExactBoundaryCapSectorInput]) -> Vec<FaceHandle> {
    sectors.iter().map(|sector| sector.cap_face).collect()
}

fn pcurve_pair_for_cap_support_diagnostic<'a>(
    pcurves: &'a PreparedExactBoundaryPcurves,
    diagnostic: &ApexConstructionDiagnostic,
) -> Option<&'a PreparedExactBoundaryPcurvePair> {
    let ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) = diagnostic else {
        return None;
    };
    pcurves
        .pairs
        .iter()
        .find(|pair| pair.edge == diagnostic.edge)
}

fn diagnostic_continuity_residual(diagnostic: &ApexConstructionDiagnostic) -> f64 {
    match diagnostic {
        ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => diagnostic
            .continuity
            .worst_g1_normal_deficit
            .unwrap_or(diagnostic.worst_g1_coplanarity),
        _ => f64::INFINITY,
    }
}

fn dense_support_samples(
    boundary: &PreparedExactBoundaryCurve,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbon: &super::super::ribbon::PreparedExactBoundaryRibbon,
) -> Option<Vec<ExactBoundaryRibbonSample>> {
    let pair = pcurves.pairs.iter().find(|pair| {
        pair.sector_name == ribbon.sector_name
            && pair.side == ribbon.side
            && pair.edge == ribbon.edge
    })?;
    let mut samples = ribbon.samples.clone();
    for index in 0..64 {
        let parameter = index as f64 / 63.0;
        if samples
            .iter()
            .any(|sample| (sample.parameter - parameter).abs() <= 1.0e-12)
        {
            continue;
        }
        if let Some(sample) = dense_support_sample(boundary, pair, ribbon, parameter) {
            samples.push(sample);
        }
    }
    samples.sort_by(|a, b| {
        a.parameter
            .partial_cmp(&b.parameter)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Some(samples)
}

fn dense_support_sample(
    boundary: &PreparedExactBoundaryCurve,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &super::super::ribbon::PreparedExactBoundaryRibbon,
    parameter: f64,
) -> Option<ExactBoundaryRibbonSample> {
    let support_t = pair.support_pcurve.t_start
        + (pair.support_pcurve.t_end - pair.support_pcurve.t_start) * parameter;
    let support_uv = pair.support_pcurve.point_at(support_t);
    let support_point = pair.support_surface.point_at(support_uv.x, support_uv.y);
    let support_normal = pair
        .support_surface
        .normal_at(support_uv.x, support_uv.y)
        .normalized()?;
    let boundary_point = boundary.curve.point_at(parameter);
    let boundary_tangent = boundary.curve.derivative_at(parameter, 1).normalized()?;
    let nearest = nearest_support_sample(&ribbon.samples, parameter)?;
    Some(ExactBoundaryRibbonSample {
        parameter,
        support_uv: (support_uv.x, support_uv.y),
        boundary_point,
        support_point,
        boundary_tangent,
        support_normal,
        cross_derivative: nearest.cross_derivative,
    })
}

fn nearest_support_sample(
    samples: &[ExactBoundaryRibbonSample],
    parameter: f64,
) -> Option<ExactBoundaryRibbonSample> {
    samples.iter().copied().min_by(|a, b| {
        (a.parameter - parameter)
            .abs()
            .partial_cmp(&(b.parameter - parameter).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    })
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
