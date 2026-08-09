//! Candidate-only multipatch route bridge for exact-boundary reports.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryRouteDecisionStep};

#[path = "facade_multipatch_support.rs"]
mod support;

use super::super::boundary::PreparedExactBoundaryCurves;
use super::super::commit::commit_exact_boundary_sector_surfaces;
use super::super::corner::PreparedExactBoundaryInnerSpokes;
use super::super::multipatch::budget::ExactBoundaryMultiPatchBudget;
use super::super::multipatch::inventory::ExactBoundaryMultiPatchInventory;
use super::super::multipatch::report::{
    ExactBoundaryMultiPatchCandidateOutcome, ExactBoundaryMultiPatchDiagnostic,
};
use super::super::multipatch::strip::ExactBoundaryStripBoundaryDirection;
use super::super::multipatch::{
    ExactBoundaryMultiPatchCandidate, ExactBoundaryMultiPatchCandidateRequest,
    prepare_multipatch_candidate,
};
use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::ribbon::PreparedExactBoundaryRibbons;
use super::super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::super::solve::ExactBoundarySolvedSurfaces;
use super::super::validation::{
    ExactBoundaryValidationDeclineReason, ExactBoundaryValidationDiagnostic,
};
use support::{
    dense_support_samples, optimize_selected_strip_live_support_seam,
    selected_strip_live_equivalent_diagnostic,
    support_authority_boundary_attempt_with_control_limit,
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
    if cap_support_continuity_diagnostic(&probe, &sector_faces(sectors), 5.0e-3, 5.0e-2).is_none() {
        return Vec::new();
    }
    pcurves
        .pairs
        .iter()
        .filter_map(|pair| {
            let cap_surface = solved_surfaces
                .surfaces
                .iter()
                .find(|surface| surface.sector_name == pair.sector_name)?;
            let report = pcurve_pair_cap_support_g1_report(cap_surface, pair)?;
            (report.residual > 5.0e-2).then(|| {
                ExactBoundaryValidationDiagnostic::constrained_strip_band(
                    pair,
                    report.residual,
                    report.sample_parameter,
                )
            })
        })
        .collect()
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
    let support_pair = pcurve_pair_for_ribbon(pcurves, ribbon)?;
    let budget = ExactBoundaryMultiPatchBudget::default();
    let reauthor_attempt = support_authority_boundary_attempt_with_control_limit(
        boundary,
        support_pair,
        budget.max_width_direction_controls,
    );
    if let Some(decline) = reauthor_attempt.decline {
        let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(sector, side);
        let strip = inventory.outer_strips.first()?;
        diagnostics.push(
            ExactBoundaryMultiPatchDiagnostic::support_boundary_reauthor_decline(
                decline.reason,
                sector,
                strip,
                decline.residual,
                decline.sample_parameter,
                decline.tolerance,
                decline.detail_code,
            )
            .into_apex_diagnostic(),
        );
    }
    let candidate_boundary = reauthor_attempt
        .candidate
        .unwrap_or_else(|| boundary.curve.clone());
    let support_samples = dense_support_samples(&candidate_boundary, support_pair, ribbon);
    let selected_sides = [side];
    let attempt = prepare_multipatch_candidate(ExactBoundaryMultiPatchCandidateRequest {
        sector,
        selected_sides: &selected_sides,
        boundary: &candidate_boundary,
        boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
        support_samples: &support_samples,
        support_authority_boundary: (candidate_boundary != boundary.curve)
            .then_some(&candidate_boundary),
        pole: inner_spokes.pole,
        budget,
    });
    match attempt.report.outcome {
        ExactBoundaryMultiPatchCandidateOutcome::CommitEligible => {
            let mut candidate = attempt.candidate?;
            optimize_selected_strip_live_support_seam(
                &mut candidate,
                &candidate_boundary,
                support_pair,
                budget,
            );
            if let Some(diagnostic) = selected_strip_live_equivalent_diagnostic(
                sector,
                &candidate,
                &candidate_boundary,
                support_pair,
                budget,
            ) {
                let diagnostic_code = diagnostic.code.clone();
                diagnostics.push(diagnostic.into_apex_diagnostic());
                return Some(ExactBoundaryMultipatchRouteAttempt {
                    step: ExactBoundaryRouteDecisionStep::declined(
                        "exact_boundary.multipatch",
                        Some(diagnostic_code),
                    ),
                    sector: sector.clone(),
                    candidate: None,
                });
            }
            Some(ExactBoundaryMultipatchRouteAttempt {
                step: ExactBoundaryRouteDecisionStep::commit_eligible("exact_boundary.multipatch"),
                sector: sector.clone(),
                candidate: Some(candidate),
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

#[derive(Debug, Clone, Copy)]
struct PairG1Report {
    residual: f64,
    sample_parameter: f64,
}

fn pcurve_pair_cap_support_g1_report(
    cap_surface: &super::super::solve::ExactBoundarySectorSurface,
    pair: &PreparedExactBoundaryPcurvePair,
) -> Option<PairG1Report> {
    let (cap_start_t, cap_end_t) = pair.cap_pcurve.parameter_range();
    let (support_start_t, support_end_t) = pair.support_pcurve.parameter_range();
    let cap_start = pcurve_surface_point(&cap_surface.surface, &pair.cap_pcurve, cap_start_t);
    let cap_end = pcurve_surface_point(&cap_surface.surface, &pair.cap_pcurve, cap_end_t);
    let support_start =
        pcurve_surface_point(&pair.support_surface, &pair.support_pcurve, support_start_t);
    let support_end =
        pcurve_surface_point(&pair.support_surface, &pair.support_pcurve, support_end_t);
    let support_forward = cap_start.distance_to(&support_start) + cap_end.distance_to(&support_end)
        <= cap_start.distance_to(&support_end) + cap_end.distance_to(&support_start);
    let mut report = PairG1Report {
        residual: 0.0,
        sample_parameter: 0.0,
    };
    for index in 0..100 {
        let parameter = 0.05 + 0.90 * index as f64 / 99.0;
        let cap_uv = pair
            .cap_pcurve
            .point_at(cap_start_t + (cap_end_t - cap_start_t) * parameter);
        let support_parameter = if support_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let support_uv = pair
            .support_pcurve
            .point_at(support_start_t + (support_end_t - support_start_t) * support_parameter);
        let residual = interior_normal_deficit_between_surfaces(
            &cap_surface.surface,
            cap_uv.x,
            cap_uv.y,
            &pair.support_surface,
            support_uv.x,
            support_uv.y,
        )?;
        if residual > report.residual {
            report = PairG1Report {
                residual,
                sample_parameter: parameter,
            };
        }
    }
    Some(report)
}

fn pcurve_surface_point(
    surface: &impl SurfaceEval,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

fn interior_normal_deficit_between_surfaces(
    cap_surface: &impl SurfaceEval,
    cap_u: f64,
    cap_v: f64,
    support_surface: &impl SurfaceEval,
    support_u: f64,
    support_v: f64,
) -> Option<f64> {
    let ((cap_u_min, cap_u_max), (cap_v_min, cap_v_max)) = cap_surface.parameter_range();
    let ((support_u_min, support_u_max), (support_v_min, support_v_max)) =
        support_surface.parameter_range();
    let cap_u = interior_parameter(cap_u, cap_u_min, cap_u_max);
    let cap_v = interior_parameter(cap_v, cap_v_min, cap_v_max);
    let support_u = interior_parameter(support_u, support_u_min, support_u_max);
    let support_v = interior_parameter(support_v, support_v_min, support_v_max);
    let cap_normal = cap_surface
        .partial_u(cap_u, cap_v)
        .cross(&cap_surface.partial_v(cap_u, cap_v));
    let support_normal = support_surface
        .partial_u(support_u, support_v)
        .cross(&support_surface.partial_v(support_u, support_v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length <= 1.0e-300 || support_length <= 1.0e-300 {
        return None;
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}

fn pcurve_pair_for_ribbon<'a>(
    pcurves: &'a PreparedExactBoundaryPcurves,
    ribbon: &super::super::ribbon::PreparedExactBoundaryRibbon,
) -> Option<&'a PreparedExactBoundaryPcurvePair> {
    pcurves.pairs.iter().find(|pair| {
        pair.sector_name == ribbon.sector_name
            && pair.side == ribbon.side
            && pair.edge == ribbon.edge
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

/// Trace-gated observability for the multipatch repair ladder: the
/// per-sector/side repair diagnostics and each attempt's decision
/// step (measured 2026-08-09: this is how the (1,1,1) two-sided
/// constrained-band class was localized).
pub(super) fn trace_multipatch_attempts(
    repair_diagnostics: &[crate::apex_construction::exact_boundary_cap::validation::ExactBoundaryValidationDiagnostic],
    attempts: &[ExactBoundaryMultipatchRouteAttempt],
) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    for repair_diagnostic in repair_diagnostics {
        eprintln!(
            "[direct-apex-trace] repair diagnostic: sector={:?} side={:?} code={} residual={:?}",
            repair_diagnostic.sector_name,
            repair_diagnostic.outer_side,
            repair_diagnostic.code,
            repair_diagnostic.residual,
        );
    }
    for attempt in attempts {
        eprintln!(
            "[direct-apex-trace] multipatch attempt: sector={} candidate={} step={:?}",
            attempt.sector.sector_name,
            attempt.candidate.is_some(),
            attempt.step,
        );
    }
}
