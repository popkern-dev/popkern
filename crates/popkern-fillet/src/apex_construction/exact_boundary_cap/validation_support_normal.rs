//! Support-normal continuity validation for exact-boundary cap candidates.
#![allow(clippy::result_large_err)]

use popkern_core::geometry::surface::SurfaceEval;

use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::ribbon::{ExactBoundaryRibbonSample, PreparedExactBoundaryRibbons};
use super::super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};
use super::{
    ExactBoundaryValidationBudget, ExactBoundaryValidationDeclineReason,
    ExactBoundaryValidationDiagnostic,
};

pub(super) fn validate_support_normal_continuity(
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
    budget: ExactBoundaryValidationBudget,
) -> Result<(), ExactBoundaryValidationDiagnostic> {
    if let Some(diagnostic) =
        collect_constrained_strip_band_diagnostics(solved, pcurves, ribbons, budget)
            .into_iter()
            .next()
    {
        return Err(diagnostic);
    }
    Ok(())
}

pub(crate) fn collect_constrained_strip_band_diagnostics(
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
    budget: ExactBoundaryValidationBudget,
) -> Vec<ExactBoundaryValidationDiagnostic> {
    let mut diagnostics = Vec::new();
    for pair in &pcurves.pairs {
        let Some(solved_surface) = solved.surfaces.iter().find(|surface| {
            surface.sector_index == pair.sector_index && surface.sector_name == pair.sector_name
        }) else {
            diagnostics.push(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::SurfaceInventoryMismatch,
                Some(pair.sector_name.clone()),
                None,
            ));
            continue;
        };
        let Some(ribbon) = ribbons.ribbons.iter().find(|ribbon| {
            ribbon.sector_index == pair.sector_index
                && ribbon.side == pair.side
                && ribbon.edge == pair.edge
                && ribbon.support_half_edge == pair.support_half_edge
        }) else {
            diagnostics.push(ExactBoundaryValidationDiagnostic::new(
                ExactBoundaryValidationDeclineReason::PcurveInventoryMismatch,
                Some(pair.sector_name.clone()),
                None,
            ));
            continue;
        };
        let measured = ribbon
            .samples
            .iter()
            .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
            .map(|sample| {
                support_normal_deficit(solved_surface, pair, sample, budget)
                    .map(|deficit| (deficit, sample.parameter))
            })
            .try_fold((0.0_f64, 0.0_f64), |worst, deficit| {
                deficit.map(|deficit| if deficit.0 > worst.0 { deficit } else { worst })
            });
        let (worst_deficit, worst_parameter) = match measured {
            Ok(measured) => measured,
            Err(diagnostic) => {
                diagnostics.push(diagnostic);
                continue;
            }
        };
        if worst_deficit > budget.max_support_normal_deficit {
            diagnostics.push(ExactBoundaryValidationDiagnostic::constrained_strip_band(
                pair,
                worst_deficit,
                worst_parameter,
            ));
        }
    }
    diagnostics
}

fn support_normal_deficit(
    solved_surface: &ExactBoundarySectorSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    sample: &ExactBoundaryRibbonSample,
    budget: ExactBoundaryValidationBudget,
) -> Result<f64, ExactBoundaryValidationDiagnostic> {
    let cap_uv = pair.cap_pcurve.point_at(
        pair.cap_pcurve.t_start
            + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
    );
    let (normal_u, normal_v) =
        interior_normal_sample_parameters(&solved_surface.surface, cap_uv.x, cap_uv.y);
    let cap_normal = solved_surface
        .surface
        .partial_u(normal_u, normal_v)
        .cross(&solved_surface.surface.partial_v(normal_u, normal_v));
    let support_normal = sample.support_normal;
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length < budget.min_normal_length
        || support_length < budget.min_normal_length
        || !cap_length.is_finite()
        || !support_length.is_finite()
    {
        return Err(ExactBoundaryValidationDiagnostic::constrained_strip_band(
            pair,
            f64::INFINITY,
            sample.parameter,
        ));
    }
    Ok(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

fn interior_normal_sample_parameters(
    surface: &popkern_core::geometry::NurbsSurface,
    u: f64,
    v: f64,
) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    (
        interior_parameter(u, u_min, u_max),
        interior_parameter(v, v_min, v_max),
    )
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}
