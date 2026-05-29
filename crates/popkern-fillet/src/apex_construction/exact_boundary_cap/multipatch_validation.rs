//! Validation gates for multi-patch exact-boundary cap candidates.
#![allow(dead_code, clippy::result_large_err)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Vec3};

use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::budget::ExactBoundaryMultiPatchBudget;
use super::inner::{ExactBoundaryInnerPolePatchPrototype, ExactBoundaryStripInnerInterfaceCurve};
use super::inventory::{
    ExactBoundaryMultiPatchInventoryDiagnostic, ExactBoundaryOuterStripPatchInventory,
    ExactBoundaryPatchInterfaceInventory,
};
use super::report::{ExactBoundaryMultiPatchDiagnostic, ExactBoundaryMultiPatchValidationOrder};
use super::strip::{ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection};

const UNIFORM_SAMPLE_COUNT: usize = 32;
const THRESHOLD_REFINEMENT_RATIO: f64 = 0.8;

/// Validated strip-only candidate evidence before inner-patch construction.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryOuterStripValidationReport {
    /// Candidate strip patch name.
    pub patch_name: String,
    /// Deterministic gate order evaluated by the validator.
    pub validation_order: Vec<ExactBoundaryMultiPatchValidationOrder>,
    /// Worst exact outer-boundary G0 residual.
    pub max_outer_boundary_residual: f64,
    /// Worst support-normal deficit.
    pub max_support_normal_deficit: f64,
    /// Minimum sampled normal length across the strip.
    pub min_normal_length: f64,
    /// Number of deterministic boundary samples evaluated.
    pub boundary_sample_count: usize,
    /// Number of support-ribbon samples evaluated.
    pub support_strip_sample_count: usize,
    /// Number of near-threshold midpoint refinements evaluated.
    pub midpoint_refinement_count: usize,
}

/// Validated inner pole G0 evidence before commit eligibility.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryInnerPoleValidationReport {
    /// Candidate inner patch name.
    pub patch_name: String,
    /// Adjacent interface name.
    pub interface_name: String,
    /// Worst interface G0 residual.
    pub max_interface_g0: f64,
    /// Worst pole/spoke G0 residual.
    pub max_pole_g0: f64,
    /// Worst measured interface G1 deficit.
    pub max_interface_g1_deficit: f64,
    /// Whether interface G1 has been validated by a later slice.
    pub interface_g1_validated: bool,
}

/// Validate one candidate strip without mutating topology or routing.
pub(crate) fn validate_outer_strip_candidate(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    inventory_diagnostic: Option<ExactBoundaryMultiPatchInventoryDiagnostic>,
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
    support_samples: &[ExactBoundaryRibbonSample],
    budget: ExactBoundaryMultiPatchBudget,
) -> Result<ExactBoundaryOuterStripValidationReport, ExactBoundaryMultiPatchDiagnostic> {
    let mut validation_order = Vec::new();
    validation_order.push(ExactBoundaryMultiPatchValidationOrder::Inventory);
    if let Some(diagnostic) = inventory_diagnostic {
        return Err(ExactBoundaryMultiPatchDiagnostic::inventory_mismatch(
            diagnostic,
        ));
    }

    validation_order.push(ExactBoundaryMultiPatchValidationOrder::SourceParameterization);
    if let Some((residual, parameter)) = pcurve_parameterization_defect(candidate, boundary) {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::pcurve_parameterization_invalid(
                sector,
                strip,
                residual,
                parameter,
                budget.outer_boundary_residual,
            ),
        );
    }

    validation_order.push(ExactBoundaryMultiPatchValidationOrder::OuterBoundaryPosition);
    let boundary_report = exact_outer_boundary_report(candidate, boundary, support_samples, budget);
    if boundary_report.max_residual > budget.outer_boundary_residual {
        return Err(ExactBoundaryMultiPatchDiagnostic::outer_boundary_position(
            sector,
            strip,
            boundary_report.max_residual,
            boundary_report.sample_parameter,
            budget.outer_boundary_residual,
        ));
    }

    validation_order.push(ExactBoundaryMultiPatchValidationOrder::SupportStripBand);
    let support_report = support_strip_band_report(candidate, support_samples);
    if support_report.max_deficit > budget.support_normal_deficit {
        return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
            sector,
            strip,
            support_report.max_deficit,
            support_report.sample_parameter,
            budget.support_normal_deficit,
            Vec::new(),
        ));
    }

    validation_order.push(ExactBoundaryMultiPatchValidationOrder::Foldover);
    let foldover_report = foldover_report(candidate, boundary);
    if foldover_report.min_normal_length < budget.min_normal_length {
        return Err(ExactBoundaryMultiPatchDiagnostic::patch_foldover(
            sector,
            strip,
            foldover_report.min_normal_length,
            foldover_report.sample_parameter,
            budget.min_normal_length,
        ));
    }

    validation_order.push(ExactBoundaryMultiPatchValidationOrder::ReportConsistency);
    Ok(ExactBoundaryOuterStripValidationReport {
        patch_name: candidate.patch_name.clone(),
        validation_order,
        max_outer_boundary_residual: boundary_report.max_residual,
        max_support_normal_deficit: support_report.max_deficit,
        min_normal_length: foldover_report.min_normal_length,
        boundary_sample_count: boundary_report.sample_count,
        support_strip_sample_count: support_report.sample_count,
        midpoint_refinement_count: boundary_report.midpoint_refinement_count,
    })
}

/// Validate inner pole prototype G0 and decline G0-only candidates before commit.
pub(crate) fn validate_inner_pole_patch_prototype(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    interface_inventory: &ExactBoundaryPatchInterfaceInventory,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    budget: ExactBoundaryMultiPatchBudget,
) -> Result<ExactBoundaryInnerPoleValidationReport, ExactBoundaryMultiPatchDiagnostic> {
    if prototype.max_interface_g0 > budget.interface_g0 {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_position(
            sector,
            strip,
            interface_inventory,
            prototype.max_interface_g0,
            interface.sample_parameter,
            budget.interface_g0,
        ));
    }
    if prototype.max_pole_g0 > budget.pole_g0 {
        return Err(ExactBoundaryMultiPatchDiagnostic::pole_position(
            sector,
            strip,
            interface_inventory,
            prototype.max_pole_g0,
            0.0,
            budget.pole_g0,
        ));
    }
    let interface_g1_report = prototype.interface_g1_report;
    if let Some(report) = interface_g1_report
        && report.max_interface_g1_deficit > budget.interface_g1_deficit
    {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_tangent(
            sector,
            strip,
            interface_inventory,
            report.max_interface_g1_deficit,
            report.sample_parameter,
            budget.interface_g1_deficit,
        ));
    }
    if !prototype.interface_g1_validated || interface_g1_report.is_none() {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_tangent(
            sector,
            strip,
            interface_inventory,
            interface_g1_report
                .map(|report| report.max_interface_g1_deficit)
                .unwrap_or(f64::INFINITY),
            interface_g1_report
                .map(|report| report.sample_parameter)
                .unwrap_or(interface.sample_parameter),
            budget.interface_g1_deficit,
        ));
    }
    Ok(ExactBoundaryInnerPoleValidationReport {
        patch_name: prototype.patch_name.clone(),
        interface_name: prototype.interface_name.clone(),
        max_interface_g0: prototype.max_interface_g0,
        max_pole_g0: prototype.max_pole_g0,
        max_interface_g1_deficit: interface_g1_report
            .map(|report| report.max_interface_g1_deficit)
            .unwrap_or(f64::INFINITY),
        interface_g1_validated: prototype.interface_g1_validated,
    })
}

fn pcurve_parameterization_defect(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
) -> Option<(f64, f64)> {
    let (t_start, t_end) = boundary.parameter_range();
    if !candidate.outer_pcurve.t_start.is_finite()
        || !candidate.outer_pcurve.t_end.is_finite()
        || candidate.outer_pcurve.t_end <= candidate.outer_pcurve.t_start
    {
        return Some((f64::INFINITY, t_start));
    }
    let start_u = candidate.outer_pcurve.point_at(t_start).x;
    let mid = 0.5 * (t_start + t_end);
    let mid_u = candidate.outer_pcurve.point_at(mid).x;
    let end_u = candidate.outer_pcurve.point_at(t_end).x;
    let range_tolerance = (t_end - t_start).abs().max(1.0) * 1.0e-9;
    if (end_u - start_u).abs() <= range_tolerance {
        return Some((0.0, mid));
    }
    let monotone = match candidate.boundary_direction {
        ExactBoundaryStripBoundaryDirection::Forward => start_u <= mid_u && mid_u <= end_u,
        ExactBoundaryStripBoundaryDirection::Reversed => start_u >= mid_u && mid_u >= end_u,
    };
    if !monotone
        || ![start_u, mid_u, end_u]
            .iter()
            .all(|value| value.is_finite())
    {
        return Some(((end_u - start_u).abs(), mid));
    }
    if start_u.min(end_u) < t_start - range_tolerance
        || start_u.max(end_u) > t_end + range_tolerance
    {
        return Some(((end_u - start_u).abs(), mid));
    }
    None
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct BoundaryValidationReport {
    max_residual: f64,
    sample_parameter: f64,
    sample_count: usize,
    midpoint_refinement_count: usize,
}

fn exact_outer_boundary_report(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
    support_samples: &[ExactBoundaryRibbonSample],
    budget: ExactBoundaryMultiPatchBudget,
) -> BoundaryValidationReport {
    let mut parameters = deterministic_boundary_parameters(boundary, support_samples);
    let mut report = boundary_residual_report(candidate, boundary, &parameters, 0);
    if report.max_residual > budget.outer_boundary_residual * THRESHOLD_REFINEMENT_RATIO
        && let Some(midpoint) = midpoint_near_sample(&parameters, report.sample_parameter)
    {
        parameters.push(midpoint);
        parameters.sort_by(|a, b| a.total_cmp(b));
        parameters.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-12);
        report = boundary_residual_report(candidate, boundary, &parameters, 1);
    }
    report
}

fn deterministic_boundary_parameters(
    boundary: &NurbsCurve,
    support_samples: &[ExactBoundaryRibbonSample],
) -> Vec<f64> {
    let (t_start, t_end) = boundary.parameter_range();
    let mut parameters = vec![t_start, t_end];
    parameters.extend(
        boundary
            .knots
            .iter()
            .copied()
            .filter(|knot| *knot >= t_start && *knot <= t_end),
    );
    parameters.extend(support_samples.iter().map(|sample| sample.parameter));
    for index in 0..UNIFORM_SAMPLE_COUNT {
        let alpha = index as f64 / (UNIFORM_SAMPLE_COUNT - 1) as f64;
        parameters.push(t_start + (t_end - t_start) * alpha);
    }
    parameters.sort_by(|a, b| a.total_cmp(b));
    parameters.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-12);
    parameters
}

fn boundary_residual_report(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
    parameters: &[f64],
    midpoint_refinement_count: usize,
) -> BoundaryValidationReport {
    let mut max_residual = 0.0;
    let mut sample_parameter = parameters.first().copied().unwrap_or(0.0);
    for parameter in parameters {
        let uv = candidate.outer_pcurve.point_at(*parameter);
        let strip_point = candidate.surface.point_at(uv.x, uv.y);
        let boundary_point = boundary.point_at(uv.x);
        let residual = strip_point.distance_to(&boundary_point);
        if residual > max_residual {
            max_residual = residual;
            sample_parameter = *parameter;
        }
    }
    BoundaryValidationReport {
        max_residual,
        sample_parameter,
        sample_count: parameters.len(),
        midpoint_refinement_count,
    }
}

fn midpoint_near_sample(parameters: &[f64], sample_parameter: f64) -> Option<f64> {
    let index = parameters
        .iter()
        .position(|parameter| (*parameter - sample_parameter).abs() <= 1.0e-12)?;
    let neighbor = if index + 1 < parameters.len() {
        parameters[index + 1]
    } else if index > 0 {
        parameters[index - 1]
    } else {
        return None;
    };
    Some(0.5 * (sample_parameter + neighbor))
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct SupportStripBandReport {
    max_deficit: f64,
    sample_parameter: f64,
    sample_count: usize,
}

fn support_strip_band_report(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    support_samples: &[ExactBoundaryRibbonSample],
) -> SupportStripBandReport {
    let mut max_deficit = 0.0;
    let mut sample_parameter = support_samples
        .first()
        .map(|sample| sample.parameter)
        .unwrap_or(0.0);
    for sample in support_samples {
        let uv = candidate.outer_pcurve.point_at(sample.parameter);
        let Some(deficit) =
            support_normal_deficit_at(&candidate.surface, uv.x, uv.y, sample.support_normal)
        else {
            if is_endpoint_sample(sample.parameter) {
                continue;
            }
            let deficit = f64::INFINITY;
            if deficit > max_deficit {
                max_deficit = deficit;
                sample_parameter = sample.parameter;
            }
            continue;
        };
        if deficit > max_deficit {
            max_deficit = deficit;
            sample_parameter = sample.parameter;
        }
    }
    SupportStripBandReport {
        max_deficit,
        sample_parameter,
        sample_count: support_samples.len(),
    }
}

fn is_endpoint_sample(parameter: f64) -> bool {
    parameter <= 1.0e-12 || (1.0 - parameter).abs() <= 1.0e-12
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct FoldoverReport {
    min_normal_length: f64,
    sample_parameter: f64,
    orientation_flip: bool,
}

fn foldover_report(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
) -> FoldoverReport {
    let parameters = deterministic_boundary_parameters(boundary, &[]);
    let mut min_normal_length = f64::INFINITY;
    let mut sample_parameter = parameters.first().copied().unwrap_or(0.0);
    let mut reference_normal = None;
    let mut orientation_flip = false;
    for parameter in parameters {
        if is_endpoint_sample(parameter) {
            continue;
        }
        for v in [0.25, 0.5, 0.75] {
            let uv = candidate.outer_pcurve.point_at(parameter);
            let normal = candidate
                .surface
                .partial_u(uv.x, v)
                .cross(&candidate.surface.partial_v(uv.x, v));
            let length = normal.length();
            if length < min_normal_length {
                min_normal_length = length;
                sample_parameter = parameter;
            }
            if length.is_finite() && length > 1.0e-12 {
                let unit = normal * (1.0 / length);
                if let Some(reference) = reference_normal {
                    if unit.dot(&reference) < -1.0e-6 {
                        orientation_flip = true;
                    }
                } else {
                    reference_normal = Some(unit);
                }
            }
        }
    }
    FoldoverReport {
        min_normal_length,
        sample_parameter,
        orientation_flip,
    }
}

fn support_normal_deficit_at(
    surface: &NurbsSurface,
    u: f64,
    v: f64,
    support_normal: Vec3,
) -> Option<f64> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u = interior_parameter(u, u_min, u_max);
    let v = interior_parameter(v, v_min, v_max);
    let cap_normal = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length < 1.0e-12
        || support_length < 1.0e-12
        || !cap_length.is_finite()
        || !support_length.is_finite()
    {
        return None;
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

fn interior_parameter(parameter: f64, min: f64, max: f64) -> f64 {
    let span = max - min;
    if span.abs() <= 1.0e-12 {
        return parameter;
    }
    let inset = span.abs() * 1.0e-6;
    parameter.clamp(min + inset, max - inset)
}

#[cfg(test)]
#[path = "multipatch_validation_tests.rs"]
mod tests;
