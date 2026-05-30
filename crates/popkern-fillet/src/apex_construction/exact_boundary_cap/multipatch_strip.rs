//! Outer support strip construction for multi-patch exact-boundary candidates.
#![allow(dead_code, clippy::result_large_err)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Vec3};

use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::budget::ExactBoundaryMultiPatchBudget;
use super::inventory::ExactBoundaryOuterStripPatchInventory;
use super::report::ExactBoundaryMultiPatchDiagnostic;
#[path = "multipatch_strip_support.rs"]
mod strip_support;
use strip_support::{
    boundary_tangent_compatibility, control_row_parameter, endpoint_fallback_width_direction,
    is_endpoint_sample, sample_at_parameter, strip_width_from_samples, support_refined_boundary,
    width_direction, width_direction_from_tangent,
};

const OUTER_BOUNDARY_SAMPLE_COUNT: usize = 9;
const MIN_STRIP_NORMAL_LENGTH: f64 = 1.0e-12;
const STRIP_ROW_COUNT: usize = 3;
const SUPPORT_NORMAL_ROW: usize = 1;
const INNER_INTERFACE_ROW: usize = 2;

/// Candidate strip direction along the copied exact-boundary curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryStripBoundaryDirection {
    /// Candidate strip `u` follows the source boundary parameter.
    Forward,
    /// Candidate strip `u` reverses the source boundary parameter.
    Reversed,
}

/// Candidate-only outer strip scaffold before support-normal row solving.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryOuterStripBoundaryScaffold {
    /// Stable semantic patch name.
    pub patch_name: String,
    /// Candidate strip surface with exact copied outer boundary row.
    pub surface: NurbsSurface,
    /// Candidate cap-side pcurve on the strip surface.
    pub outer_pcurve: Pcurve,
    /// Direction from source boundary parameter to candidate strip parameter.
    pub boundary_direction: ExactBoundaryStripBoundaryDirection,
    /// Maximum sampled residual between copied boundary row and source curve.
    pub max_outer_boundary_residual: f64,
    /// Number of samples used for the residual gate.
    pub outer_boundary_sample_count: usize,
    /// Local explicit residual tolerance.
    pub tolerance: f64,
    /// Whether the source boundary used non-unit rational weights.
    pub rational_boundary: bool,
    /// Whether the second row is intentionally a zero-width placeholder.
    pub support_normal_row_pending: bool,
    /// Support-normal solve report after MP-006 row construction.
    pub support_normal_report: Option<ExactBoundaryStripSupportNormalReport>,
}

/// Bounded evidence from the support-normal row solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryStripSupportNormalReport {
    /// Width used for the first strip row.
    pub strip_width: f64,
    /// Worst support-normal deficit measured on support ribbon samples.
    pub max_support_normal_deficit: f64,
    /// Sample parameter where the worst deficit was measured.
    pub sample_parameter: f64,
    /// Number of support strip samples evaluated.
    pub support_strip_sample_count: usize,
    /// Bounded control-row iterations consumed by this solve.
    pub bounded_iteration_count: usize,
    /// Local explicit support-normal tolerance.
    pub tolerance: f64,
}

pub(crate) fn scaffold_outer_strip_boundary(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    boundary: &NurbsCurve,
    boundary_direction: ExactBoundaryStripBoundaryDirection,
    tolerance: f64,
) -> Result<ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryMultiPatchDiagnostic> {
    if !boundary_layout_is_valid(boundary) {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::pcurve_parameterization_invalid(
                sector, strip, 0.0, 0.0, tolerance,
            ),
        );
    }
    if rational_weights_are_unsupported(boundary) {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::rational_boundary_unsupported(sector, strip),
        );
    }

    let surface = zero_width_strip_surface(boundary);
    let outer_pcurve = outer_boundary_pcurve(boundary, boundary_direction);
    let residual_report = outer_boundary_residual(&surface, boundary, &outer_pcurve);
    if residual_report.max_residual > tolerance {
        return Err(ExactBoundaryMultiPatchDiagnostic::outer_boundary_position(
            sector,
            strip,
            residual_report.max_residual,
            residual_report.sample_parameter,
            tolerance,
        ));
    }

    Ok(ExactBoundaryOuterStripBoundaryScaffold {
        patch_name: strip.patch_name.clone(),
        surface,
        outer_pcurve,
        boundary_direction,
        max_outer_boundary_residual: residual_report.max_residual,
        outer_boundary_sample_count: OUTER_BOUNDARY_SAMPLE_COUNT,
        tolerance,
        rational_boundary: boundary
            .weights
            .iter()
            .any(|weight| (*weight - 1.0).abs() > 1.0e-12),
        support_normal_row_pending: true,
        support_normal_report: None,
    })
}

pub(crate) fn solve_strip_support_normal_row(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    scaffold: &ExactBoundaryOuterStripBoundaryScaffold,
    boundary: &NurbsCurve,
    samples: &[ExactBoundaryRibbonSample],
    budget: ExactBoundaryMultiPatchBudget,
) -> Result<ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryMultiPatchDiagnostic> {
    let Some(strip_width) = strip_width_from_samples(samples) else {
        return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
            sector,
            strip,
            f64::INFINITY,
            0.0,
            budget.support_normal_deficit,
            Vec::new(),
        ));
    };
    let compatibility = boundary_tangent_compatibility(boundary, samples);
    if compatibility.min_reachable_deficit > budget.support_normal_deficit {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::boundary_tangent_incompatible(
                sector,
                strip,
                compatibility.min_reachable_deficit,
                compatibility.sample_parameter,
                budget.support_normal_deficit,
            ),
        );
    }
    let refined_boundary = match support_refined_boundary(boundary, samples) {
        Some(refined_boundary) => refined_boundary,
        None => {
            return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
                sector,
                strip,
                f64::INFINITY,
                0.0,
                budget.support_normal_deficit,
                Vec::new(),
            ));
        }
    };
    let mut surface = zero_width_strip_surface(&refined_boundary);
    if surface.num_u > budget.max_width_direction_controls {
        return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
            sector,
            strip,
            f64::INFINITY,
            0.0,
            budget.support_normal_deficit,
            Vec::new(),
        ));
    }
    let bounded_iteration_count = surface.num_u * surface.num_v.saturating_sub(1);
    for index in 0..surface.num_u {
        let parameter = control_row_parameter(&refined_boundary, index);
        let sample = sample_at_parameter(samples, parameter);
        let Some(direction) = width_direction(boundary, parameter, sample.support_normal)
            .or_else(|| {
                width_direction_from_tangent(sample.boundary_tangent, sample.support_normal)
            })
            .or_else(|| sample.cross_derivative.normalized())
            .or_else(|| endpoint_fallback_width_direction(samples, sample.parameter))
        else {
            return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
                sector,
                strip,
                f64::INFINITY,
                sample.parameter,
                budget.support_normal_deficit,
                Vec::new(),
            ));
        };
        let row0 = surface.control_points[control_index(index, 0, surface.num_v)];
        surface.control_points[control_index(index, SUPPORT_NORMAL_ROW, surface.num_v)] =
            row0 + direction * (strip_width * 0.5);
        surface.control_points[control_index(index, INNER_INTERFACE_ROW, surface.num_v)] =
            row0 + direction * strip_width;
    }

    let report = support_normal_report(&surface, scaffold, samples, budget.support_normal_deficit);
    if report.max_support_normal_deficit > budget.support_normal_deficit {
        return Err(ExactBoundaryMultiPatchDiagnostic::support_strip_band(
            sector,
            strip,
            report.max_support_normal_deficit,
            report.sample_parameter,
            budget.support_normal_deficit,
            Vec::new(),
        ));
    }

    Ok(ExactBoundaryOuterStripBoundaryScaffold {
        surface,
        support_normal_row_pending: false,
        support_normal_report: Some(ExactBoundaryStripSupportNormalReport {
            strip_width,
            bounded_iteration_count,
            ..report
        }),
        ..scaffold.clone()
    })
}

fn boundary_layout_is_valid(boundary: &NurbsCurve) -> bool {
    let control_count = boundary.control_points.len();
    if control_count < 2
        || boundary.weights.len() != control_count
        || boundary.knots.len() != control_count + boundary.degree + 1
        || boundary.degree + 1 > control_count
    {
        return false;
    }
    if !boundary
        .knots
        .windows(2)
        .all(|pair| pair[0].is_finite() && pair[1].is_finite() && pair[0] <= pair[1])
    {
        return false;
    }
    boundary
        .weights
        .iter()
        .all(|weight| *weight > 0.0 && weight.is_finite())
}

fn rational_weights_are_unsupported(boundary: &NurbsCurve) -> bool {
    boundary
        .weights
        .iter()
        .any(|weight| (*weight - 1.0).abs() > 1.0e-12)
        && boundary.weights.len() != boundary.control_points.len()
}

fn zero_width_strip_surface(boundary: &NurbsCurve) -> NurbsSurface {
    let mut control_points = Vec::with_capacity(boundary.control_points.len() * STRIP_ROW_COUNT);
    let mut weights = Vec::with_capacity(boundary.weights.len() * STRIP_ROW_COUNT);
    for (point, weight) in boundary.control_points.iter().zip(&boundary.weights) {
        for _ in 0..STRIP_ROW_COUNT {
            control_points.push(*point);
            weights.push(*weight);
        }
    }
    NurbsSurface {
        degree_u: boundary.degree,
        degree_v: 2,
        knots_u: boundary.knots.clone(),
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        num_u: boundary.control_points.len(),
        num_v: STRIP_ROW_COUNT,
        weights,
    }
}

fn control_index(u_index: usize, v_index: usize, num_v: usize) -> usize {
    u_index * num_v + v_index
}

fn outer_boundary_pcurve(
    boundary: &NurbsCurve,
    boundary_direction: ExactBoundaryStripBoundaryDirection,
) -> Pcurve {
    let (t_start, t_end) = boundary.parameter_range();
    match boundary_direction {
        ExactBoundaryStripBoundaryDirection::Forward => {
            Pcurve::linear((t_start, 0.0), (t_end, 0.0), t_start, t_end)
        }
        ExactBoundaryStripBoundaryDirection::Reversed => {
            Pcurve::linear((t_end, 0.0), (t_start, 0.0), t_start, t_end)
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct OuterBoundaryResidualReport {
    max_residual: f64,
    sample_parameter: f64,
}

fn outer_boundary_residual(
    surface: &NurbsSurface,
    boundary: &NurbsCurve,
    pcurve: &Pcurve,
) -> OuterBoundaryResidualReport {
    let (t_start, t_end) = boundary.parameter_range();
    let mut report = OuterBoundaryResidualReport {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..OUTER_BOUNDARY_SAMPLE_COUNT {
        let alpha = index as f64 / (OUTER_BOUNDARY_SAMPLE_COUNT - 1) as f64;
        let source_t = t_start + (t_end - t_start) * alpha;
        let uv = pcurve.point_at(source_t);
        let strip_point = surface.point_at(uv.x, uv.y);
        let boundary_point = boundary.point_at(uv.x);
        let residual = strip_point.distance_to(&boundary_point);
        if residual > report.max_residual {
            report = OuterBoundaryResidualReport {
                max_residual: residual,
                sample_parameter: source_t,
            };
        }
    }
    report
}

fn support_normal_report(
    surface: &NurbsSurface,
    scaffold: &ExactBoundaryOuterStripBoundaryScaffold,
    samples: &[ExactBoundaryRibbonSample],
    tolerance: f64,
) -> ExactBoundaryStripSupportNormalReport {
    let mut max_support_normal_deficit = 0.0;
    let mut sample_parameter = samples
        .first()
        .map(|sample| sample.parameter)
        .unwrap_or(0.0);
    for sample in samples {
        let uv = scaffold.outer_pcurve.point_at(sample.parameter);
        let Some(deficit) = support_normal_deficit_at(surface, uv.x, uv.y, sample.support_normal)
        else {
            if is_endpoint_sample(sample.parameter) {
                continue;
            }
            let deficit = f64::INFINITY;
            if deficit > max_support_normal_deficit {
                max_support_normal_deficit = deficit;
                sample_parameter = sample.parameter;
            }
            continue;
        };
        if deficit > max_support_normal_deficit {
            max_support_normal_deficit = deficit;
            sample_parameter = sample.parameter;
        }
    }
    ExactBoundaryStripSupportNormalReport {
        strip_width: 0.0,
        max_support_normal_deficit,
        sample_parameter,
        support_strip_sample_count: samples.len(),
        bounded_iteration_count: 0,
        tolerance,
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
    if cap_length < MIN_STRIP_NORMAL_LENGTH
        || support_length < MIN_STRIP_NORMAL_LENGTH
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
#[path = "multipatch_strip_tests.rs"]
mod tests;
