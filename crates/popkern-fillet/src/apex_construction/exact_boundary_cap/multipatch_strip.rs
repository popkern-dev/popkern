//! Outer support strip construction for multi-patch exact-boundary candidates.
#![allow(dead_code)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface};

use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::inventory::ExactBoundaryOuterStripPatchInventory;
use super::report::ExactBoundaryMultiPatchDiagnostic;

const OUTER_BOUNDARY_SAMPLE_COUNT: usize = 9;

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
}

/// Build a candidate strip surface whose outer row exactly reuses the source
/// NURBS boundary layout. The second row is a zero-width placeholder for MP-006.
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
    let mut control_points = Vec::with_capacity(boundary.control_points.len() * 2);
    let mut weights = Vec::with_capacity(boundary.weights.len() * 2);
    for (point, weight) in boundary.control_points.iter().zip(&boundary.weights) {
        control_points.push(*point);
        control_points.push(*point);
        weights.push(*weight);
        weights.push(*weight);
    }
    NurbsSurface {
        degree_u: boundary.degree,
        degree_v: 1,
        knots_u: boundary.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: boundary.control_points.len(),
        num_v: 2,
        weights,
    }
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

#[cfg(test)]
#[path = "multipatch_strip_tests.rs"]
mod tests;
