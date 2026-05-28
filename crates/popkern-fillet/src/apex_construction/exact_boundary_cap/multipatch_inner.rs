//! Inner pole patch construction for multi-patch exact-boundary candidates.
#![allow(dead_code)]

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;

use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::inventory::{
    ExactBoundaryOuterStripPatchInventory, ExactBoundaryPatchInterfaceInventory,
};
use super::report::ExactBoundaryMultiPatchDiagnostic;
use super::strip::{ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection};

const INTERFACE_SAMPLE_COUNT: usize = 17;

/// Candidate-only strip/inner interface curve from a strip inner boundary.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryStripInnerInterfaceCurve {
    /// Stable semantic interface name.
    pub interface_name: String,
    /// Adjacent strip patch name.
    pub strip_patch_name: String,
    /// Adjacent inner pole patch name reserved for MP-009.
    pub inner_patch_name: String,
    /// Interface 3-D curve copied from strip `v=1`.
    pub curve: NurbsCurve,
    /// Strip-side pcurve for this interface.
    pub strip_pcurve: Pcurve,
    /// Whether the interface parameter reverses source strip direction.
    pub strip_parameter_reversed: bool,
    /// Maximum sampled gap between interface curve and strip inner boundary.
    pub max_strip_boundary_residual: f64,
    /// Sample parameter where the maximum residual was measured.
    pub sample_parameter: f64,
    /// Number of deterministic samples evaluated.
    pub sample_count: usize,
    /// Local explicit G0 tolerance.
    pub tolerance: f64,
}

/// Generate the strip/inner interface from the strip surface inner boundary.
pub(crate) fn construct_strip_inner_interface_curve(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    interface: &ExactBoundaryPatchInterfaceInventory,
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    tolerance: f64,
) -> Result<ExactBoundaryStripInnerInterfaceCurve, ExactBoundaryMultiPatchDiagnostic> {
    if interface.strip_patch_name != strip.patch_name || interface.side != strip.side {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_position(
            sector,
            strip,
            interface,
            f64::INFINITY,
            0.0,
            tolerance,
        ));
    }
    let curve = strip_inner_boundary_curve(candidate);
    let strip_pcurve = strip_inner_boundary_pcurve(candidate, &curve);
    let residual = interface_residual(candidate, &curve, &strip_pcurve);
    if residual.max_residual > tolerance {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_position(
            sector,
            strip,
            interface,
            residual.max_residual,
            residual.sample_parameter,
            tolerance,
        ));
    }

    Ok(ExactBoundaryStripInnerInterfaceCurve {
        interface_name: interface.interface_name.clone(),
        strip_patch_name: interface.strip_patch_name.clone(),
        inner_patch_name: interface.inner_patch_name.clone(),
        curve,
        strip_pcurve,
        strip_parameter_reversed: candidate.boundary_direction
            == ExactBoundaryStripBoundaryDirection::Reversed,
        max_strip_boundary_residual: residual.max_residual,
        sample_parameter: residual.sample_parameter,
        sample_count: INTERFACE_SAMPLE_COUNT,
        tolerance,
    })
}

fn strip_inner_boundary_curve(candidate: &ExactBoundaryOuterStripBoundaryScaffold) -> NurbsCurve {
    let control_points = (0..candidate.surface.num_u)
        .map(|index| candidate.surface.control_points[index * candidate.surface.num_v + 1])
        .collect::<Vec<_>>();
    let weights = (0..candidate.surface.num_u)
        .map(|index| candidate.surface.weights[index * candidate.surface.num_v + 1])
        .collect::<Vec<_>>();
    NurbsCurve {
        degree: candidate.surface.degree_u,
        knots: candidate.surface.knots_u.clone(),
        control_points,
        weights,
    }
}

fn strip_inner_boundary_pcurve(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    curve: &NurbsCurve,
) -> Pcurve {
    let (t_start, t_end) = curve.parameter_range();
    match candidate.boundary_direction {
        ExactBoundaryStripBoundaryDirection::Forward => {
            Pcurve::linear((t_start, 1.0), (t_end, 1.0), t_start, t_end)
        }
        ExactBoundaryStripBoundaryDirection::Reversed => {
            Pcurve::linear((t_end, 1.0), (t_start, 1.0), t_start, t_end)
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct InterfaceResidual {
    max_residual: f64,
    sample_parameter: f64,
}

fn interface_residual(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    curve: &NurbsCurve,
    strip_pcurve: &Pcurve,
) -> InterfaceResidual {
    let (t_start, t_end) = curve.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let uv = strip_pcurve.point_at(parameter);
        let strip_point = candidate.surface.point_at(uv.x, uv.y);
        let curve_point = curve.point_at(uv.x);
        let residual = strip_point.distance_to(&curve_point);
        if residual > report.max_residual {
            report = InterfaceResidual {
                max_residual: residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}

#[cfg(test)]
#[path = "multipatch_inner_tests.rs"]
mod tests;
