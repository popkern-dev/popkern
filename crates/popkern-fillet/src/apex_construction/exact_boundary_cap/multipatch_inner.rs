//! Inner pole patch construction for multi-patch exact-boundary candidates.
#![allow(dead_code, clippy::result_large_err)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::inventory::{
    ExactBoundaryOuterStripPatchInventory, ExactBoundaryPatchInterfaceInventory,
};
use super::report::ExactBoundaryMultiPatchDiagnostic;
use super::strip::{ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection};

const INTERFACE_SAMPLE_COUNT: usize = 17;
const POLE_SAMPLE_COUNT: usize = 9;

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

/// Candidate-only inner pole patch prototype with G0 validation evidence.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryInnerPolePatchPrototype {
    /// Stable semantic inner patch name.
    pub patch_name: String,
    /// Adjacent interface name.
    pub interface_name: String,
    /// Candidate inner patch surface.
    pub surface: NurbsSurface,
    /// Inner patch pcurve along the strip/inner interface.
    pub interface_pcurve: Pcurve,
    /// Maximum sampled gap between inner patch and interface curve.
    pub max_interface_g0: f64,
    /// Maximum sampled gap between inner patch pole row and pole point.
    pub max_pole_g0: f64,
    /// Local explicit G0 tolerance.
    pub tolerance: f64,
    /// Interface tangent-plane compatibility evidence, when measured.
    pub interface_g1_report: Option<ExactBoundaryInterfaceG1Report>,
    /// G0 prototype has not established interface G1 yet.
    pub interface_g1_validated: bool,
}

/// Sampled interface G1 evidence between strip and inner patch surfaces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactBoundaryInterfaceG1Report {
    /// Worst tangent-plane deficit along the interface.
    pub max_interface_g1_deficit: f64,
    /// Sample parameter where the worst deficit was measured.
    pub sample_parameter: f64,
    /// Number of deterministic samples evaluated.
    pub sample_count: usize,
    /// Local explicit tangent-plane tolerance.
    pub tolerance: f64,
}

/// Build a G0-only inner pole patch prototype from interface curve to pole.
pub(crate) fn construct_inner_pole_patch_prototype(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    interface_inventory: &ExactBoundaryPatchInterfaceInventory,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    pole: Point3,
    tolerance: f64,
) -> Result<ExactBoundaryInnerPolePatchPrototype, ExactBoundaryMultiPatchDiagnostic> {
    if !curve_geometry_is_finite(&interface.curve) {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_position(
            sector,
            strip,
            interface_inventory,
            f64::INFINITY,
            0.0,
            tolerance,
        ));
    }
    if !point_is_finite(pole) {
        return Err(ExactBoundaryMultiPatchDiagnostic::pole_position(
            sector,
            strip,
            interface_inventory,
            f64::INFINITY,
            0.0,
            tolerance,
        ));
    }
    let surface = inner_pole_surface(interface, pole);
    let interface_pcurve = inner_interface_pcurve(interface);
    let interface_g0 = inner_interface_g0(&surface, interface, &interface_pcurve);
    if interface_g0.max_residual > tolerance {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_position(
            sector,
            strip,
            interface_inventory,
            interface_g0.max_residual,
            interface_g0.sample_parameter,
            tolerance,
        ));
    }
    let pole_g0 = pole_g0(&surface, pole);
    if pole_g0.max_residual > tolerance {
        return Err(ExactBoundaryMultiPatchDiagnostic::pole_position(
            sector,
            strip,
            interface_inventory,
            pole_g0.max_residual,
            pole_g0.sample_parameter,
            tolerance,
        ));
    }

    Ok(ExactBoundaryInnerPolePatchPrototype {
        patch_name: interface.inner_patch_name.clone(),
        interface_name: interface.interface_name.clone(),
        surface,
        interface_pcurve,
        max_interface_g0: interface_g0.max_residual,
        max_pole_g0: pole_g0.max_residual,
        tolerance,
        interface_g1_report: None,
        interface_g1_validated: false,
    })
}

/// Measure strip/inner interface G1 and return an updated prototype.
pub(crate) fn measure_inner_interface_g1(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    interface_inventory: &ExactBoundaryPatchInterfaceInventory,
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    tolerance: f64,
) -> Result<ExactBoundaryInnerPolePatchPrototype, ExactBoundaryMultiPatchDiagnostic> {
    let report = interface_g1_report(strip_candidate, interface, prototype, tolerance);
    if report.max_interface_g1_deficit > tolerance {
        return Err(ExactBoundaryMultiPatchDiagnostic::interface_tangent(
            sector,
            strip,
            interface_inventory,
            report.max_interface_g1_deficit,
            report.sample_parameter,
            tolerance,
        ));
    }
    let mut accepted = prototype.clone();
    accepted.interface_g1_report = Some(report);
    accepted.interface_g1_validated = true;
    Ok(accepted)
}

/// Decline a G0-only prototype before commit until MP-010 validates G1.
pub(crate) fn decline_g0_only_inner_prototype_before_commit(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    interface_inventory: &ExactBoundaryPatchInterfaceInventory,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    interface_g1_tolerance: f64,
) -> ExactBoundaryMultiPatchDiagnostic {
    let report = prototype.interface_g1_report;
    ExactBoundaryMultiPatchDiagnostic::interface_tangent(
        sector,
        strip,
        interface_inventory,
        report
            .map(|g1| g1.max_interface_g1_deficit)
            .unwrap_or(f64::INFINITY),
        report
            .map(|g1| g1.sample_parameter)
            .unwrap_or(interface.sample_parameter),
        interface_g1_tolerance,
    )
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

fn curve_geometry_is_finite(curve: &NurbsCurve) -> bool {
    curve.knots.iter().all(|value| value.is_finite())
        && curve
            .weights
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        && curve
            .control_points
            .iter()
            .all(|point| point_is_finite(*point))
}

fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
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

fn inner_pole_surface(
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    pole: Point3,
) -> NurbsSurface {
    let mut control_points = Vec::with_capacity(interface.curve.control_points.len() * 2);
    let mut weights = Vec::with_capacity(interface.curve.weights.len() * 2);
    for (point, weight) in interface
        .curve
        .control_points
        .iter()
        .zip(&interface.curve.weights)
    {
        control_points.push(*point);
        control_points.push(pole);
        weights.push(*weight);
        weights.push(1.0);
    }
    NurbsSurface {
        degree_u: interface.curve.degree,
        degree_v: 1,
        knots_u: interface.curve.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: interface.curve.control_points.len(),
        num_v: 2,
        weights,
    }
}

fn inner_interface_pcurve(interface: &ExactBoundaryStripInnerInterfaceCurve) -> Pcurve {
    let (t_start, t_end) = interface.curve.parameter_range();
    if interface.strip_parameter_reversed {
        Pcurve::linear((t_end, 0.0), (t_start, 0.0), t_start, t_end)
    } else {
        Pcurve::linear((t_start, 0.0), (t_end, 0.0), t_start, t_end)
    }
}

fn inner_interface_g0(
    surface: &NurbsSurface,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    pcurve: &Pcurve,
) -> InterfaceResidual {
    let (t_start, t_end) = interface.curve.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let uv = pcurve.point_at(parameter);
        let patch_point = surface.point_at(uv.x, uv.y);
        let interface_point = interface.curve.point_at(uv.x);
        let residual = patch_point.distance_to(&interface_point);
        if residual > report.max_residual {
            report = InterfaceResidual {
                max_residual: residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}

fn interface_g1_report(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    tolerance: f64,
) -> ExactBoundaryInterfaceG1Report {
    let (t_start, t_end) = interface.curve.parameter_range();
    let mut report = ExactBoundaryInterfaceG1Report {
        max_interface_g1_deficit: 0.0,
        sample_parameter: t_start,
        sample_count: INTERFACE_SAMPLE_COUNT,
        tolerance,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let strip_uv = interface.strip_pcurve.point_at(parameter);
        let inner_uv = prototype.interface_pcurve.point_at(parameter);
        let deficit = tangent_plane_deficit(
            &strip_candidate.surface,
            strip_uv.x,
            strip_uv.y,
            &prototype.surface,
            inner_uv.x,
            inner_uv.y,
        )
        .unwrap_or(f64::INFINITY);
        if deficit > report.max_interface_g1_deficit {
            report.max_interface_g1_deficit = deficit;
            report.sample_parameter = parameter;
        }
    }
    report
}

fn tangent_plane_deficit(
    strip_surface: &NurbsSurface,
    strip_u: f64,
    strip_v: f64,
    inner_surface: &NurbsSurface,
    inner_u: f64,
    inner_v: f64,
) -> Option<f64> {
    let strip_normal = safe_surface_normal(strip_surface, strip_u, strip_v)?;
    let inner_normal = safe_surface_normal(inner_surface, inner_u, inner_v)?;
    Some(1.0 - strip_normal.dot(&inner_normal).abs())
}

fn safe_surface_normal(surface: &NurbsSurface, u: f64, v: f64) -> Option<Vec3> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u = interior_parameter(u, u_min, u_max);
    let v = interior_parameter(v, v_min, v_max);
    surface.normal_at(u, v).normalized()
}

fn interior_parameter(parameter: f64, min: f64, max: f64) -> f64 {
    let span = max - min;
    if span.abs() <= 1.0e-12 {
        return parameter;
    }
    let inset = span.abs() * 1.0e-6;
    parameter.clamp(min + inset, max - inset)
}

fn pole_g0(surface: &NurbsSurface, pole: Point3) -> InterfaceResidual {
    let ((t_start, t_end), _) = surface.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..POLE_SAMPLE_COUNT {
        let parameter = t_start + (t_end - t_start) * index as f64 / (POLE_SAMPLE_COUNT - 1) as f64;
        let patch_point = surface.point_at(parameter, 1.0);
        let residual = patch_point.distance_to(&pole);
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
