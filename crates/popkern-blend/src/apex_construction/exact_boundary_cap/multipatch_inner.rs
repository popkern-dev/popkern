//! Inner pole patch construction for multi-patch exact-boundary candidates.
#![allow(dead_code, clippy::result_large_err)]

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};

use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

use super::inner_geometry::{
    INTERFACE_SAMPLE_COUNT, curve_geometry_is_finite, inner_interface_g0, inner_interface_pcurve,
    inner_pole_surface, interface_g1_report, interface_residual, point_is_finite, pole_g0,
    strip_inner_boundary_curve, strip_inner_boundary_pcurve,
};
use super::inner_optimizer::optimize_inner_interface_g1;
use super::inventory::{
    ExactBoundaryOuterStripPatchInventory, ExactBoundaryPatchInterfaceInventory,
};
use super::report::ExactBoundaryMultiPatchDiagnostic;
use super::strip::{ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection};

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
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
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
    let surface = inner_pole_surface(strip_candidate, interface, pole);
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
    let optimized = optimize_inner_interface_g1(strip_candidate, interface, prototype, tolerance);
    let report = interface_g1_report(strip_candidate, interface, &optimized, tolerance);
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
    let mut accepted = optimized;
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

#[cfg(test)]
#[path = "multipatch_inner_tests.rs"]
mod tests;
