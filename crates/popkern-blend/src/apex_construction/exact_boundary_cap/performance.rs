//! Candidate-only performance counters for exact-boundary apex caps.

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCandidatePerformanceReport};

use super::pcurves::PreparedExactBoundaryPcurves;
use super::ribbon::PreparedExactBoundaryRibbons;
use super::solve::ExactBoundarySolvedSurfaces;
use super::support::ExactBoundarySupportSurfaceFamily;
use super::validation::ExactBoundaryValidationReport;

/// Count support samples consumed by prepared ribbons.
pub(crate) fn support_sample_count(ribbons: &PreparedExactBoundaryRibbons) -> usize {
    ribbons
        .ribbons
        .iter()
        .map(|ribbon| ribbon.samples.len())
        .sum()
}

/// Copy support-family evidence into the durable performance report.
pub(crate) fn record_support_family_counts(
    performance: &mut ExactBoundaryCandidatePerformanceReport,
    ribbons: &PreparedExactBoundaryRibbons,
) {
    performance.plane_support_ribbon_count = 0;
    performance.cylinder_support_ribbon_count = 0;
    performance.cone_support_ribbon_count = 0;
    performance.sphere_support_ribbon_count = 0;
    performance.nurbs_support_ribbon_count = 0;
    for ribbon in &ribbons.ribbons {
        match ribbon.support_family {
            ExactBoundarySupportSurfaceFamily::Plane => {
                performance.plane_support_ribbon_count += 1;
            }
            ExactBoundarySupportSurfaceFamily::Cylinder => {
                performance.cylinder_support_ribbon_count += 1;
            }
            ExactBoundarySupportSurfaceFamily::Cone => {
                performance.cone_support_ribbon_count += 1;
            }
            ExactBoundarySupportSurfaceFamily::Sphere => {
                performance.sphere_support_ribbon_count += 1;
            }
            ExactBoundarySupportSurfaceFamily::Nurbs => {
                performance.nurbs_support_ribbon_count += 1;
            }
            ExactBoundarySupportSurfaceFamily::Unsupported => {}
        }
    }
}

/// Count solved surface control points allocated for the candidate.
pub(crate) fn surface_control_point_count(solved: &ExactBoundarySolvedSurfaces) -> usize {
    solved
        .surfaces
        .iter()
        .map(|surface| surface.surface.control_points.len())
        .sum()
}

/// Copy validation evidence into the durable performance report.
pub(crate) fn record_validation_report(
    performance: &mut ExactBoundaryCandidatePerformanceReport,
    report: ExactBoundaryValidationReport,
) {
    performance.validation_sector_count = report.validated_sector_count;
    performance.validation_pcurve_pair_count = report.validated_pcurve_pair_count;
    performance.validation_spoke_edge_count = report.validated_spoke_edge_count;
}

/// Copy pcurve preparation evidence into the durable performance report.
pub(crate) fn record_pcurve_report(
    performance: &mut ExactBoundaryCandidatePerformanceReport,
    pcurves: &PreparedExactBoundaryPcurves,
) {
    performance.pcurve_pair_count = pcurves.validated_pair_count();
}

/// Finalize candidate timing and deterministic decline evidence.
pub(crate) fn finalize_performance_report(
    performance: &mut ExactBoundaryCandidatePerformanceReport,
    started: popkern_core::time::Stopwatch,
    diagnostics: &[ApexConstructionDiagnostic],
) {
    performance.elapsed_micros = started.elapsed().map(|d| d.as_micros());
    performance.deterministic_skip_code = deterministic_skip_code(diagnostics);
}

/// Return the first known deterministic exact-boundary decline code.
///
/// These codes are pure functions of candidate inventory, tolerances, and BRep
/// topology/geometry state. Recording them here gives future runtime policy a
/// durable skip key without suppressing diagnostics from the normal report path.
pub(crate) fn deterministic_skip_code(
    diagnostics: &[ApexConstructionDiagnostic],
) -> Option<String> {
    diagnostics.iter().find_map(|diagnostic| {
        let ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) = diagnostic else {
            return None;
        };
        if is_deterministic_exact_boundary_decline(&diagnostic.code) {
            Some(diagnostic.code.clone())
        } else {
            None
        }
    })
}

fn is_deterministic_exact_boundary_decline(code: &str) -> bool {
    code.starts_with("exact_boundary_cap.")
}
