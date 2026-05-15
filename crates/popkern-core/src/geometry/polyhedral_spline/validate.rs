use super::atlas::{SeamContinuityReport, validate_spline_atlas};
use super::readiness::{PolyhedralSplineBrepReadinessReport, validate_brep_readiness};
use super::types::{PatchEdge, PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError};
use crate::geometry::Point3;

#[derive(Debug, Clone, PartialEq)]
pub struct ContinuityReport {
    pub samples: usize,
    pub max_position_gap: f64,
    pub tolerance: f64,
    pub accepted: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PatchRegularityReport {
    pub patch_index: usize,
    pub samples_per_axis: usize,
    pub min_normal_length: f64,
    pub tolerance: f64,
    pub accepted: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralSplineValidationReport {
    pub patch_count: usize,
    pub seam_reports: Vec<SeamContinuityReport>,
    pub regularity_reports: Vec<PatchRegularityReport>,
    pub brep_readiness: PolyhedralSplineBrepReadinessReport,
}

pub fn validate_finite_patch(patch: &PolyhedralPatch) -> Result<(), PolyhedralSplineError> {
    for (index, point) in patch.control_points.iter().enumerate() {
        if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
            return Err(PolyhedralSplineError::NonFinitePoint {
                family: patch.family,
                index,
            });
        }
    }
    Ok(())
}

pub fn validate_finite_spline(spline: &PolyhedralSpline) -> Result<(), PolyhedralSplineError> {
    for patch in &spline.patches {
        validate_finite_patch(patch)?;
    }
    Ok(())
}

pub fn validate_c0_edge_match(
    a: &PolyhedralPatch,
    edge_a: PatchEdge,
    b: &PolyhedralPatch,
    edge_b: PatchEdge,
    samples: usize,
    tolerance: f64,
) -> ContinuityReport {
    let samples = samples.max(2);
    let mut max_position_gap = 0.0_f64;
    for i in 0..samples {
        let t = i as f64 / (samples - 1) as f64;
        let pa = eval_edge(a, edge_a, t);
        let pb = eval_edge(b, edge_b, t);
        max_position_gap = max_position_gap.max(distance(pa, pb));
    }
    ContinuityReport {
        samples,
        max_position_gap,
        tolerance,
        accepted: max_position_gap <= tolerance,
    }
}

pub fn validate_polyhedral_spline(
    spline: &PolyhedralSpline,
    samples: usize,
    c0_tolerance: f64,
    angle_tolerance: f64,
    jacobian_tolerance: f64,
) -> Result<PolyhedralSplineValidationReport, PolyhedralSplineError> {
    validate_finite_spline(spline)?;
    let seam_reports = validate_spline_atlas(spline, samples, c0_tolerance, angle_tolerance)?;
    let brep_readiness = validate_brep_readiness(spline, samples, c0_tolerance)?;
    let mut regularity_reports = Vec::with_capacity(spline.patches.len());
    for (patch_index, patch) in spline.patches.iter().enumerate() {
        let report =
            validate_patch_interior_regularity(patch_index, patch, samples, jacobian_tolerance);
        if !report.accepted {
            return Err(PolyhedralSplineError::DegeneratePatchInterior {
                family: spline.family,
                patch_index,
                min_normal_length: report.min_normal_length,
            });
        }
        regularity_reports.push(report);
    }
    Ok(PolyhedralSplineValidationReport {
        patch_count: spline.patches.len(),
        seam_reports,
        regularity_reports,
        brep_readiness,
    })
}

pub fn validate_patch_interior_regularity(
    patch_index: usize,
    patch: &PolyhedralPatch,
    samples_per_axis: usize,
    tolerance: f64,
) -> PatchRegularityReport {
    let samples_per_axis = samples_per_axis.max(3);
    let mut min_normal_length = f64::INFINITY;
    for i in 1..(samples_per_axis - 1) {
        let u = i as f64 / (samples_per_axis - 1) as f64;
        for j in 1..(samples_per_axis - 1) {
            let v = j as f64 / (samples_per_axis - 1) as f64;
            let normal_length = patch.partial_u(u, v).cross(&patch.partial_v(u, v)).length();
            min_normal_length = min_normal_length.min(normal_length);
        }
    }
    PatchRegularityReport {
        patch_index,
        samples_per_axis,
        min_normal_length,
        tolerance,
        accepted: min_normal_length >= tolerance,
    }
}

fn eval_edge(patch: &PolyhedralPatch, edge: PatchEdge, t: f64) -> Point3 {
    match edge {
        PatchEdge::UMin => patch.evaluate(0.0, t),
        PatchEdge::UMax => patch.evaluate(1.0, t),
        PatchEdge::VMin => patch.evaluate(t, 0.0),
        PatchEdge::VMax => patch.evaluate(t, 1.0),
    }
}

fn distance(a: Point3, b: Point3) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}
