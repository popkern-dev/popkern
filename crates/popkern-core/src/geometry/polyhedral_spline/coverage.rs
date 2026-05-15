use std::collections::BTreeSet;

use super::adapter::{ControlNetSplineOrigin, construct_from_control_net_with_origins};
use super::brep_export::{PolyhedralBrepExportOptions, build_polyhedral_spline_brep};
use super::net::ControlNet;
use super::types::{PolyhedralSplineError, PolyhedralSplineFamily};
use super::validate::validate_polyhedral_spline;
use crate::diagnostics::Diagnostics;
use crate::naming::OperationId;
use crate::tolerance::Tolerance;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolyhedralCoverageOptions {
    pub tolerance: Tolerance,
    pub samples: usize,
    pub c0_tolerance: f64,
    pub angle_tolerance: f64,
    pub jacobian_tolerance: f64,
}

impl Default for PolyhedralCoverageOptions {
    fn default() -> Self {
        let tolerance = Tolerance::default();
        Self {
            tolerance,
            samples: 11,
            c0_tolerance: tolerance.linear(),
            angle_tolerance: 1.0e-6,
            jacobian_tolerance: 1.0e-12,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolyhedralCoverageStatus {
    Accepted {
        brep_edges: usize,
        boundary_edges: usize,
        seam_edges: usize,
        collapsed_singularities: usize,
    },
    Rejected {
        reason: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyhedralCandidateCoverage {
    pub index: usize,
    pub origin: ControlNetSplineOrigin,
    pub family: PolyhedralSplineFamily,
    pub patches: usize,
    pub status: PolyhedralCoverageStatus,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyhedralControlNetCoverageReport {
    pub source_vertices: usize,
    pub source_faces: usize,
    pub candidate_splines: usize,
    pub accepted_splines: usize,
    pub rejected_splines: usize,
    pub accepted_patches: usize,
    pub accepted_faces: usize,
    pub accepted_edges: usize,
    pub boundary_edges: usize,
    pub covered_anchor_vertices: Vec<usize>,
    pub covered_anchor_faces: Vec<usize>,
    pub candidates: Vec<PolyhedralCandidateCoverage>,
}

pub fn analyze_control_net_coverage(
    net: &ControlNet,
    options: PolyhedralCoverageOptions,
) -> Result<PolyhedralControlNetCoverageReport, PolyhedralSplineError> {
    let splines = construct_from_control_net_with_origins(net)?;
    let mut candidates = Vec::with_capacity(splines.len());
    let mut accepted_patches = 0usize;
    let mut accepted_faces = 0usize;
    let mut accepted_edges = 0usize;
    let mut boundary_edges = 0usize;
    let mut covered_vertices = BTreeSet::new();
    let mut covered_faces = BTreeSet::new();

    for (index, item) in splines.iter().enumerate() {
        let status = candidate_status(index, &item.spline, options);
        if let PolyhedralCoverageStatus::Accepted {
            brep_edges,
            boundary_edges: b,
            ..
        } = status
        {
            accepted_patches += item.spline.patches.len();
            accepted_faces += item.spline.patches.len();
            accepted_edges += brep_edges;
            boundary_edges += b;
            match item.origin {
                ControlNetSplineOrigin::Face { face } => {
                    covered_faces.insert(face);
                }
                ControlNetSplineOrigin::Vertex { vertex } => {
                    covered_vertices.insert(vertex);
                }
            }
        }
        candidates.push(PolyhedralCandidateCoverage {
            index,
            origin: item.origin,
            family: item.spline.family,
            patches: item.spline.patches.len(),
            status,
        });
    }

    let accepted_splines = candidates
        .iter()
        .filter(|candidate| matches!(candidate.status, PolyhedralCoverageStatus::Accepted { .. }))
        .count();
    Ok(PolyhedralControlNetCoverageReport {
        source_vertices: net.vertices.len(),
        source_faces: net.faces.len(),
        candidate_splines: candidates.len(),
        accepted_splines,
        rejected_splines: candidates.len() - accepted_splines,
        accepted_patches,
        accepted_faces,
        accepted_edges,
        boundary_edges,
        covered_anchor_vertices: covered_vertices.into_iter().collect(),
        covered_anchor_faces: covered_faces.into_iter().collect(),
        candidates,
    })
}

fn candidate_status(
    index: usize,
    spline: &super::types::PolyhedralSpline,
    options: PolyhedralCoverageOptions,
) -> PolyhedralCoverageStatus {
    if let Err(err) = validate_polyhedral_spline(
        spline,
        options.samples,
        options.c0_tolerance,
        options.angle_tolerance,
        options.jacobian_tolerance,
    ) {
        return rejection(format!("validation failed: {err}"));
    }
    let export = match build_polyhedral_spline_brep(spline, export_options(index, options)) {
        Ok(export) => export,
        Err(err) => return rejection(format!("BRep export failed: {err}")),
    };
    let diagnostics = Diagnostics::run(&export.brep);
    if !diagnostics.defects.pcurve_escapes.is_empty() {
        return rejection(format!(
            "{} pcurve escape diagnostics",
            diagnostics.defects.pcurve_escapes.len()
        ));
    }
    PolyhedralCoverageStatus::Accepted {
        brep_edges: export.brep.edges.len(),
        boundary_edges: export.boundary_edge_count,
        seam_edges: export.seam_edge_count,
        collapsed_singularities: export.collapsed_singularity_count,
    }
}

fn export_options(index: usize, options: PolyhedralCoverageOptions) -> PolyhedralBrepExportOptions {
    PolyhedralBrepExportOptions {
        tolerance: options.tolerance,
        samples: options.samples,
        c0_tolerance: options.c0_tolerance,
        operation: OperationId(index as u64),
    }
}

fn rejection(reason: String) -> PolyhedralCoverageStatus {
    PolyhedralCoverageStatus::Rejected { reason }
}
