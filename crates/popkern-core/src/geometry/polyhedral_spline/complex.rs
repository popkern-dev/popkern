use super::adapter::{ControlNetSplineOrigin, construct_from_control_net_with_origins};
use super::brep_export::{PolyhedralBrepExportOptions, build_polyhedral_spline_brep};
use super::coverage::PolyhedralCoverageOptions;
use super::net::ControlNet;
use super::types::{PolyhedralSplineError, PolyhedralSplineFamily};
use super::validate::validate_polyhedral_spline;
use crate::diagnostics::Diagnostics;
use crate::geometry::Point3;
use crate::naming::OperationId;
use crate::topology::{BRep, EdgeHandle};

#[derive(Debug, Clone)]
pub struct PolyhedralSplineSurfaceComplex {
    pub source_vertices: usize,
    pub source_faces: usize,
    pub candidate_splines: usize,
    pub components: Vec<PolyhedralSplineSurfaceComponent>,
    pub rejections: Vec<PolyhedralSplineRejection>,
}

#[derive(Debug, Clone)]
pub struct PolyhedralSplineSurfaceComponent {
    pub index: usize,
    pub origin: ControlNetSplineOrigin,
    pub family: PolyhedralSplineFamily,
    pub patches: usize,
    pub boundary_edges: usize,
    pub seam_edges: usize,
    pub collapsed_singularities: usize,
    pub brep: BRep,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolyhedralSplineRejection {
    pub index: usize,
    pub origin: ControlNetSplineOrigin,
    pub family: PolyhedralSplineFamily,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralComplexSewingReport {
    pub components: usize,
    pub boundary_edges: usize,
    pub matched_boundary_pairs: usize,
    pub unmatched_boundary_edges: usize,
    pub max_matched_endpoint_gap: f64,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundarySegment {
    component: usize,
    start: Point3,
    end: Point3,
}

pub fn build_control_net_surface_complex(
    net: &ControlNet,
    options: PolyhedralCoverageOptions,
) -> Result<PolyhedralSplineSurfaceComplex, PolyhedralSplineError> {
    let splines = construct_from_control_net_with_origins(net)?;
    let mut components = Vec::new();
    let mut rejections = Vec::new();

    for (index, item) in splines.iter().enumerate() {
        match build_component(index, item.origin, &item.spline, options) {
            Ok(component) => components.push(component),
            Err(reason) => rejections.push(PolyhedralSplineRejection {
                index,
                origin: item.origin,
                family: item.spline.family,
                reason,
            }),
        }
    }

    Ok(PolyhedralSplineSurfaceComplex {
        source_vertices: net.vertices.len(),
        source_faces: net.faces.len(),
        candidate_splines: splines.len(),
        components,
        rejections,
    })
}

pub fn analyze_surface_complex_sewing(
    complex: &PolyhedralSplineSurfaceComplex,
    tolerance: f64,
) -> PolyhedralComplexSewingReport {
    let boundaries = boundary_segments(complex);
    let mut matched = vec![false; boundaries.len()];
    let mut matched_pairs = 0usize;
    let mut max_gap = 0.0_f64;

    for a in 0..boundaries.len() {
        if matched[a] {
            continue;
        }
        let mut best: Option<(usize, f64)> = None;
        for b in (a + 1)..boundaries.len() {
            if matched[b] || boundaries[a].component == boundaries[b].component {
                continue;
            }
            let gap = reversed_endpoint_gap(boundaries[a], boundaries[b]);
            if gap <= tolerance && best.is_none_or(|(_, current)| gap < current) {
                best = Some((b, gap));
            }
        }
        if let Some((b, gap)) = best {
            matched[a] = true;
            matched[b] = true;
            matched_pairs += 1;
            max_gap = max_gap.max(gap);
        }
    }

    let matched_edges = matched.iter().filter(|&&value| value).count();
    PolyhedralComplexSewingReport {
        components: complex.components.len(),
        boundary_edges: boundaries.len(),
        matched_boundary_pairs: matched_pairs,
        unmatched_boundary_edges: boundaries.len() - matched_edges,
        max_matched_endpoint_gap: max_gap,
    }
}

fn build_component(
    index: usize,
    origin: ControlNetSplineOrigin,
    spline: &super::types::PolyhedralSpline,
    options: PolyhedralCoverageOptions,
) -> Result<PolyhedralSplineSurfaceComponent, String> {
    validate_polyhedral_spline(
        spline,
        options.samples,
        options.c0_tolerance,
        options.angle_tolerance,
        options.jacobian_tolerance,
    )
    .map_err(|err| format!("validation failed: {err}"))?;
    let export = build_polyhedral_spline_brep(
        spline,
        PolyhedralBrepExportOptions {
            tolerance: options.tolerance,
            samples: options.samples,
            c0_tolerance: options.c0_tolerance,
            operation: OperationId(index as u64),
        },
    )
    .map_err(|err| format!("BRep export failed: {err}"))?;
    let diagnostics = Diagnostics::run(&export.brep);
    if !diagnostics.defects.pcurve_escapes.is_empty() {
        return Err(format!(
            "{} pcurve escape diagnostics",
            diagnostics.defects.pcurve_escapes.len()
        ));
    }
    Ok(PolyhedralSplineSurfaceComponent {
        index,
        origin,
        family: spline.family,
        patches: spline.patches.len(),
        boundary_edges: export.boundary_edge_count,
        seam_edges: export.seam_edge_count,
        collapsed_singularities: export.collapsed_singularity_count,
        brep: export.brep,
    })
}

fn boundary_segments(complex: &PolyhedralSplineSurfaceComplex) -> Vec<BoundarySegment> {
    let mut segments = Vec::new();
    for (component_index, component) in complex.components.iter().enumerate() {
        for &edge in &component.brep.boundary_edges {
            if let Some((start, end)) = edge_segment(&component.brep, edge) {
                segments.push(BoundarySegment {
                    component: component_index,
                    start,
                    end,
                });
            }
        }
    }
    segments
}

pub(super) fn edge_segment(brep: &BRep, edge: EdgeHandle) -> Option<(Point3, Point3)> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    let start = brep.vertices.get(half_edge.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    Some((start, end))
}

pub(super) fn reversed_endpoint_gap(a: BoundarySegment, b: BoundarySegment) -> f64 {
    a.start.distance_to(&b.end).max(a.end.distance_to(&b.start))
}
