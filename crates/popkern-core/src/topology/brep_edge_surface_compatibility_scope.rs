//! Scoped edge-curve to face-surface compatibility diagnostics.

use crate::healing::deviation::edge_surface_deviation;
use crate::topology::{
    BRep, BrepEdgeSurfaceCompatibilityFinding, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityReport, FaceHandle, HalfEdgeHandle,
};

/// Analyze whether scoped edge uses lie on their adjacent face surfaces.
pub fn analyze_scoped_brep_edge_surface_compatibility(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> BrepEdgeSurfaceCompatibilityReport {
    let mut report = BrepEdgeSurfaceCompatibilityReport {
        inspected_half_edges: 0,
        incompatible_half_edges: 0,
        max_deviation: 0.0,
        findings: Vec::new(),
    };
    let mut inspected = Vec::new();
    for &half_edge in half_edges {
        if inspected.contains(&half_edge) {
            continue;
        }
        inspected.push(half_edge);
        let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
            continue;
        };
        let Some(edge_data) = brep.edges.get(half_edge_data.edge) else {
            continue;
        };
        let Some(curve) = &edge_data.curve else {
            continue;
        };
        let Some((face, surface)) = face_surface_for_half_edge(brep, half_edge) else {
            continue;
        };
        let deviation = edge_surface_deviation(
            &curve.geometry,
            (curve.t_start, curve.t_end),
            surface,
            options.samples,
        );
        report.inspected_half_edges += 1;
        report.max_deviation = report.max_deviation.max(deviation.max_deviation);
        if deviation.max_deviation > options.tolerance {
            report.incompatible_half_edges += 1;
            report.findings.push(BrepEdgeSurfaceCompatibilityFinding {
                edge: half_edge_data.edge,
                half_edge,
                face,
                max_deviation: deviation.max_deviation,
                max_curve_parameter: deviation.max_location,
            });
        }
    }
    report
}

fn face_surface_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, &crate::geometry::surface::SurfaceGeometry)> {
    let loop_ = brep.half_edges.get(half_edge)?.loop_;
    let face = brep.loops.get(loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    Some((face, surface))
}
