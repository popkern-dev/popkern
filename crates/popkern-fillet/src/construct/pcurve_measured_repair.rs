//! Measured pcurve repair for late export-readiness gates.

use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

pub(crate) fn repair_measured_plane_cylinder_pcurves(brep: &mut BRep, tolerance: f64) -> usize {
    let edges = brep.edges.keys().collect::<Vec<_>>();
    let mut repaired = 0;
    for edge in edges {
        let Some(report) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if report.max_deviation <= tolerance {
            continue;
        }
        if repair_plane_cylinder_edge_from_triple(brep, edge, tolerance) {
            repaired += 1;
        }
    }
    repaired
}

fn repair_plane_cylinder_edge_from_triple(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: f64,
) -> bool {
    let half_edges = edge_half_edges(brep, edge);
    if half_edges.len() != 2 {
        return false;
    }
    let Some((_face_a, surface_a)) = face_surface_for_half_edge(brep, half_edges[0]) else {
        return false;
    };
    let Some((_face_b, surface_b)) = face_surface_for_half_edge(brep, half_edges[1]) else {
        return false;
    };
    if !is_plane_cylinder_pair(&surface_a, &surface_b) {
        return false;
    }
    let start = brep.vertices[brep.half_edges[half_edges[0]].origin].position;
    let end = brep.vertices[brep.half_edges[brep.half_edges[half_edges[0]].next].origin].position;
    let Ok(triple) = compute_edge_geometry(
        &surface_a,
        &surface_b,
        EdgeAnchors::Open { start, end },
        None,
        &popkern_core::tolerance::Tolerance::default(),
        &EdgeTripleOptions::default(),
    ) else {
        return false;
    };
    let previous_pcurves = half_edges
        .iter()
        .map(|half_edge| (*half_edge, brep.half_edges[*half_edge].pcurve.clone()))
        .collect::<Vec<_>>();
    brep.half_edges[half_edges[0]].pcurve = Some(triple.pcurve_a.clone());
    brep.half_edges[half_edges[1]].pcurve = Some(triple.pcurve_b.clone());
    let accepted = brep
        .measure_edge_pcurve_deviation(edge)
        .is_some_and(|updated| updated.max_deviation <= tolerance)
        && brep.validate().is_ok();
    if accepted {
        true
    } else {
        for (half_edge, pcurve) in previous_pcurves {
            brep.half_edges[half_edge].pcurve = pcurve;
        }
        false
    }
}

fn is_plane_cylinder_pair(surface_a: &SurfaceGeometry, surface_b: &SurfaceGeometry) -> bool {
    matches!(surface_a, SurfaceGeometry::Plane { .. })
        && matches!(surface_b, SurfaceGeometry::Cylinder { .. })
        || matches!(surface_a, SurfaceGeometry::Cylinder { .. })
            && matches!(surface_b, SurfaceGeometry::Plane { .. })
}

fn face_surface_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, SurfaceGeometry)> {
    let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
    let surface = brep.faces.get(face)?.surface.as_ref()?.geometry.clone();
    Some((face, surface))
}

fn edge_half_edges(brep: &BRep, edge: EdgeHandle) -> Vec<HalfEdgeHandle> {
    let Some(primary) = brep.edges.get(edge).map(|edge| edge.half_edge) else {
        return Vec::new();
    };
    let mut out = vec![primary];
    if let Some(twin) = brep.half_edges.get(primary).map(|half_edge| half_edge.twin)
        && twin != primary
    {
        out.push(twin);
    }
    out
}
