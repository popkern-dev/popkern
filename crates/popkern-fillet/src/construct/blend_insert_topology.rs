//! Topology helpers for vertex-blend insertion.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    find_face_trim_anchors, find_halfedge_from_vertex_in_face, identify_corner_face,
};

/// Phase 1: Split off corner regions from each face around V using
/// the per-face fillet trim endpoints as mef anchors.
pub(super) fn split_corner_triangles(
    brep: &mut BRep,
    vertex: VertexHandle,
    face_info: &[(FaceHandle, HalfEdgeHandle)],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, String> {
    let mut corner_faces = Vec::new();

    for (face, _he_arriving_at_v) in face_info {
        let anchors = find_face_trim_anchors(brep, *face, vertex, edge_geometries, tolerance)?;
        if anchors.len() < 2 {
            return Err(format!(
                "face {face:?} has {} trim-endpoint anchors around vertex {vertex:?}, expected >=2",
                anchors.len()
            ));
        }
        let (anchor_a, anchor_b) = (anchors[0], anchors[1]);

        let he_at_a = find_halfedge_from_vertex_in_face(brep, *face, anchor_a)?;
        let he_at_b = find_halfedge_from_vertex_in_face(brep, *face, anchor_b)?;

        let loop_a = brep.half_edges[he_at_a].loop_;
        let loop_b = brep.half_edges[he_at_b].loop_;
        if loop_a != loop_b {
            return Err(format!(
                "anchor half-edges on face {face:?} are in different loops"
            ));
        }

        let mef_result = brep
            .mef(he_at_a, he_at_b)
            .map_err(|e| format!("mef failed on face {face:?}: {e}"))?;

        let corner = identify_corner_face(brep, *face, mef_result.face, vertex)?;
        corner_faces.push(corner);
    }

    Ok(corner_faces)
}

/// Normalize an under-resolved high-valence blend loop before the
/// N-sector Gregory/polar poke. When Euler merging collapses several
/// boundary intervals into fewer than `n_sides` topological edges,
/// split the longest live loop edge until the loop has `n_sides`
/// edges. The existing poke then performs its ordinary N-to-2N
/// midpoint split.
pub(super) fn normalize_underresolved_face_loop_to_n(
    brep: &mut BRep,
    face: FaceHandle,
    n_sides: usize,
) -> Result<usize, String> {
    if n_sides < 5 {
        return Ok(0);
    }
    let mut splits = 0usize;
    loop {
        let loop_hes = outer_loop_half_edges(brep, face)?;
        let count = loop_hes.len();
        if count >= n_sides {
            return Ok(splits);
        }
        let Some(edge) = longest_loop_edge(brep, &loop_hes) else {
            return Err(format!(
                "underresolved face-loop normalizer found no splittable edge at count {count}"
            ));
        };
        brep.split_edge(edge, 0.5)
            .map_err(|e| format!("underresolved face-loop split failed: {e}"))?;
        splits += 1;
        if splits > n_sides {
            return Err(format!(
                "underresolved face-loop normalizer exceeded split budget for N={n_sides}"
            ));
        }
    }
}

fn outer_loop_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| "face handle is invalid".to_string())?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| "face outer loop is invalid".to_string())?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Err("outer loop start half-edge is invalid".to_string());
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..64 {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err("outer loop walk exceeded safety cap".to_string())
}

fn longest_loop_edge(brep: &BRep, loop_hes: &[HalfEdgeHandle]) -> Option<EdgeHandle> {
    loop_hes
        .iter()
        .filter_map(|&he| {
            let he_data = brep.half_edges.get(he)?;
            let next = brep.half_edges.get(he_data.next)?;
            let a = brep.vertices.get(he_data.origin)?.position;
            let b = brep.vertices.get(next.origin)?.position;
            Some((he_data.edge, a.distance_to(&b)))
        })
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(edge, _)| edge)
}
