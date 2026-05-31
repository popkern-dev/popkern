//! Topology helpers for vertex-blend insertion.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

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
