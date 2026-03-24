//! Vertex blend face insertion into BRep topology.
//!
//! At a vertex where multiple fillets meet, carves a blend face from the
//! surrounding faces using Euler operators (`mef` / `kef`). The blend face
//! replaces the small triangular region around the vertex.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend;
use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::{
    build_fillet_surface_list, collect_anchor_vertices, collect_edges_through_vertex,
    collect_faces_around_vertex, find_face_anchors, find_halfedge_from_vertex_in_face,
    identify_corner_face, repair_vertex_outgoing,
};
use crate::types::VertexMode;

/// Insert a smooth G1 vertex blend face at a multi-fillet vertex.
///
/// For each face touching the vertex, splits off the corner triangle using
/// `mef()`, then merges all corner triangles into one N-gon by removing
/// edges through the vertex using `kef()`. Finally assigns the blend surface
/// geometry to the merged face.
///
/// # Errors
///
/// Returns a descriptive string if any topology operation fails.
pub(crate) fn insert_vertex_blend(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<(), String> {
    // Repair the vertex's outgoing half-edge if it has become stale
    // due to prior topology operations (edge splitting, face splitting).
    repair_vertex_outgoing(brep, vertex)?;

    // Collect the set of anchor vertices adjacent to V in faces around V.
    let anchor_vertices = collect_anchor_vertices(brep, vertex, edge_indices, edge_geometries)?;

    if anchor_vertices.len() < 2 {
        return Err(format!(
            "vertex {vertex:?} has fewer than 2 anchor vertices ({}) — cannot form blend face",
            anchor_vertices.len()
        ));
    }

    // Collect faces touching V (via half-edges originating from V).
    let face_info = collect_faces_around_vertex(brep, vertex)?;

    if face_info.is_empty() {
        return Err(format!("no faces found around vertex {vertex:?}"));
    }

    // Phase 1: For each face, find the two anchor vertices adjacent to V
    // in that face's loop, and call mef() to split off the corner triangle.
    let _corner_faces = split_corner_triangles(brep, vertex, &face_info, &anchor_vertices)?;

    // Phase 2: Merge all corner triangles by removing edges through V.
    let blend_face = merge_corner_triangles(brep, vertex)?;

    // Phase 3: Remove the orphaned vertex V if it has no remaining edges.
    remove_orphaned_vertex(brep, vertex);

    // Phase 4: Assign blend surface geometry to the merged face.
    assign_blend_surface(
        brep,
        blend_face,
        vertex,
        edge_indices,
        edge_geometries,
        tolerance,
    )
}

/// Phase 1: Split off corner triangles from each face around V.
fn split_corner_triangles(
    brep: &mut BRep,
    vertex: VertexHandle,
    face_info: &[(FaceHandle, popkern_core::topology::HalfEdgeHandle)],
    anchor_vertices: &std::collections::HashSet<VertexHandle>,
) -> Result<Vec<FaceHandle>, String> {
    let mut corner_faces: Vec<FaceHandle> = Vec::new();

    for (face, he_arriving_at_v) in face_info {
        let (anchor_before, anchor_after) =
            find_face_anchors(brep, *face, vertex, *he_arriving_at_v, anchor_vertices)?;

        let he_at_anchor_after = find_halfedge_from_vertex_in_face(brep, *face, anchor_after)?;
        let he_at_anchor_before = find_halfedge_from_vertex_in_face(brep, *face, anchor_before)?;

        // Verify both half-edges are in the same loop before calling mef.
        let loop_a = brep.half_edges[he_at_anchor_after].loop_;
        let loop_b = brep.half_edges[he_at_anchor_before].loop_;
        if loop_a != loop_b {
            return Err(format!(
                "anchor half-edges on face {face:?} are in different loops"
            ));
        }

        let mef_result = brep
            .mef(he_at_anchor_after, he_at_anchor_before)
            .map_err(|e| format!("mef failed on face {face:?}: {e}"))?;

        let corner = identify_corner_face(brep, *face, mef_result.face, vertex)?;
        corner_faces.push(corner);
    }

    Ok(corner_faces)
}

/// Phase 2: Merge corner triangles by removing edges through V via `kef`.
///
/// Returns the surviving blend face handle.
fn merge_corner_triangles(brep: &mut BRep, vertex: VertexHandle) -> Result<FaceHandle, String> {
    let edges_through_v = collect_edges_through_vertex(brep, vertex)?;

    let mut blend_face: Option<FaceHandle> = None;

    for edge in edges_through_v {
        if !brep.edges.contains_key(edge) {
            continue; // Edge was already removed by a prior kef.
        }

        let he = brep.edges[edge].half_edge;
        let he_loop = brep.half_edges[he].loop_;
        let he_face = brep.loops[he_loop].face;
        let twin = brep.half_edges[he].twin;
        let twin_loop = brep.half_edges[twin].loop_;
        let twin_face = brep.loops[twin_loop].face;

        if he_face == twin_face {
            continue; // Already merged.
        }

        brep.kef(he)
            .map_err(|e| format!("kef failed removing edge {edge:?} through vertex: {e}"))?;

        blend_face = Some(he_face);
    }

    blend_face
        .ok_or_else(|| "no edges were removed through vertex — blend face not created".to_string())
}

/// Phase 3: Remove vertex V if it is no longer connected to any edges.
fn remove_orphaned_vertex(brep: &mut BRep, vertex: VertexHandle) {
    let v_has_edges = brep.half_edges.values().any(|he| he.origin == vertex);
    if !v_has_edges {
        brep.vertices.remove(vertex);
    }
}

/// Phase 4: Compute and assign the blend surface to the merged face.
fn assign_blend_surface(
    brep: &mut BRep,
    blend_face: FaceHandle,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<(), String> {
    let fillet_surfaces = build_fillet_surface_list(edge_indices, edge_geometries);
    let blend_boundary =
        blend::compute_blend_boundary(vertex, edge_indices, &fillet_surfaces, brep, tolerance)
            .map_err(|e| format!("blend boundary computation failed: {e}"))?;

    let blend_result = blend::build_blend_surface(
        &blend_boundary,
        &fillet_surfaces,
        brep,
        VertexMode::Smooth,
        tolerance,
    )
    .map_err(|e| format!("blend surface computation failed: {e}"))?;

    brep.faces[blend_face].surface = Some(SurfaceOnFace {
        geometry: blend_result.surface,
        orientation: true,
    });

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_vertex_blend_too_few_anchors() {
        // A cube vertex has 3 neighbors, but if we pass no edge indices
        // the anchor set will be non-empty (all neighbors). To trigger the
        // "fewer than 2 anchors" path we need a degenerate vertex.
        // Instead, test with a default (non-existent) vertex handle.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = popkern_core::topology::VertexHandle::default();
        let result = insert_vertex_blend(&mut brep, bogus, &[], &[], &Tolerance::default());
        assert!(result.is_err(), "should fail for invalid vertex");
    }

    #[test]
    fn insert_vertex_blend_on_cube_vertex() {
        // Calling insert_vertex_blend on a real cube vertex with no
        // edge geometries should fail during blend surface computation
        // (since there are no fillet surfaces to blend), but it should
        // get past the anchor check (3 anchors >= 2).
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = brep.vertices.keys().next().unwrap();
        let result = insert_vertex_blend(&mut brep, vertex, &[], &[], &Tolerance::default());
        // This should fail (no edge geometries -> empty fillet surfaces
        // -> blend boundary computation fails), but NOT at the anchor check.
        assert!(result.is_err());
        let err = result.unwrap_err();
        // The error should be from blend computation, not "fewer than 2 anchors".
        assert!(
            !err.contains("fewer than 2"),
            "should not fail at anchor check: {err}"
        );
    }
}
