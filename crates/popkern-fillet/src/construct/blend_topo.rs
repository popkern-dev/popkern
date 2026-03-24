//! Topology walking helpers for vertex blend insertion.
//!
//! Provides functions to walk around a vertex in the BRep topology,
//! collect adjacent faces, anchor vertices, and incident edges.
//! These are used by [`super::blend_insert`] during blend face insertion.

use std::collections::HashSet;

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::fillet_surface::FilletSurfaceResult;

/// Repair a vertex's outgoing half-edge if it has become stale.
///
/// After edge splitting and face splitting, a vertex's `outgoing` field
/// may point to a half-edge that was removed. This scans all half-edges
/// to find one originating from the vertex and updates the field.
pub(crate) fn repair_vertex_outgoing(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
    let v = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?;

    if brep.half_edges.contains_key(v.outgoing) {
        return Ok(()); // Already valid.
    }

    // Scan for any half-edge originating from this vertex.
    let mut found = None;
    for (he_handle, he) in &brep.half_edges {
        if he.origin == vertex {
            found = Some(he_handle);
            break;
        }
    }

    let new_outgoing = found.ok_or_else(|| {
        let count = brep
            .half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .count();
        format!("vertex {vertex:?} has no half-edges originating from it (scan found {count})")
    })?;
    brep.vertices[vertex].outgoing = new_outgoing;
    Ok(())
}

/// Collect anchor vertices around a vertex V.
///
/// The anchor vertices are the non-V endpoints of edges that connect
/// directly to V in the current BRep topology.
pub(crate) fn collect_anchor_vertices(
    brep: &BRep,
    vertex: VertexHandle,
    _edge_indices: &[usize],
    _edge_geometries: &[EdgeGeometry],
) -> Result<HashSet<VertexHandle>, String> {
    let mut anchors = HashSet::new();

    let outgoing = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .outgoing;

    if !brep.half_edges.contains_key(outgoing) {
        return Err(format!("vertex {vertex:?} has no valid outgoing half-edge"));
    }

    // Walk around V using twin.next to visit all outgoing half-edges.
    let mut cursor = outgoing;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found during anchor walk"))?;

        let twin = brep
            .half_edges
            .get(he.twin)
            .ok_or_else(|| format!("twin {0:?} not found", he.twin))?;
        anchors.insert(twin.origin);

        cursor = brep.half_edges[he.twin].next;

        iters += 1;
        if cursor == outgoing || iters >= max_iters {
            break;
        }
    }

    Ok(anchors)
}

/// Collect faces around a vertex with the half-edge arriving at V.
///
/// Returns `(face, he_arriving_at_v)` for each face incident to V.
pub(crate) fn collect_faces_around_vertex(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<(FaceHandle, HalfEdgeHandle)>, String> {
    let outgoing = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .outgoing;

    if !brep.half_edges.contains_key(outgoing) {
        return Err(format!("vertex {vertex:?} has no valid outgoing half-edge"));
    }

    let mut faces = Vec::new();
    let mut seen_faces = HashSet::new();
    let mut cursor = outgoing;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;

    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;

        let face = brep.loops[he.loop_].face;
        if seen_faces.insert(face) {
            faces.push((face, he.prev));
        }

        cursor = brep.half_edges[he.twin].next;
        iters += 1;
        if cursor == outgoing || iters >= max_iters {
            break;
        }
    }

    Ok(faces)
}

/// Find the two anchor vertices adjacent to V in a face's loop.
///
/// Returns `(anchor_before_V, anchor_after_V)`.
pub(crate) fn find_face_anchors(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
    _he_arriving_at_v: HalfEdgeHandle,
    anchor_set: &HashSet<VertexHandle>,
) -> Result<(VertexHandle, VertexHandle), String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep.loops[outer_loop].half_edge;

    // Walk the loop to find the half-edge originating from V.
    let mut he_from_v = None;
    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        if brep.half_edges[cursor].origin == vertex {
            he_from_v = Some(cursor);
            break;
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }

    let he_from_v = he_from_v
        .ok_or_else(|| format!("vertex {vertex:?} not found in face {face:?}'s outer loop"))?;

    let he_arriving = brep.half_edges[he_from_v].prev;
    let anchor_before = brep.half_edges[he_arriving].origin;

    let he_leaving_next = brep.half_edges[he_from_v].next;
    let anchor_after = brep.half_edges[he_leaving_next].origin;

    if !anchor_set.contains(&anchor_before) {
        return Err(format!(
            "vertex before V ({anchor_before:?}) on face {face:?} is not an anchor vertex"
        ));
    }
    if !anchor_set.contains(&anchor_after) {
        return Err(format!(
            "vertex after V ({anchor_after:?}) on face {face:?} is not an anchor vertex"
        ));
    }

    Ok((anchor_before, anchor_after))
}

/// Find a half-edge originating from the given vertex in a face's loop.
pub(crate) fn find_halfedge_from_vertex_in_face(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<HalfEdgeHandle, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep.loops[outer_loop].half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        if brep.half_edges[cursor].origin == vertex {
            return Ok(cursor);
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }

    Err(format!(
        "no half-edge from vertex {vertex:?} found in face {face:?}"
    ))
}

/// Determine which of two faces is the corner triangle containing the vertex.
pub(crate) fn identify_corner_face(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    vertex: VertexHandle,
) -> Result<FaceHandle, String> {
    if face_loop_contains_vertex(brep, face_b, vertex)? {
        Ok(face_b)
    } else if face_loop_contains_vertex(brep, face_a, vertex)? {
        Ok(face_a)
    } else {
        Err(format!(
            "neither face {face_a:?} nor {face_b:?} contains vertex {vertex:?}"
        ))
    }
}

/// Check if a face's outer loop contains a vertex.
pub(crate) fn face_loop_contains_vertex(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<bool, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let start_he = brep.loops[face_data.outer_loop].half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        if brep.half_edges[cursor].origin == vertex {
            return Ok(true);
        }
        cursor = brep.half_edges[cursor].next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Ok(false)
}

/// Collect all edges incident to a vertex.
pub(crate) fn collect_edges_through_vertex(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<EdgeHandle>, String> {
    let outgoing = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .outgoing;

    if !brep.half_edges.contains_key(outgoing) {
        return Err(format!("vertex {vertex:?} has no valid outgoing half-edge"));
    }

    let mut edges = Vec::new();
    let mut seen = HashSet::new();
    let mut cursor = outgoing;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;

    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;

        if seen.insert(he.edge) {
            edges.push(he.edge);
        }

        cursor = brep.half_edges[he.twin].next;
        iters += 1;
        if cursor == outgoing || iters >= max_iters {
            break;
        }
    }

    Ok(edges)
}

/// Build `FilletSurfaceResult` list from edge geometries for blend boundary
/// computation.
pub(crate) fn build_fillet_surface_list(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
) -> Vec<FilletSurfaceResult> {
    edge_indices
        .iter()
        .filter_map(|&idx| {
            let eg = edge_geometries.get(idx)?;
            Some(FilletSurfaceResult {
                surface: eg.surface.clone(),
                trim_curve_1: eg.trim_3d_1.clone(),
                trim_curve_2: eg.trim_3d_2.clone(),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn default_tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn collect_faces_around_cube_vertex() {
        let brep = BRep::unit_cube(default_tol()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let faces = collect_faces_around_vertex(&brep, vertex);
        assert!(faces.is_ok(), "should succeed: {:?}", faces.err());
        let faces = faces.expect("already checked");
        assert_eq!(faces.len(), 3, "cube vertex touches 3 faces");
    }

    #[test]
    fn collect_anchors_around_cube_vertex() {
        let brep = BRep::unit_cube(default_tol()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let anchors = collect_anchor_vertices(&brep, vertex, &[], &[]);
        assert!(anchors.is_ok(), "should succeed: {:?}", anchors.err());
        let anchors = anchors.expect("already checked");
        assert_eq!(anchors.len(), 3, "cube vertex has 3 neighbors");
    }

    #[test]
    fn collect_edges_through_cube_vertex() {
        let brep = BRep::unit_cube(default_tol()).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let edges = collect_edges_through_vertex(&brep, vertex);
        assert!(edges.is_ok(), "should succeed: {:?}", edges.err());
        let edges = edges.expect("already checked");
        assert_eq!(edges.len(), 3, "cube vertex has 3 incident edges");
    }

    #[test]
    fn collect_faces_invalid_vertex() {
        let brep = BRep::unit_cube(default_tol()).expect("cube should build");
        let bogus = VertexHandle::default();
        let result = collect_faces_around_vertex(&brep, bogus);
        assert!(result.is_err(), "should fail for invalid vertex");
    }

    /// Helper to collect vertices of a face's outer loop.
    fn face_vertex_set(brep: &BRep, face: FaceHandle) -> HashSet<VertexHandle> {
        let start_he = brep.loops[brep.faces[face].outer_loop].half_edge;
        let mut verts = HashSet::new();
        let mut cursor = start_he;
        loop {
            verts.insert(brep.half_edges[cursor].origin);
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
        verts
    }

    /// Find a vertex NOT in the given face's outer loop.
    fn vertex_not_in_face(brep: &BRep, face: FaceHandle) -> Option<VertexHandle> {
        let fv = face_vertex_set(brep, face);
        brep.vertices.keys().find(|v| !fv.contains(v))
    }

    #[test]
    fn find_face_anchors_vertex_not_in_face() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        if let Some(v) = vertex_not_in_face(&brep, face) {
            let he = brep.half_edges.keys().next().unwrap();
            let anchors: HashSet<VertexHandle> = brep.vertices.keys().collect();
            let result = find_face_anchors(&brep, face, v, he, &anchors);
            assert!(result.is_err(), "should fail when vertex not in face loop");
        }
    }

    #[test]
    fn invalid_vertex_error_paths() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let bogus = VertexHandle::default();
        assert!(collect_faces_around_vertex(&brep, bogus).is_err());
        assert!(collect_anchor_vertices(&brep, bogus, &[], &[]).is_err());
        assert!(collect_edges_through_vertex(&brep, bogus).is_err());
    }

    #[test]
    fn find_halfedge_from_vertex_in_face_not_found() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        if let Some(v) = vertex_not_in_face(&brep, face) {
            assert!(find_halfedge_from_vertex_in_face(&brep, face, v).is_err());
        }
    }

    #[test]
    fn identify_corner_face_neither_contains() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let faces: Vec<FaceHandle> = brep.faces.keys().collect();
        let (fa, fb) = (faces[0], faces[1]);
        let va = face_vertex_set(&brep, fa);
        let vb = face_vertex_set(&brep, fb);
        if let Some(v) = brep
            .vertices
            .keys()
            .find(|v| !va.contains(v) && !vb.contains(v))
        {
            assert!(identify_corner_face(&brep, fa, fb, v).is_err());
        }
    }

    #[test]
    fn face_loop_contains_vertex_positive() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let start_he = brep.loops[brep.faces[face].outer_loop].half_edge;
        let v = brep.half_edges[start_he].origin;
        assert!(face_loop_contains_vertex(&brep, face, v).unwrap());
    }

    #[test]
    fn repair_vertex_outgoing_already_valid() {
        let mut brep = BRep::unit_cube(default_tol()).unwrap();
        let vertex = brep.vertices.keys().next().unwrap();
        assert!(repair_vertex_outgoing(&mut brep, vertex).is_ok());
    }
}
