//! Topology walking helpers for vertex blend insertion.
//!
//! Provides functions to walk around a vertex in the BRep topology,
//! collect adjacent faces, anchor vertices, and incident edges.
//! These are used by [`super::blend_insert`] during blend face insertion.

use std::collections::HashSet;

use popkern_core::tolerance::Tolerance;
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

/// Collect the fillet radii from `edge_geometries` — one per
/// edge, used by [`find_face_trim_anchors`] to match face vertices
/// at those distances from V.
pub(crate) fn collect_fillet_radii(edge_geometries: &[EdgeGeometry]) -> Vec<f64> {
    edge_geometries
        .iter()
        .map(|eg| eg.effective_radius)
        .collect()
}

/// Find the per-face fillet-trim-endpoint anchors for a face
/// incident to `vertex`. Scans the face's outer loop and returns
/// vertices whose distance-from-V matches any of the supplied
/// `radii` within `tolerance.linear() * 10` (the factor covers
/// the mef-chord intersection point `IP` at `r·√2` distance
/// staying outside the match band: IP's distance is `≥
/// r·√(n-1)/n` while a trim vertex is at exactly `r`).
///
/// **Uniqueness by distance** — a fillet trim endpoint on a cube
/// edge sits exactly at `r` from V; chord-crossing points like
/// `IP_A` at `(r_a, r_b, 0)` sit at `√(r_a² + r_b²) ≥ r · √2`
/// which the loose band excludes.
///
/// For a trihedral vertex blend the return is exactly two
/// vertices: the two trim endpoints of the two fillets that
/// bound this face.
pub(crate) fn find_face_trim_anchors(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
    radii: &[f64],
    tolerance: &Tolerance,
) -> Result<Vec<VertexHandle>, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?
        .half_edge;
    let v_pos = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .position;

    // Match band: wide enough to absorb BRep tolerance jitter, but
    // tight enough to exclude IP-crossing points at `r·√(n-1)/n`.
    // Scale with the smallest radius so absolute tolerance tracks
    // the geometry.
    let min_r = radii
        .iter()
        .copied()
        .filter(|r| r.is_finite() && *r > 0.0)
        .fold(f64::INFINITY, f64::min);
    let band = (tolerance.linear() * 100.0).max(min_r * 1e-3);

    let mut anchors: Vec<VertexHandle> = Vec::new();
    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    for _ in 0..max_iters {
        let vh = brep.half_edges[cursor].origin;
        if vh != vertex {
            let vpos = brep.vertices[vh].position;
            let d = v_pos.distance_to(&vpos);
            if radii.iter().any(|&r| (d - r).abs() < band) && !anchors.contains(&vh) {
                anchors.push(vh);
            }
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }

    Ok(anchors)
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
                continuity_downgrade: None,
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
    fn find_face_trim_anchors_vertex_not_in_face() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        if let Some(v) = vertex_not_in_face(&brep, face) {
            // Vertex not in face: `find_face_trim_anchors` walks the
            // face loop against the given `radii` and returns an
            // empty vec (no error path for missing V — that's the
            // caller's responsibility).
            let result = find_face_trim_anchors(&brep, face, v, &[0.1], &default_tol());
            assert!(result.is_ok(), "should return Ok for valid face handle");
            assert!(
                result.unwrap().is_empty(),
                "should return no anchors when V isn't on the face loop"
            );
        }
    }

    #[test]
    fn invalid_vertex_error_paths() {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let bogus = VertexHandle::default();
        assert!(collect_faces_around_vertex(&brep, bogus).is_err());
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
