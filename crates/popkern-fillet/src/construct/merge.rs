//! Face identification and merge operations for blend construction.
//!
//! Given the two halves produced by [`split_face`], determines which half
//! is "consumed" by the fillet/chamfer, and merges consumed faces via `kef`
//! to create the final blend face.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Determine which of the two split faces is the "consumed" face.
///
/// The consumed face is the one whose outer loop contains the original
/// fillet edge. Returns `(kept_face, consumed_face)`.
pub fn identify_consumed_face(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    fillet_edge: EdgeHandle,
) -> Result<(FaceHandle, FaceHandle), String> {
    if face_loop_contains_edge(brep, face_b, fillet_edge)? {
        Ok((face_a, face_b))
    } else if face_loop_contains_edge(brep, face_a, fillet_edge)? {
        Ok((face_b, face_a))
    } else {
        Err(format!(
            "neither split face contains fillet edge {fillet_edge:?}"
        ))
    }
}

/// Check if a face contains the given edge in its outer loop.
///
/// Used to re-identify consumed faces after Phase B splits may have
/// fragmented them.
pub fn identify_consumed_face_contains(brep: &BRep, face: FaceHandle, edge: EdgeHandle) -> bool {
    face_loop_contains_edge(brep, face, edge).unwrap_or(false)
}

/// Check if a face's outer loop contains a half-edge belonging to the given edge.
fn face_loop_contains_edge(
    brep: &BRep,
    face: FaceHandle,
    edge: EdgeHandle,
) -> Result<bool, String> {
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

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
        if he.edge == edge {
            return Ok(true);
        }
        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Ok(false)
}

/// Find the half-edge of the given edge that lies in the given face's loop.
pub fn find_edge_halfedge_on_face(
    brep: &BRep,
    edge: EdgeHandle,
    face: FaceHandle,
) -> Result<HalfEdgeHandle, String> {
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

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
        if he.edge == edge {
            return Ok(cursor);
        }
        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }
    Err(format!(
        "edge {edge:?} not found in face {face:?}'s outer loop"
    ))
}

/// Merge two consumed faces by removing the original fillet edge via `kef`.
///
/// After merging, the surviving face is assigned the fillet surface geometry.
/// Returns the handle of the merged (fillet) face.
pub fn merge_consumed_faces(
    brep: &mut BRep,
    consumed_1: FaceHandle,
    consumed_2: FaceHandle,
    fillet_edge: EdgeHandle,
    surface: &SurfaceGeometry,
) -> Result<FaceHandle, String> {
    // Find the half-edge of the fillet edge that lies on consumed_1's loop.
    let he_on_consumed_1 = find_edge_halfedge_on_face(brep, fillet_edge, consumed_1)?;

    // Verify the twin is on consumed_2 (sanity check before kef).
    let twin_handle = brep.half_edges[he_on_consumed_1].twin;
    let twin_face = brep.loops[brep.half_edges[twin_handle].loop_].face;
    if twin_face != consumed_2 {
        return Err(format!(
            "fillet edge twin is on face {twin_face:?}, expected {consumed_2:?}"
        ));
    }

    // kef removes the edge and merges the two faces.
    // It kills the twin's loop/face (consumed_2) and keeps consumed_1.
    brep.kef(he_on_consumed_1)
        .map_err(|e| format!("kef failed: {e}"))?;

    // After kef, consumed_2 is removed and consumed_1 survives with the merged loop.
    brep.faces[consumed_1].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    Ok(consumed_1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn identify_consumed_face_finds_correct_side() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_of_he = cube.loops[he.loop_].face;

        assert!(face_loop_contains_edge(&cube, face_of_he, edge).unwrap());
    }

    #[test]
    fn find_edge_halfedge_on_face_works() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_of_he = cube.loops[he.loop_].face;

        let found = find_edge_halfedge_on_face(&cube, edge, face_of_he).unwrap();
        assert_eq!(cube.half_edges[found].edge, edge);
    }

    #[test]
    fn identify_consumed_face_error_when_neither_contains_edge() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        // Collect two faces and an edge that belongs to neither.
        let mut faces_iter = cube.faces.keys();
        let face_a = faces_iter.next().unwrap();
        let face_b = faces_iter.next().unwrap();

        // Find an edge that is NOT on face_a or face_b.
        let edge_not_on_either = cube.edges.keys().find(|&e| {
            let he_handle = cube.edges[e].half_edge;
            let face_1 = cube.loops[cube.half_edges[he_handle].loop_].face;
            let twin_handle = cube.half_edges[he_handle].twin;
            let face_2 = cube.loops[cube.half_edges[twin_handle].loop_].face;
            face_1 != face_a && face_1 != face_b && face_2 != face_a && face_2 != face_b
        });

        if let Some(unrelated_edge) = edge_not_on_either {
            let result = identify_consumed_face(&cube, face_a, face_b, unrelated_edge);
            assert!(
                result.is_err(),
                "should fail when neither face contains the fillet edge"
            );
            let err = result.unwrap_err();
            assert!(
                err.contains("neither split face"),
                "error should mention 'neither split face': {err}"
            );
        }
        // If all edges touch one of the two faces (unlikely for a cube), skip.
    }

    #[test]
    fn find_edge_halfedge_on_face_error_when_edge_not_on_face() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();

        // Pick an edge and a face that does NOT contain it.
        let edge = cube.edges.keys().next().unwrap();
        let edge_he = cube.edges[edge].half_edge;
        let face_of_edge = cube.loops[cube.half_edges[edge_he].loop_].face;

        // Find a face that is NOT face_of_edge.
        let other_face = cube.faces.keys().find(|&f| f != face_of_edge).unwrap();

        // Also verify the twin face is not the other_face.
        let twin_handle = cube.half_edges[edge_he].twin;
        let twin_face = cube.loops[cube.half_edges[twin_handle].loop_].face;
        if twin_face != other_face {
            let result = find_edge_halfedge_on_face(&cube, edge, other_face);
            assert!(
                result.is_err(),
                "should fail when edge is not on the face's outer loop"
            );
            let err = result.unwrap_err();
            assert!(
                err.contains("not found in face"),
                "error should mention 'not found in face': {err}"
            );
        }
    }

    #[test]
    fn face_loop_contains_edge_returns_false_for_unrelated_edge() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();
        let edge_he = cube.edges[edge].half_edge;
        let face_of_edge = cube.loops[cube.half_edges[edge_he].loop_].face;

        // Check a face that does NOT contain this edge.
        if let Some(other_face) = cube.faces.keys().find(|&f| f != face_of_edge) {
            let twin_he = cube.half_edges[edge_he].twin;
            let twin_face = cube.loops[cube.half_edges[twin_he].loop_].face;
            if twin_face != other_face {
                let contains = face_loop_contains_edge(&cube, other_face, edge).unwrap();
                assert!(!contains, "unrelated face should not contain the edge");
            }
        }
    }
}
