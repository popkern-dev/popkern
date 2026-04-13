//! Face identification and merge operations for blend construction.
//!
//! Given the two halves produced by [`split_face`], determines which half
//! is "consumed" by the fillet/chamfer, and merges consumed faces via `kef`
//! to create the final blend face.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Round 5r+: descendants-aware consumed-face identifier. When
/// Phase A's perp splits fragment a fillet edge handle (because
/// that edge was the perpendicular of another fillet edge in the
/// same batch), the canonical handle now points to only one of the
/// two descendant sub-edges. This helper accepts the full set of
/// descendant handles and matches any of them — any descendant
/// sub-edge on a face's outer loop proves the face is adjacent to
/// the (original) fillet edge.
///
/// Returns `(kept_face, consumed_face)`.
pub fn identify_consumed_face_any(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    descendants: &[EdgeHandle],
) -> Result<(FaceHandle, FaceHandle), String> {
    if face_loop_contains_any_edge(brep, face_b, descendants)? {
        Ok((face_a, face_b))
    } else if face_loop_contains_any_edge(brep, face_a, descendants)? {
        Ok((face_b, face_a))
    } else {
        Err(format!(
            "neither split face contains any descendant of fillet edge {descendants:?}"
        ))
    }
}

fn face_loop_contains_any_edge(
    brep: &BRep,
    face: FaceHandle,
    edges: &[EdgeHandle],
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
        if edges.contains(&he.edge) {
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

/// Merge two consumed faces by removing a descendant fillet edge
/// via `kef`. Walks `consumed_1`'s outer loop and finds the
/// *first* half-edge whose `edge` is in `descendants` AND whose
/// twin is on `consumed_2` — the key difference from a naive
/// search is that on a cube corner, `consumed_1` may contain
/// several descendant sub-edges, only one of which is adjacent
/// across the fillet edge to `consumed_2`.
///
/// After merging, the surviving face is assigned the fillet
/// surface geometry. Returns the handle of the merged (fillet)
/// face.
pub fn merge_consumed_faces_any(
    brep: &mut BRep,
    consumed_1: FaceHandle,
    consumed_2: FaceHandle,
    descendants: &[EdgeHandle],
    surface: &SurfaceGeometry,
) -> Result<FaceHandle, String> {
    let he_on_consumed_1 =
        find_adjacent_descendant_halfedge(brep, descendants, consumed_1, consumed_2)?;

    brep.kef(he_on_consumed_1)
        .map_err(|e| format!("kef failed: {e}"))?;

    brep.faces[consumed_1].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    Ok(consumed_1)
}

/// Find a half-edge on `consumed_1`'s outer loop whose `edge` is
/// one of `descendants` AND whose twin's loop face equals
/// `consumed_2`. Returns an error if no such half-edge exists.
fn find_adjacent_descendant_halfedge(
    brep: &BRep,
    descendants: &[EdgeHandle],
    consumed_1: FaceHandle,
    consumed_2: FaceHandle,
) -> Result<HalfEdgeHandle, String> {
    let face_data = brep
        .faces
        .get(consumed_1)
        .ok_or_else(|| format!("face {consumed_1:?} not found"))?;
    let outer_loop = face_data.outer_loop;
    let start_he = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?
        .half_edge;

    let mut cursor = start_he;
    let max_iters = brep.half_edges.len() + 1;
    for _ in 0..max_iters {
        let he = brep
            .half_edges
            .get(cursor)
            .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
        if descendants.contains(&he.edge) {
            let twin_handle = he.twin;
            let twin_face = brep.loops[brep.half_edges[twin_handle].loop_].face;
            if twin_face == consumed_2 {
                return Ok(cursor);
            }
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
    Err(format!(
        "no descendant half-edge on {consumed_1:?} has its twin on {consumed_2:?}"
    ))
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

        assert!(
            face_loop_contains_any_edge(&cube, face_of_he, std::slice::from_ref(&edge)).unwrap()
        );
    }

    #[test]
    fn find_adjacent_descendant_halfedge_works() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let edge_data = &cube.edges[edge];
        let he = &cube.half_edges[edge_data.half_edge];
        let face_of_he = cube.loops[he.loop_].face;
        let twin = &cube.half_edges[he.twin];
        let twin_face = cube.loops[twin.loop_].face;

        let found = find_adjacent_descendant_halfedge(
            &cube,
            std::slice::from_ref(&edge),
            face_of_he,
            twin_face,
        )
        .unwrap();
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
            let result = identify_consumed_face_any(
                &cube,
                face_a,
                face_b,
                std::slice::from_ref(&unrelated_edge),
            );
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
    fn find_adjacent_descendant_halfedge_error_when_edge_not_on_face() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();

        // Pick an edge and a face that does NOT contain it.
        let edge = cube.edges.keys().next().unwrap();
        let edge_he = cube.edges[edge].half_edge;
        let face_of_edge = cube.loops[cube.half_edges[edge_he].loop_].face;
        let twin_handle = cube.half_edges[edge_he].twin;
        let twin_face = cube.loops[cube.half_edges[twin_handle].loop_].face;

        // Find a face that is NOT face_of_edge or twin_face.
        let other_face = cube
            .faces
            .keys()
            .find(|&f| f != face_of_edge && f != twin_face)
            .unwrap();

        let result = find_adjacent_descendant_halfedge(
            &cube,
            std::slice::from_ref(&edge),
            other_face,
            twin_face,
        );
        assert!(
            result.is_err(),
            "should fail when edge is not on the face's outer loop"
        );
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
                let contains =
                    face_loop_contains_any_edge(&cube, other_face, std::slice::from_ref(&edge))
                        .unwrap();
                assert!(!contains, "unrelated face should not contain the edge");
            }
        }
    }
}
