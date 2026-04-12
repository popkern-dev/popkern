//! Degeneracy repair: remove zero-length edges and collapse small faces.
//!
//! Handles degenerate edges at surface singularities (poles, collapsed
//! parameter lines) by detecting edges whose endpoints are coincident
//! and merging the vertex positions.

use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

/// Result of a degeneracy repair operation.
#[derive(Debug, Clone)]
pub struct DegeneracyRepairResult {
    /// Number of zero-length edges whose vertex positions were merged.
    pub zero_length_edges_removed: usize,
    /// Number of small faces collapsed (currently always 0; reserved for
    /// future iterations).
    pub small_faces_collapsed: usize,
}

/// Remove zero-length edges from a BRep.
///
/// For each edge whose two endpoint vertices are coincident (distance <
/// tolerance), removes the edge geometrically by merging the two vertex
/// positions into one.  This handles degenerate edges at surface
/// singularities (poles, collapsed parameter lines).
///
/// Edges where both half-edges originate from the same vertex (seam edges
/// with identical origin) are skipped -- they are topological artefacts,
/// not degeneracies.
///
/// **Note:** this merges vertex *positions* only.  Full topological removal
/// of the edge (via Euler operators) is deferred to a future iteration
/// because it requires careful half-edge relinking.
///
/// Returns the number of edges whose vertex positions were merged.
pub fn remove_zero_length_edges(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    // Collect candidate edges first to avoid borrowing conflicts.
    let candidates: Vec<EdgeHandle> = brep.edges.keys().collect();
    let linear_tol = tolerance.linear();
    let mut count = 0usize;

    for edge_handle in candidates {
        let Some(edge_data) = brep.edges.get(edge_handle) else {
            continue;
        };
        let he_handle = edge_data.half_edge;

        let Some(he) = brep.half_edges.get(he_handle) else {
            continue;
        };
        let start_vertex = he.origin;
        let twin_handle = he.twin;

        let Some(twin) = brep.half_edges.get(twin_handle) else {
            continue;
        };
        let end_vertex = twin.origin;

        // Skip seam edges where both half-edges originate from the same vertex.
        if start_vertex == end_vertex {
            continue;
        }

        let Some(start_data) = brep.vertices.get(start_vertex) else {
            continue;
        };
        let start_pos = start_data.position;

        let Some(end_data) = brep.vertices.get(end_vertex) else {
            continue;
        };
        let end_pos = end_data.position;

        // Check if the two endpoints are within linear tolerance.
        if start_pos.distance_to(&end_pos) >= linear_tol {
            continue;
        }

        // Merge: move the end vertex to the start vertex's position.
        if let Some(v) = brep.vertices.get_mut(end_vertex) {
            v.position = start_pos;
        }
        count += 1;
    }

    count
}

/// Repair all degeneracy defects in a BRep.
///
/// Removes zero-length edges (by merging vertex positions) and collapses
/// small faces.  Small-face collapse is reserved for a future iteration;
/// currently only zero-length edge repair is performed.
pub fn repair_degeneracies(brep: &mut BRep, tolerance: &Tolerance) -> DegeneracyRepairResult {
    let zero_length_edges_removed = remove_zero_length_edges(brep, tolerance);

    DegeneracyRepairResult {
        zero_length_edges_removed,
        small_faces_collapsed: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn valid_cube_has_no_zero_length_edges() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let removed = remove_zero_length_edges(&mut brep, &tol);
        assert_eq!(removed, 0, "valid cube should have no zero-length edges");
    }

    #[test]
    fn repair_degeneracies_on_valid_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();
        let result = repair_degeneracies(&mut brep, &tol);
        assert_eq!(result.zero_length_edges_removed, 0);
        assert_eq!(result.small_faces_collapsed, 0);
    }

    #[test]
    fn repair_degeneracies_on_empty_brep() {
        let mut brep = BRep::new(Tolerance::default());
        let tol = Tolerance::default();
        let result = repair_degeneracies(&mut brep, &tol);
        assert_eq!(result.zero_length_edges_removed, 0);
        assert_eq!(result.small_faces_collapsed, 0);
    }

    #[test]
    fn zero_length_edge_is_detected_and_merged() {
        // Build a cube and then move one edge's end vertex to match
        // its start vertex, creating a zero-length edge.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let tol = Tolerance::default();

        // Pick any edge and collapse it by moving its end vertex.
        let edge_handle = brep.edges.keys().next().unwrap();
        let edge_data = brep.edges.get(edge_handle).unwrap();
        let he = brep.half_edges.get(edge_data.half_edge).unwrap();
        let start_vertex = he.origin;
        let twin = brep.half_edges.get(he.twin).unwrap();
        let end_vertex = twin.origin;

        // Move end vertex to start vertex position (create a zero-length edge).
        let start_pos = brep.vertices.get(start_vertex).unwrap().position;
        if let Some(v) = brep.vertices.get_mut(end_vertex) {
            v.position = start_pos;
        }

        let removed = remove_zero_length_edges(&mut brep, &tol);
        // The edge we collapsed plus any other edges that now have both
        // endpoints at the same position (the moved vertex is shared by
        // multiple edges on the cube).
        assert!(
            removed >= 1,
            "should detect at least the collapsed edge, got {removed}"
        );
    }
}
