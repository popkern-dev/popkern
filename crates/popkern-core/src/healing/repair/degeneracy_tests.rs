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

#[test]
fn detect_zero_area_faces_on_valid_cube_returns_empty() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let zero_area = detect_zero_area_faces(&brep, &tol);
    assert!(
        zero_area.is_empty(),
        "valid cube should have no zero-area faces"
    );
}

#[test]
fn detect_cusps_on_valid_cube_returns_empty() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let cusps = detect_cusps(&brep, &tol);
    assert!(cusps.is_empty(), "valid cube should have no cusps");
}

#[test]
fn repair_degeneracies_full_on_valid_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = repair_degeneracies_full(&mut brep, &tol);
    assert_eq!(result.zero_length_edges_removed, 0);
    assert_eq!(result.small_faces_collapsed, 0);
    assert_eq!(result.zero_area_faces_detected, 0);
    assert_eq!(result.cusps_detected, 0);
}

#[test]
fn repair_degeneracies_full_on_empty_brep() {
    let mut brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let result = repair_degeneracies_full(&mut brep, &tol);
    assert_eq!(result.zero_length_edges_removed, 0);
    assert_eq!(result.zero_area_faces_detected, 0);
    assert_eq!(result.cusps_detected, 0);
}
