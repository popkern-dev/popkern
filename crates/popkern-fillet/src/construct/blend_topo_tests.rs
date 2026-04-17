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
        // Empty edge_geometries → no trim targets → empty anchors.
        let result = find_face_trim_anchors(&brep, face, v, &[], &default_tol());
        assert!(result.is_ok(), "should return Ok for valid face handle");
        assert!(
            result.unwrap().is_empty(),
            "should return no anchors without trim data"
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
