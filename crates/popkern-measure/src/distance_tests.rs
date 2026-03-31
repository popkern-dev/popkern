use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

#[test]
fn same_vertex_distance_is_zero() {
    let cube = make_cube();
    let v = cube.vertices.keys().next().unwrap();
    let d = min_distance(
        &cube,
        TopoHandle::Vertex(v),
        &cube,
        TopoHandle::Vertex(v),
        &Tolerance::default(),
    )
    .unwrap();
    assert!(d < 1e-10, "self-distance should be 0, got {d}");
}

#[test]
fn vertex_to_vertex_distance() {
    let cube = make_cube();
    let verts: Vec<_> = cube.vertices.keys().collect();
    let d = min_distance(
        &cube,
        TopoHandle::Vertex(verts[0]),
        &cube,
        TopoHandle::Vertex(verts[1]),
        &Tolerance::default(),
    )
    .unwrap();
    assert!(d > 0.0 && d <= 1.74, "distance = {d}");
}

#[test]
fn vertex_to_edge_distance() {
    let cube = make_cube();
    let v = cube.vertices.keys().next().unwrap();
    let e = cube.edges.keys().next().unwrap();
    let d = min_distance(
        &cube,
        TopoHandle::Vertex(v),
        &cube,
        TopoHandle::Edge(e),
        &Tolerance::default(),
    )
    .unwrap();
    assert!(d >= 0.0 && d <= 1.74, "distance = {d}");
}

#[test]
fn edge_to_vertex_symmetry() {
    let cube = make_cube();
    let v = cube.vertices.keys().next().unwrap();
    let e = cube.edges.keys().next().unwrap();
    let tol = Tolerance::default();
    let d_ve = min_distance(
        &cube,
        TopoHandle::Vertex(v),
        &cube,
        TopoHandle::Edge(e),
        &tol,
    )
    .unwrap();
    let d_ev = min_distance(
        &cube,
        TopoHandle::Edge(e),
        &cube,
        TopoHandle::Vertex(v),
        &tol,
    )
    .unwrap();
    assert!(
        (d_ve - d_ev).abs() < 1e-10,
        "symmetry violated: {d_ve} vs {d_ev}"
    );
}

#[test]
fn vertex_to_face_distance() {
    let cube = make_cube();
    let v = cube.vertices.keys().next().unwrap();
    let f = cube.faces.keys().next().unwrap();
    let tol = Tolerance::default();
    let d = min_distance(
        &cube,
        TopoHandle::Vertex(v),
        &cube,
        TopoHandle::Face(f),
        &tol,
    )
    .unwrap();
    assert!(d >= 0.0 && d <= 1.74, "distance = {d}");
}

#[test]
fn edge_to_edge_distance() {
    let cube = make_cube();
    let edges: Vec<_> = cube.edges.keys().collect();
    let tol = Tolerance::default();
    let d = min_distance(
        &cube,
        TopoHandle::Edge(edges[0]),
        &cube,
        TopoHandle::Edge(edges[0]),
        &tol,
    )
    .unwrap();
    assert!(d < 1e-10, "self edge distance should be ~0, got {d}");
}

#[test]
fn face_to_face_distance() {
    let cube = make_cube();
    let faces: Vec<_> = cube.faces.keys().collect();
    let tol = Tolerance::default();
    let d = min_distance(
        &cube,
        TopoHandle::Face(faces[0]),
        &cube,
        TopoHandle::Face(faces[0]),
        &tol,
    )
    .unwrap();
    assert!(d < 1e-10, "self face distance should be ~0, got {d}");
}

#[test]
fn invalid_vertex_handle_returns_error() {
    let cube = make_cube();
    let bad_v = VertexHandle::default();
    let v = cube.vertices.keys().next().unwrap();
    let result = min_distance(
        &cube,
        TopoHandle::Vertex(bad_v),
        &cube,
        TopoHandle::Vertex(v),
        &Tolerance::default(),
    );
    assert!(result.is_err());
}
