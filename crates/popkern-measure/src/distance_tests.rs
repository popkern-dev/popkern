use super::*;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

/// Build a unit cube translated by (dx, dy, dz) using the prism constructor.
fn make_translated_cube(dx: f64, dy: f64, dz: f64) -> BRep {
    let profile = [
        Point3::new(dx, dy, dz),
        Point3::new(dx + 1.0, dy, dz),
        Point3::new(dx + 1.0, dy + 1.0, dz),
        Point3::new(dx, dy + 1.0, dz),
    ];
    BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("translated cube")
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

// ---------------------------------------------------------------------------
// face-face distance tests
// ---------------------------------------------------------------------------

/// Find the minimum face-face distance across all pairs between two BReps.
fn min_face_face_distance(brep_a: &BRep, brep_b: &BRep, tol: &Tolerance) -> f64 {
    let mut best = f64::MAX;
    for fa in brep_a.faces.keys() {
        for fb in brep_b.faces.keys() {
            let d = min_distance(
                brep_a,
                TopoHandle::Face(fa),
                brep_b,
                TopoHandle::Face(fb),
                tol,
            )
            .unwrap();
            if d < best {
                best = d;
            }
        }
    }
    best
}

#[test]
fn face_face_parallel_cubes_gap_distance() {
    // Two unit cubes: one at origin, one at x=3. Gap is 2.0.
    // The minimum face-face distance across all pairs should be 2.0
    // (from the right face of cube_a at x=1 to the left face of cube_b at x=3).
    let cube_a = make_cube();
    let cube_b = make_translated_cube(3.0, 0.0, 0.0);
    let tol = Tolerance::default();

    let d = min_face_face_distance(&cube_a, &cube_b, &tol);
    assert!(
        (d - 2.0).abs() < 0.15,
        "parallel cube gap should be ~2.0, got {d}"
    );
}

#[test]
fn face_face_symmetry() {
    // min_distance(face_a, face_b) == min_distance(face_b, face_a) for any pair.
    let cube_a = make_cube();
    let cube_b = make_translated_cube(3.0, 0.0, 0.0);
    let tol = Tolerance::default();
    let fa = cube_a.faces.keys().next().unwrap();
    let fb = cube_b.faces.keys().next().unwrap();

    let d_ab = min_distance(
        &cube_a,
        TopoHandle::Face(fa),
        &cube_b,
        TopoHandle::Face(fb),
        &tol,
    )
    .unwrap();
    let d_ba = min_distance(
        &cube_b,
        TopoHandle::Face(fb),
        &cube_a,
        TopoHandle::Face(fa),
        &tol,
    )
    .unwrap();
    assert!(
        (d_ab - d_ba).abs() < 1e-10,
        "face-face symmetry violated: {d_ab} vs {d_ba}"
    );
}

#[test]
fn face_face_sphere_to_cube_distance() {
    // Unit sphere at origin (radius 1), cube at [3,4]x[0,1]x[0,1].
    // Closest point on sphere is at (1,0.5,0.5) approx,
    // closest point on cube is at (3,0.5,0.5) approx.
    // Min distance across all face pairs should be ~2.0.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let cube_b = make_translated_cube(3.0, 0.0, 0.0);
    let tol = Tolerance::default();

    let d = min_face_face_distance(&sphere, &cube_b, &tol);
    assert!(
        (d - 2.0).abs() < 0.2,
        "sphere-to-cube distance should be ~2.0, got {d}"
    );
}

#[test]
fn face_face_touching_cubes_zero_distance() {
    // Two unit cubes placed side by side (cube_a at origin, cube_b at x=1).
    // The closest face pair should have distance ~0.
    let cube_a = make_cube();
    let cube_b = make_translated_cube(1.0, 0.0, 0.0);
    let tol = Tolerance::default();

    let d = min_face_face_distance(&cube_a, &cube_b, &tol);
    assert!(
        d < 0.05,
        "touching cubes should have ~0 min face-face distance, got {d}"
    );
}

// ---------------------------------------------------------------------------
// edge-face distance tests
// ---------------------------------------------------------------------------

/// Find the minimum edge-face distance across all edge-face pairs from brep_a
/// edges to brep_b faces.
fn min_edge_face_distance(brep_a: &BRep, brep_b: &BRep, tol: &Tolerance) -> f64 {
    let mut best = f64::MAX;
    for ea in brep_a.edges.keys() {
        for fb in brep_b.faces.keys() {
            let d = min_distance(
                brep_a,
                TopoHandle::Edge(ea),
                brep_b,
                TopoHandle::Face(fb),
                tol,
            )
            .unwrap();
            if d < best {
                best = d;
            }
        }
    }
    best
}

#[test]
fn edge_face_same_body_touching() {
    // An edge of a cube lies on an adjacent face, so distance should be 0.
    let cube = make_cube();
    let tol = Tolerance::default();
    let edge = cube.edges.keys().next().unwrap();
    let face = cube.faces.keys().next().unwrap();
    let d = min_distance(
        &cube,
        TopoHandle::Edge(edge),
        &cube,
        TopoHandle::Face(face),
        &tol,
    )
    .unwrap();
    // Edge is on or adjacent to the face, distance <= sqrt(2) for any edge-face
    // pair on a unit cube.
    assert!(d >= 0.0 && d <= 1.5, "distance = {d}");
}

#[test]
fn edge_face_gap_between_cubes() {
    // Cube_a at origin, cube_b at x=5. The minimum edge-face distance
    // should be 4.0 (closest edge of cube_a at x=1, closest face of cube_b
    // at x=5).
    let cube_a = make_cube();
    let cube_b = make_translated_cube(5.0, 0.0, 0.0);
    let tol = Tolerance::default();

    let d = min_edge_face_distance(&cube_a, &cube_b, &tol);
    assert!(
        (d - 4.0).abs() < 0.1,
        "closest edge-to-face distance should be ~4.0, got {d}"
    );
}

#[test]
fn edge_face_symmetry() {
    // min_distance(edge, face) == min_distance(face, edge)
    let cube_a = make_cube();
    let cube_b = make_translated_cube(3.0, 0.0, 0.0);
    let tol = Tolerance::default();
    let edge_a = cube_a.edges.keys().next().unwrap();
    let face_b = cube_b.faces.keys().next().unwrap();

    let d_ef = min_distance(
        &cube_a,
        TopoHandle::Edge(edge_a),
        &cube_b,
        TopoHandle::Face(face_b),
        &tol,
    )
    .unwrap();
    let d_fe = min_distance(
        &cube_b,
        TopoHandle::Face(face_b),
        &cube_a,
        TopoHandle::Edge(edge_a),
        &tol,
    )
    .unwrap();
    assert!(
        (d_ef - d_fe).abs() < 1e-10,
        "edge-face symmetry violated: {d_ef} vs {d_fe}"
    );
}

#[test]
fn edge_face_cylinder_edge_to_cube() {
    // Cylinder at origin (r=1, h=2), cube at [5,6]x[0,1]x[0,1].
    // Min distance from any cylinder edge to any cube face should be ~4.0
    // (closest point on cylinder circle at x=1, closest cube face at x=5).
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let cube_b = make_translated_cube(5.0, 0.0, 0.0);
    let tol = Tolerance::default();

    let d = min_edge_face_distance(&cyl, &cube_b, &tol);
    assert!(
        (d - 4.0).abs() < 0.2,
        "cylinder edge to cube face should be ~4.0, got {d}"
    );
}
