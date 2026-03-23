#[cfg(test)]
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Create a closed square curve in 3D on the z=1 plane (top face of unit cube).
/// Vertices at (0.25, 0.25, 1), (0.75, 0.25, 1), (0.75, 0.75, 1), (0.25, 0.75, 1).
fn make_closed_square_3d() -> NurbsCurve {
    let pts = vec![
        Point3::new(0.25, 0.25, 1.0),
        Point3::new(0.75, 0.25, 1.0),
        Point3::new(0.75, 0.75, 1.0),
        Point3::new(0.25, 0.75, 1.0),
        Point3::new(0.25, 0.25, 1.0), // close
    ];
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: pts,
        weights: vec![1.0; 5],
    }
}

/// Create a corresponding closed square curve in UV space.
fn make_closed_square_uv() -> NurbsCurve {
    let pts = vec![
        Point3::new(0.25, 0.25, 0.0),
        Point3::new(0.75, 0.25, 0.0),
        Point3::new(0.75, 0.75, 0.0),
        Point3::new(0.25, 0.75, 0.0),
        Point3::new(0.25, 0.25, 0.0), // close
    ];
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: pts,
        weights: vec![1.0; 5],
    }
}

#[test]
fn add_inner_loop_creates_hole() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &Tolerance::default())
        .unwrap();

    assert!(!cube.faces[face].inner_loops.is_empty());
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
}

#[test]
fn add_inner_loop_half_edges_form_closed_chain() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &Tolerance::default())
        .unwrap();

    // Walk the inner loop and verify it forms a closed cycle of length 4.
    let loop_data = &cube.loops[result.loop_];
    let start_he = loop_data.half_edge;
    let mut cursor = start_he;
    let mut count = 0;
    loop {
        count += 1;
        cursor = cube.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
        assert!(count < 100, "cycle did not close");
    }
    assert_eq!(count, 4);

    // Verify prev chain also forms a closed cycle.
    cursor = start_he;
    count = 0;
    loop {
        count += 1;
        cursor = cube.half_edges[cursor].prev;
        if cursor == start_he {
            break;
        }
        assert!(count < 100, "prev cycle did not close");
    }
    assert_eq!(count, 4);
}

#[test]
fn add_inner_loop_too_few_segments_fails() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube.add_inner_loop(face, &curve_3d, &curve_uv, 2, &Tolerance::default());
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("num_segments"));
}

#[test]
fn add_inner_loop_open_curve_fails() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    // An open curve (not closed).
    let open_curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let curve_uv = make_closed_square_uv();

    let result = cube.add_inner_loop(face, &open_curve, &curve_uv, 4, &Tolerance::default());
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("not closed"));
}

#[test]
fn add_inner_loop_invalid_face_fails() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    // Use a default (invalid) face handle.
    let bad_face = Default::default();
    let result = cube.add_inner_loop(bad_face, &curve_3d, &curve_uv, 4, &Tolerance::default());
    assert!(result.is_err());
}

#[test]
fn add_inner_loop_twins_are_paired() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &Tolerance::default())
        .unwrap();

    // Each inner half-edge should have a valid twin, and twin's twin should point back.
    let start_he = cube.loops[result.loop_].half_edge;
    let mut cursor = start_he;
    loop {
        let he = &cube.half_edges[cursor];
        let twin = he.twin;
        assert_ne!(twin, cursor, "half-edge should not be self-twinned");
        assert_eq!(
            cube.half_edges[twin].twin, cursor,
            "twin's twin should point back"
        );
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
}

#[test]
fn add_inner_loop_pcurves_are_set() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &Tolerance::default())
        .unwrap();

    // Every inner half-edge should have a pcurve.
    let start_he = cube.loops[result.loop_].half_edge;
    let mut cursor = start_he;
    loop {
        let he = &cube.half_edges[cursor];
        assert!(he.pcurve.is_some(), "inner half-edge should have a pcurve");
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
}

#[test]
fn add_inner_loop_edges_have_curves() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &Tolerance::default())
        .unwrap();

    for &edge_h in &result.edges {
        assert!(
            cube.edges[edge_h].curve.is_some(),
            "each edge should have curve geometry"
        );
    }
}

#[test]
fn add_inner_loop_vertices_at_sample_points() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let curve_3d = make_closed_square_3d();
    let curve_uv = make_closed_square_uv();
    let tol = Tolerance::default();

    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 4, &tol)
        .unwrap();

    // Vertices should be at the sampled 3D points.
    let (t_start, t_end) = curve_3d.parameter_range();
    for (i, &vh) in result.vertices.iter().enumerate() {
        let t = t_start + (t_end - t_start) * (i as f64) / 4.0;
        let expected = curve_3d.point_at(t);
        let actual = cube.vertices[vh].position;
        assert!(
            tol.points_coincident(&actual, &expected),
            "vertex {i} position mismatch: expected {:?}, got {:?}",
            expected,
            actual
        );
    }
}
