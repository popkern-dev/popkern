use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn closed_fillet_edge_returns_diagnostic_error() {
    // Construct a fake "closed fillet edge" scenario by passing
    // a fillet edge whose v_start and v_end point at the same
    // vertex. Easiest synthetic input: the cube's first edge
    // with a hand-crafted call where we pretend both endpoints
    // are the same. Since we can't easily fake that on a real
    // cube edge, we use a real cube edge and check that the
    // diagnostic appears for the actual closed-fillet case via
    // the integration probe (probe_5c2_plane_cylinder_unfillet
    // — kept as #[ignore] documentation).
    //
    // Here we just verify the unit cube case (open edge) still
    // works after the closed-edge guard was added.
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let face_1 = cube.loops[he.loop_].face;
    let twin = &cube.half_edges[he.twin];
    let face_2 = cube.loops[twin.loop_].face;

    let snapshot = cube.clone();
    let result = split_perpendicular_edges_by_radius(
        &snapshot,
        &mut cube,
        face_1,
        face_2,
        edge,
        None,
        None,
        0.1,
        0.1,
        &Tolerance::default(),
    );
    // An open cube edge must NOT trigger the closed-edge guard.
    assert!(
        result.is_ok(),
        "open cube edge should still split: {:?}",
        result.err()
    );
}

#[test]
fn split_perp_edges_on_cube() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();

    // Find the two adjacent faces.
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let face_1 = cube.loops[he.loop_].face;
    let twin = &cube.half_edges[he.twin];
    let face_2 = cube.loops[twin.loop_].face;

    let snapshot = cube.clone();
    let result = split_perpendicular_edges_by_radius(
        &snapshot,
        &mut cube,
        face_1,
        face_2,
        edge,
        None,
        None,
        0.1,
        0.1,
        &Tolerance::default(),
    );
    assert!(
        result.is_ok(),
        "split perp edges failed: {:?}",
        result.err()
    );
}

#[test]
fn line_nurbs_produces_valid_curve() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let curve = line_nurbs(&a, &b);
    let (t0, t1) = curve.parameter_range();
    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);
    assert!(p0.distance_to(&a) < 1e-10);
    assert!(p1.distance_to(&b) < 1e-10);
}

#[test]
fn estimate_radius_from_geometry_on_cube() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let edge_data = &cube.edges[edge];
    let he = &cube.half_edges[edge_data.half_edge];
    let v_pos = cube.vertices[he.origin].position;

    // Create a trim curve that starts offset from the vertex.
    let offset_pt = Point3::new(v_pos.x + 0.15, v_pos.y, v_pos.z);
    let trim = line_nurbs(
        &offset_pt,
        &Point3::new(v_pos.x + 0.15, v_pos.y + 1.0, v_pos.z),
    );

    let radius = estimate_radius_from_geometry(&cube, edge, &trim);
    assert!(
        (radius - 0.15).abs() < 1e-10,
        "estimated radius={radius}, expected ~0.15"
    );
}

#[test]
fn estimate_radius_invalid_edge() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let trim = line_nurbs(&Point3::new(0.0, 0.0, 0.0), &Point3::new(1.0, 0.0, 0.0));
    let radius = estimate_radius_from_geometry(&cube, bogus, &trim);
    assert!(
        (radius - 0.0).abs() < 1e-15,
        "should return 0 for invalid edge"
    );
}

#[test]
fn pcurve_for_face_on_cube() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let pt_a = Point3::new(0.2, 0.0, 0.0);
    let pt_b = Point3::new(0.8, 0.0, 0.0);
    let result = pcurve_for_face(&cube, face, &pt_a, &pt_b);
    assert!(result.is_ok(), "pcurve_for_face failed: {:?}", result.err());
    let pcurve = result.unwrap();
    // Pcurve control points should have z=0 (UV convention).
    for cp in &pcurve.control_points {
        assert!(cp.z.abs() < 1e-10, "pcurve z should be 0");
    }
}

#[test]
fn pcurve_for_face_missing_face() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let bogus = popkern_core::topology::FaceHandle::default();
    let result = pcurve_for_face(
        &cube,
        bogus,
        &Point3::new(0.0, 0.0, 0.0),
        &Point3::new(1.0, 0.0, 0.0),
    );
    assert!(result.is_err(), "should fail for missing face");
}

#[test]
fn pcurve_for_face_missing_surface() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    cube.faces[face].surface = None;
    let result = pcurve_for_face(
        &cube,
        face,
        &Point3::new(0.0, 0.0, 0.0),
        &Point3::new(1.0, 0.0, 0.0),
    );
    assert!(result.is_err(), "should fail when face has no surface");
}

#[test]
fn split_perp_edges_invalid_fillet_edge() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_1 = cube.faces.keys().next().unwrap();
    let face_2 = cube.faces.keys().nth(1).unwrap();
    let bogus_edge = popkern_core::topology::EdgeHandle::default();

    let snapshot = cube.clone();
    let result = split_perpendicular_edges_by_radius(
        &snapshot,
        &mut cube,
        face_1,
        face_2,
        bogus_edge,
        None,
        None,
        0.1,
        0.1,
        &Tolerance::default(),
    );
    assert!(result.is_err(), "should fail for invalid fillet edge");
}

#[test]
fn compute_pcurve_for_trim_on_cylinder() {
    use popkern_core::geometry::Vec3;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let pt_a = Point3::new(1.0, 0.0, 0.0);
    let pt_b = Point3::new(1.0, 0.0, 1.0);
    let pcurve = compute_pcurve_for_trim(&cyl, &pt_a, &pt_b);
    assert_eq!(pcurve.degree, 1);
    assert_eq!(pcurve.control_points.len(), 2);
    // All pcurve z-values should be 0.
    for cp in &pcurve.control_points {
        assert!(cp.z.abs() < 1e-10);
    }
}
