//! End-to-end integration tests for the sketch-to-BRep pipeline.
//!
//! These tests exercise the full workflow: build a sketch, solve
//! constraints, convert to BRep, and verify the result.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{Constraint, SketchOnPlane, SketchPlane};

/// Test 1: Rectangle sketch with constraints -> to_brep -> valid face.
#[test]
fn sketch_rectangle_to_face() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p1 = sop.sketch.add_point(0.0, 0.0);
    let p2 = sop.sketch.add_point(5.0, 0.0);
    let p3 = sop.sketch.add_point(5.0, 3.0);
    let p4 = sop.sketch.add_point(0.0, 3.0);
    let l1 = sop.sketch.add_line(p1, p2);
    let l2 = sop.sketch.add_line(p2, p3);
    let l3 = sop.sketch.add_line(p3, p4);
    let l4 = sop.sketch.add_line(p4, p1);

    sop.sketch.add_constraint(Constraint::Fixed(p1));
    sop.sketch.add_constraint(Constraint::Horizontal(l1));
    sop.sketch.add_constraint(Constraint::Vertical(l2));
    sop.sketch.add_constraint(Constraint::Horizontal(l3));
    sop.sketch.add_constraint(Constraint::Vertical(l4));
    sop.sketch.add_constraint(Constraint::Distance(p1, p2, 5.0));
    sop.sketch.add_constraint(Constraint::Distance(p2, p3, 3.0));

    let solve_result = sop.sketch.solve().unwrap();
    assert!(
        solve_result.residual_norm < 1e-8,
        "solver should converge, residual={}",
        solve_result.residual_norm
    );

    let tol = Tolerance::default();
    let result = sop.to_brep(tol).unwrap();
    assert_eq!(result.regions.len(), 1);
    assert!(result.brep.faces.len() >= 1);
    assert!(result.brep.edges.len() >= 4);
}

/// Test 2: Sketch on a cube face with imported reference geometry.
#[test]
fn sketch_on_cube_face() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let plane = SketchPlane::from_face(&cube, face).unwrap();
    let mut sop = SketchOnPlane::new(plane);

    // Import face edges as construction geometry.
    let imported = sop.import_face_edges(&cube, face).unwrap();
    assert!(!imported.is_empty());

    // Add a small rectangle inside the face.
    let p1 = sop.sketch.add_point(0.25, 0.25);
    let p2 = sop.sketch.add_point(0.75, 0.25);
    let p3 = sop.sketch.add_point(0.75, 0.75);
    let p4 = sop.sketch.add_point(0.25, 0.75);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p4);
    let _l4 = sop.sketch.add_line(p4, p1);

    sop.sketch.add_constraint(Constraint::Fixed(p1));
    sop.sketch.add_constraint(Constraint::Fixed(p2));
    sop.sketch.add_constraint(Constraint::Fixed(p3));
    sop.sketch.add_constraint(Constraint::Fixed(p4));

    sop.sketch.solve().unwrap();
    let result = sop.to_brep(tol).unwrap();
    assert_eq!(result.regions.len(), 1);
}

/// Test 3: Two semicircular arcs forming a full circle -> to_brep -> 1 face region.
#[test]
fn circle_sketch_to_face() {
    use std::f64::consts::PI;

    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Two semicircular arcs forming a closed circle at origin, radius 2.
    let _arc1 = sop.sketch.add_arc(0.0, 0.0, 2.0, 0.0, PI);
    let _arc2 = sop.sketch.add_arc(0.0, 0.0, 2.0, PI, 2.0 * PI);

    let tol = Tolerance::default();
    let result = sop.to_brep(tol).unwrap();

    assert_eq!(
        result.regions.len(),
        1,
        "two semicircles should produce 1 face region"
    );
}

/// Test 4: Outer rectangle + inner circle hole -> 1 region with 1 hole.
#[test]
fn rectangle_with_circle_hole() {
    use std::f64::consts::PI;

    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Outer rectangle: (0,0) - (10,0) - (10,10) - (0,10)
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(10.0, 0.0);
    let p2 = sop.sketch.add_point(10.0, 10.0);
    let p3 = sop.sketch.add_point(0.0, 10.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    // Inner circle via two semicircular arcs: center (5,5), radius 2
    let _arc1 = sop.sketch.add_arc(5.0, 5.0, 2.0, 0.0, PI);
    let _arc2 = sop.sketch.add_arc(5.0, 5.0, 2.0, PI, 2.0 * PI);

    let tol = Tolerance::default();
    let result = sop.to_brep(tol).unwrap();

    assert_eq!(result.regions.len(), 1, "should have 1 face region");
    assert_eq!(
        result.regions[0].holes.len(),
        1,
        "should have 1 hole (the circle)"
    );
}

/// Test 5: Two disjoint rectangles -> 2 separate face regions.
#[test]
fn two_disjoint_rectangles() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Rectangle A: (0,0) - (3,0) - (3,2) - (0,2)
    let a0 = sop.sketch.add_point(0.0, 0.0);
    let a1 = sop.sketch.add_point(3.0, 0.0);
    let a2 = sop.sketch.add_point(3.0, 2.0);
    let a3 = sop.sketch.add_point(0.0, 2.0);
    let _la0 = sop.sketch.add_line(a0, a1);
    let _la1 = sop.sketch.add_line(a1, a2);
    let _la2 = sop.sketch.add_line(a2, a3);
    let _la3 = sop.sketch.add_line(a3, a0);

    // Rectangle B: (5,0) - (8,0) - (8,2) - (5,2)
    let b0 = sop.sketch.add_point(5.0, 0.0);
    let b1 = sop.sketch.add_point(8.0, 0.0);
    let b2 = sop.sketch.add_point(8.0, 2.0);
    let b3 = sop.sketch.add_point(5.0, 2.0);
    let _lb0 = sop.sketch.add_line(b0, b1);
    let _lb1 = sop.sketch.add_line(b1, b2);
    let _lb2 = sop.sketch.add_line(b2, b3);
    let _lb3 = sop.sketch.add_line(b3, b0);

    let tol = Tolerance::default();
    let result = sop.to_brep(tol).unwrap();

    assert_eq!(
        result.regions.len(),
        2,
        "two disjoint rectangles should produce 2 face regions"
    );
}

/// Test 6: Sketch on a tilted plane -> verify 3D geometry lies on that plane.
#[test]
fn tilted_plane_geometry() {
    let normal = Vec3::new(1.0, 1.0, 1.0);
    let origin = Point3::new(0.0, 0.0, 0.0);
    let plane = SketchPlane::from_normal(origin, normal).unwrap();

    let mut sop = SketchOnPlane::new(plane);

    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(2.0, 0.0);
    let p2 = sop.sketch.add_point(2.0, 1.0);
    let p3 = sop.sketch.add_point(0.0, 1.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    sop.sketch.add_constraint(Constraint::Fixed(p0));
    sop.sketch.add_constraint(Constraint::Fixed(p1));
    sop.sketch.add_constraint(Constraint::Fixed(p2));
    sop.sketch.add_constraint(Constraint::Fixed(p3));
    sop.sketch.solve().unwrap();

    let tol = Tolerance::default();
    let result = sop.to_brep(tol).unwrap();

    assert_eq!(result.regions.len(), 1);

    // Verify all vertices satisfy the plane equation:
    // dot(position - origin, normal) = 0
    let n = plane.normal;
    let o = plane.origin;
    for (_, v) in &result.brep.vertices {
        let d = v.position - o;
        let dist = d.dot(&n).abs();
        assert!(
            dist < 1e-8,
            "vertex {:?} should lie on tilted plane, distance={dist}",
            v.position
        );
    }
}
