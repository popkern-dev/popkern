//! End-to-end pipeline tests for sketch-to-BRep conversion.

use popkern_sketch::{SketchOnPlane, SketchPlane, Tolerance};

#[test]
fn rectangle_to_brep_produces_valid_region() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(2.0, 0.0);
    let p2 = sop.sketch.add_point(2.0, 1.0);
    let p3 = sop.sketch.add_point(0.0, 1.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    // 1 face region.
    assert_eq!(result.regions.len(), 1);
    // Front face exists.
    assert!(result.brep.faces.contains_key(result.regions[0].face));
    // 4 vertices, 4 edges.
    assert_eq!(result.brep.vertices.len(), 4);
    assert_eq!(result.brep.edges.len(), 4);
    // All edges have curves.
    for (_, edge) in &result.brep.edges {
        assert!(edge.curve.is_some());
    }
    // All half-edges have pcurves.
    for (_, he) in &result.brep.half_edges {
        assert!(he.pcurve.is_some());
    }
}

#[test]
fn triangle_to_brep_euler_passes() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(3.0, 0.0);
    let p2 = sop.sketch.add_point(1.5, 2.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p0);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert_eq!(result.regions.len(), 1);
    assert!(result.brep.validate_euler_poincare().is_ok());
}

#[test]
fn no_edges_returns_error() {
    let sop = SketchOnPlane::new(SketchPlane::xy());
    let result = sop.to_brep(Tolerance::default());
    assert!(result.is_err());
}

#[test]
fn open_wire_returns_no_closed_loops() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let _l = sop.sketch.add_line(p0, p1);

    let result = sop.to_brep(Tolerance::default());
    assert!(result.is_err());
}

#[test]
fn rectangle_with_arc_edge() {
    use std::f64::consts::FRAC_PI_2;
    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Rectangle-like shape: two lines + two arcs
    // Line: (0,0) -> (2,0)
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(2.0, 0.0);
    let _l0 = sop.sketch.add_line(p0, p1);

    // Arc: center (2,1), radius 1, from -PI/2 to PI/2 -> (2,0) to (2,2)
    let _arc1 = sop.sketch.add_arc(2.0, 1.0, 1.0, -FRAC_PI_2, FRAC_PI_2);

    // Line: (2,2) -> (0,2)
    let p2 = sop.sketch.add_point(2.0, 2.0);
    let p3 = sop.sketch.add_point(0.0, 2.0);
    let _l1 = sop.sketch.add_line(p2, p3);

    // Arc: center (0,1), radius 1, from PI/2 to 3*PI/2 -> (0,2) to (0,0)
    let _arc2 = sop
        .sketch
        .add_arc(0.0, 1.0, 1.0, FRAC_PI_2, FRAC_PI_2 + std::f64::consts::PI);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.brep.edges.len(), 4);
}

#[test]
fn rectangle_on_xz_plane() {
    let mut sop = SketchOnPlane::new(SketchPlane::xz());
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let p2 = sop.sketch.add_point(1.0, 1.0);
    let p3 = sop.sketch.add_point(0.0, 1.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert_eq!(result.regions.len(), 1);

    // Verify that vertex positions are in the XZ plane (Y = 0).
    for (_, v) in &result.brep.vertices {
        assert!(
            v.position.y.abs() < 1e-10,
            "vertex should be in XZ plane, got y={}",
            v.position.y
        );
    }
}

#[test]
fn rectangle_with_hole_pipeline() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Outer rectangle: (0,0) - (6,0) - (6,6) - (0,6)
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(6.0, 0.0);
    let p2 = sop.sketch.add_point(6.0, 6.0);
    let p3 = sop.sketch.add_point(0.0, 6.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    // Inner rectangle: (2,2) - (4,2) - (4,4) - (2,4)
    let q0 = sop.sketch.add_point(2.0, 2.0);
    let q1 = sop.sketch.add_point(4.0, 2.0);
    let q2 = sop.sketch.add_point(4.0, 4.0);
    let q3 = sop.sketch.add_point(2.0, 4.0);
    let _m0 = sop.sketch.add_line(q0, q1);
    let _m1 = sop.sketch.add_line(q1, q2);
    let _m2 = sop.sketch.add_line(q2, q3);
    let _m3 = sop.sketch.add_line(q3, q0);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert_eq!(result.regions.len(), 1, "should have 1 face region");
    assert_eq!(
        result.regions[0].holes.len(),
        1,
        "should have 1 hole in the region"
    );
}

#[test]
fn warnings_from_degenerate_edge() {
    let mut sop = SketchOnPlane::new(SketchPlane::xy());

    // Rectangle + a degenerate (zero-length) line.
    let p0 = sop.sketch.add_point(0.0, 0.0);
    let p1 = sop.sketch.add_point(1.0, 0.0);
    let p2 = sop.sketch.add_point(1.0, 1.0);
    let p3 = sop.sketch.add_point(0.0, 1.0);
    let _l0 = sop.sketch.add_line(p0, p1);
    let _l1 = sop.sketch.add_line(p1, p2);
    let _l2 = sop.sketch.add_line(p2, p3);
    let _l3 = sop.sketch.add_line(p3, p0);

    // Degenerate line at same point.
    let d0 = sop.sketch.add_point(5.0, 5.0);
    let d1 = sop.sketch.add_point(5.0, 5.0);
    let _ld = sop.sketch.add_line(d0, d1);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert!(
        !result.warnings.is_empty(),
        "should have a degenerate edge warning"
    );
    assert_eq!(result.regions.len(), 1);
}

#[test]
fn full_circle_to_brep() {
    use std::f64::consts::PI;

    let mut sop = SketchOnPlane::new(SketchPlane::xy());
    let _arc = sop.sketch.add_arc(0.0, 0.0, 2.0, 0.0, 2.0 * PI);

    let result = sop.to_brep(Tolerance::default()).unwrap();

    assert_eq!(
        result.regions.len(),
        1,
        "full circle should produce 1 face region"
    );
    // Full circle split into 2 semicircular edges, 2 vertices.
    assert_eq!(result.brep.edges.len(), 2);
    assert_eq!(result.brep.vertices.len(), 2);
    assert!(result.warnings.is_empty(), "no warnings expected");
}
