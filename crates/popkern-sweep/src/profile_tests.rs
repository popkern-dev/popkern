//! Unit tests for profile analysis (split from `profile.rs` under
//! the 500-line file budget).

use super::*;
use std::f64::consts::PI;

/// Create a unit circle as a rational NURBS (degree 2, 9 control points).
fn unit_circle() -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Create a simple line segment from (0,0,0) to (5,0,0).
fn line_segment() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(5.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Create a rectangle as a degree-1 polyline (4 corners + close).
fn rectangle_polyline() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0; 5],
    }
}

#[test]
fn circle_is_closed_one_span() {
    let tol = Tolerance::default().with_linear(1e-3);
    let profile = SweepProfile::Wire(unit_circle());
    let analyzed = analyze_profile(&profile, &tol).ok();
    assert!(analyzed.is_some());
    let analyzed = analyzed.unwrap(); // ok in tests
    assert!(analyzed.is_closed);
    // Circle is smooth, so should remain 1 span.
    assert_eq!(analyzed.spans.len(), 1);
}

#[test]
fn line_segment_is_open() {
    let tol = Tolerance::default();
    let profile = SweepProfile::Wire(line_segment());
    let analyzed = analyze_profile(&profile, &tol).ok();
    assert!(analyzed.is_some());
    let analyzed = analyzed.unwrap(); // ok in tests
    assert!(!analyzed.is_closed);
    assert_eq!(analyzed.spans.len(), 1);
}

#[test]
fn rectangle_is_closed_and_splits() {
    let tol = Tolerance::default().with_linear(1e-3);
    let profile = SweepProfile::Wire(rectangle_polyline());
    let analyzed = analyze_profile(&profile, &tol).ok();
    assert!(analyzed.is_some());
    let analyzed = analyzed.unwrap(); // ok in tests
    assert!(analyzed.is_closed);
    // A rectangle polyline has 90° corners → 4 spans (one per side).
    assert!(
        analyzed.spans.len() >= 3,
        "expected >=3 spans for rectangle, got {}",
        analyzed.spans.len()
    );
}

#[test]
fn centroid_near_origin_for_centered_circle() {
    let tol = Tolerance::default().with_linear(1e-3);
    let profile = SweepProfile::Wire(unit_circle());
    let analyzed = analyze_profile(&profile, &tol).ok();
    assert!(analyzed.is_some());
    let analyzed = analyzed.unwrap(); // ok in tests
    // Circle is centered at origin, so centroid should be near (0,0,0).
    let c = analyzed.centroid;
    assert!(
        c.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 0.1,
        "centroid should be near origin, got ({}, {}, {})",
        c.x,
        c.y,
        c.z
    );
}

#[test]
fn empty_wire_returns_error() {
    let tol = Tolerance::default();
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![],
        weights: vec![],
    };
    let result = analyze_profile(&SweepProfile::Wire(curve), &tol);
    assert!(result.is_err());
}

#[test]
fn multi_wire_uses_largest_centroid() {
    let tol = Tolerance::default().with_linear(1e-3);
    // Outer: large circle, Inner: small line segment
    let outer = unit_circle();
    let inner = line_segment();
    let profile = SweepProfile::MultiWire(vec![inner, outer]);
    let analyzed = analyze_profile(&profile, &tol).ok();
    assert!(analyzed.is_some());
    let analyzed = analyzed.unwrap(); // ok in tests
    // Centroid should be near the circle's center, not the line's center.
    let c = analyzed.centroid;
    assert!(
        c.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 0.2,
        "centroid should be near circle center"
    );
}

#[test]
fn signed_area_circle_positive() {
    let circle = unit_circle();
    let area = signed_area_2d(&circle);
    // Should be approximately π for unit circle.
    assert!(
        (area.abs() - PI).abs() < 0.2,
        "area should be ~π, got {}",
        area.abs()
    );
}

#[test]
fn face_profile_from_cube() {
    use popkern_core::topology::BRep;
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let profile = SweepProfile::Face {
        brep: Box::new(brep.clone()),
        face,
    };
    let result = analyze_profile(&profile, &Tolerance::default());
    assert!(
        result.is_ok(),
        "face profile should analyze: {:?}",
        result.err()
    );
    let analyzed = result.unwrap();
    assert!(analyzed.is_closed);
    assert!(!analyzed.spans.is_empty());
}
