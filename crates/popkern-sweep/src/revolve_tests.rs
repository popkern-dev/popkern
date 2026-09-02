// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

/// A line segment in the XZ plane, offset from the Z axis by 1 unit.
fn line_profile_xz() -> SweepProfile {
    SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 2.0)],
        weights: vec![1.0, 1.0],
    })
}

/// A closed rectangular profile in the XZ plane for solid revolution.
fn closed_rect_profile() -> SweepProfile {
    SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 5],
    })
}

// ---- basic revolve (line profile) ----

#[test]
fn revolve_line_full() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(!brep.faces.is_empty(), "expected >= 1 face");
    assert!(!brep.edges.is_empty(), "expected edges");
    assert!(!brep.vertices.is_empty(), "expected vertices");
}

#[test]
fn revolve_line_partial() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(!brep.faces.is_empty(), "expected >= 1 face");
}

// ---- closed profile revolution ----

#[test]
fn revolve_closed_partial_has_caps() {
    let result = revolve(
        &closed_rect_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Partial revolution of closed profile: side faces + 2 caps.
    assert!(brep.faces.len() >= 3, "expected >= 3 faces");
}

#[test]
fn revolve_closed_full_no_caps() {
    let result = revolve(
        &closed_rect_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Full revolution wraps: 4 spans -> >= 4 side faces, no caps.
    assert!(brep.faces.len() >= 4, "expected >= 4 faces");
}

// ---- error cases ----

#[test]
fn revolve_zero_angle_fails() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        0.0,
        Tolerance::default(),
        None,
    );
    assert!(result.is_err(), "zero angle should fail");
}

#[test]
fn revolve_negative_angle_fails() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        -1.0,
        Tolerance::default(),
        None,
    );
    assert!(result.is_err(), "negative angle should fail");
}

// ---- various angles ----

#[test]
fn revolve_90_degrees_closed() {
    let result = revolve(
        &closed_rect_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "90deg failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // 4 side + 2 caps = 6.
    assert!(brep.faces.len() >= 6, "90deg: expected >= 6 faces");
}

#[test]
fn revolve_180_degrees_closed() {
    let result = revolve(
        &closed_rect_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "180deg failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(brep.faces.len() >= 3, "180deg: expected >= 3 faces");
}

#[test]
fn revolve_45_degrees_line() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_4,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "45deg failed: {:?}", result.err());
    assert!(!result.unwrap().brep.faces.is_empty()); // ok in tests
}

#[test]
fn revolve_radius_2_surface() {
    // Revolve a line at x=2 around Z axis.
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 3.0)],
        weights: vec![1.0, 1.0],
    });
    let result = revolve(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
}

// ---- different axis positions ----

#[test]
fn revolve_offset_axis() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(-3.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "offset axis failed: {:?}", result.err());
    assert!(!result.unwrap().brep.faces.is_empty()); // ok in tests
}

#[test]
fn revolve_x_axis() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 1.0), Point3::new(3.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    });
    let result = revolve(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        TAU,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "X-axis failed: {:?}", result.err());
    assert!(!result.unwrap().brep.faces.is_empty()); // ok in tests
}

#[test]
fn revolve_y_axis() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 3.0)],
        weights: vec![1.0, 1.0],
    });
    let result = revolve(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 1.0, 0.0),
        PI,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "Y-axis failed: {:?}", result.err());
    assert!(!result.unwrap().brep.faces.is_empty()); // ok in tests
}

// ---- warnings ----

#[test]
fn revolve_no_warnings() {
    let result = revolve(
        &line_profile_xz(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
        Tolerance::default(),
        None,
    );
    assert!(result.is_ok(), "revolve failed: {:?}", result.err());
    let res = result.unwrap(); // ok in tests
    assert!(
        res.warnings.is_empty(),
        "unexpected warnings: {:?}",
        res.warnings
    );
}

// ---- internal helpers ----

#[test]
fn build_cap_curves_returns_span_copies() {
    let spans = vec![NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }];
    let caps = build_cap_curves(&spans);
    assert_eq!(caps.len(), 1);
    assert_eq!(caps[0].control_points.len(), 2);
}

#[test]
fn rotated_caps_90deg_z() {
    let spans = vec![NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }];
    let rotated = build_rotated_cap_curves(
        &spans,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    assert_eq!(rotated.len(), 1);
    // (1,0,0) -> (0,1,0) after 90deg around Z.
    let p = &rotated[0].control_points[0];
    assert!(p.x.abs() < 1e-10, "expected x~0, got {}", p.x);
    assert!((p.y - 1.0).abs() < 1e-10, "expected y~1, got {}", p.y);
}

#[test]
fn rotated_caps_180deg_flips_x() {
    let spans = vec![NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }];
    let rotated = build_rotated_cap_curves(
        &spans,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
    );
    let p = &rotated[0].control_points[0];
    assert!((p.x + 1.0).abs() < 1e-10, "expected x~-1, got {}", p.x);
    assert!(p.y.abs() < 1e-10, "expected y~0, got {}", p.y);
}

#[test]
fn rotated_caps_preserves_degree_and_knots() {
    let spans = vec![NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.5),
            Point3::new(2.0, 0.0, 1.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    }];
    let rotated = build_rotated_cap_curves(
        &spans,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
    );
    assert_eq!(rotated[0].degree, 2);
    assert_eq!(rotated[0].knots, spans[0].knots);
    assert_eq!(rotated[0].weights, spans[0].weights);
    assert_eq!(rotated[0].control_points.len(), 3);
}

#[test]
fn rotated_caps_multiple_spans() {
    let spans = vec![
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        },
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        },
    ];
    let rotated = build_rotated_cap_curves(
        &spans,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    assert_eq!(rotated.len(), 2);
    // Both spans rotated independently.
    assert!(rotated[0].control_points[0].x.abs() < 1e-10);
    assert!(rotated[1].control_points[0].x.abs() < 1e-10);
}

// ---- PK-408: closed profiles serve valid solids ----

fn counts(brep: &popkern_core::topology::BRep) -> (usize, usize, usize) {
    (brep.vertices.len(), brep.edges.len(), brep.faces.len())
}

#[test]
fn full_revolution_of_a_rectangle_is_a_welded_ring() {
    let result = revolve(
        &closed_rect_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
        Tolerance::default(),
        None,
    )
    .expect("ring");
    // Two cylinders with seams, two planar annuli (outer loop + hole).
    assert_eq!(counts(&result.brep), (4, 6, 4), "V E F of a ring");
    result.brep.validate().expect("ring validates");
    // Geometry stays where the profile was: radii 1 and 2, not centred.
    let radii: Vec<f64> = result
        .brep
        .vertices
        .values()
        .map(|v| (v.position.x.powi(2) + v.position.y.powi(2)).sqrt())
        .collect();
    assert!(
        radii
            .iter()
            .all(|r| (r - 1.0).abs() < 1e-9 || (r - 2.0).abs() < 1e-9),
        "{radii:?}"
    );
}

#[test]
fn partial_revolution_of_a_rectangle_validates_with_caps() {
    for angle in [FRAC_PI_2, PI, 3.0 * FRAC_PI_2] {
        let result = revolve(
            &closed_rect_profile(),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            angle,
            Tolerance::default(),
            None,
        )
        .expect("partial");
        assert_eq!(counts(&result.brep), (8, 12, 6), "angle {angle}");
        result
            .brep
            .validate()
            .unwrap_or_else(|e| panic!("angle {angle}: {e}"));
    }
}

#[test]
fn full_revolution_of_a_profile_touching_the_axis_is_a_disc() {
    let profile = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.5),
            Point3::new(0.0, 0.0, 0.5),
            Point3::new(0.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 5],
    });
    let result = revolve(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
        Tolerance::default(),
        None,
    )
    .expect("disc");
    // The span on the axis sweeps nothing and is dropped; the planar
    // discs are single circle loops, so no pole vertex is needed: a
    // cylinder of three faces.
    assert_eq!(counts(&result.brep), (2, 3, 3));
    result.brep.validate().expect("disc validates");
}

#[test]
fn clockwise_and_counter_clockwise_profiles_both_serve_outward() {
    let ccw = closed_rect_profile();
    let SweepProfile::Wire(curve) = &ccw else {
        unreachable!()
    };
    let cw = SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots: curve.knots.clone(),
        control_points: curve.control_points.iter().rev().copied().collect(),
        weights: curve.weights.clone(),
    });
    for (label, profile) in [("ccw", &ccw), ("cw", &cw)] {
        let result = revolve(
            profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            PI,
            Tolerance::default(),
            None,
        )
        .expect(label);
        result
            .brep
            .validate()
            .unwrap_or_else(|e| panic!("{label}: {e}"));
    }
}
