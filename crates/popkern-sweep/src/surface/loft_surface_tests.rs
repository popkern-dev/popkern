//! Tests for loft surface generation (split per the file-size rule).

use super::*;
use popkern_core::geometry::surface::SurfaceEval;

/// Create a unit circle as a rational NURBS (degree 2, 9 control points).
fn unit_circle_at_z(z: f64, radius: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, z),
            Point3::new(r, r, z),
            Point3::new(0.0, r, z),
            Point3::new(-r, r, z),
            Point3::new(-r, 0.0, z),
            Point3::new(-r, -r, z),
            Point3::new(0.0, -r, z),
            Point3::new(r, -r, z),
            Point3::new(r, 0.0, z),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

fn make_section(curve: NurbsCurve) -> LoftSection {
    LoftSection {
        curve,
        position: None,
    }
}

#[test]
fn loft_two_identical_circles() {
    let sections = vec![
        make_section(unit_circle_at_z(0.0, 1.0)),
        make_section(unit_circle_at_z(1.0, 1.0)),
    ];
    let options = LoftOptions::default();
    let tol = Tolerance::default();
    let result = generate_loft_surface(&sections, &options, &tol);
    assert!(result.is_ok(), "loft failed: {:?}", result.err());

    let result = result.unwrap(); // ok in tests
    let surf = &result.surface;

    // Should have 9 control points in u (from circle) and 2 in v.
    assert_eq!(surf.num_u, 9);
    assert_eq!(surf.num_v, 2);

    // Evaluate at the midpoint: should be approximately on a unit circle at z=0.5.
    let mid = surf.point_at(0.0, 0.5);
    assert!((mid.z - 0.5).abs() < 0.3, "expected z~0.5, got {}", mid.z);
    let r = (mid.x * mid.x + mid.y * mid.y).sqrt();
    assert!((r - 1.0).abs() < 0.3, "expected radius~1.0, got {}", r);
}

#[test]
fn loft_circle_to_larger_circle() {
    let sections = vec![
        make_section(unit_circle_at_z(0.0, 1.0)),
        make_section(unit_circle_at_z(1.0, 2.0)),
    ];
    let options = LoftOptions::default();
    let tol = Tolerance::default();
    let result = generate_loft_surface(&sections, &options, &tol);
    assert!(result.is_ok(), "loft failed: {:?}", result.err());

    let result = result.unwrap(); // ok in tests
    let surf = &result.surface;

    // At v=0, should be on the smaller circle.
    let p0 = surf.point_at(0.0, 0.0);
    let r0 = (p0.x * p0.x + p0.y * p0.y).sqrt();
    assert!((r0 - 1.0).abs() < 0.3, "expected r~1.0 at v=0, got {}", r0);

    // At v=1, should be on the larger circle.
    let p1 = surf.point_at(0.0, 1.0);
    let r1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
    assert!((r1 - 2.0).abs() < 0.3, "expected r~2.0 at v=1, got {}", r1);
}

#[test]
fn loft_three_sections() {
    let sections = vec![
        make_section(unit_circle_at_z(0.0, 1.0)),
        make_section(unit_circle_at_z(1.0, 1.5)),
        make_section(unit_circle_at_z(2.0, 1.0)),
    ];
    let result = generate_loft_surface(&sections, &LoftOptions::default(), &Tolerance::default());
    let surf = &result.unwrap().surface; // ok in tests
    assert_eq!(surf.num_u, 9);
    assert_eq!(surf.num_v, 3);
    assert!(surf.point_at(0.0, 0.0).z.abs() < 0.3);
    assert!((surf.point_at(0.0, 1.0).z - 2.0).abs() < 0.3);
}

#[test]
fn loft_of_circles_is_exactly_rational() {
    // Homogeneous-space column interpolation: a loft of rational
    // circles must reproduce the circle EXACTLY at every station
    // and stay rational in between — fitting Euclidean control
    // points with weights stamped to 1.0 gave a parabolic sheet
    // 6.1% of the radius off (PK-443).
    let sections = vec![
        make_section(unit_circle_at_z(0.0, 1.0)),
        make_section(unit_circle_at_z(1.0, 1.0)),
    ];
    let result = generate_loft_surface(&sections, &LoftOptions::default(), &Tolerance::default())
        .expect("loft");
    let surf = &result.surface;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    assert!(
        surf.weights.iter().any(|&x| (x - w).abs() < 1e-12),
        "rational weights preserved through the loft"
    );
    for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for i in 0..32 {
            let u = i as f64 / 32.0;
            let p = surf.point_at(u, v);
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-12,
                "radius {r} at (u={u}, v={v}) — the loft left the circle"
            );
        }
    }
}

#[test]
fn loft_too_few_sections_errors() {
    let opts = LoftOptions::default();
    let tol = Tolerance::default();
    // Zero sections.
    let r0 = generate_loft_surface(&[], &opts, &tol);
    assert!(r0.unwrap_err().to_string().contains("at least 2"));
    // One section.
    let r1 = generate_loft_surface(&[make_section(unit_circle_at_z(0.0, 1.0))], &opts, &tol);
    assert!(r1.unwrap_err().to_string().contains("at least 2"));
}

#[test]
fn loft_incompatible_sections_resampled() {
    // Circle (9 CPs, degree 2) vs line segment (2 CPs, degree 1).
    let circle = unit_circle_at_z(0.0, 1.0);
    let line = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(-1.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    let sections = vec![make_section(circle), make_section(line)];
    let options = LoftOptions::default();
    let tol = Tolerance::default();
    let result = generate_loft_surface(&sections, &options, &tol);
    // Should succeed after resampling.
    assert!(result.is_ok(), "loft failed: {:?}", result.err());

    let result = result.unwrap(); // ok in tests
    // Both should have been resampled to RESAMPLE_COUNT points.
    assert_eq!(result.surface.num_u, RESAMPLE_COUNT);
}

#[test]
fn loft_caps_match_sections() {
    let sections = vec![
        make_section(unit_circle_at_z(0.0, 1.0)),
        make_section(unit_circle_at_z(3.0, 1.0)),
    ];
    let result = generate_loft_surface(&sections, &LoftOptions::default(), &Tolerance::default());
    let result = result.unwrap(); // ok in tests
    assert!(result.start_cap.point_at(0.0).z.abs() < 0.01);
    assert!((result.end_cap.point_at(0.0).z - 3.0).abs() < 0.01);
}

#[test]
fn v_params_uniform_for_coincident_centroids() {
    // Two circles at the same location.
    let sections = vec![
        unit_circle_at_z(0.0, 1.0),
        unit_circle_at_z(0.0, 1.0),
        unit_circle_at_z(0.0, 1.0),
    ];
    let options = LoftOptions::default();
    let params = compute_v_params(&sections, &options);
    assert_eq!(params.len(), 3);
    assert!((params[0]).abs() < 1e-10);
    assert!((params[1] - 0.5).abs() < 1e-10);
    assert!((params[2] - 1.0).abs() < 1e-10);
}

#[test]
fn build_clamped_knots_basic() {
    let knots = build_clamped_knots(5, 2);
    // Length = 5 + 2 + 1 = 8
    assert_eq!(knots.len(), 8);
    assert_eq!(knots[0], 0.0);
    assert_eq!(knots[1], 0.0);
    assert_eq!(knots[2], 0.0);
    assert_eq!(knots[7], 1.0);
    assert_eq!(knots[6], 1.0);
    assert_eq!(knots[5], 1.0);
}
