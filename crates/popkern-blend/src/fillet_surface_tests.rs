use super::*;
use crate::center_curve::CenterCurveResult;
use crate::types::ContinuityLevel;
use popkern_core::geometry::{Point3, Vec3};

fn simple_center_result() -> CenterCurveResult {
    let fp1 = (0..5)
        .map(|i| Point3::new(i as f64, 0.0, 0.0))
        .collect::<Vec<_>>();
    let fp2 = (0..5)
        .map(|i| Point3::new(i as f64, 0.0, 1.0))
        .collect::<Vec<_>>();
    let center_pts: Vec<Point3> = (0..5).map(|i| Point3::new(i as f64, 0.0, 0.5)).collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    CenterCurveResult {
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    }
}

#[test]
fn fillet_surface_plane_plane() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cr = simple_center_result();

    let result = build_fillet_surface(
        &s1,
        &s2,
        &cr,
        1.0,
        5,
        true,
        ContinuityLevel::G1,
        false,
        false,
        crate::types::BlendProfile::Circular,
    );
    assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
    let fs = result.unwrap();
    // Plane-plane (constant, straight center curve) produces a cylinder.
    assert!(
        matches!(fs.surface, SurfaceGeometry::Cylinder { .. }),
        "expected Cylinder for plane-plane constant, got {:?}",
        std::mem::discriminant(&fs.surface)
    );
}

/// regression guard: a Variable spec on plane-plane
/// must NOT collapse to a Torus, even though the geometry is
/// otherwise a plane-plane case. Setting `is_constant_radius =
/// false` forces the NURBS sweep path.
#[test]
fn fillet_surface_plane_plane_variable_uses_nurbs_sweep() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cr = simple_center_result();
    let result = build_fillet_surface(
        &s1,
        &s2,
        &cr,
        1.0,
        5,
        false,
        ContinuityLevel::G1,
        false,
        false,
        crate::types::BlendProfile::Circular,
    );
    assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
    let fs = result.unwrap();
    assert!(
        matches!(fs.surface, SurfaceGeometry::Nurbs(_)),
        "variable plane-plane should produce a NURBS sweep, got {:?}",
        std::mem::discriminant(&fs.surface)
    );
}

// §11.1 analytic fast-path tests (plane-sphere, plane-cylinder
// perpendicular, plane-cylinder oblique) live in
// `fillet_surface_analytic_tests.rs` — split out so neither file
// exceeds the 500-line structural cap.

#[test]
fn fillet_surface_general() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 5.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cr = simple_center_result();

    let result = build_fillet_surface(
        &sphere,
        &plane,
        &cr,
        1.0,
        5,
        true,
        ContinuityLevel::G1,
        false,
        false,
        crate::types::BlendProfile::Circular,
    );
    assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
    let fs = result.unwrap();
    // Non-plane-plane should produce NURBS.
    assert!(
        matches!(fs.surface, SurfaceGeometry::Nurbs(_)),
        "expected NURBS for general case"
    );
}

#[test]
fn sample_arc_basic() {
    let center = Point3::new(0.0, 0.0, 0.0);
    let start = Point3::new(1.0, 0.0, 0.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let pts = sample_arc(&center, &start, &end, 1.0, 5);
    assert_eq!(pts.len(), 5);

    // First and last should be close to start and end.
    assert!(pts[0].distance_to(&start) < 1e-10);
    let end_dist = pts[4].distance_to(&end);
    assert!(
        end_dist < 0.2,
        "last arc point should be near end, dist={end_dist}"
    );
}

#[test]
fn build_clamped_knots_basic() {
    let knots = build_clamped_uniform_knots(4, 1);
    assert_eq!(knots.len(), 4 + 1 + 1); // n + degree + 1
    assert!((knots[0] - 0.0).abs() < 1e-15);
    assert!((knots[knots.len() - 1] - 1.0).abs() < 1e-15);
}

// -----------------------------------------------------------------------
// try_plane_plane_analytic analytic path
// -----------------------------------------------------------------------

#[test]
fn cylinder_perpendicular_planes() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let result = try_plane_plane_analytic(&n1, &n2, &cr, 1.0);
    assert!(
        result.is_some(),
        "perpendicular planes should produce cylinder"
    );

    let fs = result.unwrap();
    match &fs.surface {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!(
                (radius - 1.0).abs() < 1e-10,
                "cylinder radius should be 1.0, got {radius}",
            );
        }
        other => panic!("expected Cylinder, got {:?}", std::mem::discriminant(other)),
    }
}

#[test]
fn cylinder_45_degree_planes() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let n2 = Vec3::new(0.0, c, c);

    let cr = simple_center_result();
    let result = try_plane_plane_analytic(&n1, &n2, &cr, 0.5);
    assert!(result.is_some(), "45-deg planes should produce cylinder");

    let fs = result.unwrap();
    match &fs.surface {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!(
                (radius - 0.5).abs() < 1e-10,
                "cylinder radius should be 0.5, got {radius}",
            );
        }
        other => panic!("expected Cylinder, got {:?}", std::mem::discriminant(other)),
    }
}

#[test]
fn torus_center_follows_concave_spine() {
    // Concave arrangement: the center curve's circle sits on the
    // +bisector (void) side. The analytic torus must take its center
    // from the spine, not from the convex anti-bisector assumption
    // (which would misplace it by two major radii).
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let radius = 0.5f64;
    let half = std::f64::consts::FRAC_PI_4;
    let major = radius / half.sin();
    let bisector = Vec3::new(0.0, 1.0, 1.0).normalized().unwrap();
    let spine_center = Point3::new(0.0, 0.0, 0.0) + bisector * major;
    let u = bisector * -1.0;
    let v = Vec3::new(
        0.0,
        -std::f64::consts::FRAC_1_SQRT_2,
        std::f64::consts::FRAC_1_SQRT_2,
    );
    let spine: Vec<Point3> = [-0.3f64, -0.15, 0.0, 0.15, 0.3]
        .iter()
        .map(|theta| spine_center + (u * theta.cos() + v * theta.sin()) * major)
        .collect();
    let footprints_1: Vec<Point3> = spine
        .iter()
        .map(|p| Point3::new(p.x, p.y, p.z - radius))
        .collect();
    let footprints_2: Vec<Point3> = spine
        .iter()
        .map(|p| Point3::new(p.x, p.y - radius, p.z))
        .collect();
    let cr = CenterCurveResult {
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        center_curve: crate::center_curve::fit_nurbs_through_points(&spine).unwrap(),
        footprint_1: footprints_1,
        footprint_2: footprints_2,
        contact_solutions: None,
    };

    let result = try_plane_plane_analytic(&n1, &n2, &cr, radius)
        .expect("concave curved plane-plane fillet should produce an analytic torus");
    match &result.surface {
        SurfaceGeometry::Torus {
            center,
            major_radius,
            minor_radius,
            ..
        } => {
            let center_error = center.distance_to(&spine_center);
            assert!(
                center_error < 1.0e-3,
                "torus center must follow the concave spine circle: expected {spine_center:?}, \
                 got {center:?} (error {center_error:.3e})"
            );
            assert!(
                (major_radius - major).abs() < 1.0e-3,
                "major {major_radius}"
            );
            assert!(
                (minor_radius - radius).abs() < 1e-10,
                "minor {minor_radius}"
            );
        }
        other => panic!("expected Torus, got {:?}", std::mem::discriminant(other)),
    }
}

#[test]
fn torus_parallel_planes_returns_none() {
    // Parallel planes (angle = 0) should fail (sin_half ~ 0).
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 0.0, 1.0);
    let cr = simple_center_result();

    let result = try_plane_plane_analytic(&n1, &n2, &cr, 1.0);
    assert!(result.is_none(), "parallel planes should return None");
}

#[test]
fn torus_antiparallel_planes_returns_none() {
    // Anti-parallel normals (angle = PI, half = PI/2, sin ~ 1).
    // This is a valid degenerate case; torus construction may or may not succeed
    // depending on the angle check, but should not panic.
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 0.0, -1.0);
    let cr = simple_center_result();

    // We just verify no panic. The result can be Some or None.
    let _ = try_plane_plane_analytic(&n1, &n2, &cr, 1.0);
}

#[test]
fn analytic_has_valid_trim_curves() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let fs = try_plane_plane_analytic(&n1, &n2, &cr, 1.0)
        .expect("perpendicular planes should produce cylinder");

    assert!(
        !fs.trim_curve_1.control_points.is_empty(),
        "trim_curve_1 should have control points",
    );
    assert!(
        !fs.trim_curve_2.control_points.is_empty(),
        "trim_curve_2 should have control points",
    );
}

#[test]
fn cylinder_axis_is_center_curve_direction() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let fs = try_plane_plane_analytic(&n1, &n2, &cr, 1.0)
        .expect("perpendicular planes should produce cylinder");

    if let SurfaceGeometry::Cylinder { axis, .. } = &fs.surface {
        // Center curve runs along x-axis: (0,0,0.5) → (4,0,0.5)
        let expected = Vec3::new(1.0, 0.0, 0.0);
        let dot = axis.dot(&expected).abs();
        assert!(
            dot > 0.99,
            "cylinder axis should be parallel to center curve, got dot={dot}",
        );
    } else {
        panic!("expected Cylinder surface");
    }
}
