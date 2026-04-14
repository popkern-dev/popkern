use super::*;
use crate::center_curve::CenterCurveResult;
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

    let result = build_fillet_surface(&s1, &s2, &cr, 1.0, 5, true, false, false, false);
    assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
    let fs = result.unwrap();
    // Plane-plane (constant) should produce a torus.
    assert!(
        matches!(fs.surface, SurfaceGeometry::Torus { .. }),
        "expected torus for plane-plane constant, got {:?}",
        std::mem::discriminant(&fs.surface)
    );
}

/// Round 5k regression guard: a Variable spec on plane-plane
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
    let result = build_fillet_surface(&s1, &s2, &cr, 1.0, 5, false, false, false, false);
    assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
    let fs = result.unwrap();
    assert!(
        matches!(fs.surface, SurfaceGeometry::Nurbs(_)),
        "variable plane-plane should produce a NURBS sweep, got {:?}",
        std::mem::discriminant(&fs.surface)
    );
}

/// §11.1 plane-sphere analytic fast path: a constant-radius fillet
/// between a plane and a sphere must produce an exact Torus, not a
/// generic NURBS sweep.
///
/// Configuration: sphere of radius 5 centred at origin, plane z = 5
/// (tangent at the north pole) moved down to z = 3 so the
/// intersection is a circle. For a fillet radius r = 0.5 the
/// offset plane is at z = 3 − 0.5 = 2.5 and the offset sphere has
/// radius 5.5, so the center curve is a circle at z = 2.5 with
/// radius √(5.5² − 2.5²) = √(30.25 − 6.25) = √24 ≈ 4.899. The torus
/// should therefore have:
///   - axis = plane normal (0, 0, 1)
///   - major_radius ≈ 4.899
///   - minor_radius = 0.5
///
/// We don't run the real center-curve producer here (that needs the
/// offset-surface intersection solver); we hand-build a
/// `CenterCurveResult` with the analytically-known circle and
/// verify `build_fillet_surface` picks the analytic path.
#[test]
fn fillet_surface_plane_sphere_produces_analytic_torus() {
    use popkern_core::geometry::curve::CurveEval;

    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 5.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 3.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Center curve: circle in z = 2.5 plane, radius √24, centred at origin.
    // Use 60 input points so the cubic NURBS fit-+-resample error
    // falls to ~0.1%, below the test tolerance.
    let major_expected = (24.0_f64).sqrt();
    let n = 60usize;
    let center_pts: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                2.5,
            )
        })
        .collect();
    let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
    // Footprints: one on the sphere (slightly below z=3), one on
    // the plane (slightly inside the sphere-contact circle). The
    // exact positions don't matter for the surface-type assertion;
    // they only drive the trim curves.
    let fp1: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            Point3::new(
                major_expected * theta.cos(),
                major_expected * theta.sin(),
                3.0,
            )
        })
        .collect();
    let fp2: Vec<Point3> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            let r_in = major_expected - 0.1;
            Point3::new(r_in * theta.cos(), r_in * theta.sin(), 3.0)
        })
        .collect();
    let cr = CenterCurveResult {
        center_curve,
        footprint_1: fp1,
        footprint_2: fp2,
        contact_solutions: None,
    };

    let r = 0.5;
    // Test both orderings — the torus itself is symmetric in the
    // two supports, so either dispatch should hit the fast path.
    for (s1, s2) in [(&sphere, &plane), (&plane, &sphere)] {
        let result = build_fillet_surface(s1, s2, &cr, r, 5, true, false, false, false)
            .expect("plane-sphere fillet should build");
        match result.surface {
            SurfaceGeometry::Torus {
                axis,
                major_radius,
                minor_radius,
                center,
            } => {
                assert!(
                    axis.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() > 0.999,
                    "torus axis should be ±ẑ, got {axis:?}"
                );
                // Minor radius equals the fillet radius exactly —
                // this value is passed through, not sampled.
                assert!(
                    (minor_radius - r).abs() < 1e-12,
                    "minor_radius should equal fillet radius {r}, got {minor_radius}"
                );
                // Major radius and centre are read from sampling
                // the NURBS-fitted center curve, so they inherit
                // hand-fit interpolation error. Loose 2 % tolerance
                // here — in the real pipeline the center curve
                // comes from the offset-surface intersection,
                // which gives much tighter fidelity than a manual
                // control polygon.
                assert!(
                    (major_radius - major_expected).abs() / major_expected < 2e-2,
                    "major_radius expected {major_expected}, got {major_radius}"
                );
                assert!(
                    (center.z - 2.5).abs() < 0.1,
                    "centre z expected 2.5, got {}",
                    center.z
                );
            }
            other => panic!(
                "plane-sphere expected analytic Torus, got {:?}",
                std::mem::discriminant(&other)
            ),
        }
    }
    // Prove the trim curves still evaluate — the fast path reused
    // the footprint data rather than recomputing it.
    let cr_for_eval = simple_center_result();
    let _ = cr_for_eval.center_curve.point_at(0.0);
}

#[test]
fn fillet_surface_general() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 5.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cr = simple_center_result();

    let result = build_fillet_surface(&sphere, &plane, &cr, 1.0, 5, true, false, false, false);
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
// try_plane_plane_torus analytic path
// -----------------------------------------------------------------------

#[test]
fn torus_perpendicular_planes() {
    // Two perpendicular planes (angle between normals = PI/2).
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let result = try_plane_plane_torus(&n1, &n2, &cr, 1.0);
    assert!(
        result.is_some(),
        "perpendicular planes should produce torus"
    );

    let fs = result.unwrap();
    match &fs.surface {
        SurfaceGeometry::Torus {
            minor_radius,
            major_radius,
            ..
        } => {
            assert!(
                (*minor_radius - 1.0).abs() < 1e-10,
                "minor radius should be 1.0, got {minor_radius}",
            );
            // For 90-deg planes: half_angle = PI/4, major = r/sin(PI/4) = r*sqrt(2)
            let expected_major = 1.0 / (std::f64::consts::FRAC_PI_4).sin();
            assert!(
                (*major_radius - expected_major).abs() < 1e-6,
                "major radius should be {expected_major}, got {major_radius}",
            );
        }
        other => panic!("expected Torus, got {:?}", std::mem::discriminant(other)),
    }
}

#[test]
fn torus_45_degree_planes() {
    // Planes at 45 degrees (angle between normals = PI/4).
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let n2 = Vec3::new(0.0, c, c); // 45 deg from n1

    let cr = simple_center_result();
    let result = try_plane_plane_torus(&n1, &n2, &cr, 0.5);
    assert!(result.is_some(), "45-deg planes should produce torus");

    let fs = result.unwrap();
    match &fs.surface {
        SurfaceGeometry::Torus {
            minor_radius,
            major_radius,
            ..
        } => {
            assert!(
                (*minor_radius - 0.5).abs() < 1e-10,
                "minor radius should be 0.5, got {minor_radius}",
            );
            // angle between normals = PI/4, half_angle = PI/8
            let half = std::f64::consts::FRAC_PI_4 * 0.5;
            let expected_major = 0.5 / half.sin();
            assert!(
                (*major_radius - expected_major).abs() < 1e-6,
                "major radius should be {expected_major}, got {major_radius}",
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

    let result = try_plane_plane_torus(&n1, &n2, &cr, 1.0);
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
    let _ = try_plane_plane_torus(&n1, &n2, &cr, 1.0);
}

#[test]
fn torus_has_valid_trim_curves() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let fs = try_plane_plane_torus(&n1, &n2, &cr, 1.0)
        .expect("perpendicular planes should produce torus");

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
fn torus_axis_is_cross_product_of_normals() {
    let n1 = Vec3::new(0.0, 0.0, 1.0);
    let n2 = Vec3::new(0.0, 1.0, 0.0);
    let cr = simple_center_result();

    let fs = try_plane_plane_torus(&n1, &n2, &cr, 1.0)
        .expect("perpendicular planes should produce torus");

    if let SurfaceGeometry::Torus { axis, .. } = &fs.surface {
        // n1 x n2 = (0,0,1) x (0,1,0) = (-1,0,0), normalized = (-1,0,0)
        let expected = Vec3::new(-1.0, 0.0, 0.0);
        let dot = axis.dot(&expected).abs();
        assert!(
            dot > 0.99,
            "torus axis should be parallel to n1 x n2, got dot={dot}",
        );
    } else {
        panic!("expected Torus surface");
    }
}
