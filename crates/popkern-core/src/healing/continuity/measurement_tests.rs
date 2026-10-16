use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// Two identical planes should achieve G3 with all deviations zero.
#[test]
fn identical_planes_g3() {
    let plane = xy_plane();
    let samples: Vec<BoundarySample> = (0..10)
        .map(|i| {
            let t = i as f64 / 9.0;
            BoundarySample {
                u_a: t,
                v_a: 0.0,
                u_b: t,
                v_b: 0.0,
            }
        })
        .collect();

    let tol = Tolerance::default();
    let result = cross_boundary_continuity(&plane, &plane, &samples, &tol);

    assert_eq!(result.achieved_level, ContinuityLevel::G3);
    assert!(result.g0_max < 1e-15, "g0_max: {}", result.g0_max);
    assert!(result.g0_mean < 1e-15, "g0_mean: {}", result.g0_mean);
    assert!(result.g1_max_angle < 1e-15, "g1: {}", result.g1_max_angle);
    assert!(
        result.g2_max_curvature_diff < 1e-15,
        "g2: {}",
        result.g2_max_curvature_diff
    );
    assert!(
        result.g3_max_curvature_rate_diff < 1e-10,
        "g3: {}",
        result.g3_max_curvature_rate_diff
    );
    assert_eq!(result.sample_count, 10);
}

/// Two planes with a 45-degree angular offset should be G0 but not G1.
#[test]
fn angled_planes_g0_not_g1() {
    let plane_a = xy_plane();
    let angle = std::f64::consts::FRAC_PI_4;
    let plane_b = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, -angle.sin(), angle.cos()),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, angle.cos(), angle.sin()),
    };

    let samples: Vec<BoundarySample> = (0..10)
        .map(|i| {
            let t = i as f64 / 9.0;
            BoundarySample {
                u_a: t,
                v_a: 0.0,
                u_b: t,
                v_b: 0.0,
            }
        })
        .collect();

    let tol = Tolerance::default();
    let result = cross_boundary_continuity(&plane_a, &plane_b, &samples, &tol);

    assert!(
        result.g0_max < 1e-12,
        "planes share boundary: g0_max={}",
        result.g0_max
    );
    assert!(
        result.g1_max_angle > 0.1,
        "should detect angular deviation: {}",
        result.g1_max_angle
    );
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn empty_samples_returns_zero_stats() {
    let plane = xy_plane();
    let tol = Tolerance::default();
    let result = cross_boundary_continuity(&plane, &plane, &[], &tol);
    assert_eq!(result.sample_count, 0);
    assert!(result.g0_max < 1e-15);
}

#[test]
fn percentile_95_basic() {
    let mut values: Vec<f64> = (0..100).map(|i| i as f64).collect();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p95 = super::percentile_95(&values);
    assert!((p95 - 95.0).abs() < 1e-10, "p95={p95}");
}

#[test]
fn percentile_95_empty() {
    assert!((super::percentile_95(&[]) - 0.0).abs() < 1e-15);
}

#[test]
fn plane_has_no_internal_discontinuities() {
    let plane = xy_plane();
    let tol = Tolerance::default();
    let records = internal_discontinuities(&plane, ContinuityLevel::G2, &tol, 8);
    assert!(
        records.is_empty(),
        "Plane should have no internal discontinuities, got {}",
        records.len()
    );
}

#[test]
fn sphere_has_no_internal_discontinuities_with_wide_tolerance() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let tol = Tolerance::default().with_angular(1.0);
    let records = internal_discontinuities(&sphere, ContinuityLevel::G1, &tol, 20);
    assert!(
        records.is_empty(),
        "Sphere should have no G1 discontinuities at wide tolerance, got {}",
        records.len()
    );
}
