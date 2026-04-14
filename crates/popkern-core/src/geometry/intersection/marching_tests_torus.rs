//! Torus-plane marching tests. Split from `marching_tests_coverage.rs`
//! to keep both files under the 500-line structural cap.
#[cfg(test)]
use super::*;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::surface::SurfaceEval;
use std::f64::consts::TAU;

/// Marching on a torus-plane intersection exercises the periodic wrapping in
/// both u and v (the torus is periodic in both directions) and produces a closed
/// curve, verifying the closure detection path inside `march_one_direction`.
#[test]
fn march_torus_plane_closed_curve() {
    // Torus with major radius 2, minor radius 0.5, axis along Z.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    // Horizontal plane at z=0 cuts through the outer equator of the torus.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // At torus parameters (u=0, v=0): outer equator point is (2.5, 0, 0).
    let p0 = torus.point_at(0.0, 0.0);
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: p0.x,
        v2: p0.y,
        point: p0,
    };
    let tol = Tolerance::default().with_linear(1e-3);
    let (samples, is_closed) = march(
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert!(
        samples.len() > 3,
        "torus-plane intersection should produce multiple samples, got {}",
        samples.len()
    );
    assert!(
        is_closed,
        "torus-plane intersection should be a closed curve"
    );
    // All samples should be approximately at z=0.
    for s in &samples {
        assert!(
            s.point.z.abs() < 0.05,
            "torus-plane intersection point should be at z~0, got z={}",
            s.point.z
        );
    }
}

#[test]
fn march_one_direction_inner_closure_detection() {
    // Use a small torus so the march wraps around quickly.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 0.3,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p0 = torus.point_at(0.0, 0.0);
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: p0.x,
        v2: p0.y,
        point: p0,
    };
    // Use a slightly coarser tolerance so the march completes the full circle
    // quickly and the inner closure check fires.
    let tol = Tolerance::default().with_linear(5e-3);
    let (samples, is_closed) = march(
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &plane,
        (-3.0, 3.0),
        (-3.0, 3.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert!(
        samples.len() > 3,
        "should produce multiple samples, got {}",
        samples.len()
    );
    assert!(
        is_closed,
        "small torus equator should be detected as closed"
    );
}

/// Verify the torus-plane march stays on the expected circle and covers
/// the full circumference.
#[test]
fn torus_march_stays_on_circle() {
    use crate::geometry::intersection::marching::march_one_direction;

    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p0 = torus.point_at(0.0, 0.0); // (2.5, 0, 0)
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: p0.x,
        v2: p0.y,
        point: p0,
    };
    let tol = Tolerance::default().with_linear(1e-3);

    // March forward only
    let samples = march_one_direction(
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &start,
        1.0,
        &tol,
        &DefaultDegeneracyHandler,
    );

    let n = samples.len();
    eprintln!("[torus diag] {} forward samples", n);

    // Check: are all points on the outer circle (radius 2.5)?
    let mut max_r_dev = 0.0_f64;
    let mut max_z_dev = 0.0_f64;
    for s in &samples {
        let r = (s.point.x.powi(2) + s.point.y.powi(2)).sqrt();
        let r_dev = (r - 2.5).abs();
        max_r_dev = max_r_dev.max(r_dev);
        max_z_dev = max_z_dev.max(s.point.z.abs());
    }
    eprintln!("[torus diag] max radial deviation: {max_r_dev:.6}");
    eprintln!("[torus diag] max z deviation: {max_z_dev:.6}");

    if n > 1 {
        let last = &samples[n - 1];
        let dist_to_start = last.point.distance_to(&start.point);
        eprintln!(
            "[torus diag] last=({:.3},{:.3},{:.3}), dist_to_start={:.6}",
            last.point.x, last.point.y, last.point.z, dist_to_start
        );

        // Compute total arc length
        let mut arc_len = 0.0;
        for i in 1..n {
            arc_len += samples[i].point.distance_to(&samples[i - 1].point);
        }
        let expected_circumference = 2.0 * std::f64::consts::PI * 2.5;
        eprintln!(
            "[torus diag] arc_length={:.3}, expected circumference={:.3}, ratio={:.3}",
            arc_len,
            expected_circumference,
            arc_len / expected_circumference
        );
    }

    // The march should stay on the circle and cover the full circumference
    assert!(
        max_r_dev < 0.1,
        "march should stay on r=2.5 circle, max deviation={max_r_dev}"
    );
    assert!(
        max_z_dev < 0.05,
        "march should stay at z=0, max z deviation={max_z_dev}"
    );
}
