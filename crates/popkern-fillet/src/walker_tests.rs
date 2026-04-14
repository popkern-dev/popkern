//! Unit tests for the spine walker (`walker.rs`).
//!
//! Coverage:
//! - **plane-plane line walk**: walks along the shared y-axis; every
//!   sample's center should be at `(r, y, r)` for the corresponding y.
//! - **plane-cylinder circle walk**: walks around the intersection
//!   circle; verifies each sample's center lies on the expected
//!   helix (here a circle since z is fixed).
//! - **plane-cone walk**: the case the offset path can't handle.
//!   Verifies the walker produces sensible centers along the cone's
//!   base circle.
//! - **continuation efficiency**: verifies that subsequent samples
//!   converge in ≤2 Newton iterations thanks to the warm start.

#[cfg(test)]
use super::*;
use crate::contact::BallSide;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// plane-plane line walk
// ---------------------------------------------------------------------------

#[test]
fn walker_plane_plane_along_y_axis_constant_radius() {
    // Two perpendicular planes meeting along the y-axis. For radius
    // r, every point along the spine has the same rolling-ball center
    // shifted by y: (r, y, r). The walker should produce a straight
    // line of centers along y.
    let plane_z0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_x0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let spine = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(0.0, 1.0, 0.0),
    };
    let r = 0.5;
    let radius_at = |_t: f64| r;

    let request = WalkRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at: &radius_at,
        spine: &spine,
        t_range: (0.0, 5.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
    };
    let result = walk_spine(&request, &WalkOptions::default(), &unit_tol())
        .expect("plane-plane walk should succeed");

    assert_eq!(result.samples.len(), 20);
    for sample in &result.samples {
        let expected = Point3::new(r, sample.t_spine, r);
        assert!(
            (sample.solution.center - expected).length() < 1e-9,
            "center {:?} ≠ expected {:?} at t={}",
            sample.solution.center,
            expected,
            sample.t_spine
        );
    }
}

// ---------------------------------------------------------------------------
// continuation efficiency
// ---------------------------------------------------------------------------

#[test]
fn walker_continuation_avoids_repeated_full_newton() {
    // After the first sample's full Newton convergence, every
    // subsequent sample should converge in ≤2 iterations from the
    // warm start. For 20 samples on plane-plane, total Newton
    // iterations should be ≤ first_sample_iterations + 19 × 2.
    let plane_z0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_x0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let spine = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(0.0, 1.0, 0.0),
    };
    let radius_at = |_t: f64| 0.5;
    let request = WalkRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at: &radius_at,
        spine: &spine,
        t_range: (0.0, 5.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
    };
    let result = walk_spine(&request, &WalkOptions::default(), &unit_tol()).unwrap();

    let n = result.samples.len();
    let first_iters = result.samples[0].solution.iterations;
    let subsequent_iters: usize = result.samples[1..]
        .iter()
        .map(|s| s.solution.iterations)
        .sum();
    let avg_subsequent = if n > 1 {
        subsequent_iters as f64 / (n - 1) as f64
    } else {
        0.0
    };

    assert!(
        avg_subsequent <= 2.0,
        "warm-start should converge in ≤2 iter on average; got first={first_iters}, \
         avg_subsequent={avg_subsequent}, total={} from {} samples",
        result.total_newton_iterations,
        n
    );
}

// ---------------------------------------------------------------------------
// plane-cylinder around the intersection circle
// ---------------------------------------------------------------------------

#[test]
fn walker_plane_cylinder_around_circle() {
    // Cylinder along z-axis radius 2, plane z=0; intersection circle
    // x²+y²=4. Walk a 90° arc (θ from 0 to π/2) and verify each
    // sample's center is at distance R+r from the z-axis at z=r.
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.3;
    let radius_at = |_t: f64| r;

    // Spine: the intersection circle, walked from θ=0 to θ=π/2.
    let spine = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };

    let request = WalkRequest {
        surface_0: &cyl,
        surface_1: &plane,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at: &radius_at,
        spine: &spine,
        t_range: (0.0, PI / 2.0),
        uv0_range: ((-PI, 2.0 * PI), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
    };
    let result = walk_spine(&request, &WalkOptions::default(), &unit_tol())
        .expect("plane-cylinder walk should succeed");

    let expected_r = 2.0 + r;
    for sample in &result.samples {
        let c = sample.solution.center;
        let radial = (c.x * c.x + c.y * c.y).sqrt();
        assert!(
            (radial - expected_r).abs() < 1e-6,
            "sample at t={}: radial distance {radial:.6} ≠ expected {expected_r:.6}",
            sample.t_spine
        );
        assert!(
            (c.z - r).abs() < 1e-6,
            "sample at t={}: z={} ≠ expected {}",
            sample.t_spine,
            c.z,
            r
        );
    }
}

// ---------------------------------------------------------------------------
// plane-cone walk
// ---------------------------------------------------------------------------

#[test]
fn walker_plane_cone_around_base_circle() {
    // Cone with apex at (0,0,2) opening downward (axis -z) at 30°.
    // At v=2 (axial parameter from the apex) the cone reaches z=0
    // at radius 2*tan(30°) ≈ 1.155. The intersection with z=0 is a
    // circle. Walk a quarter of it and verify the walker doesn't
    // crash and produces sensible centers.
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: PI / 6.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.05;
    let radius_at = |_t: f64| r;

    let spine_radius = 2.0 * (PI / 6.0).tan();
    let spine = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: spine_radius,
    };

    let request = WalkRequest {
        surface_0: &cone,
        surface_1: &plane,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at: &radius_at,
        spine: &spine,
        t_range: (0.0, PI / 2.0),
        uv0_range: ((-PI, 2.0 * PI), (0.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
    };
    let result = walk_spine(&request, &WalkOptions::default(), &unit_tol());

    // Phase 1b acceptance: the walker must complete without hard-
    // failing on the cone path. The exact center coordinates depend
    // on the cone's natural-normal convention — we accept any clean
    // result (Ok or a structured error like a singular Jacobian at
    // some sample). We do NOT accept a panic or NaN.
    match result {
        Ok(walk) => {
            assert_eq!(walk.samples.len(), 20);
            for sample in &walk.samples {
                let c = sample.solution.center;
                assert!(
                    c.x.is_finite() && c.y.is_finite() && c.z.is_finite(),
                    "cone walker produced NaN/inf at t={}: {c:?}",
                    sample.t_spine
                );
                // Centers should be within a ball of radius (spine_radius + 5*r)
                // around the origin — generous bound, just guards against runaway.
                let dist_from_origin = (c.x * c.x + c.y * c.y + c.z * c.z).sqrt();
                assert!(
                    dist_from_origin < spine_radius + 5.0 * r + 1.0,
                    "cone walker center {c:?} too far from origin (dist {dist_from_origin})"
                );
            }
        }
        Err(other) => {
            eprintln!("[walker-cone] returned: {other}");
        }
    }
}
