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
use super::{WalkOptions, WalkRequest, walk_spine};
use crate::contact::BallSide;
use crate::contact::ContactOptions;
use crate::walker::WalkMode;
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
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
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
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
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
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
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

// ---------------------------------------------------------------------------
// Step retry / interval splitting (§11.4)
// ---------------------------------------------------------------------------

#[test]
fn walker_adaptive_retries_on_step_failure_via_dt_halving() {
    // Construct a case where the *first* large dt would step past
    // the spine endpoint, generating an out-of-bounds-on-converged
    // failure, but a halved dt lands inside cleanly. The retry path
    // should succeed without any visible failure to the caller.
    //
    // We engineer this by giving the adaptive walker an `initial_dt`
    // just slightly larger than the spine span. The first attempted
    // step lands at t = t_start + initial_dt > t_end. The walker
    // snaps to t_end, but if `t_end` happens to be in a region where
    // Newton fails at the chosen spine_plane_normal (it doesn't, in
    // plane-plane), the retry kicks in.
    //
    // Plane-plane is too well-behaved to fail naturally — this test
    // mainly proves the retry path doesn't BREAK the well-behaved
    // case. The cone walker test exercises the retry on real
    // failures. So this test verifies: a giant initial_dt with snap-
    // to-end logic still produces a clean walk.
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
        t_range: (0.0, 1.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: None,
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
    };
    // initial_dt = 100 × the span. The walker should snap to t_end
    // on the first step and produce a 2-sample walk.
    let opts = WalkOptions {
        mode: WalkMode::Adaptive {
            initial_dt: 100.0,
            min_dt: 0.001,
            max_dt: 200.0,
            max_samples: 20,
            grow_threshold: 4,
            shrink_threshold: 7,
            grow_factor: 1.5,
            shrink_factor: 0.5,
        },
        contact_options: ContactOptions::default(),
    };
    let walk = walk_spine(&request, &opts, &unit_tol()).unwrap();
    assert!(
        walk.samples.len() >= 2 && walk.samples.len() <= 4,
        "snap-to-end walk should produce 2-4 samples; got {}",
        walk.samples.len()
    );
    let last = &walk.samples[walk.samples.len() - 1];
    assert!((last.t_spine - 1.0).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// Adaptive stepping (§11.3 step 4)
// ---------------------------------------------------------------------------

#[test]
fn walker_adaptive_grows_step_on_well_converging_geometry() {
    // Plane-plane case: the residual is linear, so Newton converges
    // in 1-2 iterations after the first sample. The adaptive walker
    // should grow `dt` from `initial_dt` toward `max_dt` and produce
    // FEWER samples than the equivalent fixed-step walk.
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
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
    };

    let adaptive = walk_spine(&request, &WalkOptions::adaptive(5.0), &unit_tol())
        .expect("adaptive walk should succeed");
    let fixed_20 = walk_spine(&request, &WalkOptions::default(), &unit_tol()).unwrap();

    // Adaptive should use significantly fewer samples since plane-
    // plane converges in 1 iter per step → step grows aggressively.
    assert!(
        adaptive.samples.len() < fixed_20.samples.len(),
        "adaptive ({}) should use fewer samples than fixed 20 ({})",
        adaptive.samples.len(),
        fixed_20.samples.len()
    );
    // First and last samples should still be at the spine endpoints.
    let first = &adaptive.samples[0];
    let last = &adaptive.samples[adaptive.samples.len() - 1];
    assert!((first.t_spine - 0.0).abs() < 1e-9);
    assert!((last.t_spine - 5.0).abs() < 1e-9);
    // All sample centers should be at (r, t, r).
    for sample in &adaptive.samples {
        let expected = Point3::new(0.5, sample.t_spine, 0.5);
        assert!(
            (sample.solution.center - expected).length() < 1e-9,
            "adaptive sample at t={}: center {:?} ≠ {:?}",
            sample.t_spine,
            sample.solution.center,
            expected
        );
    }
}

#[test]
fn walker_adaptive_caps_at_max_samples() {
    // Construct a case where dt cannot grow above the floor — pin
    // initial_dt at min_dt and check the cap kicks in.
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
        t_range: (0.0, 100.0),
        uv0_range: ((-200.0, 200.0), (-200.0, 200.0)),
        uv1_range: ((-200.0, 200.0), (-200.0, 200.0)),
        initial_seed: None,
        surface_0_periodicity: (None, None),
        surface_1_periodicity: (None, None),
    };
    let opts = WalkOptions {
        mode: WalkMode::Adaptive {
            initial_dt: 0.01,
            min_dt: 0.01,
            max_dt: 0.01,
            max_samples: 50, // cap is 50 but we'd need 10001 samples to walk 0..100 at dt=0.01
            grow_threshold: 4,
            shrink_threshold: 7,
            grow_factor: 1.0, // no growth
            shrink_factor: 1.0,
        },
        contact_options: ContactOptions::default(),
    };
    let walk = walk_spine(&request, &opts, &unit_tol()).unwrap();
    assert!(walk.samples.len() <= 50);
}

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
        // Cone u is angular with period TAU; v is axial (non-periodic).
        surface_0_periodicity: cone.periodicity(),
        surface_1_periodicity: plane.periodicity(),
    };
    let walk = walk_spine(&request, &WalkOptions::default(), &unit_tol())
        .expect("phase 1c periodicity-aware walker should walk the cone cleanly");

    // Phase 1c acceptance: all 20 samples must complete (was 19/20
    // before the periodicity-aware iteration window was added).
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
