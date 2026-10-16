//! Robustness tests for `solve_rolling_ball_contact`: preview-mode
//! iteration cap, opt-in trace, parallel-normals singularity, bounds
//! clamping (initial-guess and converged-on-boundary).

#[cfg(test)]
use crate::contact::{
    BallSide, ContactError, ContactOptions, ContactQuality, ContactRequest, ContactTrace,
    SingularityKind, solve_rolling_ball_contact,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// ContactQuality::Preview
// ---------------------------------------------------------------------------

#[test]
fn preview_quality_caps_iterations() {
    // Plane-sphere convex case — Newton takes a few iterations.
    // Preview must consume ≤3 iterations and return within ~1e-3.
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 1.0),
        radius: 2.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.5;
    let spine_x = 3.0_f64.sqrt();
    let request = ContactRequest {
        surface_0: &sphere,
        surface_1: &plane,
        side_0: BallSide::Inside,
        side_1: BallSide::Outside,
        radius: r,
        initial_uv0: (0.0, PI / 2.0),
        uv0_range: ((-PI, 2.0 * PI), (0.0, PI)),
        initial_uv1: (spine_x, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(spine_x, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let preview_opts = ContactOptions {
        quality: ContactQuality::Preview,
        ..ContactOptions::default()
    };
    let solution = solve_rolling_ball_contact(&request, &preview_opts, &unit_tol(), None)
        .expect("preview should converge to its relaxed tolerance");
    assert!(
        solution.iterations <= 3,
        "preview should use ≤3 iterations, got {}",
        solution.iterations
    );
    let expected = Point3::new(6.0_f64.sqrt(), 0.0, r);
    assert!(
        (solution.center - expected).length() < 1.0e-3,
        "preview center {:?} too far from expected {:?}",
        solution.center,
        expected
    );
}

// ---------------------------------------------------------------------------
// ContactTrace
// ---------------------------------------------------------------------------

#[test]
fn trace_records_per_iteration_state() {
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
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: 0.5,
        initial_uv0: (4.0, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 4.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let mut trace = ContactTrace::default();
    let solution = solve_rolling_ball_contact(
        &request,
        &ContactOptions::default(),
        &unit_tol(),
        Some(&mut trace),
    )
    .expect("plane-plane should converge");

    assert!(!trace.steps.is_empty());
    for (i, step) in trace.steps.iter().enumerate() {
        assert_eq!(step.iteration, i);
    }
    let first = &trace.steps[0];
    let last = &trace.steps[trace.steps.len() - 1];
    assert!(last.residual <= first.residual);
    for step in &trace.steps {
        assert!(step.condition.is_finite() && step.condition > 0.0);
    }
    if trace.steps.len() > 1 {
        assert!(trace.steps[0].step_length > 0.0);
    }
    for step in &trace.steps {
        assert_eq!(step.damping, 0.0);
    }
    assert_eq!(solution.iterations + 1, trace.steps.len());
}

#[test]
fn trace_none_has_zero_overhead_path() {
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
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: 0.5,
        initial_uv0: (2.0, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 2.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let with_none =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None)
            .unwrap();
    let mut trace = ContactTrace::default();
    let with_some = solve_rolling_ball_contact(
        &request,
        &ContactOptions::default(),
        &unit_tol(),
        Some(&mut trace),
    )
    .unwrap();
    assert_eq!(with_none.iterations, with_some.iterations);
    assert!((with_none.center - with_some.center).length() < 1e-12);
}

// ---------------------------------------------------------------------------
// Singularity detection
// ---------------------------------------------------------------------------

#[test]
fn parallel_planes_return_normals_parallel_singularity() {
    let plane_a = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_b = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let request = ContactRequest {
        surface_0: &plane_a,
        surface_1: &plane_b,
        side_0: BallSide::Outside,
        side_1: BallSide::Inside,
        radius: 1.0,
        initial_uv0: (0.0, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let result =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None);
    match result {
        Err(ContactError::SingularJacobian {
            kind: SingularityKind::NormalsParallel,
            ..
        }) => {}
        other => panic!("expected SingularJacobian/NormalsParallel, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Bounds clamping
// ---------------------------------------------------------------------------

#[test]
fn out_of_bounds_initial_guess_returns_oob_error() {
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
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: 0.5,
        initial_uv0: (10.0, 0.0),
        uv0_range: ((0.0, 5.0), (-5.0, 5.0)),
        initial_uv1: (0.0, 0.0),
        uv1_range: ((-5.0, 5.0), (-5.0, 5.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let result =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None);
    match result {
        Err(ContactError::OutOfBounds {
            surface_index: 0,
            uv: (10.0, 0.0),
            clamped_uv: (5.0, 0.0),
        }) => {}
        other => panic!("expected OutOfBounds for surface 0, got {other:?}"),
    }
}

#[test]
fn converged_solution_on_face_boundary_reports_failure() {
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
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: 2.0,
        initial_uv0: (0.5, 0.0),
        uv0_range: ((0.0, 1.0), (-5.0, 5.0)),
        initial_uv1: (0.0, 0.5),
        uv1_range: ((-5.0, 5.0), (-5.0, 5.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let result =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None);
    match result {
        Err(ContactError::OutOfBounds { .. })
        | Err(ContactError::NotConverged { .. })
        | Err(ContactError::Diverged { .. }) => {}
        other => panic!("expected OutOfBounds/NotConverged/Diverged, got {other:?}"),
    }
}
