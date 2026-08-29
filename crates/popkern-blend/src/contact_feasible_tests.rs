//! Tests for `pre_check_contact` and `solve_with_max_feasible`.

#[cfg(test)]
use crate::contact::{
    BallSide, ContactError, ContactFeasibility, ContactOptions, ContactRequest, SingularityKind,
    pre_check_contact, solve_with_max_feasible,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// solve_with_max_feasible
// ---------------------------------------------------------------------------

#[test]
fn max_feasible_returns_original_radius_when_fits() {
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
        initial_uv0: (0.5, 0.0),
        uv0_range: ((0.0, 5.0), (-5.0, 5.0)),
        initial_uv1: (0.0, 0.5),
        uv1_range: ((-5.0, 5.0), (0.0, 5.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let (_solution, achieved) =
        solve_with_max_feasible(&request, &ContactOptions::default(), &unit_tol())
            .expect("should succeed at the requested radius");
    assert!((achieved - 0.5).abs() < 1e-12);
}

#[test]
fn max_feasible_finds_largest_radius_for_constrained_face() {
    // Plane z=0 has u in [0, 1.0]. The natural contact for radius r
    // is at u = r, so the largest feasible r is 1.0.
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
        radius: 5.0, // way too big
        initial_uv0: (0.5, 0.0),
        uv0_range: ((0.0, 1.0), (-5.0, 5.0)),
        initial_uv1: (0.0, 0.5),
        uv1_range: ((-5.0, 5.0), (0.0, 1.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let (_solution, achieved) =
        solve_with_max_feasible(&request, &ContactOptions::default(), &unit_tol())
            .expect("binary search should find a feasible radius");
    assert!(
        (achieved - 1.0).abs() / 1.0 < 0.05,
        "max feasible should be ~1.0, got {achieved}"
    );
}

#[test]
fn max_feasible_propagates_non_size_errors() {
    // Parallel planes have no fillet at any radius; the size search
    // shouldn't help and the error should propagate as-is.
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
    let result = solve_with_max_feasible(&request, &ContactOptions::default(), &unit_tol());
    match result {
        Err(ContactError::SingularJacobian {
            kind: SingularityKind::NormalsParallel,
            ..
        }) => {}
        other => panic!("expected NormalsParallel, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// pre_check_contact
// ---------------------------------------------------------------------------

#[test]
fn pre_check_perpendicular_planes_is_fine() {
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
        initial_uv0: (0.5, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 0.5),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    match pre_check_contact(&request, &unit_tol()) {
        ContactFeasibility::LikelyFine => {}
        other => panic!("expected LikelyFine, got {other:?}"),
    }
}

#[test]
fn pre_check_parallel_planes_is_infeasible() {
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
        side_1: BallSide::Outside,
        radius: 1.0,
        initial_uv0: (0.0, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    match pre_check_contact(&request, &unit_tol()) {
        ContactFeasibility::LikelyInfeasible { .. } => {}
        other => panic!("expected LikelyInfeasible, got {other:?}"),
    }
}

#[test]
fn pre_check_negative_radius_is_infeasible() {
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
        radius: -0.5,
        initial_uv0: (0.5, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 0.5),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    match pre_check_contact(&request, &unit_tol()) {
        ContactFeasibility::LikelyInfeasible { .. } => {}
        other => panic!("expected LikelyInfeasible for r<0, got {other:?}"),
    }
}

#[test]
fn pre_check_radius_larger_than_face_extent_is_tight() {
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
    // Tiny face: 0.1 × 0.1. Radius 5 is way bigger.
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: 5.0,
        initial_uv0: (0.05, 0.05),
        uv0_range: ((0.0, 0.1), (0.0, 0.1)),
        initial_uv1: (0.05, 0.05),
        uv1_range: ((0.0, 0.1), (0.0, 0.1)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    match pre_check_contact(&request, &unit_tol()) {
        ContactFeasibility::LikelyTight { .. } => {}
        other => panic!("expected LikelyTight, got {other:?}"),
    }
}
