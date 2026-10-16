//! Unit tests for the math primitives in `contact_math.rs`:
//! `solve_4x4`, `condition_estimate_4x4`, `signed_unit_normal`,
//! Weingarten map (`weingarten_dn_du`/`dv`), and the
//! `build_jacobian_and_residual` builder. Plus one cross-cutting
//! sanity test on the Jacobian's condition number.

#[cfg(test)]
use crate::contact::BallSide;
use crate::contact::math::{
    build_jacobian_and_residual, condition_estimate_4x4, signed_unit_normal, solve_4x4,
    weingarten_dn_du, weingarten_dn_dv,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use std::f64::consts::PI;

const FLOOR: f64 = 1e-12;

// ---------------------------------------------------------------------------
// solve_4x4
// ---------------------------------------------------------------------------

#[test]
fn solve_4x4_identity_returns_b() {
    let a = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let b = [3.0, 5.0, 7.0, 11.0];
    let x = solve_4x4(a, b, FLOOR).expect("identity is invertible");
    for i in 0..4 {
        assert!((x[i] - b[i]).abs() < 1e-15);
    }
}

#[test]
fn solve_4x4_known_system() {
    let a = [
        [2.0, 0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0, 0.0],
        [0.0, 0.0, 5.0, 0.0],
        [0.0, 0.0, 0.0, 7.0],
    ];
    let b = [2.0, 6.0, 15.0, 28.0];
    let x = solve_4x4(a, b, FLOOR).unwrap();
    assert!((x[0] - 1.0).abs() < 1e-12);
    assert!((x[1] - 2.0).abs() < 1e-12);
    assert!((x[2] - 3.0).abs() < 1e-12);
    assert!((x[3] - 4.0).abs() < 1e-12);
}

#[test]
fn solve_4x4_dense_system() {
    let a = [
        [4.0, 3.0, 2.0, 1.0],
        [1.0, 4.0, 3.0, 2.0],
        [2.0, 1.0, 4.0, 3.0],
        [3.0, 2.0, 1.0, 4.0],
    ];
    let x_true = [1.0, 2.0, 3.0, 4.0];
    let mut b = [0.0; 4];
    for r in 0..4 {
        for c in 0..4 {
            b[r] += a[r][c] * x_true[c];
        }
    }
    let x = solve_4x4(a, b, FLOOR).unwrap();
    for i in 0..4 {
        assert!(
            (x[i] - x_true[i]).abs() < 1e-10,
            "x[{i}] = {} expected {}",
            x[i],
            x_true[i]
        );
    }
}

#[test]
fn solve_4x4_requires_pivoting() {
    let a = [
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let b = [2.0, 3.0, 5.0, 7.0];
    let x = solve_4x4(a, b, FLOOR).unwrap();
    assert!((x[0] - 3.0).abs() < 1e-15);
    assert!((x[1] - 2.0).abs() < 1e-15);
    assert!((x[2] - 5.0).abs() < 1e-15);
    assert!((x[3] - 7.0).abs() < 1e-15);
}

#[test]
fn solve_4x4_singular_returns_none() {
    let a = [
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 3.0, 4.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let b = [1.0, 1.0, 1.0, 1.0];
    assert!(solve_4x4(a, b, FLOOR).is_none());
}

// ---------------------------------------------------------------------------
// condition_estimate_4x4
// ---------------------------------------------------------------------------

#[test]
fn condition_of_identity_is_one() {
    let i = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let cond = condition_estimate_4x4(&i);
    assert!((cond - 1.0).abs() < 1e-12);
}

#[test]
fn condition_of_singular_is_infinity() {
    let m = [
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 3.0, 4.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    assert_eq!(condition_estimate_4x4(&m), f64::INFINITY);
}

#[test]
fn condition_of_diag_2_3_5_7() {
    // diag(2, 3, 5, 7) has 1-norm 7 and inverse-1-norm 1/2 → condition 3.5.
    let d = [
        [2.0, 0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0, 0.0],
        [0.0, 0.0, 5.0, 0.0],
        [0.0, 0.0, 0.0, 7.0],
    ];
    let cond = condition_estimate_4x4(&d);
    assert!((cond - 3.5).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// signed_unit_normal
// ---------------------------------------------------------------------------

#[test]
fn signed_unit_normal_outside_matches_normal_at() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let n = signed_unit_normal(&plane, 0.0, 0.0, BallSide::Outside, FLOOR).unwrap();
    assert!((n - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-12);
}

#[test]
fn signed_unit_normal_inside_flips_sign() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let n = signed_unit_normal(&plane, 0.0, 0.0, BallSide::Inside, FLOOR).unwrap();
    assert!((n - Vec3::new(0.0, 0.0, -1.0)).length() < 1e-12);
}

// ---------------------------------------------------------------------------
// weingarten_dn_du / dn_dv
// ---------------------------------------------------------------------------

#[test]
fn weingarten_on_plane_is_zero() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let n = signed_unit_normal(&plane, 0.5, 0.7, BallSide::Outside, FLOOR).unwrap();
    let dn_du = weingarten_dn_du(&plane, 0.5, 0.7, &n, FLOOR);
    let dn_dv = weingarten_dn_dv(&plane, 0.5, 0.7, &n, FLOOR);
    assert!(dn_du.length() < 1e-12);
    assert!(dn_dv.length() < 1e-12);
}

#[test]
fn weingarten_on_unit_sphere_radial() {
    // For a unit sphere centered at origin, the Weingarten map should
    // match a centered-difference approximation of dn̂/du and dn̂/dv.
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let u = 1.0;
    let v = PI / 3.0;
    let n = signed_unit_normal(&sphere, u, v, BallSide::Outside, FLOOR).unwrap();
    let dn_du_analytic = weingarten_dn_du(&sphere, u, v, &n, FLOOR);
    let dn_dv_analytic = weingarten_dn_dv(&sphere, u, v, &n, FLOOR);

    let h = 1e-6;
    let n_up_u = signed_unit_normal(&sphere, u + h, v, BallSide::Outside, FLOOR).unwrap();
    let n_dn_u = signed_unit_normal(&sphere, u - h, v, BallSide::Outside, FLOOR).unwrap();
    let dn_du_fd = (n_up_u - n_dn_u) * (1.0 / (2.0 * h));

    let n_up_v = signed_unit_normal(&sphere, u, v + h, BallSide::Outside, FLOOR).unwrap();
    let n_dn_v = signed_unit_normal(&sphere, u, v - h, BallSide::Outside, FLOOR).unwrap();
    let dn_dv_fd = (n_up_v - n_dn_v) * (1.0 / (2.0 * h));

    let err_u = (dn_du_analytic - dn_du_fd).length();
    let err_v = (dn_dv_analytic - dn_dv_fd).length();
    assert!(
        err_u < 1e-6,
        "dn̂/du on sphere: analytic {dn_du_analytic:?}, FD {dn_du_fd:?}, err {err_u}"
    );
    assert!(
        err_v < 1e-6,
        "dn̂/dv on sphere: analytic {dn_dv_analytic:?}, FD {dn_dv_fd:?}, err {err_v}"
    );
}

#[test]
fn weingarten_on_cylinder_radial_only() {
    // Cylinder along z-axis: nonzero curvature in u (around the
    // barrel), zero in v (along the axis).
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let u = 0.5;
    let v = 1.3;
    let n = signed_unit_normal(&cyl, u, v, BallSide::Outside, FLOOR).unwrap();
    let dn_du = weingarten_dn_du(&cyl, u, v, &n, FLOOR);
    let dn_dv = weingarten_dn_dv(&cyl, u, v, &n, FLOOR);
    assert!(
        dn_dv.length() < 1e-10,
        "dn̂/dv on cylinder should be zero, got {dn_dv:?}"
    );
    assert!(
        dn_du.length() > 0.1,
        "dn̂/du on cylinder should be nonzero, got {dn_du:?}"
    );

    let h = 1e-6;
    let n_up = signed_unit_normal(&cyl, u + h, v, BallSide::Outside, FLOOR).unwrap();
    let n_dn = signed_unit_normal(&cyl, u - h, v, BallSide::Outside, FLOOR).unwrap();
    let dn_du_fd = (n_up - n_dn) * (1.0 / (2.0 * h));
    let err = (dn_du - dn_du_fd).length();
    assert!(
        err < 1e-6,
        "cylinder dn̂/du: analytic {dn_du:?}, FD {dn_du_fd:?}, err {err}"
    );
}

// ---------------------------------------------------------------------------
// build_jacobian_and_residual
// ---------------------------------------------------------------------------

#[test]
fn jacobian_residual_at_known_solution_is_zero() {
    // Two perpendicular planes meeting at the y-axis. With the correct
    // contact seed (u=1 on plane z=0, v=1 on plane x=0) the residual
    // is exactly zero.
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
    let spine_point = Point3::new(0.0, 0.0, 0.0);
    let spine_plane_normal = Vec3::new(0.0, 1.0, 0.0);

    // Wrong-seed sanity: residual must be nonzero away from the truth.
    let (_j, residual, c0, c1, _n0, _n1) = build_jacobian_and_residual(
        &plane_z0,
        &plane_x0,
        BallSide::Outside,
        BallSide::Outside,
        1.0,
        (0.0, 0.0),
        (0.0, 0.0),
        spine_point,
        spine_plane_normal,
        FLOOR,
    )
    .expect("plane-plane Jacobian should build");
    assert!(c0.distance_to(&c1) > 0.5);
    assert!(residual.iter().any(|&x| x.abs() > 0.1));

    // Correct seed: ball center should be (1, 0, 1).
    let (_j, residual, c0, c1, _n0, _n1) = build_jacobian_and_residual(
        &plane_z0,
        &plane_x0,
        BallSide::Outside,
        BallSide::Outside,
        1.0,
        (1.0, 0.0),
        (0.0, 1.0),
        spine_point,
        spine_plane_normal,
        FLOOR,
    )
    .unwrap();
    assert!(
        (c0 - Point3::new(1.0, 0.0, 1.0)).length() < 1e-12,
        "c0 should be (1,0,1), got {c0:?}"
    );
    assert!(
        (c1 - Point3::new(1.0, 0.0, 1.0)).length() < 1e-12,
        "c1 should be (1,0,1), got {c1:?}"
    );
    for (i, r) in residual.iter().enumerate().take(3) {
        assert!(r.abs() < 1e-12, "residual[{i}] = {r} should be 0");
    }
    assert!(
        residual[3].abs() < 1e-12,
        "spine residual should be 0, got {}",
        residual[3]
    );
}

#[test]
fn jacobian_4x4_invertible_for_plane_plane() {
    // At the converged plane-plane solution, the Jacobian should be
    // well-conditioned (NOT singular).
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
    let (j, _residual, _c0, _c1, _n0, _n1) = build_jacobian_and_residual(
        &plane_z0,
        &plane_x0,
        BallSide::Outside,
        BallSide::Outside,
        1.0,
        (1.0, 0.0),
        (0.0, 1.0),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        FLOOR,
    )
    .unwrap();
    let cond = condition_estimate_4x4(&j);
    assert!(
        cond.is_finite() && cond < 100.0,
        "plane-plane Jacobian condition should be small, got {cond}"
    );
}
