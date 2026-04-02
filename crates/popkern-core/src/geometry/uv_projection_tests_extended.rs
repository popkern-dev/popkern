use super::*;
use crate::geometry::Vec3;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

// ── Plane projection: exact UV ─────────────────────────────────────────

#[test]
fn plane_projection_off_plane_point() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let (u, v) = project_point_to_uv(&surface, &Point3::new(7.0, -3.0, 99.0))
        .expect("plane always projects");
    assert!((u - 7.0).abs() < 1e-10, "u={u}");
    assert!((v - (-3.0)).abs() < 1e-10, "v={v}");
}

#[test]
fn plane_projection_origin_is_zero() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 2.0, 3.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(1.0, 2.0, 3.0)).expect("plane always projects");
    assert!(u.abs() < 1e-10);
    assert!(v.abs() < 1e-10);
}

// ── Cylinder projection: periodic handling ─────────────────────────────

#[test]
fn cylinder_periodic_wrap_near_seam() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let (u, _v) = project_point_to_uv(&surface, &Point3::new(1.0, -0.001, 0.0))
        .expect("cylinder should project");
    assert!(u > TAU - 0.01, "u={u} should be near TAU");
    assert!(u <= TAU, "u={u}");
}

#[test]
fn cylinder_point_at_pi() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let (u, v) = project_point_to_uv(&surface, &Point3::new(-1.0, 0.0, 3.0))
        .expect("cylinder should project");
    assert!((u - PI).abs() < 1e-10, "u={u}");
    assert!((v - 3.0).abs() < 1e-10, "v={v}");
}

// ── Sphere projection: polar singularity ───────────────────────────────

#[test]
fn sphere_south_pole_projects_to_v_pi() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let (_, v) =
        project_point_to_uv(&surface, &Point3::new(0.0, 0.0, -1.0)).expect("sphere should project");
    assert!((v - PI).abs() < 1e-10, "v={v}");
}

#[test]
fn sphere_equator_x_positive() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(1.0, 0.0, 0.0)).expect("sphere should project");
    assert!(u.abs() < 1e-10 || (u - TAU).abs() < 1e-10, "u={u}");
    assert!((v - FRAC_PI_2).abs() < 1e-10, "v={v}");
}

#[test]
fn sphere_near_pole_does_not_panic() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let (_, v) = project_point_to_uv(&surface, &Point3::new(1e-14, 1e-14, 1.0))
        .expect("near-pole should project");
    assert!(v < 0.01, "v={v}");
}

// ── project_point_to_uv_near: NURBS hint-based projection ─────────────

#[test]
fn nurbs_near_matches_global_for_on_surface_point() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    use crate::geometry::surface::SurfaceEval;
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 2.0, 2.0],
        knots_v: vec![0.0, 0.0, 3.0, 3.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let pt = surface.point_at(1.0, 1.5);
    let (u, v) = project_point_to_uv_near(&surface, &pt, 0.9, 1.4).expect("near should converge");
    assert!((u - 1.0).abs() < 0.01, "u={u}");
    assert!((v - 1.5).abs() < 0.01, "v={v}");
}

#[test]
fn nurbs_near_converges_from_exact_hint() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    use crate::geometry::surface::SurfaceEval;
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 2.0, 2.0],
        knots_v: vec![0.0, 0.0, 3.0, 3.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let pt = surface.point_at(0.5, 1.0);
    let (u, v) = project_point_to_uv_near(&surface, &pt, 0.5, 1.0)
        .expect("exact hint should converge immediately");
    assert!((u - 0.5).abs() < 1e-8, "u={u}");
    assert!((v - 1.0).abs() < 1e-8, "v={v}");
}

#[test]
fn near_analytic_ignores_hint() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let (u, v) = project_point_to_uv_near(&surface, &Point3::new(0.0, 1.0, 5.0), 999.0, 999.0)
        .expect("cylinder near should project");
    assert!((u - FRAC_PI_2).abs() < 1e-10, "u={u}");
    assert!((v - 5.0).abs() < 1e-10, "v={v}");
}

#[test]
fn nurbs_near_far_point_returns_none() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    // 100 units away exceeds MAX_DIST=1.0.
    let result = project_point_to_uv_near(&surface, &Point3::new(0.5, 0.5, 100.0), 0.5, 0.5);
    assert!(result.is_none(), "Far point should return None");
}

// ── bounded_parameter_range variants ───────────────────────────────────

#[test]
fn tight_bounded_range_cylinder() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let vertices = vec![Point3::new(1.0, 0.0, 2.0), Point3::new(1.0, 0.0, 8.0)];
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&vertices, &surface);
    assert!(u0 >= 0.0 && u1 <= TAU, "u: ({u0}, {u1})");
    let span: f64 = 8.0 - 2.0;
    let expected_pad = (span * 0.05).max(0.05);
    assert!((v0 - (2.0 - expected_pad)).abs() < 0.01, "v0={v0}");
    assert!((v1 - (8.0 + expected_pad)).abs() < 0.01, "v1={v1}");
}

#[test]
fn exact_bounded_range_cylinder() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let vertices = vec![Point3::new(1.0, 0.0, 2.0), Point3::new(1.0, 0.0, 8.0)];
    let ((u0, u1), (v0, v1)) = exact_bounded_parameter_range(&vertices, &surface);
    assert!(u0 >= 0.0 && u1 <= TAU);
    assert!((v0 - 2.0).abs() < 1e-10, "v0={v0}");
    assert!((v1 - 8.0).abs() < 1e-10, "v1={v1}");
}

#[test]
fn tight_bounded_range_sphere_tighter_than_native() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let vertices = vec![
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
    ];
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&vertices, &surface);
    assert!(u0 >= 0.0, "u0={u0}");
    assert!(u1 <= TAU, "u1={u1}");
    assert!(v0 >= 0.0 && v0 < FRAC_PI_2, "v0={v0}");
    assert!(v1 > FRAC_PI_2 && v1 <= PI, "v1={v1}");
}

#[test]
fn exact_bounded_range_plane() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let vertices = vec![Point3::new(1.0, 2.0, 0.0), Point3::new(3.0, 4.0, 0.0)];
    let ((u0, u1), (v0, v1)) = exact_bounded_parameter_range(&vertices, &surface);
    assert!((u0 - 1.0).abs() < 1e-10, "u0={u0}");
    assert!((u1 - 3.0).abs() < 1e-10, "u1={u1}");
    assert!((v0 - 2.0).abs() < 1e-10, "v0={v0}");
    assert!((v1 - 4.0).abs() < 1e-10, "v1={v1}");
}

#[test]
fn bounded_range_50_percent_padding_plane() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let vertices = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 0.0)];
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
    assert!((u0 - (-5.0)).abs() < 1e-10, "u0={u0}");
    assert!((u1 - 15.0).abs() < 1e-10, "u1={u1}");
    assert!((v0 - (-5.0)).abs() < 1e-10, "v0={v0}");
    assert!((v1 - 15.0).abs() < 1e-10, "v1={v1}");
}

#[test]
fn bounded_range_minimum_padding_applies() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let vertices = vec![Point3::new(5.0, 5.0, 0.0)];
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
    assert!((u0 - 4.0).abs() < 1e-10, "u0={u0}");
    assert!((u1 - 6.0).abs() < 1e-10, "u1={u1}");
    assert!((v0 - 4.0).abs() < 1e-10, "v0={v0}");
    assert!((v1 - 6.0).abs() < 1e-10, "v1={v1}");
}

#[test]
fn nurbs_tight_range_uses_native_knots() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 2.0, 2.0],
        knots_v: vec![0.0, 0.0, 3.0, 3.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let vertices = vec![Point3::new(1.0, 1.5, 0.0)];
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&vertices, &surface);
    assert!((u0).abs() < 1e-10, "u0={u0}");
    assert!((u1 - 2.0).abs() < 1e-10, "u1={u1}");
    assert!((v0).abs() < 1e-10, "v0={v0}");
    assert!((v1 - 3.0).abs() < 1e-10, "v1={v1}");
}
