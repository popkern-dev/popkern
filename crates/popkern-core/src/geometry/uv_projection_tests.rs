use super::*;
use crate::geometry::Vec3;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[test]
fn plane_projection_round_trips() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(3.0, 4.0, 0.0)).expect("plane should project");
    assert!((u - 3.0).abs() < 1e-10);
    assert!((v - 4.0).abs() < 1e-10);
}

#[test]
fn plane_projection_round_trips_non_orthogonal_axes() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, -2.0, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(2.0, 0.0, 0.0),
        v_axis: Vec3::new(0.75, 1.5, 0.0),
    };
    let expected_u = 0.37;
    let expected_v = -0.24;
    let point = surface.point_at(expected_u, expected_v);
    let (u, v) =
        project_point_to_uv(&surface, &point).expect("non-orthogonal plane should project");
    assert!((u - expected_u).abs() < 1e-10, "u={u}");
    assert!((v - expected_v).abs() < 1e-10, "v={v}");
}

#[test]
fn cylinder_point_at_r_0_h() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    // Point on cylinder at angle 0, height 5
    let (u, v) = project_point_to_uv(&surface, &Point3::new(2.0, 0.0, 5.0))
        .expect("cylinder should project");
    assert!(u.abs() < 1e-10, "u={u}");
    assert!((v - 5.0).abs() < 1e-10, "v={v}");
}

#[test]
fn cylinder_point_at_0_r_0() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    // Point on cylinder at angle pi/2, height 0
    let (u, v) = project_point_to_uv(&surface, &Point3::new(0.0, 2.0, 0.0))
        .expect("cylinder should project");
    assert!((u - FRAC_PI_2).abs() < 1e-10, "u={u}");
    assert!(v.abs() < 1e-10, "v={v}");
}

#[test]
fn cone_apex_projects_to_v_zero() {
    let surface = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: PI / 6.0,
    };
    let (_, v) =
        project_point_to_uv(&surface, &Point3::new(0.0, 0.0, 0.0)).expect("cone should project");
    assert!(v.abs() < 1e-10, "v={v}");
}

#[test]
fn sphere_north_pole_projects_to_v_zero() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    // North pole is (0, 0, r) -> colatitude = 0.
    let (_, v) =
        project_point_to_uv(&surface, &Point3::new(0.0, 0.0, 1.0)).expect("sphere should project");
    assert!(v.abs() < 1e-10, "v={v}");
}

#[test]
fn sphere_center_projects_safely() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: 5.0,
        axis: crate::geometry::Vec3::Z,
    };
    let result = project_point_to_uv(&surface, &Point3::new(1.0, 2.0, 3.0));
    assert!(result.is_some(), "center should not panic");
}

#[test]
fn torus_projection_outer_equator() {
    let surface = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    // Point on the outer equator at angle u=0, v=0: (R+r, 0, 0) = (4, 0, 0)
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(4.0, 0.0, 0.0)).expect("torus should project");
    assert!(u.abs() < 1e-10, "u={u}");
    assert!(v.abs() < 1e-10, "v={v}");
}

#[test]
fn torus_projection_top() {
    let surface = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    // Point on top of the torus at angle u=0, v=pi/2: (3, 0, 1)
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(3.0, 0.0, 1.0)).expect("torus should project");
    assert!(u.abs() < 1e-10, "u={u}");
    assert!((v - FRAC_PI_2).abs() < 1e-10, "v={v}");
}

#[test]
fn bounded_range_cylinder_vertices() {
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    // Cylinder u range is [0, TAU] (finite), v is infinite.
    // Provide vertices at heights 0 and 5.
    let vertices = vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 5.0)];
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
    // u should remain [0, TAU] because it's already finite.
    assert!((u0).abs() < 1e-10);
    assert!((u1 - TAU).abs() < 1e-10);
    // v should be bounded around [0, 5] with padding.
    assert!(v0 < 0.0);
    assert!(v1 > 5.0);
    assert!(v0.is_finite());
    assert!(v1.is_finite());
}

#[test]
fn bounded_range_plane_uses_vertex_projection() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let vertices = vec![Point3::new(1.0, 2.0, 0.0), Point3::new(3.0, 4.0, 0.0)];
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
    assert!(u0 < 1.0 && u1 > 3.0, "u range: ({u0}, {u1})");
    assert!(v0 < 2.0 && v1 > 4.0, "v range: ({v0}, {v1})");
    assert!(u0.is_finite() && u1.is_finite());
    assert!(v0.is_finite() && v1.is_finite());
}

#[test]
fn bounded_range_sphere_returns_own_range() {
    let surface = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let vertices = vec![Point3::new(1.0, 0.0, 0.0)];
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
    // Sphere already has finite range -> returned as-is.
    assert!((u0).abs() < 1e-10);
    assert!((u1 - TAU).abs() < 1e-10);
    assert!((v0).abs() < 1e-10);
    assert!((v1 - PI).abs() < 1e-10);
}

#[test]
fn bounded_range_empty_vertices_returns_raw() {
    let surface = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let ((u0, u1), (v0, v1)) = bounded_parameter_range(&[], &surface);
    assert!(u0.is_infinite());
    assert!(u1.is_infinite());
    assert!(v0.is_infinite());
    assert!(v1.is_infinite());
}

#[test]
fn nurbs_flat_surface_projects_correctly() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    // A flat NURBS surface in the XY plane, [0,2]x[0,3]
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
    // Project the point (1, 1.5, 0) -- should be at u=1, v=1.5
    let (u, v) =
        project_point_to_uv(&surface, &Point3::new(1.0, 1.5, 0.0)).expect("NURBS should project");
    assert!((u - 1.0).abs() < 0.01, "u={u}, expected ~1.0");
    assert!((v - 1.5).abs() < 0.01, "v={v}, expected ~1.5");
}

#[test]
fn nurbs_curved_surface_projects() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    // A curved NURBS surface (quarter cylinder approximation)
    let surface = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    // Project a point near the surface -- should get Some
    let result = project_point_to_uv(&surface, &Point3::new(0.5, 0.5, 0.5));
    assert!(result.is_some(), "curved NURBS should project");
    let (u, v) = result.expect("already checked is_some");
    assert!((0.0..=1.0).contains(&u), "u={u} out of range");
    assert!((0.0..=1.0).contains(&v), "v={v} out of range");
}
