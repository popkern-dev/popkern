use super::*;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn detect_plane_from_planar_points() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = detect_analytical_form(&plane, ((-1.0, 1.0), (-1.0, 1.0)), 1e-6);
    assert!(result.is_some(), "should detect plane");
    let (kind, _, dev) = result.unwrap();
    assert_eq!(kind, AnalyticalKind::Plane);
    assert!(dev < 1e-10, "deviation should be near zero for exact plane");
}

#[test]
fn detect_cylinder_from_analytical() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    // Sample a partial cylinder arc with finite height.
    let result = detect_analytical_form(&cylinder, ((0.2, 2.5), (-1.0, 1.0)), 1e-6);
    assert!(result.is_some(), "should detect cylinder");
    let (kind, geom, dev) = result.unwrap();
    assert_eq!(kind, AnalyticalKind::Cylinder);
    assert!(dev < 1e-6, "deviation {dev} should be within tolerance");
    match geom {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!((radius - 1.0).abs() < 1e-3, "radius should be ~1.0");
        }
        _ => panic!("expected Cylinder geometry"),
    }
}

#[test]
fn detect_cylinder_from_nurbs() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 2.0, 3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.5,
    };
    let nurbs = SurfaceGeometry::Nurbs(cylinder.to_nurbs_in_range((0.2, 2.5), (-1.0, 1.0)));
    let param_range = nurbs.parameter_range();
    let result = detect_analytical_form(&nurbs, param_range, 0.1);
    assert!(result.is_some(), "should detect cylinder from NURBS");
    let (kind, geom, dev) = result.unwrap();
    assert_eq!(kind, AnalyticalKind::Cylinder);
    assert!(dev < 0.1, "deviation {dev} should be within tolerance");
    match geom {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!(
                (radius - 2.5).abs() < 0.15,
                "radius should be ~2.5, got {radius}"
            );
        }
        _ => panic!("expected Cylinder geometry"),
    }
}

#[test]
fn detect_sphere_from_analytical() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 3.0,
    };
    // Sample a patch (not full sphere to avoid degenerate poles).
    let result = detect_analytical_form(&sphere, ((0.2, 2.0), (0.5, 2.5)), 1e-6);
    assert!(result.is_some(), "should detect sphere");
    let (kind, geom, dev) = result.unwrap();
    assert_eq!(kind, AnalyticalKind::Sphere);
    assert!(dev < 1e-3, "deviation {dev} should be small");
    match geom {
        SurfaceGeometry::Sphere { radius, .. } => {
            assert!(
                (radius - 3.0).abs() < 0.01,
                "radius should be ~3.0, got {radius}"
            );
        }
        _ => panic!("expected Sphere geometry"),
    }
}

#[test]
fn detect_sphere_from_nurbs() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 0.0),
        radius: 2.0,
    };
    let nurbs = SurfaceGeometry::Nurbs(sphere.to_nurbs_in_range((0.2, 2.0), (0.5, 2.5)));
    let param_range = nurbs.parameter_range();
    // NURBS sampling introduces approximation error; use wider tolerance.
    let result = detect_analytical_form(&nurbs, param_range, 0.2);
    assert!(result.is_some(), "should detect sphere from NURBS");
    let (kind, geom, dev) = result.unwrap();
    assert_eq!(kind, AnalyticalKind::Sphere);
    assert!(dev < 0.2, "deviation {dev} should be within tolerance");
    match geom {
        SurfaceGeometry::Sphere { radius, .. } => {
            assert!(
                (radius - 2.0).abs() < 0.2,
                "radius should be ~2.0, got {radius}"
            );
        }
        _ => panic!("expected Sphere geometry"),
    }
}

#[test]
fn non_analytical_surface_returns_none() {
    // A saddle surface: z = x*y, which is not plane, cylinder, or sphere.
    use crate::geometry::nurbs_surface::NurbsSurface;
    let n = 5;
    let mut control_points = Vec::with_capacity(n * n);
    for i in 0..n {
        let x = -1.0 + 2.0 * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let y = -1.0 + 2.0 * (j as f64) / ((n - 1) as f64);
            control_points.push(Point3::new(x, y, x * y));
        }
    }
    let nurbs = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![-1.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.0],
        knots_v: vec![-1.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.0],
        control_points,
        num_u: n,
        num_v: n,
        weights: vec![1.0; n * n],
    });
    let param_range = nurbs.parameter_range();
    let result = detect_analytical_form(&nurbs, param_range, 1e-6);
    assert!(
        result.is_none(),
        "saddle surface should not match any analytical form"
    );
}

#[test]
fn simplify_face_on_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    // Cube faces are already Plane geometry, so simplification should be a no-op.
    let result = simplify_face(&mut brep, face, 1e-6);
    assert!(
        result.is_none(),
        "already-analytical face should not be simplified"
    );
}

#[test]
fn simplify_all_faces_on_cube_returns_empty() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let results = simplify_all_faces(&mut brep, 1e-6);
    assert!(results.is_empty(), "all cube faces are already planar");
}

#[test]
fn analytical_kind_values() {
    assert_ne!(AnalyticalKind::Plane, AnalyticalKind::Cylinder);
    assert_ne!(AnalyticalKind::Sphere, AnalyticalKind::Torus);
}
