//! Tests for the analytic (exact) surface offset arms.

use super::*;
use popkern_core::geometry::{Point3, Vec3};

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

#[test]
fn plane_outward() {
    let result = offset_analytic(&xy_plane(), 1.0, false).unwrap();
    match result {
        SurfaceGeometry::Plane { origin, .. } => {
            assert!((origin.z - 1.0).abs() < 1e-12);
        }
        _ => panic!("expected Plane"),
    }
}

#[test]
fn plane_inward() {
    let result = offset_analytic(&xy_plane(), 1.0, true).unwrap();
    match result {
        SurfaceGeometry::Plane { origin, .. } => {
            assert!((origin.z - (-1.0)).abs() < 1e-12);
        }
        _ => panic!("expected Plane"),
    }
}

#[test]
fn cylinder_outward() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let result = offset_analytic(&cyl, 2.0, false).unwrap();
    match result {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!((radius - 7.0).abs() < 1e-12);
        }
        _ => panic!("expected Cylinder"),
    }
}

#[test]
fn cylinder_inward() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let result = offset_analytic(&cyl, 2.0, true).unwrap();
    match result {
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!((radius - 3.0).abs() < 1e-12);
        }
        _ => panic!("expected Cylinder"),
    }
}

#[test]
fn cylinder_radius_too_small() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let result = offset_analytic(&cyl, 2.0, true);
    assert!(result.is_err());
}

#[test]
fn sphere_offset_follows_stored_normal() {
    // The standalone Sphere parametrization's normal points TOWARD
    // the center (measured via radial_normal_sign), so under the
    // along--n contract the flag meanings are: inward=false moves
    // along +n (radius shrinks), inward=true along -n (radius
    // grows). The old per-arm convention ("outward increases radii")
    // ignored the stored normal and disagreed with the plane arm.
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 3.0,
        axis: popkern_core::geometry::Vec3::Z,
    };
    assert_eq!(radial_normal_sign(&sphere), Some(-1.0));
    match offset_analytic(&sphere, 1.0, false).unwrap() {
        SurfaceGeometry::Sphere { radius, .. } => {
            assert!((radius - 2.0).abs() < 1e-12, "along +n = toward center");
        }
        _ => panic!("expected Sphere"),
    }
    match offset_analytic(&sphere, 1.0, true).unwrap() {
        SurfaceGeometry::Sphere { radius, .. } => {
            assert!((radius - 4.0).abs() < 1e-12, "along -n = away from center");
        }
        _ => panic!("expected Sphere"),
    }
}

#[test]
fn torus_inward() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 2.0,
    };
    let result = offset_analytic(&torus, 1.0, true).unwrap();
    match result {
        SurfaceGeometry::Torus { minor_radius, .. } => {
            assert!((minor_radius - 1.0).abs() < 1e-12);
        }
        _ => panic!("expected Torus"),
    }
}

#[test]
fn cone_offset_translates_apex_exactly() {
    // Offset along +n̂ (outward from the axis) moves the apex by
    // −d/sinα along the unit axis; along −n̂ by +d/sinα. Verified
    // against direct point sampling: every offset-surface point
    // must sit exactly `d` along the original's unit normal.
    use popkern_core::geometry::surface::SurfaceEval;
    let half_angle = 0.5_f64;
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle,
    };
    let d = 0.25;
    let offset = match offset_analytic(&cone, d, false) {
        Ok(SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle: ha,
        }) => {
            // axis = (0,0,-1); apex' = apex − â·d/sinα = z + d/sinα.
            assert!((ha - half_angle).abs() < 1e-15);
            assert!((axis.z + 1.0).abs() < 1e-15);
            assert!(
                (apex.z - (2.0 + d / half_angle.sin())).abs() < 1e-12,
                "apex.z {}",
                apex.z
            );
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle: ha,
            }
        }
        other => panic!("expected exact Cone offset, got {other:?}"),
    };
    // Sample check: original point + d·n̂ lies on the offset cone
    // (radial residual against the translated apex is zero).
    for &(u, v) in &[(0.3, 0.7), (2.0, 1.5), (5.0, 0.2)] {
        let p = cone.point_at(u, v);
        let n = cone.normal_at(u, v).normalized().unwrap();
        let q = p + n * d;
        if let SurfaceGeometry::Cone { apex, axis, .. } = &offset {
            let rel = q - *apex;
            let t = rel.dot(axis);
            let radial = (rel - *axis * t).length();
            assert!(
                (radial - t * half_angle.tan()).abs() < 1e-12,
                "offset sample off the translated cone at ({u},{v})"
            );
        }
    }
}

#[test]
fn torus_spindle_offset_declines() {
    // Growing the minor radius to/past the major radius produces a
    // self-intersecting spindle torus — a structured decline with
    // the safe bound, never served or approximated.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 2.0,
    };
    // Stored torus normal points away from the tube center
    // (radial_normal_sign = +1), so the growing direction is +n
    // (inward=false).
    let result = offset_analytic(&torus, 3.5, false);
    match result {
        Err(crate::offset::OffsetError::SpindleTorus {
            new_minor,
            major,
            max_safe,
        }) => {
            assert!((new_minor - 5.5).abs() < 1e-12);
            assert!((major - 5.0).abs() < 1e-12);
            assert!(max_safe > 2.9 && max_safe < 3.0, "max_safe {max_safe}");
        }
        other => panic!("expected SpindleTorus decline, got {other:?}"),
    }
}

#[test]
fn cone_corresponding_window_shifts_axial_parameter() {
    let half_angle = 0.5_f64;
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle,
    };
    let d = 0.4;
    let (u_w, v_w) = corresponding_offset_window(&cone, d, false, (0.0, 6.0), (1.0, 3.0));
    let shift = d * half_angle.cos().powi(2) / half_angle.sin();
    assert_eq!(u_w, (0.0, 6.0));
    assert!((v_w.0 - (1.0 + shift)).abs() < 1e-12);
    assert!((v_w.1 - (3.0 + shift)).abs() < 1e-12);
    // Non-cone classes keep their windows verbatim.
    let plane = xy_plane();
    assert_eq!(
        corresponding_offset_window(&plane, d, true, (0.0, 1.0), (0.0, 1.0)),
        ((0.0, 1.0), (0.0, 1.0))
    );
}
