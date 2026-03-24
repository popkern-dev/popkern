//! Analytic (exact) surface offset for primitive surface types.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Compute an exact offset of an analytic surface.
///
/// For a plane, the origin is translated by `±distance * normal`. For
/// cylinder, sphere, and torus the radius is adjusted. The `inward` flag
/// reduces radii (or translates against the normal); outward increases them.
///
/// Returns `Err` for surface types that have no closed-form offset (Cone,
/// NURBS).
pub fn offset_analytic(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
) -> Result<SurfaceGeometry, String> {
    let sign = if inward { -1.0 } else { 1.0 };

    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => {
            let offset_vec: Vec3 = *normal * (sign * distance);
            Ok(SurfaceGeometry::Plane {
                origin: *origin + offset_vec,
                normal: *normal,
                u_axis: *u_axis,
                v_axis: *v_axis,
            })
        }

        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let new_radius = radius + sign * distance;
            if new_radius <= 0.0 {
                return Err(format!(
                    "cylinder offset radius would be {new_radius:.6} (≤ 0)"
                ));
            }
            Ok(SurfaceGeometry::Cylinder {
                origin: *origin,
                axis: *axis,
                radius: new_radius,
            })
        }

        SurfaceGeometry::Sphere { center, radius } => {
            let new_radius = radius + sign * distance;
            if new_radius <= 0.0 {
                return Err(format!(
                    "sphere offset radius would be {new_radius:.6} (≤ 0)"
                ));
            }
            Ok(SurfaceGeometry::Sphere {
                center: *center,
                radius: new_radius,
            })
        }

        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let new_minor = minor_radius + sign * distance;
            if new_minor <= 0.0 {
                return Err(format!(
                    "torus offset minor radius would be {new_minor:.6} (≤ 0)"
                ));
            }
            Ok(SurfaceGeometry::Torus {
                center: *center,
                axis: *axis,
                major_radius: *major_radius,
                minor_radius: new_minor,
            })
        }

        SurfaceGeometry::Cone { .. } => Err("analytic offset not supported for cones".to_string()),

        SurfaceGeometry::Nurbs(_) => {
            Err("analytic offset not supported for NURBS surfaces".to_string())
        }
    }
}

#[cfg(test)]
mod tests {
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
    fn sphere_outward() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 3.0,
        };
        let result = offset_analytic(&sphere, 1.0, false).unwrap();
        match result {
            SurfaceGeometry::Sphere { radius, .. } => {
                assert!((radius - 4.0).abs() < 1e-12);
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
    fn cone_not_supported() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: 0.5,
        };
        assert!(offset_analytic(&cone, 1.0, false).is_err());
    }
}
