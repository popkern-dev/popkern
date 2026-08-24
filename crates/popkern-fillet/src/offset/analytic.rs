//! Analytic (exact) surface offset for primitive surface types.
//!
//! Direction contract: `inward = true` means "offset along the stored
//! geometric normal's NEGATIVE direction" (-n), uniformly across all
//! arms. Radius-form surfaces (cylinder/sphere/torus) translate that
//! into a radius delta through the measured radial sign of their stored
//! normal — a parametrization whose normal points toward the axis must
//! GROW its radius under an along--n offset.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Sign of the stored normal along the radius-growing direction: `+1.0`
/// when the surface normal points away from its axis/center (radius
/// grows along `+n`), `-1.0` when it points toward it. `None` for
/// surface classes without a radius form.
pub fn radial_normal_sign(surface: &SurfaceGeometry) -> Option<f64> {
    let radial = match surface {
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let p = surface.point_at(0.0, 0.0);
            let axis_n = axis.normalized()?;
            let rel = p - *origin;
            rel - axis_n * rel.dot(&axis_n)
        }
        SurfaceGeometry::Sphere { center, .. } => {
            let p = surface.point_at(0.0, 0.5);
            p - *center
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            ..
        } => {
            let p = surface.point_at(0.0, 0.0);
            let axis_n = axis.normalized()?;
            let rel = p - *center;
            let ring_dir = (rel - axis_n * rel.dot(&axis_n)).normalized()?;
            let ring_center = *center + ring_dir * *major_radius;
            p - ring_center
        }
        _ => return None,
    };
    let n = match surface {
        SurfaceGeometry::Sphere { .. } => surface.normal_at(0.0, 0.5),
        _ => surface.normal_at(0.0, 0.0),
    };
    let dot = n.dot(&radial);
    if dot == 0.0 {
        return None;
    }
    Some(if dot > 0.0 { 1.0 } else { -1.0 })
}

/// Outcome of attempting a closed-form analytic offset.
#[derive(Debug, Clone)]
pub enum AnalyticOffset {
    /// The surface class has a closed-form offset and it exists at the
    /// requested distance.
    Exact(SurfaceGeometry),
    /// No closed-form arm exists for this surface class (cone, NURBS,
    /// Gregory, polyhedral spline) — approximation is the honest route.
    Unsupported,
}

/// Compute an exact offset of an analytic surface, distinguishing
/// "no closed-form arm for this class" from "the offset does not exist".
///
/// For a plane, the origin is translated by `±distance * normal`. For
/// cylinder, sphere, and torus the radius is adjusted by the offset
/// distance folded with the stored normal's measured radial sign, so
/// `inward` uniformly means "offset along −n" (see the module doc).
///
/// # Errors
///
/// `Err` is a genuine geometric decline: the offset surface does not
/// exist at the requested distance (a radius driven through zero).
/// Callers must NOT fall back to approximation on `Err` — an
/// approximation of a nonexistent offset is inverted geometry.
/// Surface classes without a closed form return
/// [`AnalyticOffset::Unsupported`] instead.
pub fn offset_analytic_classified(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
) -> Result<AnalyticOffset, String> {
    let sign = if inward { -1.0 } else { 1.0 };

    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => {
            let offset_vec: Vec3 = *normal * (sign * distance);
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Plane {
                origin: *origin + offset_vec,
                normal: *normal,
                u_axis: *u_axis,
                v_axis: *v_axis,
            }))
        }

        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let radial = radial_normal_sign(surface)
                .ok_or_else(|| "cylinder with degenerate axis".to_string())?;
            let new_radius = radius + sign * radial * distance;
            if new_radius <= 0.0 {
                return Err(format!(
                    "cylinder offset radius would be {new_radius:.6} (≤ 0)"
                ));
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Cylinder {
                origin: *origin,
                axis: *axis,
                radius: new_radius,
            }))
        }

        SurfaceGeometry::Sphere { center, radius } => {
            let radial = radial_normal_sign(surface)
                .ok_or_else(|| "sphere with degenerate normal".to_string())?;
            let new_radius = radius + sign * radial * distance;
            if new_radius <= 0.0 {
                return Err(format!(
                    "sphere offset radius would be {new_radius:.6} (≤ 0)"
                ));
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Sphere {
                center: *center,
                radius: new_radius,
            }))
        }

        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let radial = radial_normal_sign(surface)
                .ok_or_else(|| "torus with degenerate axis".to_string())?;
            let new_minor = minor_radius + sign * radial * distance;
            if new_minor <= 0.0 {
                return Err(format!(
                    "torus offset minor radius would be {new_minor:.6} (≤ 0)"
                ));
            }
            Ok(AnalyticOffset::Exact(SurfaceGeometry::Torus {
                center: *center,
                axis: *axis,
                major_radius: *major_radius,
                minor_radius: new_minor,
            }))
        }

        SurfaceGeometry::Cone { .. }
        | SurfaceGeometry::Nurbs(_)
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => Ok(AnalyticOffset::Unsupported),
    }
}

/// Compute an exact offset of an analytic surface.
///
/// Compatibility wrapper over [`offset_analytic_classified`]: surface
/// classes without a closed-form arm are reported as `Err` alongside
/// genuine geometric declines. Dispatch code that must distinguish the
/// two (to avoid approximating a nonexistent offset) uses the
/// classified form.
pub fn offset_analytic(
    surface: &SurfaceGeometry,
    distance: f64,
    inward: bool,
) -> Result<SurfaceGeometry, String> {
    match offset_analytic_classified(surface, distance, inward)? {
        AnalyticOffset::Exact(s) => Ok(s),
        AnalyticOffset::Unsupported => {
            let class = match surface {
                SurfaceGeometry::Cone { .. } => "cone",
                SurfaceGeometry::Nurbs(_) => "NURBS",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::PolyhedralSpline(_) => "polyhedral spline",
                _ => "analytic",
            };
            Err(format!(
                "analytic offset not supported for {class} surfaces"
            ))
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
    fn cone_not_supported() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: 0.5,
        };
        assert!(offset_analytic(&cone, 1.0, false).is_err());
    }
}
