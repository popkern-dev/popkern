//! `SurfaceEval` implementations for analytic surface types.

use std::f64::consts::TAU;

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

/// Given a unit vector `axis`, return two unit vectors `(perp1, perp2)` such that
/// `(perp1, perp2, axis)` forms a right-handed orthonormal frame.
///
/// `perp1` is the Gram–Schmidt projection of a reference cardinal direction onto the
/// plane perpendicular to `axis`, so the frame is stable and deterministic.
fn orthonormal_basis(axis: Vec3) -> (Vec3, Vec3) {
    // Choose a reference cardinal that is least parallel to axis.
    let reference = if axis.x.abs() <= axis.y.abs() && axis.x.abs() <= axis.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if axis.y.abs() <= axis.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };

    // Gram–Schmidt: project `reference` onto the plane perpendicular to `axis`.
    let proj = reference - axis * reference.dot(&axis);
    let perp1 = proj.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    // perp2 completes the right-handed frame: perp2 = perp1 × axis... but we want
    // (perp1, perp2, axis) right-handed, so perp2 = axis × perp1 gives left-hand order.
    // Use perp2 = perp1.cross(axis) to keep right-handed (perp1 × perp2 = axis).
    // Actually for right-handed: perp1 × perp2 = axis, so perp2 = axis × perp1 negated,
    // i.e. perp2 = perp1 × axis... Let's just do perp2 = axis.cross(perp1) and negate
    // if needed. Simplest: perp2 = perp1 rotated 90° around axis = axis × perp1 (but this
    // gives perp2 such that perp2 × axis = perp1, i.e. right-handed is axis,perp1,perp2).
    // For our surfaces we just need two perpendicular vectors; use cross product directly.
    let perp2 = axis.cross(&perp1);
    (perp1, perp2)
}

impl SurfaceEval for SurfaceGeometry {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        match self {
            SurfaceGeometry::Plane {
                origin,
                u_axis,
                v_axis,
                ..
            } => *origin + *u_axis * u + *v_axis * v,
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let radial = perp1 * (radius * u.cos()) + perp2 * (radius * u.sin());
                *origin + radial + *axis * v
            }
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let r = v * half_angle.tan();
                let radial = perp1 * (r * u.cos()) + perp2 * (r * u.sin());
                *apex + *axis * v + radial
            }
            SurfaceGeometry::Sphere { center, radius } => {
                // Standard spherical: u = longitude [0, 2π), v = colatitude [0, π]
                // (0,0,1) is north pole at v=0
                let x = radius * v.sin() * u.cos();
                let y = radius * v.sin() * u.sin();
                let z = radius * v.cos();
                *center + Vec3::new(x, y, z)
            }
            SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let tube_center_offset =
                    perp1 * (major_radius * u.cos()) + perp2 * (major_radius * u.sin());
                let tube_radial = (perp1 * u.cos() + perp2 * u.sin()) * (minor_radius * v.cos());
                let tube_axial = *axis * (minor_radius * v.sin());
                *center + tube_center_offset + tube_radial + tube_axial
            }
        }
    }

    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { u_axis, .. } => *u_axis,
            SurfaceGeometry::Cylinder { axis, radius, .. } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                perp1 * (-radius * u.sin()) + perp2 * (radius * u.cos())
            }
            SurfaceGeometry::Cone {
                axis, half_angle, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let r = v * half_angle.tan();
                perp1 * (-r * u.sin()) + perp2 * (r * u.cos())
            }
            SurfaceGeometry::Sphere { radius, .. } => {
                // d/du: (-sin(u)*sin(v), cos(u)*sin(v), 0)
                Vec3::new(-radius * v.sin() * u.sin(), radius * v.sin() * u.cos(), 0.0)
            }
            SurfaceGeometry::Torus {
                axis,
                major_radius,
                minor_radius,
                ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                // d/du of tube_center_offset + tube_radial
                let d_center = perp1 * (-major_radius * u.sin()) + perp2 * (major_radius * u.cos());
                let d_radial = (perp1 * (-u.sin()) + perp2 * u.cos()) * (minor_radius * v.cos());
                d_center + d_radial
            }
        }
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { v_axis, .. } => *v_axis,
            SurfaceGeometry::Cylinder { axis, .. } => *axis,
            SurfaceGeometry::Cone {
                axis, half_angle, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let tan_ha = half_angle.tan();
                perp1 * (tan_ha * u.cos()) + perp2 * (tan_ha * u.sin()) + *axis
            }
            SurfaceGeometry::Sphere { radius, .. } => {
                // d/dv: (cos(u)*cos(v), sin(u)*cos(v), -sin(v))
                Vec3::new(
                    radius * u.cos() * v.cos(),
                    radius * u.sin() * v.cos(),
                    -radius * v.sin(),
                )
            }
            SurfaceGeometry::Torus {
                axis, minor_radius, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                // d/dv of tube_radial + tube_axial
                let d_radial = (perp1 * u.cos() + perp2 * u.sin()) * (-minor_radius * v.sin());
                let d_axial = *axis * (minor_radius * v.cos());
                d_radial + d_axial
            }
        }
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { normal, .. } => *normal,
            _ => {
                let du = self.partial_u(u, v);
                let dv = self.partial_v(u, v);
                du.cross(&dv)
            }
        }
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        match self {
            SurfaceGeometry::Plane { .. } => (
                (f64::NEG_INFINITY, f64::INFINITY),
                (f64::NEG_INFINITY, f64::INFINITY),
            ),
            SurfaceGeometry::Cylinder { .. } => ((0.0, TAU), (f64::NEG_INFINITY, f64::INFINITY)),
            SurfaceGeometry::Cone { .. } => ((0.0, TAU), (0.0, f64::INFINITY)),
            SurfaceGeometry::Sphere { .. } => ((0.0, TAU), (0.0, std::f64::consts::PI)),
            SurfaceGeometry::Torus { .. } => ((0.0, TAU), (0.0, TAU)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn plane_point_at_origin() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p = plane.point_at(0.0, 0.0);
        assert!((p.x).abs() < 1e-10);
        assert!((p.y).abs() < 1e-10);
    }

    #[test]
    fn plane_point_at_offset() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p = plane.point_at(3.0, 4.0);
        assert!((p.x - 3.0).abs() < 1e-10);
        assert!((p.y - 4.0).abs() < 1e-10);
    }

    #[test]
    fn plane_normal_is_constant() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let n = plane.normal_at(5.0, 10.0);
        assert!((n.x).abs() < 1e-10);
        assert!((n.y).abs() < 1e-10);
        assert!(n.z > 0.0);
    }

    #[test]
    fn sphere_point_at_north_pole() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        // v=0 should be north pole (0,0,radius) in standard spherical
        let p = sphere.point_at(0.0, 0.0);
        let dist = ((p.x).powi(2) + (p.y).powi(2) + (p.z).powi(2)).sqrt();
        assert!((dist - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_all_points_on_sphere() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(1.0, 2.0, 3.0),
            radius: 5.0,
        };
        for u in [0.0, 1.0, PI, 5.0] {
            for v in [0.0, 0.5, FRAC_PI_2, PI] {
                let p = sphere.point_at(u, v);
                let dist = ((p.x - 1.0).powi(2) + (p.y - 2.0).powi(2) + (p.z - 3.0).powi(2)).sqrt();
                assert!(
                    (dist - 5.0).abs() < 1e-10,
                    "Point at u={u}, v={v} not on sphere: dist={dist}"
                );
            }
        }
    }

    #[test]
    fn cylinder_point_on_surface() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let p = cyl.point_at(0.0, 5.0);
        let radial_dist = ((p.x).powi(2) + (p.y).powi(2)).sqrt();
        assert!((radial_dist - 2.0).abs() < 1e-10);
        assert!((p.z - 5.0).abs() < 1e-10);
    }

    #[test]
    fn torus_point_on_surface() {
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 3.0,
            minor_radius: 1.0,
        };
        // At u=0, v=0: should be at (major+minor, 0, 0) = (4, 0, 0)
        let p = torus.point_at(0.0, 0.0);
        assert!((p.x - 4.0).abs() < 1e-10);
        assert!((p.y).abs() < 1e-10);
        assert!((p.z).abs() < 1e-10);
    }
}
