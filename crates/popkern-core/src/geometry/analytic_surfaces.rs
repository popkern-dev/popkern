//! `SurfaceEval` implementations for analytic surface types.

use std::f64::consts::TAU;

use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::{
    SurfaceEval, SurfaceGeometry, bounding_box_by_sampling, closest_point_gauss_newton_bounded,
};
use crate::geometry::{Point3, Vec3};

/// Given a unit vector `axis`, return two unit vectors `(perp1, perp2)` such that
/// `(perp1, perp2, axis)` forms a right-handed orthonormal frame.
///
/// `perp1` is the Gram–Schmidt projection of a reference cardinal direction onto the
/// plane perpendicular to `axis`, so the frame is stable and deterministic.
/// Build an orthonormal basis perpendicular to `axis`.
///
/// Returns `(perp1, perp2)` where both are unit vectors perpendicular
/// to `axis` and to each other. The STEP exporter uses `perp1` as the
/// `ref_direction` for `AXIS2_PLACEMENT_3D` to ensure UV frame alignment.
pub fn orthonormal_basis(axis: Vec3) -> (Vec3, Vec3) {
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
            SurfaceGeometry::Nurbs(nurbs) => nurbs.point_at(u, v),
            SurfaceGeometry::Gregory(g) => g.point_at(u, v),
            SurfaceGeometry::PolyhedralSpline(p) => p.point_at(u, v),
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
            SurfaceGeometry::Nurbs(nurbs) => nurbs.partial_u(u, v),
            SurfaceGeometry::Gregory(g) => g.partial_u(u, v),
            SurfaceGeometry::PolyhedralSpline(p) => p.partial_u(u, v),
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
            SurfaceGeometry::Nurbs(nurbs) => nurbs.partial_v(u, v),
            SurfaceGeometry::Gregory(g) => g.partial_v(u, v),
            SurfaceGeometry::PolyhedralSpline(p) => p.partial_v(u, v),
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

    fn parameter_periodicity(&self) -> (Option<f64>, Option<f64>) {
        match self {
            SurfaceGeometry::Plane { .. } => (None, None),
            SurfaceGeometry::Cylinder { .. } => (Some(TAU), None),
            SurfaceGeometry::Cone { .. } => (Some(TAU), None),
            SurfaceGeometry::Sphere { .. } => (Some(TAU), None),
            SurfaceGeometry::Torus { .. } => (Some(TAU), Some(TAU)),
            SurfaceGeometry::Nurbs(_) => (None, None),
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => (None, None),
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
            SurfaceGeometry::Nurbs(nurbs) => nurbs.parameter_range(),
            // Gregory: both variants evaluate on the unit square
            // `[0, 1] × [0, 1]`. The triangular form restricts to
            // the simplex `u + v ≤ 1`; callers that need to clip
            // outside-domain samples should check the domain
            // directly. The bounding box in parameter space is
            // the full unit square for uniformity with tensor-
            // product surfaces.
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
                ((0.0, 1.0), (0.0, 1.0))
            }
        }
    }

    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { .. } => Vec3::zero(),
            SurfaceGeometry::Cylinder { axis, radius, .. } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                perp1 * (-radius * u.cos()) + perp2 * (-radius * u.sin())
            }
            SurfaceGeometry::Cone {
                axis, half_angle, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let r = v * half_angle.tan();
                perp1 * (-r * u.cos()) + perp2 * (-r * u.sin())
            }
            SurfaceGeometry::Sphere { radius, .. } => Vec3::new(
                -radius * v.sin() * u.cos(),
                -radius * v.sin() * u.sin(),
                0.0,
            ),
            SurfaceGeometry::Torus {
                axis,
                major_radius,
                minor_radius,
                ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let d2_center =
                    perp1 * (-major_radius * u.cos()) + perp2 * (-major_radius * u.sin());
                let d2_radial =
                    (perp1 * (-u.cos()) + perp2 * (-u.sin())) * (minor_radius * v.cos());
                d2_center + d2_radial
            }
            SurfaceGeometry::Nurbs(nurbs) => nurbs.second_partial_uu(u, v),
            // Gregory: finite-difference on the closed-form
            // first partial. The rational inner-CP blend has a
            // closed-form 2nd derivative (a ratio with `D^3` in
            // the denominator) but we defer the analytic form to
            // a later commit; FD here is already 1e-7-accurate
            // on unit-scale patches and matches the trait's
            // default implementation exactly.
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
                let h = 1e-7;
                let fwd = self.partial_u(u + h, v);
                let bwd = self.partial_u(u - h, v);
                (fwd - bwd) * (1.0 / (2.0 * h))
            }
        }
    }

    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { .. } => Vec3::zero(),
            SurfaceGeometry::Cylinder { .. } => Vec3::zero(),
            SurfaceGeometry::Cone {
                axis, half_angle, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let tan_ha = half_angle.tan();
                perp1 * (-tan_ha * u.sin()) + perp2 * (tan_ha * u.cos())
            }
            SurfaceGeometry::Sphere { radius, .. } => {
                Vec3::new(-radius * u.sin() * v.cos(), radius * u.cos() * v.cos(), 0.0)
            }
            SurfaceGeometry::Torus {
                axis, minor_radius, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                (perp1 * (-u.sin()) + perp2 * u.cos()) * (-minor_radius * v.sin())
            }
            SurfaceGeometry::Nurbs(nurbs) => nurbs.second_partial_uv(u, v),
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
                let h = 1e-7;
                let fwd = self.partial_u(u, v + h);
                let bwd = self.partial_u(u, v - h);
                (fwd - bwd) * (1.0 / (2.0 * h))
            }
        }
    }

    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        match self {
            SurfaceGeometry::Plane { .. } => Vec3::zero(),
            SurfaceGeometry::Cylinder { .. } => Vec3::zero(),
            SurfaceGeometry::Cone { .. } => Vec3::zero(),
            SurfaceGeometry::Sphere { radius, .. } => Vec3::new(
                -radius * u.cos() * v.sin(),
                -radius * u.sin() * v.sin(),
                -radius * v.cos(),
            ),
            SurfaceGeometry::Torus {
                axis, minor_radius, ..
            } => {
                let (perp1, perp2) = orthonormal_basis(*axis);
                let d2_radial = (perp1 * u.cos() + perp2 * u.sin()) * (-minor_radius * v.cos());
                let d2_axial = *axis * (-minor_radius * v.sin());
                d2_radial + d2_axial
            }
            SurfaceGeometry::Nurbs(nurbs) => nurbs.second_partial_vv(u, v),
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
                let h = 1e-7;
                let fwd = self.partial_v(u, v + h);
                let bwd = self.partial_v(u, v - h);
                (fwd - bwd) * (1.0 / (2.0 * h))
            }
        }
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        // Honor the domain: non-periodic axes must not return
        // extrapolated parameters the surface never takes (the
        // periodic axes keep their wrap-around freedom).
        closest_point_gauss_newton_bounded(self, point, self.is_periodic())
    }

    fn bounding_box(&self, u_range: (f64, f64), v_range: (f64, f64)) -> Aabb {
        bounding_box_by_sampling(self, u_range, v_range)
    }
}

#[cfg(test)]
#[path = "analytic_surfaces_tests.rs"]
mod tests;

/// O(1) signed distance from a point to an analytic surface's
/// implicit form: positive on the side the surface normal's
/// parameterization treats as "outside the primitive" (outside the
/// cylinder/sphere/torus tube, the normal side of a plane, outside
/// the cone). `None` for free-form surfaces, which have no implicit
/// form.
pub fn analytic_signed_distance(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*p - *origin).dot(normal)),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let d = *p - *origin;
            let along = d.dot(axis);
            let radial = d - *axis * along;
            Some(radial.length() - radius)
        }
        SurfaceGeometry::Sphere { center, radius } => Some(p.distance_to(center) - radius),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let d = *p - *apex;
            let along = d.dot(axis);
            let radial = (d - *axis * along).length();
            Some(radial * half_angle.cos() - along * half_angle.sin())
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let d = *p - *center;
            let along = d.dot(axis);
            let radial = (d - *axis * along).length();
            let in_plane = (radial - major_radius).hypot(along);
            Some(in_plane - minor_radius)
        }
        _ => None,
    }
}
