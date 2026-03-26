//! UV parameter-space projection utilities.
//!
//! Projects 3D points onto surface UV spaces and computes bounded parameter
//! ranges for faces with infinite parameter domains.

use std::f64::consts::TAU;

use crate::geometry::Point3;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Project a 3D point onto a surface's UV parameter space.
///
/// Returns `Some((u, v))` for analytic surface types where an inverse
/// parametrisation is available, or `None` for NURBS and other types
/// that require iterative solving.
pub fn project_point_to_uv(surface: &SurfaceGeometry, point: &Point3) -> Option<(f64, f64)> {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => {
            let d = *point - *origin;
            let u_len_sq = u_axis.dot(u_axis);
            let v_len_sq = v_axis.dot(v_axis);
            let u = if u_len_sq > 1e-30 {
                d.dot(u_axis) / u_len_sq
            } else {
                d.dot(u_axis)
            };
            let v = if v_len_sq > 1e-30 {
                d.dot(v_axis) / v_len_sq
            } else {
                d.dot(v_axis)
            };
            Some((u, v))
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *origin;
            let mut u = f64::atan2(d.dot(&perp2), d.dot(&perp1));
            if u < 0.0 {
                u += TAU;
            }
            let v = d.dot(axis);
            Some((u, v))
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *apex;
            let mut u = f64::atan2(d.dot(&perp2), d.dot(&perp1));
            if u < 0.0 {
                u += TAU;
            }
            let v = d.dot(axis);
            Some((u, v))
        }
        SurfaceGeometry::Sphere { center, radius } => {
            let d = *point - *center;
            let r = d.length();
            // Avoid division by zero at the centre.
            if r < 1e-15 {
                return Some((0.0, 0.0));
            }
            let dz = d.z / r;
            // Clamp to avoid NaN from acos outside [-1, 1].
            let v = dz.clamp(-1.0, 1.0).acos();
            let mut u = f64::atan2(d.y, d.x);
            if u < 0.0 {
                u += TAU;
            }
            let _ = radius; // used only for the forward map
            Some((u, v))
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius: _,
        } => {
            let (perp1, perp2) = orthonormal_basis(*axis);
            let d = *point - *center;
            // Project onto the equatorial plane to find the major angle u.
            let radial_1 = d.dot(&perp1);
            let radial_2 = d.dot(&perp2);
            let mut u = f64::atan2(radial_2, radial_1);
            if u < 0.0 {
                u += TAU;
            }
            // Find the minor angle v: direction from the major-circle point
            // to the actual point, decomposed into (radial, axial).
            let h = d.dot(axis);
            let radial_dist = (radial_1 * radial_1 + radial_2 * radial_2).sqrt();
            let mut v = f64::atan2(h, radial_dist - major_radius);
            if v < 0.0 {
                v += TAU;
            }
            Some((u, v))
        }
        SurfaceGeometry::Nurbs(_) => None,
    }
}

/// Compute bounded parameter ranges for a face's surface.
///
/// If the surface's native parameter range is already finite in both u and v,
/// returns it unchanged. Otherwise, projects the given vertex positions into UV
/// space, computes the axis-aligned bounding box, and adds 50 % padding
/// (minimum 1.0 in each direction) for intersection-solver headroom.
pub fn bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // If both ranges are finite, return them as-is.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        return ((u_min, u_max), (v_min, v_max));
    }

    if vertices.is_empty() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Project each vertex to UV and accumulate the bounding box.
    let mut uv_lo = (f64::INFINITY, f64::INFINITY);
    let mut uv_hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any_projected = false;

    for pos in vertices {
        if let Some((u, v)) = project_point_to_uv(surface, pos) {
            uv_lo.0 = uv_lo.0.min(u);
            uv_lo.1 = uv_lo.1.min(v);
            uv_hi.0 = uv_hi.0.max(u);
            uv_hi.1 = uv_hi.1.max(v);
            any_projected = true;
        }
    }

    if !any_projected {
        return ((u_min, u_max), (v_min, v_max));
    }

    // 50 % padding (minimum 1.0) so downstream intersection solvers have
    // enough room; precise 3-D trimming happens later.
    let u_pad = ((uv_hi.0 - uv_lo.0) * 0.5).max(1.0);
    let v_pad = ((uv_hi.1 - uv_lo.1) * 0.5).max(1.0);

    // Only override dimensions that were infinite.
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        (u_min, u_max)
    } else {
        (uv_lo.0 - u_pad, uv_hi.0 + u_pad)
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        (v_min, v_max)
    } else {
        (uv_lo.1 - v_pad, uv_hi.1 + v_pad)
    };

    (u_range, v_range)
}

/// Compute tight bounded parameter ranges for a face's surface.
///
/// Same projection algorithm as [`bounded_parameter_range`] but with
/// 5 % padding (minimum 0.05) instead of 50 %. Suitable for
/// tessellation and AABB sampling where a small margin is needed.
pub fn tight_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    vertex_bounded_parameter_range(vertices, surface, 0.05, 0.05)
}

/// Compute exact bounded parameter ranges from face vertices with no
/// padding.
///
/// Returns the exact UV bounding box of the given vertex projections.
/// Suitable for self-intersection validation where any padding would
/// cause false positives at shared edges.
pub fn exact_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    vertex_bounded_parameter_range(vertices, surface, 0.0, 0.0)
}

/// Internal: compute bounded parameter ranges with configurable padding.
///
/// Always projects face vertices into UV space when possible. For
/// dimensions where projection succeeds, uses the tighter of the
/// vertex-projected range (with padding) and the native surface range.
/// This ensures that faces using a subset of a torus or sphere are
/// tessellated only over the portion they occupy.
fn vertex_bounded_parameter_range(
    vertices: &[Point3],
    surface: &SurfaceGeometry,
    padding_fraction: f64,
    padding_minimum: f64,
) -> ((f64, f64), (f64, f64)) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    if vertices.is_empty() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Project each vertex to UV and accumulate the bounding box.
    let mut uv_lo = (f64::INFINITY, f64::INFINITY);
    let mut uv_hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any_projected = false;

    for pos in vertices {
        if let Some((u, v)) = project_point_to_uv(surface, pos) {
            uv_lo.0 = uv_lo.0.min(u);
            uv_lo.1 = uv_lo.1.min(v);
            uv_hi.0 = uv_hi.0.max(u);
            uv_hi.1 = uv_hi.1.max(v);
            any_projected = true;
        }
    }

    if !any_projected {
        // No vertices could be projected (e.g., NURBS surface).
        // Fall back to the native range.
        return ((u_min, u_max), (v_min, v_max));
    }

    let u_pad = ((uv_hi.0 - uv_lo.0) * padding_fraction).max(padding_minimum);
    let v_pad = ((uv_hi.1 - uv_lo.1) * padding_fraction).max(padding_minimum);

    let projected_u = (uv_lo.0 - u_pad, uv_hi.0 + u_pad);
    let projected_v = (uv_lo.1 - v_pad, uv_hi.1 + v_pad);

    // For each dimension: use the tighter of (native range, projected range).
    // If native range is infinite, use projected. If native is finite but
    // wider than projected, use projected (clamped to native bounds).
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        (projected_u.0.max(u_min), projected_u.1.min(u_max))
    } else {
        projected_u
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        (projected_v.0.max(v_min), projected_v.1.min(v_max))
    } else {
        projected_v
    };

    (u_range, v_range)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Vec3;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn plane_projection_round_trips() {
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let (u, v) = project_point_to_uv(&surface, &Point3::new(3.0, 4.0, 0.0))
            .expect("plane should project");
        assert!((u - 3.0).abs() < 1e-10);
        assert!((v - 4.0).abs() < 1e-10);
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
        // Point on cylinder at angle π/2, height 0
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
        let (_, v) = project_point_to_uv(&surface, &Point3::new(0.0, 0.0, 0.0))
            .expect("cone should project");
        assert!(v.abs() < 1e-10, "v={v}");
    }

    #[test]
    fn sphere_north_pole_projects_to_v_zero() {
        let surface = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        // North pole is (0, 0, r) → colatitude = 0.
        let (_, v) = project_point_to_uv(&surface, &Point3::new(0.0, 0.0, 1.0))
            .expect("sphere should project");
        assert!(v.abs() < 1e-10, "v={v}");
    }

    #[test]
    fn sphere_center_projects_safely() {
        let surface = SurfaceGeometry::Sphere {
            center: Point3::new(1.0, 2.0, 3.0),
            radius: 5.0,
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
        let (u, v) = project_point_to_uv(&surface, &Point3::new(4.0, 0.0, 0.0))
            .expect("torus should project");
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
        // Point on top of the torus at angle u=0, v=π/2: (3, 0, 1)
        let (u, v) = project_point_to_uv(&surface, &Point3::new(3.0, 0.0, 1.0))
            .expect("torus should project");
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
        };
        let vertices = vec![Point3::new(1.0, 0.0, 0.0)];
        let ((u0, u1), (v0, v1)) = bounded_parameter_range(&vertices, &surface);
        // Sphere already has finite range → returned as-is.
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
}
