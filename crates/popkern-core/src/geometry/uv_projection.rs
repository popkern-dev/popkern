//! UV parameter-space projection utilities.
//!
//! Projects 3D points onto surface UV spaces and computes bounded parameter
//! ranges for faces with infinite parameter domains.

use std::f64::consts::TAU;

use crate::geometry::Point3;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Project a 3D point onto a surface's UV parameter space.
///
/// Returns `Some((u, v))` for analytic surface types (closed-form inverse)
/// and for NURBS surfaces (iterative Gauss-Newton projection).
/// Returns `None` only when no starting point converges close enough.
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
        SurfaceGeometry::Nurbs(nurbs) => nurbs_gauss_newton_project(surface, nurbs, point),
    }
}

/// Multi-start Gauss-Newton projection of a 3D point onto a NURBS surface.
///
/// Tries five starting points (parameter domain centre plus four corners) and
/// returns the `(u, v)` with the smallest final Euclidean distance, or `None`
/// if no starting point converges to within 1.0 units.
fn nurbs_gauss_newton_project(
    surface: &SurfaceGeometry,
    nurbs: &NurbsSurface,
    point: &Point3,
) -> Option<(f64, f64)> {
    const MAX_ITER: usize = 20;
    const STEP_TOL: f64 = 1e-12;
    const MAX_DIST: f64 = 1.0;

    let ((u_lo, u_hi), (v_lo, v_hi)) = nurbs.parameter_range();
    let u_mid = 0.5 * (u_lo + u_hi);
    let v_mid = 0.5 * (v_lo + v_hi);

    let starts = [
        (u_mid, v_mid),
        (u_lo, v_lo),
        (u_lo, v_hi),
        (u_hi, v_lo),
        (u_hi, v_hi),
    ];

    let mut best_uv: Option<(f64, f64)> = None;
    let mut best_dist = f64::INFINITY;

    for (u0, v0) in starts {
        let mut u = u0;
        let mut v = v0;

        for _ in 0..MAX_ITER {
            let s = surface.point_at(u, v);
            let diff = *point - s;

            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);

            let rhs_u = diff.dot(&su);
            let rhs_v = diff.dot(&sv);

            let a00 = su.dot(&su);
            let a01 = su.dot(&sv);
            let a11 = sv.dot(&sv);
            let det = a00 * a11 - a01 * a01;

            if det.abs() < f64::EPSILON {
                break;
            }

            let du = (a11 * rhs_u - a01 * rhs_v) / det;
            let dv = (a00 * rhs_v - a01 * rhs_u) / det;

            u = (u + du).clamp(u_lo, u_hi);
            v = (v + dv).clamp(v_lo, v_hi);

            if du * du + dv * dv < STEP_TOL {
                break;
            }
        }

        let dist = point.distance_to(&surface.point_at(u, v));
        if dist < best_dist {
            best_dist = dist;
            best_uv = Some((u, v));
        }
    }

    if best_dist <= MAX_DIST { best_uv } else { None }
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
#[path = "uv_projection_tests.rs"]
mod tests;
