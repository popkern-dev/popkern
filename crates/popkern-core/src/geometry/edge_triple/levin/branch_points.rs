//! Stage 2b.50: branch-point geometry helpers for multi-branch
//! enumeration on Levin quadric pairs.
//!
//! When the rank-2 fibre of the pencil decomposes into two real
//! planes, the two intersection branches each lie in one plane.
//! Both branches pass through the **branch points** — the points
//! where the kernel line of the rank-2 fibre meets the quadric
//! surfaces. These branch points serve as valid anchors for either
//! branch's `EdgeTriple` construction (they lie on both planes by
//! definition, and on the surface by construction).
//!
//! For Steinmetz (cyl_z radius 1, cyl_x radius 1):
//! - Plane normals: `(1,0,1,0)` and `(1,0,-1,0)`.
//! - Kernel line: the y-axis.
//! - Y-axis ∩ cyl_z = `(0, ±1, 0)` — the two branch points.
//!
//! These two points are the natural endpoints of either Steinmetz
//! ellipse arc.

#![allow(dead_code)]

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Solve `a·t² + b·t + c = 0` over the reals. Returns the roots in
/// arbitrary order, or empty if there are no real roots. Handles the
/// degenerate `a = 0` (linear) case.
pub(crate) fn solve_quadratic_real_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    const TINY: f64 = 1e-15;
    if a.abs() < TINY {
        if b.abs() < TINY {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    let sqrt_d = disc.sqrt();
    let two_a = 2.0 * a;
    vec![(-b + sqrt_d) / two_a, (-b - sqrt_d) / two_a]
}

/// Line of intersection of two planes (each in 4-vector form
/// `n_x·x + n_y·y + n_z·z + d = 0`). Returns `(origin, dir)` with
/// `dir` unit-length. Returns `None` for parallel planes.
pub(crate) fn line_of_two_planes(p1: &[f64; 4], p2: &[f64; 4]) -> Option<(Point3, Vec3)> {
    let n1 = Vec3::new(p1[0], p1[1], p1[2]);
    let n2 = Vec3::new(p2[0], p2[1], p2[2]);
    let dir = n1.cross(&n2);
    let len = dir.length();
    if len < 1e-12 {
        return None;
    }
    let dir = dir * (1.0 / len);

    // Find a point on both planes by zeroing the coordinate axis most
    // parallel to `dir` and solving the resulting 2×2 linear system.
    let abs = [dir.x.abs(), dir.y.abs(), dir.z.abs()];
    let zero_idx = if abs[0] >= abs[1] && abs[0] >= abs[2] {
        0
    } else if abs[1] >= abs[2] {
        1
    } else {
        2
    };
    let (i, j) = match zero_idx {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let n1_arr = [p1[0], p1[1], p1[2]];
    let n2_arr = [p2[0], p2[1], p2[2]];
    let det = n1_arr[i] * n2_arr[j] - n1_arr[j] * n2_arr[i];
    if det.abs() < 1e-12 {
        return None;
    }
    let r1 = -p1[3];
    let r2 = -p2[3];
    let xi = (r1 * n2_arr[j] - r2 * n1_arr[j]) / det;
    let xj = (n1_arr[i] * r2 - n2_arr[i] * r1) / det;
    let mut origin = [0.0; 3];
    origin[i] = xi;
    origin[j] = xj;
    // origin[zero_idx] stays 0
    Some((Point3::new(origin[0], origin[1], origin[2]), dir))
}

/// Find where the parametric line `p(t) = origin + t·dir` meets a
/// quadric surface. Returns the real `t` roots (0, 1, or 2).
/// Supported quadrics: `Cylinder`, `Cone`, `Sphere`. Returns empty
/// vec for unsupported variants (Plane, Torus, NURBS).
pub(crate) fn line_meets_quadric_surface(
    origin: Point3,
    dir: Vec3,
    surf: &SurfaceGeometry,
) -> Vec<f64> {
    match surf {
        SurfaceGeometry::Cylinder {
            origin: o_c,
            axis,
            radius,
        } => {
            // |(p − o_c)_⊥|² = r² where (·)_⊥ removes the axial component.
            let q = origin - *o_c;
            let q_dot_a = q.dot(axis);
            let q_perp = q - (*axis * q_dot_a);
            let dir_dot_a = dir.dot(axis);
            let dir_perp = dir - (*axis * dir_dot_a);
            // |q_perp + t·dir_perp|² = r²
            let aa = dir_perp.dot(&dir_perp);
            let bb = 2.0 * dir_perp.dot(&q_perp);
            let cc = q_perp.dot(&q_perp) - radius * radius;
            solve_quadratic_real_roots(aa, bb, cc)
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            // ((p − apex)·a)² = cos²θ · |p − apex|²
            let cos2 = half_angle.pcos().powi(2);
            let q = origin - *apex;
            let q_a = q.dot(axis);
            let d_a = dir.dot(axis);
            let q_q = q.dot(&q);
            let d_d = dir.dot(&dir);
            let q_d = q.dot(&dir);
            let aa = d_a * d_a - cos2 * d_d;
            let bb = 2.0 * (q_a * d_a - cos2 * q_d);
            let cc = q_a * q_a - cos2 * q_q;
            solve_quadratic_real_roots(aa, bb, cc)
        }
        SurfaceGeometry::Sphere { center, radius } => {
            // |p − c|² = r²
            let q = origin - *center;
            let aa = dir.dot(&dir);
            let bb = 2.0 * q.dot(&dir);
            let cc = q.dot(&q) - radius * radius;
            solve_quadratic_real_roots(aa, bb, cc)
        }
        _ => Vec::new(),
    }
}

/// Synthesise `(start, end)` anchor points on the secondary plane
/// for Levin multi-branch enumeration. Returns the two branch points
/// — the intersections of the kernel line (= line common to both
/// fibre planes) with `surf`. These are guaranteed to lie on both
/// planes (by definition of the kernel) and on `surf` (by
/// construction), so they're valid anchors for the curve in either
/// plane.
///
/// Returns `None` if:
/// - the two planes are parallel (no kernel line),
/// - the kernel line meets `surf` in fewer than 2 points,
/// - or `surf` is a non-quadric variant.
pub(crate) fn branch_points_on_kernel_line(
    plane_a: &[f64; 4],
    plane_b: &[f64; 4],
    surf: &SurfaceGeometry,
) -> Option<(Point3, Point3)> {
    let (origin, dir) = line_of_two_planes(plane_a, plane_b)?;
    let ts = line_meets_quadric_surface(origin, dir, surf);
    if ts.len() < 2 {
        return None;
    }
    let p1 = origin + (dir * ts[0]);
    let p2 = origin + (dir * ts[1]);
    Some((p1, p2))
}
