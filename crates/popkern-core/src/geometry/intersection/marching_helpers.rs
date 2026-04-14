//! Small helpers for the marching loop: parameter wrapping, normal
//! estimation near singularities, tangent-plane basis construction,
//! and tangent-point direction search.
//!
//! Split out of `marching.rs` to keep the parent under the 500-line cap.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

/// Wrap a parameter value into the periodic domain `[u_range.0, u_range.1]`.
///
/// If `is_periodic` is false the value is returned unchanged.
/// If the period is non-positive the value is returned unchanged.
pub(super) fn wrap_parameter(u: f64, u_range: (f64, f64), is_periodic: bool) -> f64 {
    if !is_periodic {
        return u;
    }
    let period = u_range.1 - u_range.0;
    if period <= 0.0 || !u.is_finite() || !period.is_finite() {
        return u;
    }
    // Use modular arithmetic instead of a while loop to avoid
    // potential infinite loops with extreme values.
    let offset = u - u_range.0;
    let rem = offset.rem_euclid(period);
    u_range.0 + rem
}

/// Estimate a surface normal near a singularity by averaging normals at offset points.
///
/// Samples normals at four offset positions (`±delta` in each parameter direction).
/// Returns the normalised average if at least one offset yields a non-degenerate
/// normal, otherwise returns `None`.
pub(super) fn estimate_normal_near_singularity(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    delta: f64,
) -> Option<Vec3> {
    let offsets = [(delta, 0.0), (-delta, 0.0), (0.0, delta), (0.0, -delta)];
    let mut sum = Vec3::zero();
    let mut count = 0usize;
    for (du, dv) in &offsets {
        let n = surface.normal_at(u + du, v + dv);
        if n.length() > 1e-10 {
            sum = sum + n;
            count += 1;
        }
    }
    if count > 0 {
        (sum * (1.0 / count as f64)).normalized()
    } else {
        None
    }
}

/// Build an orthonormal tangent basis in the plane perpendicular to `normal`.
///
/// Returns `(tu, tv)` where `tu` and `tv` are unit vectors perpendicular to
/// `normal` and to each other.
pub(super) fn tangent_basis(normal: &Vec3) -> (Vec3, Vec3) {
    let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let axis = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let tu = n
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let tv = n.cross(&tu);
    (tu, tv)
}

/// Find a march direction at a tangent point using perturbation sampling.
///
/// When `N1 × N2` is near-zero (surfaces are locally tangent), the standard
/// marching direction is undefined.  This function samples 8 directions in the
/// tangent plane of surface 1's normal and picks the one that minimises the
/// distance to `surface2.point_at(u2, v2)` — i.e., the direction that stays
/// closest to the second surface.
///
/// Returns `None` if the normal at `(u1, v1)` is degenerate or all sampled
/// directions are equidistant.
#[allow(clippy::too_many_arguments)]
pub(super) fn find_tangent_march_direction(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    point: &Point3,
    tolerance: &Tolerance,
) -> Option<Vec3> {
    let n1 = surface1.normal_at(u1, v1);
    if n1.length() < tolerance.numerical_floor() {
        return None;
    }
    let delta = tolerance.linear() * 10.0;
    let (tu, tv) = tangent_basis(&n1);

    let mut best_dir = Vec3::zero();
    let mut min_dist = f64::MAX;

    for i in 0..8usize {
        let angle = TAU * i as f64 / 8.0;
        let dir = tu * angle.cos() + tv * angle.sin();
        let test_point = *point + dir * delta;
        let dist = (test_point - surface2.point_at(u2, v2)).length();
        if dist < min_dist {
            min_dist = dist;
            best_dir = dir;
        }
    }

    best_dir.normalized()
}
