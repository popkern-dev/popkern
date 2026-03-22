//! Differential marching for surface-surface intersection curves.

#[path = "marching_linalg.rs"]
mod marching_linalg;
use marching_linalg::solve_4x4;

use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// A sample point on a surface-surface intersection curve.
#[derive(Debug, Clone)]
pub struct MarchSample {
    /// Parameter on surface 1, u direction.
    pub u1: f64,
    /// Parameter on surface 1, v direction.
    pub v1: f64,
    /// Parameter on surface 2, u direction.
    pub u2: f64,
    /// Parameter on surface 2, v direction.
    pub v2: f64,
    /// The 3D point (equal to surface1.point_at(u1, v1) at convergence).
    pub point: Point3,
}

/// Project onto both surfaces using Newton 4×4.
///
/// Solves the system:
/// - `S1(u1,v1) - S2(u2,v2) = 0`  (3 scalar equations)
/// - `(S1(u1,v1) - target) · march_dir = 0`  (advancing front constraint)
///
/// Returns `None` if the solver fails to converge or the Jacobian is singular.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march_project(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    for _ in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);

        // Rows 0-2 residual: S1 - S2
        let f = p1 - p2;
        // Row 3 residual: (S1 - target) · march_dir
        let f4 = (p1 - *target).dot(march_dir);

        if f.length() < tolerance.linear && f4.abs() < tolerance.linear {
            return Some(MarchSample {
                u1,
                v1,
                u2,
                v2,
                point: p1,
            });
        }

        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);

        // Build 4×4 Jacobian J such that J * [du1, dv1, du2, dv2]^T = [f.x, f.y, f.z, f4]:
        //
        // Row 0: [S1_u.x,         S1_v.x,         -S2_u.x,  -S2_v.x ]
        // Row 1: [S1_u.y,         S1_v.y,         -S2_u.y,  -S2_v.y ]
        // Row 2: [S1_u.z,         S1_v.z,         -S2_u.z,  -S2_v.z ]
        // Row 3: [T·S1_u,         T·S1_v,          0,        0       ]
        //
        // Derivation of row 3:
        //   F4 = (S1(u1,v1) - target) · T
        //   ∂F4/∂u1 = S1_u · T,  ∂F4/∂v1 = S1_v · T,  ∂F4/∂u2 = 0,  ∂F4/∂v2 = 0
        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            [march_dir.dot(&s1u), march_dir.dot(&s1v), 0.0, 0.0],
        ];
        let rhs = [f.x, f.y, f.z, f4];

        let delta = solve_4x4(j, rhs)?;

        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }

    None
}

/// Like `march_project`, but also returns the number of Newton iterations used.
///
/// Returns `None` if the solver fails to converge.
#[allow(clippy::too_many_arguments)]
fn march_project_with_iters(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<(MarchSample, usize)> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    for iter in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);

        let f = p1 - p2;
        let f4 = (p1 - *target).dot(march_dir);

        if f.length() < tolerance.linear && f4.abs() < tolerance.linear {
            return Some((
                MarchSample {
                    u1,
                    v1,
                    u2,
                    v2,
                    point: p1,
                },
                iter + 1,
            ));
        }

        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);

        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            [march_dir.dot(&s1u), march_dir.dot(&s1v), 0.0, 0.0],
        ];
        let rhs = [f.x, f.y, f.z, f4];

        let delta = solve_4x4(j, rhs)?;

        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }

    None
}

/// March along a surface-surface intersection curve from a starting point.
///
/// Marches bidirectionally (+T and -T), merges results.
/// Returns `(samples, is_closed)`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    start: &MarchSample,
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> (Vec<MarchSample>, bool) {
    // March in +T direction.
    let forward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, 1.0, tolerance, handler,
    );

    // March in -T direction.
    let mut backward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, -1.0, tolerance, handler,
    );

    // Reverse backward so it goes from far-negative to start.
    backward.reverse();

    // Merge: backward (without duplicate start) + forward.
    // backward ends at start, forward starts at start, so drop last of backward.
    if !backward.is_empty() {
        backward.pop(); // remove duplicate start point
    }
    backward.extend(forward);
    let samples = backward;

    // Check if closed: last sample within tolerance of first, and enough steps.
    let is_closed = if samples.len() > 3 {
        let first = &samples[0];
        let last = &samples[samples.len() - 1];
        first.point.distance_to(&last.point) < tolerance.linear
    } else {
        false
    };

    (samples, is_closed)
}

/// Wrap a parameter value into the periodic domain `[u_range.0, u_range.1]`.
///
/// If `is_periodic` is false the value is returned unchanged.
/// If the period is non-positive the value is returned unchanged.
fn wrap_parameter(u: f64, u_range: (f64, f64), is_periodic: bool) -> f64 {
    if !is_periodic {
        return u;
    }
    let period = u_range.1 - u_range.0;
    if period <= 0.0 {
        return u;
    }
    let mut wrapped = u;
    while wrapped > u_range.1 {
        wrapped -= period;
    }
    while wrapped < u_range.0 {
        wrapped += period;
    }
    wrapped
}

/// March in one direction along the intersection curve.
#[allow(clippy::too_many_arguments)]
fn march_one_direction(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    start: &MarchSample,
    direction_sign: f64,
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> Vec<MarchSample> {
    let max_steps = 10000;
    let epsilon = 1e-12;
    let min_step = tolerance.linear;
    let max_step = tolerance.linear * 1000.0;
    let mut step_size = tolerance.linear * 50.0;

    let mut samples = vec![start.clone()];

    for _ in 0..max_steps {
        let current = &samples[samples.len() - 1];

        // Compute surface normals at current parameters.
        let n1 = surface1.normal_at(current.u1, current.v1);
        let n2 = surface2.normal_at(current.u2, current.v2);

        // Check singularity: degenerate normal.
        if n1.length() < epsilon && handler.on_singularity(&current.point) {
            break;
        }
        if n2.length() < epsilon && handler.on_singularity(&current.point) {
            break;
        }

        // Tangent direction = N1 x N2.
        let t_raw = n1.cross(&n2);
        let t_len = t_raw.length();

        // Check tangency: surfaces are tangent when cross product is near zero.
        if t_len < epsilon {
            let angle = n1.angle_to(&n2);
            if handler.on_tangent_during_march(&current.point, angle) {
                break;
            }
        }

        // Normalize tangent direction and apply direction sign.
        let t_normalized = match t_raw.normalized() {
            Some(t) => t * direction_sign,
            None => break, // cannot normalize zero vector
        };

        // Try Newton projection, with retry on failure.
        let mut succeeded = false;
        let mut retry_step = step_size;
        for _retry in 0..3 {
            let retry_target = current.point + t_normalized * retry_step;
            match march_project_with_iters(
                surface1,
                current.u1,
                current.v1,
                surface2,
                current.u2,
                current.v2,
                &retry_target,
                &t_normalized,
                tolerance,
            ) {
                Some((mut sample, iters)) => {
                    // Adaptive step based on iteration count.
                    if iters < 3 {
                        step_size = (retry_step * 1.5).min(max_step);
                    } else if iters > 8 {
                        step_size = (retry_step * 0.5).max(min_step);
                    } else {
                        step_size = retry_step;
                    }

                    // Wrap parameters for periodic surfaces before domain check.
                    let (periodic_u1, periodic_v1) = surface1.is_periodic();
                    let (periodic_u2, periodic_v2) = surface2.is_periodic();
                    sample.u1 = wrap_parameter(sample.u1, u1_range, periodic_u1);
                    sample.v1 = wrap_parameter(sample.v1, v1_range, periodic_v1);
                    sample.u2 = wrap_parameter(sample.u2, u2_range, periodic_u2);
                    sample.v2 = wrap_parameter(sample.v2, v2_range, periodic_v2);

                    // Check domain bounds — only stop for non-periodic directions.
                    let exits_u1 =
                        !periodic_u1 && (sample.u1 < u1_range.0 || sample.u1 > u1_range.1);
                    let exits_v1 =
                        !periodic_v1 && (sample.v1 < v1_range.0 || sample.v1 > v1_range.1);
                    let exits_u2 =
                        !periodic_u2 && (sample.u2 < u2_range.0 || sample.u2 > u2_range.1);
                    let exits_v2 =
                        !periodic_v2 && (sample.v2 < v2_range.0 || sample.v2 > v2_range.1);
                    if exits_u1 || exits_v1 || exits_u2 || exits_v2 {
                        // Exited domain; stop marching in this direction.
                        return samples;
                    }

                    // Check for closure: if we've marched at least a few
                    // steps and the new point is back near the start, stop.
                    let min_steps_before_closure_check = 3;
                    if samples.len() > min_steps_before_closure_check
                        && sample.point.distance_to(&start.point) < tolerance.linear
                    {
                        samples.push(sample);
                        return samples;
                    }

                    samples.push(sample);
                    succeeded = true;
                    break;
                }
                None => {
                    // Halve step and retry.
                    retry_step *= 0.5;
                    if retry_step < min_step {
                        break;
                    }
                }
            }
        }

        if !succeeded {
            break;
        }
    }

    samples
}

#[cfg(test)]
#[path = "marching_tests.rs"]
mod tests;
