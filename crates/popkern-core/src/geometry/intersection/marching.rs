//! Differential marching for surface-surface intersection curves.

#[path = "marching_linalg.rs"]
mod marching_linalg;

#[path = "marching_newton.rs"]
mod marching_newton;

pub(crate) use marching_newton::march_project;
use marching_newton::march_project_with_iters;

use crate::geometry::intersection::adaptive_step::{
    curvature_based_step, estimate_intersection_curvature,
};
use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

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

    // Check if closed: last sample within a step-sized tolerance of first.
    // The march steps in discrete increments, so exact closure is impossible.
    // Use the average step length as the closure tolerance.
    let avg_step = if samples.len() > 1 {
        let total_arc: f64 = samples
            .windows(2)
            .map(|w| w[0].point.distance_to(&w[1].point))
            .sum();
        total_arc / (samples.len() - 1) as f64
    } else {
        tolerance.linear
    };
    let closure_tol = (avg_step * 2.0).max(tolerance.linear);
    let is_closed = if samples.len() > 10 {
        let first = &samples[0];
        let last = &samples[samples.len() - 1];
        first.point.distance_to(&last.point) < closure_tol
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
fn estimate_normal_near_singularity(
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
fn tangent_basis(normal: &Vec3) -> (Vec3, Vec3) {
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
fn find_tangent_march_direction(
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
    if n1.length() < 1e-12 {
        return None;
    }
    let delta = tolerance.linear * 10.0;
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

/// March in one direction along the intersection curve.
#[allow(clippy::too_many_arguments)]
pub(crate) fn march_one_direction(
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
    let max_steps = 10_000;
    let epsilon = 1e-12;
    let min_step = tolerance.linear;
    let max_step = tolerance.linear * 1000.0;
    let mut step_size = tolerance.linear * 50.0;

    let mut samples = vec![start.clone()];

    for _ in 0..max_steps {
        let current = &samples[samples.len() - 1];

        // Compute surface normals at current parameters.
        let mut n1 = surface1.normal_at(current.u1, current.v1);
        let mut n2 = surface2.normal_at(current.u2, current.v2);

        // Check singularity: degenerate normal on surface 1.
        if n1.length() < epsilon {
            if handler.on_singularity(&current.point) {
                break;
            }
            // Handler returned false (robust): attempt normal estimation near singularity.
            let delta = tolerance.linear * 10.0;
            match estimate_normal_near_singularity(surface1, current.u1, current.v1, delta) {
                Some(estimated) => n1 = estimated,
                None => {
                    // Estimation failed: reduce step size and retry up to 3 times.
                    let mut recovered = false;
                    let mut reduced_step = step_size * 0.5;
                    for _attempt in 0..3 {
                        let test_n = estimate_normal_near_singularity(
                            surface1,
                            current.u1 + reduced_step,
                            current.v1,
                            delta,
                        );
                        if let Some(en) = test_n {
                            n1 = en;
                            step_size = reduced_step;
                            recovered = true;
                            break;
                        }
                        reduced_step *= 0.5;
                    }
                    if !recovered {
                        break;
                    }
                }
            }
        }

        // Check singularity: degenerate normal on surface 2.
        if n2.length() < epsilon {
            if handler.on_singularity(&current.point) {
                break;
            }
            // Handler returned false (robust): attempt normal estimation near singularity.
            let delta = tolerance.linear * 10.0;
            match estimate_normal_near_singularity(surface2, current.u2, current.v2, delta) {
                Some(estimated) => n2 = estimated,
                None => {
                    // Estimation failed: reduce step size and retry up to 3 times.
                    let mut recovered = false;
                    let mut reduced_step = step_size * 0.5;
                    for _attempt in 0..3 {
                        let test_n = estimate_normal_near_singularity(
                            surface2,
                            current.u2 + reduced_step,
                            current.v2,
                            delta,
                        );
                        if let Some(en) = test_n {
                            n2 = en;
                            step_size = reduced_step;
                            recovered = true;
                            break;
                        }
                        reduced_step *= 0.5;
                    }
                    if !recovered {
                        break;
                    }
                }
            }
        }

        // Tangent direction = N1 x N2.
        let t_raw = n1.cross(&n2);
        let t_len = t_raw.length();

        // Check tangency: surfaces are tangent when cross product is near zero.
        let t_normalized = if t_len < epsilon {
            let angle = n1.angle_to(&n2);
            if handler.on_tangent_during_march(&current.point, angle) {
                break;
            }
            // Handler returned false (robust): use perturbation-based direction.
            // Reduce step size to tolerance.linear for finer tracing.
            step_size = tolerance.linear;
            match find_tangent_march_direction(
                surface1,
                current.u1,
                current.v1,
                surface2,
                current.u2,
                current.v2,
                &current.point,
                tolerance,
            ) {
                Some(dir) => dir * direction_sign,
                None => break,
            }
        } else {
            // Normalize tangent direction and apply direction sign.
            match t_raw.normalized() {
                Some(t) => t * direction_sign,
                None => break, // cannot normalize zero vector
            }
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

                    // Curvature-based cap: limit step in high-curvature regions
                    // to maintain chord-deviation quality. Uses second derivatives
                    // to estimate the intersection curve curvature.
                    let curvature = estimate_intersection_curvature(
                        surface1, sample.u1, sample.v1, surface2, sample.u2, sample.v2,
                    );
                    let curvature_cap =
                        curvature_based_step(curvature, tolerance.linear, min_step, max_step);
                    if curvature_cap < step_size {
                        step_size = step_size * 0.7 + curvature_cap * 0.3;
                    }

                    // Near-tangent step reduction: when surfaces are nearly
                    // parallel, the intersection curve direction is poorly
                    // defined and smaller steps improve stability.
                    let n1_post = surface1.normal_at(sample.u1, sample.v1);
                    let n2_post = surface2.normal_at(sample.u2, sample.v2);
                    let cross_len = n1_post.cross(&n2_post).length();
                    let norms = n1_post.length() * n2_post.length();
                    if norms > 1e-15 {
                        let sin_angle = cross_len / norms;
                        // Only intervene for very near-tangent (sin < 0.01 ≈ < 0.6°)
                        if sin_angle < 0.01 {
                            let factor = (sin_angle * 100.0).max(0.2);
                            step_size = (step_size * factor).max(min_step);
                        }
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

                    // Check for closure: use step_size as the tolerance since
                    // the march can only land within one step of the start.
                    // Require enough steps to have traversed a meaningful
                    // distance (avoid false closure on first few steps).
                    let min_steps_before_closure_check = 20;
                    let closure_tol = step_size * 1.5;
                    if samples.len() > min_steps_before_closure_check
                        && sample.point.distance_to(&start.point) < closure_tol
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

#[cfg(test)]
#[path = "marching_tests_coverage.rs"]
mod coverage_tests;
