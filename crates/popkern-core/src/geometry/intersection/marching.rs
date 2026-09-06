//! Differential marching for surface-surface intersection curves.

#[path = "marching_linalg.rs"]
mod marching_linalg;

#[path = "marching_newton.rs"]
pub(crate) mod marching_newton;

#[path = "marching_helpers.rs"]
mod marching_helpers;

use marching_helpers::{
    estimate_normal_near_singularity, find_tangent_march_direction, principal_branch_seed,
    wrap_parameter,
};

pub use marching_newton::BoundaryPin;
pub(crate) use marching_newton::march_project;
use marching_newton::march_project_with_iters;
use marching_newton::newton_boundary_crossing;

use crate::geometry::intersection::adaptive_step::{
    curvature_based_step, estimate_intersection_curvature,
};
use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
#[allow(unused_imports)]
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// A detected branch point where the intersection curve splits.
#[derive(Debug, Clone)]
pub struct BranchPoint {
    /// The 3D location of the branch point.
    pub point: Point3,
    /// Parameters on surface 1.
    pub u1: f64,
    /// Parameters on surface 1.
    pub v1: f64,
    /// The Jacobian condition indicator (sin of normal angle) at this point.
    pub condition: f64,
}

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
    // The SEED must sit on the principal period branch every marched
    // step is wrapped onto — mechanism and measured rationale in
    // [`principal_branch_seed`] (PK-454).
    let start = principal_branch_seed(
        start, surface1, u1_range, v1_range, surface2, u2_range, v2_range,
    );
    let start = &start;

    // March in +T direction.
    let forward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, 1.0, tolerance, handler,
    );

    // A forward trace that CLOSED already covers the whole loop —
    // the backward march would re-trace it lap-for-lap (measured on
    // the coil-grind SSI: every closed marched loop carried 2x
    // samples as two duplicate laps that seam-splitting turned into
    // duplicate records and consolidation had to dedupe — PK-450).
    let forward_closed = marching_helpers::forward_trace_closed(&forward, tolerance);

    // March in -T direction (only when the loop is still open).
    let mut backward = if forward_closed {
        Vec::new()
    } else {
        march_one_direction(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, -1.0, tolerance,
            handler,
        )
    };

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
        tolerance.linear()
    };
    let closure_tol = (avg_step * 2.0).max(tolerance.linear());
    let is_closed = if samples.len() > 10 {
        let first = &samples[0];
        let last = &samples[samples.len() - 1];
        first.point.distance_to(&last.point) < closure_tol
    } else {
        false
    };

    (samples, is_closed)
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
    let epsilon = tolerance.solver();
    let min_step = tolerance.linear();
    let max_step = tolerance.linear() * 1000.0;
    let mut step_size = tolerance.linear() * 50.0;
    // Branch point candidates detected during march (spec §2.4).
    // Used by find_branch_points in starting_points.rs to seed
    // additional curve traces from branch locations.
    #[allow(unused_assignments)]
    let mut branch_points: Vec<BranchPoint> = Vec::new();

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
            let delta = tolerance.linear() * 10.0;
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
            let delta = tolerance.linear() * 10.0;
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
            // Reduce step size to tolerance.linear() for finer tracing.
            step_size = tolerance.linear();
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
                        curvature_based_step(curvature, tolerance.linear(), min_step, max_step);
                    if curvature_cap < step_size {
                        step_size = step_size * 0.7 + curvature_cap * 0.3;
                    }

                    // Near-tangent / branch point detection: monitor the
                    // Jacobian condition via |N1 × N2| / (|N1| * |N2|).
                    // When this drops toward zero, the surfaces are tangent
                    // and the intersection curve may branch.
                    let n1_post = surface1.normal_at(sample.u1, sample.v1);
                    let n2_post = surface2.normal_at(sample.u2, sample.v2);
                    let cross_len = n1_post.cross(&n2_post).length();
                    let norms = n1_post.length() * n2_post.length();
                    if norms > tolerance.numerical_floor() {
                        let sin_angle = cross_len / norms;

                        // Branch point detection: record location where
                        // condition drops below threshold (spec §2.4).
                        // Also performs second-order contact analysis (§2.5)
                        // using curvature comparison to determine crossing
                        // vs touching.
                        if sin_angle < 0.05 && sin_angle > tolerance.solver() {
                            // Second-order analysis: compare curvatures
                            let suu1 = surface1.second_partial_uu(sample.u1, sample.v1);
                            let svv1 = surface1.second_partial_vv(sample.u1, sample.v1);
                            let suu2 = surface2.second_partial_uu(sample.u2, sample.v2);
                            let svv2 = surface2.second_partial_vv(sample.u2, sample.v2);
                            let n_unit =
                                n1_post * (1.0 / n1_post.length().max(tolerance.numerical_floor()));
                            let curv_diff =
                                (suu1 - suu2).dot(&n_unit).abs() + (svv1 - svv2).dot(&n_unit).abs();

                            // If curvature difference is significant, it's a
                            // crossing or branch; if near zero, it's touching.
                            if curv_diff > tolerance.linear() {
                                branch_points.push(BranchPoint {
                                    point: sample.point,
                                    u1: sample.u1,
                                    v1: sample.v1,
                                    condition: sin_angle,
                                });
                            }
                        }

                        // Near-tangent stepping: the |S1-S2| < ε
                        // convergence band localizes the curve only to
                        // a width of ~ε/sin across, so steps BELOW that
                        // width trace residual noise, not the curve —
                        // the former reduction crushed the step toward
                        // min_step here and ground 10k-step fragments
                        // on grazing fillet contact. The band is the
                        // honest step floor; when it exceeds the
                        // maximum step, the intersection is not
                        // resolvable as a curve at this tolerance and
                        // the trace ends.
                        // Near-tangent stepping: the |S1-S2| < ε
                        // convergence band localizes the curve only to
                        // a width of ~ε/sin across, so steps BELOW the
                        // band width trace residual noise, not the
                        // curve — the former reduction crushed the
                        // step toward min_step exactly here, grinding
                        // 10k-step fragments on grazing fillet
                        // contact. The band (capped at max_step) is
                        // the honest step FLOOR. The march never
                        // terminates on unresolvability: intersection
                        // curves legitimately pass THROUGH ultra-
                        // tangent osculation stretches (measured: the
                        // octagon8 setback pipeline needs such a
                        // curve whole), and crossing them at max_step
                        // bounds both the lateral noise and the cost.
                        if sin_angle < 0.01 {
                            let band = tolerance.linear() / sin_angle.max(epsilon);
                            step_size = step_size.max(band.min(max_step));
                        }
                    }

                    // Wrap parameters for periodic surfaces before the
                    // domain check — including geometrically-closed
                    // NURBS charts, whose seam is not a boundary of
                    // the surface, only of the chart (the wrap keeps
                    // evaluation inside the chart; see
                    // `geometry::chart_periods`).
                    let (periodic_u1, periodic_v1) =
                        crate::geometry::chart_periods::marching_periodicity(surface1);
                    let (periodic_u2, periodic_v2) =
                        crate::geometry::chart_periods::marching_periodicity(surface2);
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
                        // Exited domain; stop marching in this
                        // direction. The exiting sample IS emitted:
                        // window trimming downstream refines the
                        // exact boundary crossing from the straddling
                        // pair, and discarding the outside sample
                        // left the trace up to one whole step short
                        // of the window when the last inside sample
                        // did not reach the pad band (measured: a
                        // blade-girdle record ended 9.1e-5 shy of its
                        // span-boundary junction, the open-chain weld
                        // bridged the gap, and the welded tail bent
                        // 3.3e-5 off the span surface — PK-443).
                        samples.push(sample);
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

/// Refine a surface×surface intersection point exactly onto a chart
/// iso-boundary (see [`BoundaryPin`]): Newton on `S1 = S2` with one
/// chart coordinate pinned to the boundary value. Used by boolean
/// junction-endpoint refinement — per-pair clip crossings evaluate
/// FITTED record curves and carry the fit's lateral slop (~1e-4–1e-3
/// measured at grazing junction clusters), while this solve lands on
/// the TRUE intersection∩boundary point. Returns the refined chart
/// coordinates of both surfaces and the 3D point, or `None` when the
/// solve is singular (curve tangent to the boundary) or does not
/// converge — callers keep the unrefined endpoint then.
#[allow(clippy::too_many_arguments)]
pub fn refine_boundary_crossing(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    pin: BoundaryPin,
    tolerance: &Tolerance,
) -> Option<(f64, f64, f64, f64, Point3)> {
    newton_boundary_crossing(surface1, u1, v1, surface2, u2, v2, pin, tolerance)
        .map(|s| (s.u1, s.v1, s.u2, s.v2, s.point))
}

#[cfg(test)]
#[path = "marching_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "marching_tests_coverage.rs"]
mod coverage_tests;

#[cfg(test)]
#[path = "marching_tests_torus.rs"]
mod torus_tests;
