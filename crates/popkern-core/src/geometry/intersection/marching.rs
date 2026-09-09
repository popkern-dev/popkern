//! Differential marching for surface-surface intersection curves.

#[path = "marching_linalg.rs"]
mod marching_linalg;

#[path = "marching_newton.rs"]
pub mod marching_newton;

#[path = "marching_helpers.rs"]
mod marching_helpers;

#[allow(unused_imports)]
use marching_helpers::estimate_normal_near_singularity;
use marching_helpers::{
    find_tangent_march_direction, principal_branch_seed, recover_singular_normal, wrap_parameter,
};

pub(crate) use marching_newton::march_project;
use marching_newton::march_project_with_iters;
pub use marching_newton::{BoundaryPin, refine_boundary_crossing};

use crate::geometry::intersection::adaptive_step::{
    curvature_based_step, estimate_intersection_curvature,
};
use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
#[allow(unused_imports)]
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[path = "marching_types.rs"]
mod marching_types;
pub use marching_types::{BranchPoint, MarchSample};

#[path = "marching_rebranch.rs"]
mod marching_rebranch;
use marching_rebranch::rebranch_out_of_window;
pub(crate) use marching_rebranch::rebranch_sample_pub;

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
    let max_steps = 200_000;
    let epsilon = tolerance.solver();
    let min_step = tolerance.linear();
    let max_step = tolerance.linear() * 1000.0;
    // The step ceiling is the curvature-quality (chord-sag) bound
    // when the local curvature justifies more than the
    // tolerance-scaled default: a tolerance-scaled ceiling caps the
    // marchable arc length at max_steps × 1000·tol regardless of the
    // geometry's scale, and an M6 thread's helical flank × shank
    // curve (arc ≈ 28 at radius 3, PK-460) silently truncated at
    // arc 10 while the sag bound allowed 3.4× the step. Flat regions
    // (curvature below resolution) keep the tolerance-scaled ceiling
    // — sag imposes no bound there and unbounded growth would jump
    // features. The ×16 headroom bounds the worst case explicitly.
    let mut step_ceiling = max_step;
    let mut step_size = tolerance.linear() * 50.0;
    // Branch point candidates detected during march (spec §2.4).
    // Used by find_branch_points in starting_points.rs to seed
    // additional curve traces from branch locations.
    #[allow(unused_assignments)]
    let mut branch_points: Vec<BranchPoint> = Vec::new();

    let mut samples = vec![start.clone()];
    // Termination cause for the march trace: every break overwrites it.
    let mut termination = "step-budget";
    // Wrong-preimage landings allow a few halve-and-retry attempts
    // (a smaller step near the phantom boundary can stay real), but
    // a march that keeps converging into the extrapolated basin must
    // STOP: unbounded retries re-projected the same phantom family
    // ~800k times on one two-start rod cut (5× wall time, measured).
    let mut phantom_retries = 0u32;
    let mut phantom_cache: Vec<Point3> = Vec::new();
    // ORBIT GUARD state (see the near-return check after the strict
    // closure test): accumulated arc, arc at the last counted
    // near-return, and how many well-separated near-returns failed
    // strict closure.
    let mut arc_acc = 0.0_f64;
    let mut last_near_arc = f64::NEG_INFINITY;
    let mut near_returns = 0u32;

    for _ in 0..max_steps {
        let current = &samples[samples.len() - 1];

        // Compute surface normals at current parameters.
        let mut n1 = surface1.normal_at(current.u1, current.v1);
        let mut n2 = surface2.normal_at(current.u2, current.v2);

        // Check singularity: degenerate normal on surface 1.
        if n1.length() < epsilon {
            if handler.on_singularity(&current.point) {
                termination = "singularity-1";
                break;
            }
            // Handler returned false (robust): estimate near the
            // singularity, walking forward with a halving step.
            match recover_singular_normal(surface1, current.u1, current.v1, step_size, tolerance) {
                Some((estimated, reduced)) => {
                    n1 = estimated;
                    step_size = reduced;
                }
                None => {
                    termination = "singular-normal-1-unrecovered";
                    break;
                }
            }
        }

        // Check singularity: degenerate normal on surface 2.
        if n2.length() < epsilon {
            if handler.on_singularity(&current.point) {
                termination = "singularity-2";
                break;
            }
            match recover_singular_normal(surface2, current.u2, current.v2, step_size, tolerance) {
                Some((estimated, reduced)) => {
                    n2 = estimated;
                    step_size = reduced;
                }
                None => {
                    termination = "singular-normal-2-unrecovered";
                    break;
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
                termination = "tangency";
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
                None => {
                    termination = "tangent-direction-lost";
                    break;
                }
            }
        } else {
            // Normalize tangent direction and apply direction sign.
            match t_raw.normalized() {
                Some(t) => t * direction_sign,
                None => {
                    termination = "zero-tangent";
                    break;
                }
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
                    // Wrong-preimage landings re-branch into the query
                    // windows before the domain check; a landing that
                    // cannot re-branch is an extrapolation phantom and
                    // counts as a failed projection (the halve-step
                    // retry below keeps the march on the real curve —
                    // see `rebranch_out_of_window`).
                    if !rebranch_out_of_window(
                        surface1,
                        u1_range,
                        v1_range,
                        surface2,
                        u2_range,
                        v2_range,
                        &mut sample,
                        tolerance,
                        &mut phantom_cache,
                    ) {
                        phantom_retries += 1;
                        retry_step *= 0.5;
                        if retry_step < min_step || phantom_retries > 8 {
                            break;
                        }
                        continue;
                    }
                    // Adaptive step based on iteration count.
                    if iters < 3 {
                        step_size = (retry_step * 1.5).min(step_ceiling);
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
                    let curvature_cap = curvature_based_step(
                        curvature,
                        tolerance.linear(),
                        min_step,
                        max_step * 16.0,
                    );
                    if curvature_cap < step_size {
                        step_size = step_size * 0.7 + curvature_cap * 0.3;
                    }
                    // Measurable curvature makes the sag bound the
                    // honest ceiling (never below the tolerance-scaled
                    // default; see the declaration).
                    step_ceiling = if curvature > 1e-12 {
                        max_step.max(curvature_cap)
                    } else {
                        max_step
                    };

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
                    let (period_u1, period_v1) =
                        crate::geometry::chart_periods::marching_periods(surface1);
                    let (period_u2, period_v2) =
                        crate::geometry::chart_periods::marching_periods(surface2);
                    let (periodic_u1, periodic_v1) = (period_u1.is_some(), period_v1.is_some());
                    let (periodic_u2, periodic_v2) = (period_u2.is_some(), period_v2.is_some());
                    sample.u1 = wrap_parameter(sample.u1, u1_range, period_u1);
                    sample.v1 = wrap_parameter(sample.v1, v1_range, period_v1);
                    sample.u2 = wrap_parameter(sample.u2, u2_range, period_u2);
                    sample.v2 = wrap_parameter(sample.v2, v2_range, period_v2);

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
                        termination = "domain-exit";
                        if crate::trace::trace_enabled("march") {
                            eprintln!(
                                "      [march-exit] u1={:.6} v1={:.6} u2={:.6} v2={:.6} flags {}{}{}{} w1={u1_range:?} {v1_range:?} w2={u2_range:?} {v2_range:?}",
                                sample.u1,
                                sample.v1,
                                sample.u2,
                                sample.v2,
                                exits_u1 as u8,
                                exits_v1 as u8,
                                exits_u2 as u8,
                                exits_v2 as u8
                            );
                        }
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
                        trace_march_termination(termination, &samples);
                        return samples;
                    }

                    // Check for closure: use step_size as the tolerance since
                    // the march can only land within one step of the start.
                    // Require enough steps to have traversed a meaningful
                    // distance (avoid false closure on first few steps).
                    let min_steps_before_closure_check = 20;
                    let closure_tol = step_size * 1.5;
                    let dist_to_start = sample.point.distance_to(&start.point);
                    if samples.len() > min_steps_before_closure_check && dist_to_start < closure_tol
                    {
                        samples.push(sample);
                        trace_march_termination("closure", &samples);
                        return samples;
                    }

                    // ORBIT GUARD: a trace that keeps returning NEAR
                    // its start without ever strictly closing is
                    // re-traversing the same near-closed contact ring
                    // — a numerical orbit, not new geometry. Measured
                    // (PK-460, conic cube corner): a near-tangent
                    // ring of circumference ~1 returned to within
                    // 3.8e-3..7.6e-3 of its start (strict closure
                    // needed 2.8e-3) and orbited ~400 laps to the
                    // step budget — 200k samples per trace, 41 GB
                    // peak RSS, an OOM SIGKILL in the sweep. Two
                    // near-returns separated by a full lap of arc
                    // certify the orbit; terminating there hands
                    // downstream the same honestly-named truncation
                    // the step budget produced, at lap ~2 instead of
                    // lap ~400. Open marches (a thread helix) never
                    // return near their start; genuinely closed
                    // loops strictly close on the first lap.
                    //
                    // A near-return counts ONLY when the march passes
                    // the start region travelling in the START'S OWN
                    // direction: an orbit re-treads the same curve, so
                    // it arrives aligned, while a legitimate open
                    // trace can pass NEAR its start transversally — a
                    // deep-grind plane cuts a coil tube in two arcs of
                    // one S-shaped trace whose far arc passes within
                    // the band travelling the other way (measured:
                    // the direction-blind guard truncated it mid-curve
                    // and the boolean failed validation at 6.9e-5).
                    let near_band = closure_tol * 8.0;
                    let prev_point = samples[samples.len() - 1].point;
                    arc_acc += sample.point.distance_to(&prev_point);
                    if samples.len() > min_steps_before_closure_check
                        && dist_to_start < near_band
                        && arc_acc - last_near_arc > 50.0 * near_band
                    {
                        let aligned = {
                            let step_dir = sample.point - prev_point;
                            let start_dir = samples[1].point - samples[0].point;
                            let lens = step_dir.length() * start_dir.length();
                            lens > 1e-30 && step_dir.dot(&start_dir) / lens > 0.9
                        };
                        if aligned {
                            last_near_arc = arc_acc;
                            near_returns += 1;
                            if near_returns >= 2 {
                                samples.push(sample);
                                trace_march_termination("orbit-no-closure", &samples);
                                return samples;
                            }
                        }
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
            termination = "newton-exhausted";
            break;
        }
    }
    trace_march_termination(termination, &samples);

    samples
}

#[path = "marching_bidirectional.rs"]
mod marching_bidirectional;
pub(crate) use marching_bidirectional::march;
use marching_bidirectional::trace_march_termination;

#[cfg(test)]
#[path = "marching_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "marching_tests_coverage.rs"]
mod coverage_tests;

#[cfg(test)]
#[path = "marching_tests_torus.rs"]
mod torus_tests;
