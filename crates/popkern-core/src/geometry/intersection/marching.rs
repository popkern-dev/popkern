//! Differential marching for surface-surface intersection curves.

#[path = "marching_helpers.rs"]
mod marching_helpers;
#[path = "marching_linalg.rs"]
mod marching_linalg;
#[path = "marching_newton.rs"]
pub mod marching_newton;
use marching_helpers::{
    arrival_alignment, exits_window, find_tangent_march_direction, principal_branch_seed,
    recover_singular_normal, wrap_parameter,
};

pub(crate) use marching_newton::march_project;
use marching_newton::march_project_with_iters;
pub use marching_newton::{BoundaryPin, refine_boundary_crossing};

use crate::geometry::Point3;
use crate::geometry::intersection::adaptive_step::{
    curvature_based_step, estimate_intersection_curvature, sag_step_ceiling,
};
use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

#[path = "marching_types.rs"]
mod marching_types;
pub(crate) use marching_types::MarchEnd;
pub use marching_types::{BranchPoint, MarchSample};

#[path = "marching_crossing.rs"]
mod marching_crossing;
pub(crate) use marching_crossing::departures_at_ends;
#[path = "marching_rebranch.rs"]
mod marching_rebranch;
use marching_rebranch::rebranch_out_of_window;
pub(crate) use marching_rebranch::rebranch_sample_pub;

/// March in one direction along the intersection curve: its samples,
/// and how it ended (`MarchEnd`).
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
) -> (Vec<MarchSample>, MarchEnd) {
    let max_steps = 400_000;
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
    // features. Past the ×16 cap: `sag_step_ceiling`.
    let mut step_ceiling = max_step;
    let mut step_size = tolerance.linear() * 50.0;
    // Branch point candidates detected during march (spec §2.4), seeds
    // for further traces (find_branch_points in starting_points.rs).
    #[allow(unused_assignments)]
    let mut branch_points: Vec<BranchPoint> = Vec::new();

    let mut samples = vec![start.clone()];
    let mut termination = "step-budget"; // the trace's cause; breaks overwrite
    // Wrong-preimage landings allow a few halve-and-retry attempts
    // (a smaller step near the phantom boundary can stay real), but
    // a march that keeps converging into the extrapolated basin must
    // STOP: unbounded retries re-projected the same phantom family
    // ~800k times on one two-start rod cut (5× wall time, measured).
    let mut phantom_retries = 0u32;
    let mut phantom_cache: Vec<Point3> = Vec::new();
    // ORBIT GUARD state (see the near-return check after the strict closure
    // test): arc so far, at the last near-return, and the near-returns.
    let mut arc_acc = 0.0_f64;
    let mut last_near_arc = f64::NEG_INFINITY;
    let mut near_returns = 0u32;
    let mut stalls = 0u32; // STALL GUARD run length
    // A tangent crossing ends the trace (`marching_crossing`).
    let crossing = marching_crossing::Crossing::new(
        (surface1, (u1_range, v1_range)),
        (surface2, (u2_range, v2_range)),
        tolerance,
    );

    'march: for _ in 0..max_steps {
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

        // Tangent where the normals' angle has a near-zero sine. Chart
        // normals carry the charts' speeds: at a thousandth scale the
        // bare |N1 x N2| read 3e-13 at a 0.17 crossing and the march
        // stalled at the tolerance step (SW-E 765).
        let t_normalized = if t_len < epsilon * n1.length() * n2.length() {
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
                    let sag_ceiling = sag_step_ceiling(
                        (surface1, (sample.u1, sample.v1)),
                        (surface2, (sample.u2, sample.v2)),
                        curvature,
                        tolerance.linear(),
                        max_step * 16.0,
                    );
                    let curvature_cap =
                        curvature_based_step(curvature, tolerance.linear(), min_step, sag_ceiling);
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

                    // Wrap periodic parameters before the domain check,
                    // closed NURBS charts too (their seam bounds only the
                    // chart; `geometry::chart_periods`).
                    let (period_u1, period_v1) =
                        crate::geometry::chart::ChartDescriptor::periods_of(surface1);
                    let (period_u2, period_v2) =
                        crate::geometry::chart::ChartDescriptor::periods_of(surface2);
                    sample.u1 = wrap_parameter(sample.u1, u1_range, period_u1);
                    sample.v1 = wrap_parameter(sample.v1, v1_range, period_v1);
                    sample.u2 = wrap_parameter(sample.u2, u2_range, period_u2);
                    sample.v2 = wrap_parameter(sample.v2, v2_range, period_v2);

                    // Check domain bounds on every axis the window
                    // bounds: non-periodic axes and sub-period windows
                    // of periodic ones (see `exits_window`).
                    let exits_u1 = exits_window(sample.u1, u1_range, period_u1);
                    let exits_v1 = exits_window(sample.v1, v1_range, period_v1);
                    let exits_u2 = exits_window(sample.u2, u2_range, period_u2);
                    let exits_v2 = exits_window(sample.v2, v2_range, period_v2);
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
                        return (samples, MarchEnd::Open);
                    }

                    // Closure: back within a step of the start, the gap under
                    // a quarter of the arc so far (a chord falls that short of
                    // its arc only where the trace turned back), the arc past
                    // the ball, arriving the start's way. A step count closed a
                    // loop in 15 on lap two (1392); the ball alone, a third of a
                    // small loop under the near-tangent floor (1043, 935); a
                    // thin loop's far branch passed its start 8.6e-4 off (375).
                    let min_steps_before_closure_check = 20;
                    let closure_tol = step_size * 1.5;
                    let dist_to_start = sample.point.distance_to(&start.point);
                    let prev_point = samples[samples.len() - 1].point;
                    arc_acc += sample.point.distance_to(&prev_point);
                    let alignment = arrival_alignment(&samples, prev_point, sample.point);
                    if dist_to_start < closure_tol.min(0.25 * arc_acc)
                        && arc_acc > 2.0 * closure_tol
                        && alignment > 0.0
                    {
                        samples.push(sample);
                        trace_march_termination("closure", &samples);
                        return (samples, MarchEnd::Closed);
                    }

                    // ORBIT GUARD: a trace returning NEAR its start without
                    // strictly closing re-traverses the same near-closed
                    // ring — a numerical orbit, not new geometry. Measured
                    // (PK-460, conic cube corner): a near-tangent
                    // ring of circumference ~1 returned to within
                    // 3.8e-3..7.6e-3 of its start (strict closure
                    // needed 2.8e-3) and orbited ~400 laps to the
                    // step budget — 200k samples per trace, 41 GB
                    // peak RSS, an OOM SIGKILL in the sweep. Two
                    // near-returns separated by a full lap of arc
                    // certify the orbit; terminating there hands
                    // downstream the step budget's truncation at lap ~2,
                    // not ~400. Open marches (a thread helix) never return
                    // near their start; closed loops close on lap one.
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
                    if samples.len() > min_steps_before_closure_check
                        && dist_to_start < near_band
                        && arc_acc - last_near_arc > 50.0 * near_band
                        && alignment > 0.9
                    {
                        last_near_arc = arc_acc;
                        near_returns += 1;
                        if near_returns >= 2 {
                            samples.push(sample);
                            trace_march_termination("orbit-no-closure", &samples);
                            return (samples, MarchEnd::Open);
                        }
                    }

                    // STALL GUARD: steps each advancing under 1% of the
                    // step asked for are a point, not a curve (a 70 um
                    // Gregory cap's isolated corner contact ground 200k
                    // steps per trace, 8.3 s per check).
                    let stalled = sample.point.distance_to(&prev_point) < 0.01 * step_size;
                    stalls = if stalled { stalls + 1 } else { 0 };
                    if stalls >= 64 {
                        samples.push(sample);
                        trace_march_termination("stalled", &samples);
                        return (samples, MarchEnd::Open);
                    }
                    if crossing.push(&mut samples, sample, t_normalized, step_size) {
                        termination = "tangent-crossing";
                        break 'march;
                    }
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
            // Newton gives out at a tangent crossing: the trace ends there.
            termination = match crossing.rescue(&mut samples, t_normalized, step_size) {
                true => "tangent-crossing",
                false => "newton-exhausted",
            };
            if termination == "newton-exhausted" {
                trace_march_exhaustion(surface1, surface2, samples.last());
            }
            break;
        }
    }
    trace_march_termination(termination, &samples);
    let end = match termination {
        "tangent-crossing" => MarchEnd::Node,
        _ => MarchEnd::Open,
    };
    (samples, end)
}

#[path = "marching_bidirectional.rs"]
mod marching_bidirectional;
pub(crate) use marching_bidirectional::march;
use marching_bidirectional::{trace_march_exhaustion, trace_march_termination};

#[cfg(test)]
#[path = "marching_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "marching_tests_coverage.rs"]
mod coverage_tests;

#[cfg(test)]
#[path = "marching_tests_torus.rs"]
mod torus_tests;
