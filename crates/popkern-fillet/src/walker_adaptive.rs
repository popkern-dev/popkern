//! Adaptive-stepping inner loop for the spine walker, covering
//! §11.3 step 4 (adaptive dt) and §11.4 interval splitting (bisect
//! dt on step failure). Split out of `walker.rs` to keep that file
//! under the 500-line structural cap.

use popkern_core::tolerance::Tolerance;

use crate::contact::ContactOptions;
use crate::walker::{WalkRequest, WalkResult, WalkSample, WalkerError, step_at};

/// Concrete adaptive parameters extracted from `WalkMode::Adaptive`
/// for passing into the inner loop.
#[derive(Debug, Clone, Copy)]
pub(crate) struct AdaptiveParams {
    pub initial_dt: f64,
    pub min_dt: f64,
    pub max_dt: f64,
    pub max_samples: usize,
    pub grow_threshold: usize,
    pub shrink_threshold: usize,
    pub grow_factor: f64,
    pub shrink_factor: f64,
}

/// Adaptive walker inner loop.
///
/// Starts at `initial_dt`, grows when contact-solver convergence is
/// fast (< `grow_threshold` Newton iters), shrinks when slow
/// (> `shrink_threshold`). Caps total samples at `max_samples` and
/// bails if `dt` falls below `min_dt`.
///
/// On step failure (the contact solver returns a `ContactError` for
/// any reason — singular Jacobian, drifted off, divergence, …) the
/// loop **bisects `dt` and retries** up to `STEP_RETRY_LIMIT` times
/// from the same `t`. This is §11.4 "interval splitting": a
/// smaller spine step often lands on the feasible side of a
/// transient singularity (tangent contact, near-degenerate normal,
/// etc.). If all retries fail the original error propagates.
pub(crate) fn walk_adaptive(
    request: &WalkRequest,
    tolerance: &Tolerance,
    contact_options: &ContactOptions,
    mut current_uv0: (f64, f64),
    mut current_uv1: (f64, f64),
    params: AdaptiveParams,
) -> Result<WalkResult, WalkerError> {
    let (t_start, t_end) = request.t_range;
    let direction = (t_end - t_start).signum();
    let mut samples: Vec<WalkSample> = Vec::with_capacity(params.max_samples.min(64));
    let mut total_newton_iterations = 0_usize;

    // Sample 0 at t_start.
    let solution = step_at(
        request,
        tolerance,
        contact_options,
        0,
        t_start,
        current_uv0,
        current_uv1,
    )?;
    current_uv0 = solution.uv0;
    current_uv1 = solution.uv1;
    total_newton_iterations += solution.iterations;
    samples.push(WalkSample {
        t_spine: t_start,
        solution,
    });

    let mut t = t_start;
    let mut dt = params.initial_dt.abs() * direction;
    const STEP_RETRY_LIMIT: usize = 6;

    while samples.len() < params.max_samples {
        let next_t = t + dt;
        let snap_to_end = (next_t - t_end) * direction >= 0.0;
        let target_t = if snap_to_end { t_end } else { next_t };

        let i = samples.len();
        let mut try_t = target_t;
        let mut try_dt = dt;
        let mut snapping = snap_to_end;
        let mut step_solution = None;
        let mut last_error = None;
        for _retry in 0..STEP_RETRY_LIMIT {
            match step_at(
                request,
                tolerance,
                contact_options,
                i,
                try_t,
                current_uv0,
                current_uv1,
            ) {
                Ok(s) => {
                    step_solution = Some(s);
                    dt = try_dt;
                    break;
                }
                Err(e) => {
                    last_error = Some(e);
                    try_dt *= 0.5;
                    if try_dt.abs() < params.min_dt {
                        break;
                    }
                    try_t = t + try_dt;
                    snapping = (try_t - t_end) * direction >= 0.0;
                    if snapping {
                        try_t = t_end;
                    }
                }
            }
        }
        let solution = match step_solution {
            Some(s) => s,
            None => return Err(last_error.expect("retry loop must record at least one error")),
        };

        let iters = solution.iterations;
        current_uv0 = solution.uv0;
        current_uv1 = solution.uv1;
        total_newton_iterations += iters;
        samples.push(WalkSample {
            t_spine: try_t,
            solution,
        });
        t = try_t;

        if snapping {
            break;
        }

        if iters < params.grow_threshold {
            dt *= params.grow_factor;
        } else if iters > params.shrink_threshold {
            dt *= params.shrink_factor;
        }
        let dt_abs = dt.abs().clamp(params.min_dt, params.max_dt);
        dt = dt_abs * direction;

        if dt.abs() < params.min_dt - 1.0e-15 {
            break;
        }
    }

    Ok(WalkResult {
        samples,
        total_newton_iterations,
    })
}
