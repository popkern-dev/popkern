//! Spine walker (§11.3) — phase 1b.
//!
//! Walks along a spine curve, calling [`solve_rolling_ball_contact`]
//! at each sample with the previous solution as the warm-start
//! initial guess. This continuation pattern is what makes the contact
//! solver practical for variable-radius fillets and for cones: only
//! the *first* sample pays the multi-start `closest_point` cost; every
//! subsequent sample converges in 1–3 Newton iterations from the
//! warm start.
//!
//! ## Phase 1b status (MVP)
//!
//! - **Fixed-step sampling** (`options.n_samples` evenly spaced).
//!   Adaptive step sizing per §11.3 step 4 is a phase-1b refinement.
//! - **Continuation seeding**: subsequent samples warm-start from the
//!   previous solution's `(uv0, uv1)`.
//! - **Initial seed**: the first sample uses `closest_point` on each
//!   surface against the spine endpoint (paid once per walk).
//! - **Failure handling**: any solver failure aborts the walk with
//!   the failing spine parameter and the underlying [`ContactError`].
//!   Phase 1c will add bisection / parametric detour recovery.
//!
//! See `docs/superpowers/plans/2026-04-14-fillets-rolling-ball-contact-solver.md`
//! for the broader phase-1 roadmap.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;

use crate::contact::{
    BallSide, ContactError, ContactOptions, ContactRequest, ContactSolution,
    solve_rolling_ball_contact,
};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Inputs to a spine walk.
///
/// The walker steps along `spine.parameter_range()` and produces one
/// [`ContactSolution`] per sample. The radius can vary along the
/// spine via the `radius_at` closure.
pub struct WalkRequest<'a> {
    /// First contact surface.
    pub surface_0: &'a dyn SurfaceEval,
    /// Second contact surface.
    pub surface_1: &'a dyn SurfaceEval,
    /// Which side of `surface_0` the rolling-ball center sits on.
    pub side_0: BallSide,
    /// Which side of `surface_1` the rolling-ball center sits on.
    pub side_1: BallSide,
    /// Radius profile along the spine. Receives the spine parameter
    /// (the actual `CurveEval` parameter, not normalized) and returns
    /// the local rolling-ball radius. For a constant-radius fillet
    /// pass `&|_| r`.
    pub radius_at: &'a dyn Fn(f64) -> f64,
    /// Spine curve. The walker uses `point_at` for the spine point
    /// and `tangent_at` as the spine plane normal at each sample.
    pub spine: &'a dyn CurveEval,
    /// Spine parameter range to walk. The walker steps from `t_range.0`
    /// to `t_range.1` in `n_samples` evenly-spaced steps. This is an
    /// explicit input (rather than `spine.parameter_range()`) because
    /// some `CurveEval` implementations report infinite ranges
    /// (e.g. `Line`) and the caller usually wants to walk only the
    /// portion corresponding to a specific edge.
    pub t_range: (f64, f64),
    /// Valid `(u, v)` range for `surface_0`.
    pub uv0_range: ((f64, f64), (f64, f64)),
    /// Valid `(u, v)` range for `surface_1`.
    pub uv1_range: ((f64, f64), (f64, f64)),
    /// Optional initial seed `(uv0, uv1)` at the start of the walk.
    /// `None` triggers a one-shot `closest_point` seed against
    /// `spine.point_at(t_start)` — this is the typical use.
    pub initial_seed: Option<((f64, f64), (f64, f64))>,
}

/// Tunable knobs for the walker. Use [`WalkOptions::default`] for a
/// commercial-grade preset (20 samples, contact-solver defaults).
#[derive(Debug, Clone, Copy)]
pub struct WalkOptions {
    /// Number of evenly-spaced spine samples to produce. Default 20.
    pub n_samples: usize,
    /// Per-sample contact-solver options. Defaults to
    /// [`ContactOptions::default`]. Subsequent samples typically
    /// converge in 1–3 iterations thanks to the warm start, so the
    /// default `max_iterations` of 20 has plenty of headroom.
    pub contact_options: ContactOptions,
}

impl Default for WalkOptions {
    fn default() -> Self {
        Self {
            n_samples: 20,
            contact_options: ContactOptions::default(),
        }
    }
}

/// One spine sample produced by [`walk_spine`].
#[derive(Debug, Clone)]
pub struct WalkSample {
    /// Spine curve parameter at this sample.
    pub t_spine: f64,
    /// The contact solution at this sample.
    pub solution: ContactSolution,
}

/// Aggregate output of [`walk_spine`].
#[derive(Debug, Clone)]
pub struct WalkResult {
    /// One sample per spine parameter, in walking order.
    pub samples: Vec<WalkSample>,
    /// Total Newton iterations across all samples (telemetry).
    pub total_newton_iterations: usize,
}

/// Failure modes for [`walk_spine`].
#[derive(Debug, Clone)]
pub enum WalkerError {
    /// The first contact solution at `t_start` failed. The spine seed
    /// was either wrong or the geometry is infeasible at the start.
    InitialContactFailed(ContactError),
    /// A subsequent step failed. `t_spine` identifies *which* sample.
    /// Phase 1c adds bisection / parametric detour to recover.
    StepFailed {
        /// Spine parameter where the step failed.
        t_spine: f64,
        /// 0-based index of the failing sample.
        sample_index: usize,
        /// The underlying contact-solver error.
        error: ContactError,
    },
    /// Could not seed the initial guess via `closest_point` — both
    /// surfaces' `closest_point` returned `None` at the spine start
    /// point.
    InitialSeedFailed {
        /// Which surface failed to project the spine start point.
        surface_index: usize,
    },
}

impl std::fmt::Display for WalkerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WalkerError::InitialContactFailed(e) => {
                write!(f, "spine walker: initial contact failed: {e}")
            }
            WalkerError::StepFailed {
                t_spine,
                sample_index,
                error,
            } => {
                write!(
                    f,
                    "spine walker: step {sample_index} at t={t_spine:.4} failed: {error}"
                )
            }
            WalkerError::InitialSeedFailed { surface_index } => {
                write!(
                    f,
                    "spine walker: surface {surface_index} could not seed initial guess"
                )
            }
        }
    }
}

impl std::error::Error for WalkerError {}

// ---------------------------------------------------------------------------
// Walker entry point
// ---------------------------------------------------------------------------

/// Walk a spine, producing a sequence of rolling-ball contact
/// solutions along it.
///
/// The first sample pays a one-shot `closest_point` cost (multi-start
/// Gauss-Newton against the spine start point); every subsequent
/// sample warm-starts from the previous sample's converged `(uv0, uv1)`,
/// so the inner loop typically converges in 1–3 Newton iterations
/// regardless of surface complexity.
///
/// **Phase 1b status**: fixed-step sampling. Adaptive step sizing
/// per §11.3 step 4 is a future refinement.
pub fn walk_spine(
    request: &WalkRequest,
    options: &WalkOptions,
    tolerance: &Tolerance,
) -> Result<WalkResult, WalkerError> {
    let n_samples = options.n_samples.max(2);
    let (t_start, t_end) = request.t_range;

    // Seed the initial guess at t = t_start.
    let p_start = request.spine.point_at(t_start);
    let (mut current_uv0, mut current_uv1) = match request.initial_seed {
        Some(seed) => seed,
        None => {
            let uv0 = request
                .surface_0
                .closest_point(&p_start)
                .ok_or(WalkerError::InitialSeedFailed { surface_index: 0 })?;
            let uv1 = request
                .surface_1
                .closest_point(&p_start)
                .ok_or(WalkerError::InitialSeedFailed { surface_index: 1 })?;
            (uv0, uv1)
        }
    };

    let mut samples = Vec::with_capacity(n_samples);
    let mut total_newton_iterations = 0_usize;

    for i in 0..n_samples {
        let t = if n_samples == 1 {
            t_start
        } else {
            t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1) as f64)
        };
        let p_spine = request.spine.point_at(t);
        let t_spine_tangent = request.spine.tangent_at(t);
        let radius = (request.radius_at)(t);

        let contact_request = ContactRequest {
            surface_0: request.surface_0,
            surface_1: request.surface_1,
            side_0: request.side_0,
            side_1: request.side_1,
            radius,
            initial_uv0: current_uv0,
            uv0_range: request.uv0_range,
            initial_uv1: current_uv1,
            uv1_range: request.uv1_range,
            spine_point: p_spine,
            spine_plane_normal: t_spine_tangent,
        };
        let solution = match solve_rolling_ball_contact(
            &contact_request,
            &options.contact_options,
            tolerance,
            None,
        ) {
            Ok(s) => s,
            Err(error) => {
                if i == 0 {
                    return Err(WalkerError::InitialContactFailed(error));
                }
                return Err(WalkerError::StepFailed {
                    t_spine: t,
                    sample_index: i,
                    error,
                });
            }
        };

        // Continuation: next sample warm-starts from this one.
        current_uv0 = solution.uv0;
        current_uv1 = solution.uv1;
        total_newton_iterations += solution.iterations;

        samples.push(WalkSample {
            t_spine: t,
            solution,
        });
    }

    Ok(WalkResult {
        samples,
        total_newton_iterations,
    })
}

#[cfg(test)]
#[path = "walker_tests.rs"]
mod tests;
