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
    /// Periodicity of `surface_0`'s `(u, v)` parameters, in the
    /// `(period_u, period_v)` form returned by
    /// `SurfaceGeometry::periodicity()`. `Some(period)` indicates a
    /// wrap-around dimension; `None` means non-periodic.
    ///
    /// The walker uses periodicity to compute a per-iteration UV
    /// window centered on the previous sample's solution that's wide
    /// enough to let Newton step across the period boundary without
    /// hitting `uv*_range` and reporting a spurious `OutOfBounds`.
    /// Without this, a walk around the base of a cone hits the
    /// boundary at u = ±π and the last sample fails.
    ///
    /// Default: `(None, None)` (no periodicity adjustment) — pass
    /// `surface.periodicity()` to enable.
    pub surface_0_periodicity: (Option<f64>, Option<f64>),
    /// Periodicity of `surface_1`'s `(u, v)` parameters. See
    /// [`Self::surface_0_periodicity`].
    pub surface_1_periodicity: (Option<f64>, Option<f64>),
}

/// How the walker advances along the spine.
#[derive(Debug, Clone, Copy)]
pub enum WalkMode {
    /// Fixed-step: produce exactly `n_samples` evenly-spaced samples.
    /// Simple and predictable; preferred when the caller wants a
    /// specific output size for surface fitting.
    Fixed {
        /// Number of evenly-spaced samples to produce.
        n_samples: usize,
    },
    /// Adaptive: start at `initial_dt`, grow when contact-solver
    /// convergence is fast (`< grow_threshold` Newton iterations),
    /// shrink when slow (`> shrink_threshold`) or when the Jacobian
    /// condition exceeds `condition_threshold`. Targets 2–5
    /// iterations per step in steady state per §11.3 step 4.
    /// Caps total samples at `max_samples` as a safety net.
    Adaptive {
        /// Starting step size in spine-parameter units.
        initial_dt: f64,
        /// Lower bound on step size; below this the walker bails.
        min_dt: f64,
        /// Upper bound on step size.
        max_dt: f64,
        /// Maximum number of samples (safety cap on the loop).
        max_samples: usize,
        /// Grow `dt` by `grow_factor` when the previous step's
        /// contact solver used fewer than this many Newton iterations.
        grow_threshold: usize,
        /// Shrink `dt` by `shrink_factor` when the previous step's
        /// contact solver used more than this many Newton iterations.
        shrink_threshold: usize,
        /// Multiplicative factor applied to `dt` on growth.
        grow_factor: f64,
        /// Multiplicative factor applied to `dt` on shrinkage.
        shrink_factor: f64,
    },
}

/// Tunable knobs for the walker. Use [`WalkOptions::default`] for a
/// commercial-grade preset (fixed 20 samples, contact-solver defaults).
#[derive(Debug, Clone, Copy)]
pub struct WalkOptions {
    /// Stepping mode — fixed or adaptive.
    pub mode: WalkMode,
    /// Per-sample contact-solver options. Defaults to
    /// [`ContactOptions::default`]. Subsequent samples typically
    /// converge in 1–3 iterations thanks to the warm start.
    pub contact_options: ContactOptions,
}

impl Default for WalkOptions {
    fn default() -> Self {
        Self {
            mode: WalkMode::Fixed { n_samples: 20 },
            contact_options: ContactOptions::default(),
        }
    }
}

impl WalkOptions {
    /// Commercial-grade adaptive preset. Initial step is `span / 10`,
    /// bounded `[span / 200, span / 5]`, capped at 200 samples,
    /// targeting 2–5 Newton iterations per step. Construct with the
    /// spine span (`t_end - t_start`) so the bounds scale to the
    /// actual walk.
    pub fn adaptive(spine_span: f64) -> Self {
        let span = spine_span.abs().max(1.0e-12);
        Self {
            mode: WalkMode::Adaptive {
                initial_dt: span / 10.0,
                min_dt: span / 200.0,
                max_dt: span / 5.0,
                max_samples: 200,
                grow_threshold: 4,
                shrink_threshold: 7,
                grow_factor: 1.5,
                shrink_factor: 0.5,
            },
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
    let (t_start, _t_end) = request.t_range;

    // Seed the initial guess at t = t_start.
    let p_start = request.spine.point_at(t_start);
    let (initial_uv0, initial_uv1) = match request.initial_seed {
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

    match options.mode {
        WalkMode::Fixed { n_samples } => walk_fixed(
            request,
            tolerance,
            &options.contact_options,
            initial_uv0,
            initial_uv1,
            n_samples.max(2),
        ),
        WalkMode::Adaptive {
            initial_dt,
            min_dt,
            max_dt,
            max_samples,
            grow_threshold,
            shrink_threshold,
            grow_factor,
            shrink_factor,
        } => crate::walker_adaptive::walk_adaptive(
            request,
            tolerance,
            &options.contact_options,
            initial_uv0,
            initial_uv1,
            crate::walker_adaptive::AdaptiveParams {
                initial_dt,
                min_dt,
                max_dt,
                max_samples,
                grow_threshold,
                shrink_threshold,
                grow_factor,
                shrink_factor,
            },
        ),
    }
}

/// Fixed-step inner loop. Walks `n_samples` evenly-spaced spine
/// parameters from `t_start` to `t_end`, warm-starting each Newton
/// solve from the previous sample's solution.
fn walk_fixed(
    request: &WalkRequest,
    tolerance: &Tolerance,
    contact_options: &ContactOptions,
    mut current_uv0: (f64, f64),
    mut current_uv1: (f64, f64),
    n_samples: usize,
) -> Result<WalkResult, WalkerError> {
    let (t_start, t_end) = request.t_range;
    let mut samples = Vec::with_capacity(n_samples);
    let mut total_newton_iterations = 0_usize;

    for i in 0..n_samples {
        let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1) as f64);
        let solution = step_at(
            request,
            tolerance,
            contact_options,
            i,
            t,
            current_uv0,
            current_uv1,
        )?;
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

// Adaptive walk inner loop lives in `walker_adaptive.rs` to keep
// `walker.rs` under the 500-line structural cap.

/// One contact solve at spine parameter `t`, warm-starting from
/// `(prev_uv0, prev_uv1)`. Shared by both fixed and adaptive paths.
/// Applies the periodicity-aware iteration window for samples after
/// the first.
pub(crate) fn step_at(
    request: &WalkRequest,
    tolerance: &Tolerance,
    contact_options: &ContactOptions,
    sample_index: usize,
    t: f64,
    prev_uv0: (f64, f64),
    prev_uv1: (f64, f64),
) -> Result<ContactSolution, WalkerError> {
    let p_spine = request.spine.point_at(t);
    let t_spine_tangent = request.spine.tangent_at(t);
    let radius = (request.radius_at)(t);

    let iter_uv0_range = if sample_index == 0 {
        request.uv0_range
    } else {
        recenter_range(request.uv0_range, prev_uv0, request.surface_0_periodicity)
    };
    let iter_uv1_range = if sample_index == 0 {
        request.uv1_range
    } else {
        recenter_range(request.uv1_range, prev_uv1, request.surface_1_periodicity)
    };

    let contact_request = ContactRequest {
        surface_0: request.surface_0,
        surface_1: request.surface_1,
        side_0: request.side_0,
        side_1: request.side_1,
        radius,
        initial_uv0: prev_uv0,
        uv0_range: iter_uv0_range,
        initial_uv1: prev_uv1,
        uv1_range: iter_uv1_range,
        spine_point: p_spine,
        spine_plane_normal: t_spine_tangent,
    };
    solve_rolling_ball_contact(&contact_request, contact_options, tolerance, None).map_err(
        |error| {
            if sample_index == 0 {
                WalkerError::InitialContactFailed(error)
            } else {
                WalkerError::StepFailed {
                    t_spine: t,
                    sample_index,
                    error,
                }
            }
        },
    )
}

/// Recenter the contact solver's UV range around `current_uv` for
/// periodic dimensions, leaving non-periodic dimensions untouched.
///
/// For a periodic dimension with period `T`, the new bounds are
/// `(current_u - 2·T, current_u + 2·T)` — four periods wide,
/// centered on the previous solution. Newton's natural step from a
/// good warm start is at most one period (typically a small fraction
/// of one), so the bounds are unreachable in normal use.
///
/// The 4×period factor (rather than the obvious 2×period centered
/// window) is deliberate: with a 2-period-wide window, Newton can
/// converge to a periodically-equivalent point that lands *exactly*
/// on one of the new boundaries (the same 3D point as `current_u`
/// reached by adding/subtracting one period). The post-convergence
/// `uv_on_boundary` check in the contact solver then fires and
/// reports a spurious `OutOfBounds`. Widening to 4×period puts the
/// boundaries 2 periods away from `current_u`, which is far beyond
/// any plausible step.
///
/// Non-periodic dimensions inherit `base_range` unchanged so the
/// walker still respects real face bounds (e.g. the cone's axial v
/// parameter).
fn recenter_range(
    base_range: ((f64, f64), (f64, f64)),
    current_uv: (f64, f64),
    periodicity: (Option<f64>, Option<f64>),
) -> ((f64, f64), (f64, f64)) {
    let (period_u, period_v) = periodicity;
    let u_range = match period_u {
        Some(p) => (current_uv.0 - 2.0 * p, current_uv.0 + 2.0 * p),
        None => base_range.0,
    };
    let v_range = match period_v {
        Some(p) => (current_uv.1 - 2.0 * p, current_uv.1 + 2.0 * p),
        None => base_range.1,
    };
    (u_range, v_range)
}

#[cfg(test)]
#[path = "walker_tests.rs"]
mod tests;
