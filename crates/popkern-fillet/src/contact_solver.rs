//! Damped Newton solver for the rolling-ball contact problem
//! ([`solve_rolling_ball_contact`]) plus its small private helpers
//! (line search, bounds clamp, residual norm, midpoint).
//!
//! Split from `contact.rs` to keep individual files under the
//! 500-line structural cap. Re-exported through `crate::contact`.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

use crate::contact::math::{
    build_jacobian_and_residual, condition_estimate_4x4, lm_step, signed_unit_normal, solve_4x4,
};
use crate::contact::{
    ContactError, ContactOptions, ContactQuality, ContactRequest, ContactSolution, ContactTrace,
    ContactTraceStep, SingularityKind,
};

/// Solve the rolling-ball contact problem at one spine point.
///
/// Returns a [`ContactSolution`] on success; see [`ContactError`] for
/// failure modes. Pass `Some(&mut trace)` to capture per-iteration
/// state for debugging.
///
/// The algorithm is damped Newton iteration on the 4×4 system
/// `(F(x), g(x)) = 0` where `F` is the contact-equality residual and
/// `g` is the spine-plane scalar. The Jacobian is computed analytically
/// via the Weingarten map. When the Jacobian condition exceeds 1e8 the
/// solver switches to Levenberg-Marquardt regularization with damping
/// λ growing 10× per failed step, capped at 1.0.
pub fn solve_rolling_ball_contact(
    request: &ContactRequest,
    options: &ContactOptions,
    tolerance: &Tolerance,
    mut trace: Option<&mut ContactTrace>,
) -> Result<ContactSolution, ContactError> {
    let floor = tolerance.numerical_floor();
    let residual_tol = options.residual_tolerance.max(floor);

    // Per-quality preset overrides on top of `options`.
    let (effective_max_iter, effective_tol) = match options.quality {
        ContactQuality::Final => (options.max_iterations, residual_tol),
        ContactQuality::Preview => (3.min(options.max_iterations), residual_tol * 10.0),
    };

    // Validate or clamp initial UV guesses. Default: reject
    // out-of-bounds (existing behaviour). With `clamp_initial_uv`
    // (§12.3 blend-on-blend), silently clamp to the boundary.
    //
    // PERIODIC charts wrap FIRST: a contact whose azimuth lands at
    // u = −ε on a cylinder chart spanning [0, 2π] is on the surface —
    // rejecting it as out-of-bounds killed every fillet whose seam
    // start sat on the chart seam (measured: flagship boss/plate
    // circular seam, "initial contact failed" at u = −0.0000, the
    // failure class behind PK-099 stage 3 and the round_5c fixtures).
    let period_0 = request.surface_0.parameter_periodicity();
    let period_1 = request.surface_1.parameter_periodicity();
    let wrapped_uv0 = wrap_uv_periodic(request.initial_uv0, request.uv0_range, period_0);
    let wrapped_uv1 = wrap_uv_periodic(request.initial_uv1, request.uv1_range, period_1);
    let clamped_uv0 = clamp_uv(wrapped_uv0, request.uv0_range);
    let clamped_uv1 = clamp_uv(wrapped_uv1, request.uv1_range);
    if !options.clamp_initial_uv {
        if clamped_uv0 != wrapped_uv0 {
            return Err(ContactError::OutOfBounds {
                surface_index: 0,
                uv: request.initial_uv0,
                clamped_uv: clamped_uv0,
            });
        }
        if clamped_uv1 != wrapped_uv1 {
            return Err(ContactError::OutOfBounds {
                surface_index: 1,
                uv: request.initial_uv1,
                clamped_uv: clamped_uv1,
            });
        }
    }
    let mut uv0 = clamped_uv0;
    let mut uv1 = clamped_uv1;

    // Track best-residual state for diagnostics on failure.
    let mut best_residual = f64::INFINITY;
    let mut best_uv0 = uv0;
    let mut best_uv1 = uv1;
    let mut diverge_streak = 0_usize;
    const DIVERGE_STREAK_LIMIT: usize = 4;

    // §18 telemetry: total line-search rejections across all
    // Newton iterations of this solve. Bubbles up into the
    // returned `ContactSolution`.
    let mut total_step_rejections = 0_usize;

    for iteration in 0..effective_max_iter {
        let built = build_jacobian_and_residual(
            request.surface_0,
            request.surface_1,
            request.side_0,
            request.side_1,
            request.radius,
            uv0,
            uv1,
            request.spine_point,
            request.spine_plane_normal,
            floor,
        );
        let (jacobian, residual_vec, c0, c1, n0, n1) = match built {
            Some(t) => t,
            None => {
                let p0 = request.surface_0.point_at(uv0.0, uv0.1);
                let p1 = request.surface_1.point_at(uv1.0, uv1.1);
                let raw_n0 = request.surface_0.normal_at(uv0.0, uv0.1);
                let (surface_index, uv_bad, location) = if raw_n0.length() < floor {
                    (0_usize, uv0, p0)
                } else {
                    (1_usize, uv1, p1)
                };
                return Err(ContactError::DegenerateNormal {
                    surface_index,
                    uv: uv_bad,
                    location,
                });
            }
        };

        let residual_norm = (residual_vec[0].powi(2)
            + residual_vec[1].powi(2)
            + residual_vec[2].powi(2)
            + residual_vec[3].powi(2))
        .sqrt();

        let condition = condition_estimate_4x4(&jacobian);
        if std::env::var_os("PK_CONTACT_TRACE").is_some() {
            eprintln!(
                "    [contact] iter {iteration}: res {residual_norm:.3e} uv0 ({:.4},{:.4}) uv1 ({:.4},{:.4})",
                uv0.0, uv0.1, uv1.0, uv1.1
            );
        }
        if let Some(t) = trace.as_deref_mut() {
            t.steps.push(ContactTraceStep {
                iteration,
                residual: residual_norm,
                step_length: 0.0,
                condition,
                damping: 0.0,
                uv0,
                uv1,
            });
        }

        // Convergence check. If we converged onto a UV bound, surface
        // it as OutOfBounds — the unconstrained Newton solution lies
        // outside the face, which means the radius is too large here.
        // A full-period axis of a periodic chart has no boundary: a
        // solution resting exactly on the seam iso-line is interior.
        if residual_norm < effective_tol {
            if uv_on_boundary_periodic(uv0, request.uv0_range, period_0) {
                return Err(ContactError::OutOfBounds {
                    surface_index: 0,
                    uv: uv0,
                    clamped_uv: uv0,
                });
            }
            if uv_on_boundary_periodic(uv1, request.uv1_range, period_1) {
                return Err(ContactError::OutOfBounds {
                    surface_index: 1,
                    uv: uv1,
                    clamped_uv: uv1,
                });
            }
            return Ok(ContactSolution {
                uv0,
                uv1,
                center: midpoint(c0, c1),
                contact0: request.surface_0.point_at(uv0.0, uv0.1),
                contact1: request.surface_1.point_at(uv1.0, uv1.1),
                iterations: iteration,
                step_rejections: total_step_rejections,
                residual: residual_norm,
                condition,
                jacobian,
            });
        }

        if residual_norm < best_residual {
            best_residual = residual_norm;
            best_uv0 = uv0;
            best_uv1 = uv1;
            diverge_streak = 0;
        } else {
            diverge_streak += 1;
            if diverge_streak >= DIVERGE_STREAK_LIMIT {
                return Err(ContactError::Diverged {
                    iterations: iteration + 1,
                    last_residual: residual_norm,
                    last_uv0: uv0,
                    last_uv1: uv1,
                });
            }
        }

        // Detect parallel normals before the linear solve — much more
        // useful diagnostic than "singular Jacobian." n0 and n1 are
        // unit length (signed_unit_normal normalises), so |n0×n1| =
        // |sin(angle between normals)|.
        let cross_len = n0.cross(&n1).length();
        if cross_len < tolerance.angular() {
            return Err(ContactError::SingularJacobian {
                iterations: iteration,
                condition,
                location: midpoint(c0, c1),
                kind: SingularityKind::NormalsParallel,
            });
        }

        // Solve J·Δx = -residual; switch to Levenberg-Marquardt if
        // the Jacobian is ill-conditioned.
        let neg_residual = [
            -residual_vec[0],
            -residual_vec[1],
            -residual_vec[2],
            -residual_vec[3],
        ];
        const LM_TRIGGER_CONDITION: f64 = 1.0e8;
        let used_damping;
        let step_result = if condition < LM_TRIGGER_CONDITION {
            used_damping = 0.0;
            solve_4x4(jacobian, neg_residual, floor)
        } else {
            let initial = options.initial_damping.max(1.0e-6);
            used_damping = initial;
            lm_step(&jacobian, &residual_vec, initial, floor)
        };
        if let Some(t) = trace.as_deref_mut()
            && let Some(last) = t.steps.last_mut()
        {
            last.damping = used_damping;
        }
        let step = match step_result {
            Some(s) => s,
            None => {
                return Err(ContactError::SingularJacobian {
                    iterations: iteration,
                    condition,
                    location: midpoint(c0, c1),
                    kind: SingularityKind::JacobianRankDeficient,
                });
            }
        };

        let (new_uv0, new_uv1, _accepted_residual, rejections_this_step) =
            line_search(request, uv0, uv1, &step, residual_norm, floor);
        total_step_rejections += rejections_this_step;
        let step_length = ((new_uv0.0 - uv0.0).powi(2)
            + (new_uv0.1 - uv0.1).powi(2)
            + (new_uv1.0 - uv1.0).powi(2)
            + (new_uv1.1 - uv1.1).powi(2))
        .sqrt();
        if let Some(t) = trace.as_deref_mut()
            && let Some(last) = t.steps.last_mut()
        {
            last.step_length = step_length;
        }

        uv0 = new_uv0;
        uv1 = new_uv1;
    }

    Err(ContactError::NotConverged {
        iterations: effective_max_iter,
        residual: best_residual,
        last_uv0: best_uv0,
        last_uv1: best_uv1,
    })
}

/// Backtracking line search over the Newton step. Each candidate
/// state is clamped to the request's `uv*_range` boxes before
/// evaluation. Returns `(new_uv0, new_uv1, residual_at_new_state,
/// rejections)`, where `rejections` counts how many candidate `α`
/// values were discarded before a usable step was found — §18
/// step-rejection telemetry. A healthy well-conditioned contact
/// converges with 0 rejections; a count near 6 means all attempts
/// failed to improve residual and the line search fell through
/// without progress.
fn line_search(
    request: &ContactRequest,
    uv0: (f64, f64),
    uv1: (f64, f64),
    step: &[f64; 4],
    current_residual: f64,
    floor: f64,
) -> ((f64, f64), (f64, f64), f64, usize) {
    let mut alpha = 1.0_f64;
    let mut accepted_uv0 = uv0;
    let mut accepted_uv1 = uv1;
    let mut accepted_res = current_residual;
    let mut rejections = 0_usize;

    for _ in 0..6 {
        let raw_uv0 = (uv0.0 + alpha * step[0], uv0.1 + alpha * step[1]);
        let raw_uv1 = (uv1.0 + alpha * step[2], uv1.1 + alpha * step[3]);
        let try_uv0 = clamp_uv(raw_uv0, request.uv0_range);
        let try_uv1 = clamp_uv(raw_uv1, request.uv1_range);
        let try_residual = residual_norm_at(request, try_uv0, try_uv1, floor);
        if let Some(r) = try_residual
            && (r < current_residual || alpha < 1.0 / 32.0)
        {
            accepted_uv0 = try_uv0;
            accepted_uv1 = try_uv1;
            accepted_res = r;
            return (accepted_uv0, accepted_uv1, accepted_res, rejections);
        }
        rejections += 1;
        alpha *= 0.5;
    }
    (accepted_uv0, accepted_uv1, accepted_res, rejections)
}

/// Clamp `(u, v)` into the given `((u_min, u_max), (v_min, v_max))` box.
pub(crate) fn clamp_uv(uv: (f64, f64), range: ((f64, f64), (f64, f64))) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = range;
    (uv.0.clamp(u_min, u_max), uv.1.clamp(v_min, v_max))
}

/// Wrap each periodic axis of `(u, v)` into its range when a whole
/// number of periods brings it inside. Values on a non-periodic axis,
/// or genuinely outside a PARTIAL periodic window (no period multiple
/// lands them inside), are returned unchanged for the caller's clamp
/// or rejection to handle.
fn wrap_uv_periodic(
    uv: (f64, f64),
    range: ((f64, f64), (f64, f64)),
    periodicity: (Option<f64>, Option<f64>),
) -> (f64, f64) {
    let wrap = |value: f64, (lo, hi): (f64, f64), period: Option<f64>| -> f64 {
        let Some(p) = period else { return value };
        if p <= 0.0 || (lo..=hi).contains(&value) {
            return value;
        }
        // Land in [lo, lo + p); accept if that is inside the window
        // (with a relative slop for values that sit exactly on the
        // seam iso-line, e.g. u = -1e-17 on a [0, 2π] chart).
        let wrapped = value - p * ((value - lo) / p).floor();
        let slop = (hi - lo).abs().max(p) * 1.0e-9;
        if wrapped <= hi + slop {
            wrapped.clamp(lo, hi)
        } else {
            value
        }
    };
    (
        wrap(uv.0, range.0, periodicity.0),
        wrap(uv.1, range.1, periodicity.1),
    )
}

/// True if `uv` is on the boundary of its range. Used to detect when
/// a converged solution actually rests on a UV bound. An axis whose
/// window spans a full period of a periodic chart has NO boundary —
/// the seam iso-line is interior by wraparound.
fn uv_on_boundary_periodic(
    uv: (f64, f64),
    range: ((f64, f64), (f64, f64)),
    periodicity: (Option<f64>, Option<f64>),
) -> bool {
    let ((u_min, u_max), (v_min, v_max)) = range;
    let u_span = (u_max - u_min).abs().max(1e-12);
    let v_span = (v_max - v_min).abs().max(1e-12);
    let tol_u = u_span * 1e-9;
    let tol_v = v_span * 1e-9;
    let axis_closed = |span: f64, period: Option<f64>| {
        period.is_some_and(|p| p > 0.0 && span >= p * (1.0 - 1.0e-9))
    };
    let u_boundary = !axis_closed(u_span, periodicity.0)
        && ((uv.0 - u_min).abs() < tol_u || (uv.0 - u_max).abs() < tol_u);
    let v_boundary = !axis_closed(v_span, periodicity.1)
        && ((uv.1 - v_min).abs() < tol_v || (uv.1 - v_max).abs() < tol_v);
    u_boundary || v_boundary
}

/// Evaluate the residual norm at a candidate state without building
/// the full Jacobian. Returns `None` if surface evaluation produces
/// a degenerate normal.
fn residual_norm_at(
    request: &ContactRequest,
    uv0: (f64, f64),
    uv1: (f64, f64),
    floor: f64,
) -> Option<f64> {
    let n0 = signed_unit_normal(request.surface_0, uv0.0, uv0.1, request.side_0, floor)?;
    let n1 = signed_unit_normal(request.surface_1, uv1.0, uv1.1, request.side_1, floor)?;
    let p0 = request.surface_0.point_at(uv0.0, uv0.1);
    let p1 = request.surface_1.point_at(uv1.0, uv1.1);
    let c0 = p0 + n0 * request.radius;
    let c1 = p1 + n1 * request.radius;
    let f = c0 - c1;
    let g = (c0 - request.spine_point).dot(&request.spine_plane_normal);
    Some((f.x * f.x + f.y * f.y + f.z * f.z + g * g).sqrt())
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}
