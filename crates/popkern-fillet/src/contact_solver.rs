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
    let clamped_uv0 = clamp_uv(request.initial_uv0, request.uv0_range);
    let clamped_uv1 = clamp_uv(request.initial_uv1, request.uv1_range);
    if !options.clamp_initial_uv {
        if clamped_uv0 != request.initial_uv0 {
            return Err(ContactError::OutOfBounds {
                surface_index: 0,
                uv: request.initial_uv0,
                clamped_uv: clamped_uv0,
            });
        }
        if clamped_uv1 != request.initial_uv1 {
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
        if residual_norm < effective_tol {
            if uv_on_boundary(uv0, request.uv0_range) {
                return Err(ContactError::OutOfBounds {
                    surface_index: 0,
                    uv: uv0,
                    clamped_uv: uv0,
                });
            }
            if uv_on_boundary(uv1, request.uv1_range) {
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

/// True if `uv` is on the boundary of its range. Used to detect when
/// a converged solution actually rests on a UV bound.
fn uv_on_boundary(uv: (f64, f64), range: ((f64, f64), (f64, f64))) -> bool {
    let ((u_min, u_max), (v_min, v_max)) = range;
    let u_span = (u_max - u_min).abs().max(1e-12);
    let v_span = (v_max - v_min).abs().max(1e-12);
    let tol_u = u_span * 1e-9;
    let tol_v = v_span * 1e-9;
    (uv.0 - u_min).abs() < tol_u
        || (uv.0 - u_max).abs() < tol_u
        || (uv.1 - v_min).abs() < tol_v
        || (uv.1 - v_max).abs() < tol_v
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
