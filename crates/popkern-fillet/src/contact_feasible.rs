//! Feasibility wrappers around [`solve_rolling_ball_contact`]:
//! cheap pre-check for UI hover feedback ([`pre_check_contact`]) and
//! the binary-search wrapper that finds the largest feasible radius
//! on size-related failures ([`solve_with_max_feasible`]).
//!
//! Split from `contact.rs` to keep individual files under the 500-line
//! structural cap. Re-exported through `crate::contact`.

use popkern_core::tolerance::Tolerance;

use crate::contact::math::signed_unit_normal;
use crate::contact::solver::{clamp_uv, solve_rolling_ball_contact};
use crate::contact::{
    ContactError, ContactFeasibility, ContactOptions, ContactRequest, ContactSolution,
};

/// Cheap pre-check before attempting [`solve_rolling_ball_contact`].
///
/// Performs ~5 surface evaluations and no nonlinear solve, suitable
/// for UI hover-feedback at sub-millisecond cost. Returns a
/// [`ContactFeasibility`] classification with a reason string when
/// the inputs are unlikely to converge.
///
/// The check inspects the dihedral angle at the initial guess and the
/// signed local curvature against the requested radius. It does NOT
/// guarantee feasibility (the full Newton solve may still fail on
/// pathological geometry); it only catches the cheap-to-detect
/// hopeless cases.
pub fn pre_check_contact(request: &ContactRequest, tolerance: &Tolerance) -> ContactFeasibility {
    let floor = tolerance.numerical_floor();

    if request.radius <= 0.0 {
        return ContactFeasibility::LikelyInfeasible {
            reason: format!("radius must be > 0, got {}", request.radius),
            max_feasible_radius: None,
        };
    }
    if clamp_uv(request.initial_uv0, request.uv0_range) != request.initial_uv0 {
        return ContactFeasibility::LikelyInfeasible {
            reason: "initial uv0 is outside the supplied uv0_range".into(),
            max_feasible_radius: None,
        };
    }
    if clamp_uv(request.initial_uv1, request.uv1_range) != request.initial_uv1 {
        return ContactFeasibility::LikelyInfeasible {
            reason: "initial uv1 is outside the supplied uv1_range".into(),
            max_feasible_radius: None,
        };
    }

    let n0 = match signed_unit_normal(
        request.surface_0,
        request.initial_uv0.0,
        request.initial_uv0.1,
        request.side_0,
        floor,
    ) {
        Some(n) => n,
        None => {
            return ContactFeasibility::LikelyInfeasible {
                reason: "surface 0 has a degenerate normal at the initial guess".into(),
                max_feasible_radius: None,
            };
        }
    };
    let n1 = match signed_unit_normal(
        request.surface_1,
        request.initial_uv1.0,
        request.initial_uv1.1,
        request.side_1,
        floor,
    ) {
        Some(n) => n,
        None => {
            return ContactFeasibility::LikelyInfeasible {
                reason: "surface 1 has a degenerate normal at the initial guess".into(),
                max_feasible_radius: None,
            };
        }
    };

    let cos_dihedral = n0.dot(&n1).clamp(-1.0, 1.0);
    let cross_len = n0.cross(&n1).length();
    let angular_tol = tolerance.angular();

    if cross_len < angular_tol {
        if cos_dihedral > 0.0 {
            return ContactFeasibility::LikelyInfeasible {
                reason: "surface normals are nearly co-directed — no rolling-ball fillet exists between parallel faces".into(),
                max_feasible_radius: None,
            };
        }
        return ContactFeasibility::LikelyInfeasible {
            reason: "surface normals are nearly anti-directed — fillet geometry is degenerate"
                .into(),
            max_feasible_radius: None,
        };
    }

    // Local face-size sanity from the surface tangents at the initial
    // guess. If the requested radius is large compared to the local
    // surface scale (||S_u||, ||S_v||) the ball is unlikely to fit
    // within the supplied uv ranges.
    let s0_u = request
        .surface_0
        .partial_u(request.initial_uv0.0, request.initial_uv0.1);
    let s0_v = request
        .surface_0
        .partial_v(request.initial_uv0.0, request.initial_uv0.1);
    let s1_u = request
        .surface_1
        .partial_u(request.initial_uv1.0, request.initial_uv1.1);
    let s1_v = request
        .surface_1
        .partial_v(request.initial_uv1.0, request.initial_uv1.1);
    let scale_0 = s0_u.length().max(s0_v.length());
    let scale_1 = s1_u.length().max(s1_v.length());

    let span_u0 = (request.uv0_range.0.1 - request.uv0_range.0.0).abs();
    let span_v0 = (request.uv0_range.1.1 - request.uv0_range.1.0).abs();
    let span_u1 = (request.uv1_range.0.1 - request.uv1_range.0.0).abs();
    let span_v1 = (request.uv1_range.1.1 - request.uv1_range.1.0).abs();
    let face_extent_0 = scale_0 * span_u0.min(span_v0);
    let face_extent_1 = scale_1 * span_u1.min(span_v1);

    if face_extent_0.is_finite() && request.radius > face_extent_0 {
        return ContactFeasibility::LikelyTight {
            reason: format!(
                "radius {:.3e} is comparable to surface 0's local extent {:.3e}; \
                 expect a tight or boundary-rest solution",
                request.radius, face_extent_0
            ),
        };
    }
    if face_extent_1.is_finite() && request.radius > face_extent_1 {
        return ContactFeasibility::LikelyTight {
            reason: format!(
                "radius {:.3e} is comparable to surface 1's local extent {:.3e}; \
                 expect a tight or boundary-rest solution",
                request.radius, face_extent_1
            ),
        };
    }

    ContactFeasibility::LikelyFine
}

/// Solve the rolling-ball contact problem and, on size-related
/// failure, binary-search downward to find the largest feasible
/// radius from `request`'s starting configuration.
///
/// Returns `(solution, achieved_radius)` on success — `achieved_radius`
/// equals `request.radius` if the original radius worked, otherwise
/// the largest radius the binary search could prove feasible. Returns
/// the original error if the failure is not radius-related (parallel
/// normals, degenerate normals).
///
/// Used by:
/// - §10.4 `feasible_radius_bounds` (calls this with progressively
///   tighter radii to bracket `[r_min, r_max]`);
/// - §13.3 "offer max feasible size on failure" — the recommended
///   error-recovery path for a UI.
pub fn solve_with_max_feasible(
    request: &ContactRequest,
    options: &ContactOptions,
    tolerance: &Tolerance,
) -> Result<(ContactSolution, f64), ContactError> {
    let first = solve_rolling_ball_contact(request, options, tolerance, None);
    let initial_err = match first {
        Ok(solution) => return Ok((solution, request.radius)),
        Err(e) => e,
    };

    let radius_fixable = matches!(
        initial_err,
        ContactError::OutOfBounds { .. }
            | ContactError::NotConverged { .. }
            | ContactError::Diverged { .. }
    );
    if !radius_fixable {
        return Err(initial_err);
    }

    let r_max_failing = request.radius;
    let mut r_min_succeeding = 0.0_f64;
    let mut best_solution: Option<ContactSolution> = None;
    let mut hi = r_max_failing;
    let mut lo = tolerance.numerical_floor();

    // First pass: walk `lo` upward by doubling until we find a
    // working radius, capped at 30 attempts.
    for _ in 0..30 {
        let try_radius = lo;
        let try_request = ContactRequest {
            radius: try_radius,
            ..*request
        };
        if let Ok(sol) = solve_rolling_ball_contact(&try_request, options, tolerance, None) {
            r_min_succeeding = try_radius;
            best_solution = Some(sol);
            break;
        }
        lo *= 2.0;
        if lo >= hi {
            break;
        }
    }

    if best_solution.is_none() {
        return Err(initial_err);
    }

    // Bisect (r_min_succeeding, hi).
    for _ in 0..40 {
        let mid = 0.5 * (r_min_succeeding + hi);
        if (hi - r_min_succeeding) / mid.max(tolerance.linear()) < 1.0e-3 {
            break;
        }
        let mid_request = ContactRequest {
            radius: mid,
            ..*request
        };
        match solve_rolling_ball_contact(&mid_request, options, tolerance, None) {
            Ok(sol) => {
                r_min_succeeding = mid;
                best_solution = Some(sol);
            }
            Err(_) => {
                hi = mid;
            }
        }
    }

    // `best_solution` is guaranteed `Some` at this point: the early
    // return above (line 209-211) bails out if the doubling pass
    // never found a working radius, and the bisection loop only
    // ever assigns `Some(sol)` — it never resets to `None`. The
    // `ok_or(...)?` fallback is dead in practice; it exists only
    // so the type system sees a Result path instead of an expect.
    let solution = best_solution.ok_or(ContactError::NotConverged {
        iterations: 0,
        residual: f64::INFINITY,
        last_uv0: request.initial_uv0,
        last_uv1: request.initial_uv1,
    })?;
    Ok((solution, r_min_succeeding))
}
