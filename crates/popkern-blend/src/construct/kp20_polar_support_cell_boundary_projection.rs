//! Project a 3D point onto a `Kp20Bi6Patch` boundary edge.
//!
//! Building block for the cap-vertex-aware support-cell split (Stage A2
//! of `docs/decisions/2026-06-21-kp20-cap-vertex-aware-support-cell-
//! split.md`). The macro piece is bi-6 (`Kp20Bi6Patch`, 7×7 grid), so
//! each boundary edge is a degree-6 Bernstein curve. Newton iteration
//! on squared distance, warm-started by a 33-sample chord search,
//! returns the parameter `t ∈ [0, 1]` of the closest on-curve point —
//! or signals that the target lies off the boundary or that Newton
//! could not converge inside the iteration cap.

use popkern_core::geometry::gregory::Kp20Bi6Patch;
use popkern_core::geometry::{Point3, Vec3};

/// Which of the four boundary edges of a macro piece (`Kp20Bi6Patch`)
/// a projection is taken against. UMin / UMax index `u` to 0 / 1 and
/// trace the boundary in `v`; VMin / VMax trace in `u`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SupportCellPatchSide {
    UMin,
    UMax,
    VMin,
    VMax,
}

impl SupportCellPatchSide {
    pub(crate) fn all() -> [Self; 4] {
        [Self::UMin, Self::UMax, Self::VMin, Self::VMax]
    }
}

/// Outcome of projecting a target onto a macro-piece boundary edge.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)] // first caller lands in Stage A3
pub(crate) enum ProjectionOutcome {
    /// Target lies on the boundary curve within `ON_BOUNDARY_TOL`;
    /// `t` is the parameter of the closest on-curve point and is in
    /// the open interval `(DEDUP_TOL, 1 - DEDUP_TOL)`.
    Converged(f64),
    /// Newton converged but either the residual exceeds
    /// `ON_BOUNDARY_TOL` (target is off this boundary) or `t` collapsed
    /// onto an endpoint (corner-coincident; would double-count across
    /// adjacent sides). `closest_boundary` is the converged `t`, clamped
    /// to `[0, 1]`.
    OutsideDomain { closest_boundary: f64 },
    /// Newton did not converge within `MAX_NEWTON_ITERATIONS`; the
    /// warm-start estimate is returned for caller fallback handling.
    DidNotConverge { warm_start: f64 },
}

const WARM_START_SAMPLES: usize = 33;
const MAX_NEWTON_ITERATIONS: usize = 16;

/// Parametric convergence criterion. A Newton step of less than this
/// is treated as converged. Order-of-magnitude tighter than the f64
/// epsilon on the unit-length parameter interval.
const PARAMETER_CONVERGENCE_TOL: f64 = 1.0e-12;

/// Residual convergence criterion. Squared distance below this is
/// treated as on-curve (≈ 3.2e-8 absolute distance on unit-extent
/// geometry, well below `ON_BOUNDARY_TOL`).
const RESIDUAL_CONVERGENCE_TOL: f64 = 1.0e-15;

/// Maximum on-curve residual (absolute distance) for the projection
/// to declare `Converged`. Targets farther than this from the boundary
/// curve return `OutsideDomain`. 1e-6 is sub-millimeter for the
/// unit-extent macro pieces in the fillet routes, well below the
/// 5e-3 user gate at `kp_gate_helpers::G0_GATE`.
const ON_BOUNDARY_TOL: f64 = 1.0e-6;

/// Endpoint deadband matching the `cap_vertex_splits_along_*` first pass.
/// A converged `t` within this distance of 0 or 1 is corner-coincident
/// and routes to `OutsideDomain` to avoid double-contribution across
/// adjacent sides.
const DEDUP_TOL: f64 = 1.0e-9;

/// Projects `target` onto the boundary edge of `base` identified by
/// `side`. See [`ProjectionOutcome`] for outcome semantics.
#[allow(dead_code)] // first caller lands in Stage A3
pub(crate) fn project_to_boundary(
    base: &Kp20Bi6Patch,
    side: SupportCellPatchSide,
    target: &Point3,
) -> ProjectionOutcome {
    let cps = boundary_cps(base, side);
    let warm_start = warm_start_parameter(&cps, target);
    let (converged_t, residual_sq, converged) = newton_iterate(&cps, target, warm_start);
    if !converged {
        return ProjectionOutcome::DidNotConverge { warm_start };
    }
    if residual_sq.sqrt() > ON_BOUNDARY_TOL {
        return ProjectionOutcome::OutsideDomain {
            closest_boundary: converged_t,
        };
    }
    if converged_t < DEDUP_TOL {
        return ProjectionOutcome::OutsideDomain {
            closest_boundary: 0.0,
        };
    }
    if converged_t > 1.0 - DEDUP_TOL {
        return ProjectionOutcome::OutsideDomain {
            closest_boundary: 1.0,
        };
    }
    ProjectionOutcome::Converged(converged_t)
}

fn boundary_cps(base: &Kp20Bi6Patch, side: SupportCellPatchSide) -> [Point3; 7] {
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 7];
    match side {
        SupportCellPatchSide::UMin => {
            for (k, slot) in cps.iter_mut().enumerate() {
                *slot = base.cps[0][k];
            }
        }
        SupportCellPatchSide::UMax => {
            for (k, slot) in cps.iter_mut().enumerate() {
                *slot = base.cps[6][k];
            }
        }
        SupportCellPatchSide::VMin => {
            for (k, slot) in cps.iter_mut().enumerate() {
                *slot = base.cps[k][0];
            }
        }
        SupportCellPatchSide::VMax => {
            for (k, slot) in cps.iter_mut().enumerate() {
                *slot = base.cps[k][6];
            }
        }
    }
    cps
}

fn eval_boundary_at(cps: &[Point3; 7], t: f64) -> Point3 {
    let mut work = *cps;
    for level in 1..=6 {
        for index in 0..=(6 - level) {
            work[index] = lerp(work[index], work[index + 1], t);
        }
    }
    work[0]
}

fn eval_boundary_derivative_at(cps: &[Point3; 7], t: f64) -> Vec3 {
    let mut work = [Vec3::new(0.0, 0.0, 0.0); 6];
    for (k, slot) in work.iter_mut().enumerate() {
        *slot = (cps[k + 1] - cps[k]) * 6.0;
    }
    for level in 1..=5 {
        for index in 0..=(5 - level) {
            work[index] = work[index] + (work[index + 1] - work[index]) * t;
        }
    }
    work[0]
}

fn warm_start_parameter(cps: &[Point3; 7], target: &Point3) -> f64 {
    let mut best_t = 0.0;
    let mut best_dist_sq = f64::INFINITY;
    for index in 0..WARM_START_SAMPLES {
        let t = index as f64 / (WARM_START_SAMPLES - 1) as f64;
        let point = eval_boundary_at(cps, t);
        let dist_sq = (point - *target).length_squared();
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            best_t = t;
        }
    }
    best_t
}

fn newton_iterate(cps: &[Point3; 7], target: &Point3, warm_start: f64) -> (f64, f64, bool) {
    let mut t = warm_start;
    let mut residual_sq = f64::INFINITY;
    for _ in 0..MAX_NEWTON_ITERATIONS {
        let point = eval_boundary_at(cps, t);
        let derivative = eval_boundary_derivative_at(cps, t);
        let diff = point - *target;
        residual_sq = diff.length_squared();
        if residual_sq < RESIDUAL_CONVERGENCE_TOL {
            return (t, residual_sq, true);
        }
        let numerator = diff.dot(&derivative);
        let denominator = derivative.length_squared();
        if denominator < 1.0e-30 {
            return (t.clamp(0.0, 1.0), residual_sq, false);
        }
        let delta = -numerator / denominator;
        let new_t = (t + delta).clamp(0.0, 1.0);
        if (new_t - t).abs() < PARAMETER_CONVERGENCE_TOL {
            let final_point = eval_boundary_at(cps, new_t);
            let final_residual_sq = (final_point - *target).length_squared();
            return (new_t, final_residual_sq, true);
        }
        t = new_t;
    }
    (t, residual_sq, false)
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A degree-6 patch whose boundary curves are linear in 3D: each
    /// boundary CP grid row/column is a uniform linear interpolation
    /// between the two corners on that side. This makes the closest-
    /// point parameter equal to the linear projection of the target
    /// onto the side's corner-to-corner segment — an analytic answer.
    fn linear_boundary_patch() -> Kp20Bi6Patch {
        let mut patch = Kp20Bi6Patch::zero();
        for row in 0..7 {
            for column in 0..7 {
                patch.cps[row][column] = Point3::new(row as f64 / 6.0, column as f64 / 6.0, 0.0);
            }
        }
        patch
    }

    fn assert_converged_near(outcome: ProjectionOutcome, expected: f64) {
        match outcome {
            ProjectionOutcome::Converged(t) => {
                assert!(
                    (t - expected).abs() < 1.0e-12,
                    "converged t={t} expected≈{expected}"
                );
            }
            other => panic!("expected Converged near {expected}, got {other:?}"),
        }
    }

    #[test]
    fn projection_converges_to_analytic_parameter_on_each_side() {
        let patch = linear_boundary_patch();
        // UMin (u=0): boundary varies in v from (0,0,0) to (0,1,0).
        // Target (0, 0.6, 0) projects to t = 0.6.
        let target_umin = Point3::new(0.0, 0.6, 0.0);
        assert_converged_near(
            project_to_boundary(&patch, SupportCellPatchSide::UMin, &target_umin),
            0.6,
        );
        // UMax (u=1): boundary varies in v from (1,0,0) to (1,1,0).
        let target_umax = Point3::new(1.0, 0.4, 0.0);
        assert_converged_near(
            project_to_boundary(&patch, SupportCellPatchSide::UMax, &target_umax),
            0.4,
        );
        // VMin (v=0): boundary varies in u from (0,0,0) to (1,0,0).
        let target_vmin = Point3::new(0.7, 0.0, 0.0);
        assert_converged_near(
            project_to_boundary(&patch, SupportCellPatchSide::VMin, &target_vmin),
            0.7,
        );
        // VMax (v=1): boundary varies in u from (0,1,0) to (1,1,0).
        let target_vmax = Point3::new(0.25, 1.0, 0.0);
        assert_converged_near(
            project_to_boundary(&patch, SupportCellPatchSide::VMax, &target_vmax),
            0.25,
        );
    }

    #[test]
    fn off_boundary_target_returns_outside_domain() {
        let patch = linear_boundary_patch();
        // Target 0.1 in z, far off the UMin (planar z=0) boundary.
        let target = Point3::new(0.0, 0.5, 0.1);
        let outcome = project_to_boundary(&patch, SupportCellPatchSide::UMin, &target);
        match outcome {
            ProjectionOutcome::OutsideDomain { closest_boundary } => {
                assert!((closest_boundary - 0.5).abs() < 1.0e-12);
            }
            other => panic!("expected OutsideDomain, got {other:?}"),
        }
    }

    #[test]
    fn corner_coincident_target_returns_outside_domain() {
        let patch = linear_boundary_patch();
        // Corner at u=0, v=0 lives on both UMin and VMin sides.
        // Projecting it onto UMin would converge to t≈0, which the
        // helper reroutes to OutsideDomain to prevent double-counting.
        let corner = Point3::new(0.0, 0.0, 0.0);
        let outcome = project_to_boundary(&patch, SupportCellPatchSide::UMin, &corner);
        match outcome {
            ProjectionOutcome::OutsideDomain { closest_boundary } => {
                assert_eq!(closest_boundary, 0.0);
            }
            other => panic!("expected OutsideDomain at endpoint, got {other:?}"),
        }
    }

    /// A patch whose UMin boundary collapses to a single point. The
    /// derivative is zero across the whole parameter range, so Newton
    /// cannot make progress. The helper must return `DidNotConverge`
    /// rather than divide-by-zero or loop until the iteration cap.
    #[test]
    fn degenerate_zero_derivative_boundary_returns_did_not_converge() {
        let mut patch = Kp20Bi6Patch::zero();
        // Every CP on UMin is the same point; eval_boundary_derivative
        // returns Vec3::ZERO at every t.
        for column in 0..7 {
            patch.cps[0][column] = Point3::new(0.0, 0.0, 0.0);
        }
        // Other rows pulled away so the warm start doesn't accidentally
        // satisfy the residual criterion on a non-zero CP.
        for row in 1..7 {
            for column in 0..7 {
                patch.cps[row][column] = Point3::new(row as f64, column as f64, 0.0);
            }
        }
        // Target away from the collapsed point.
        let target = Point3::new(0.5, 0.5, 0.0);
        let outcome = project_to_boundary(&patch, SupportCellPatchSide::UMin, &target);
        assert!(
            matches!(outcome, ProjectionOutcome::DidNotConverge { .. }),
            "expected DidNotConverge, got {outcome:?}"
        );
    }
}
