//! Direct exact-boundary K-P binonic spoke-row optimizer.

use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use popkern_core::geometry::gregory::{
    GregoryPatch, RectangularGregoryBiNonic, binonic_ring_crossing,
    binonic_ring_feasible_start_crossing,
};

/// Linear tolerance of the PK-283 ring-level crossing gates — the same
/// value every per-patch fold gate on the promotion ladder uses, so the
/// two gate families measure the same geometry class at the same scale.
pub(in crate::construct) const RING_CROSSING_TOLERANCE: f64 = 1.0e-6;

pub(super) struct ExactBoundaryBinonicResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    pub final_residual: f64,
}

pub(super) fn optimize_exact_boundary_binonic(
    initial: &[RectangularGregoryBiNonic],
) -> Option<ExactBoundaryBinonicResult> {
    if initial.is_empty() || !feasible_start(initial) {
        return None;
    }
    let baseline = measure_ring_residual(initial);
    let x0 = pack_spoke_rows(initial);
    let result = run_lm(
        DirectSpokeRowEval {
            base: initial.to_vec(),
        },
        x0,
        &LmConfig {
            fd_eps: 1e-4,
            lambda_init: 1e-3,
            lambda_max: 1e6,
            max_iter: 15,
            target_residual: 1e-6,
            require_feasible_start: true,
        },
        "exact-boundary-binonic-spoke-row",
    );
    if result.final_residual >= baseline {
        return None;
    }
    let mut binonic = initial.to_vec();
    unpack_spoke_rows(&result.x, &mut binonic);
    if !acceptable_ring(&binonic) {
        return None;
    }
    Some(ExactBoundaryBinonicResult {
        binonic,
        final_residual: result.final_residual,
    })
}

struct DirectSpokeRowEval {
    base: Vec<RectangularGregoryBiNonic>,
}

impl DirectSpokeRowEval {
    fn build_ring(&self, x: &[f64]) -> Vec<RectangularGregoryBiNonic> {
        let mut ring = self.base.clone();
        unpack_spoke_rows(x, &mut ring);
        ring
    }
}

impl LmEval for DirectSpokeRowEval {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        Some(residual_vector(&self.build_ring(x)))
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        !any_self_intersects(&self.build_ring(x))
    }

    // PK-275 measured decision: this eval deliberately KEEPS the FD
    // Jacobian. An analytic selector-chain Jacobian was implemented
    // and FD-agreement-verified here, but (a) it bought no wall time
    // — this optimizer's cost is the fold-safety Gregory
    // self-intersection check, not derivatives (gdb profile,
    // PK-282) — and (b) its more-accurate trajectory converged
    // pyramid4_uniform to a lower-residual ring whose ADJACENT
    // patches intersect each other, a class the per-patch fold gates
    // could not see. PK-283 closed (b): the acceptance is now
    // crossing-honest at ring level (`binonic_ring_crossing` gates
    // entry and accept above), so any trajectory landing in that
    // class is declined rather than served. (a) still holds, so FD
    // remains the keep-better choice on measurement; a future
    // analytic Jacobian is safe to reconsider purely on wall-time
    // evidence. The reverted implementation was never committed
    // (verified against 77894f2a and dangling objects, 2026-08-16).
}

fn pack_spoke_rows(patches: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let mut x = Vec::with_capacity(patches.len() * 16 * 3);
    for patch in patches {
        for j in 1..=8 {
            push_point(&mut x, patch.grid[8][j]);
        }
        for i in 1..=8 {
            push_point(&mut x, patch.grid[i][1]);
        }
    }
    x
}

fn unpack_spoke_rows(x: &[f64], patches: &mut [RectangularGregoryBiNonic]) {
    let mut idx = 0;
    for patch in patches {
        for j in 1..=8 {
            patch.grid[8][j].x = x[idx];
            patch.grid[8][j].y = x[idx + 1];
            patch.grid[8][j].z = x[idx + 2];
            idx += 3;
        }
        for i in 1..=8 {
            patch.grid[i][1].x = x[idx];
            patch.grid[i][1].y = x[idx + 1];
            patch.grid[i][1].z = x[idx + 2];
            idx += 3;
        }
    }
}

fn push_point(x: &mut Vec<f64>, point: popkern_core::geometry::Point3) {
    x.push(point.x);
    x.push(point.y);
    x.push(point.z);
}

/// Entry feasibility: every patch individually fold-free AND no
/// macro-scale pairwise crossing (PK-283). The coarse feasible-start
/// screen — not the full accept semantics — because raw guide rings
/// legitimately carry spoke-hover contact an order deeper than
/// converged rings (measured: pyramid3_uniform's initial ring at
/// 1.886e-3, polished to 2.38e-5 by this optimizer and served clean).
fn feasible_start(binonic: &[RectangularGregoryBiNonic]) -> bool {
    if any_self_intersects(binonic) {
        return false;
    }
    if let Some(crossing) = binonic_ring_feasible_start_crossing(binonic, RING_CROSSING_TOLERANCE) {
        if popkern_core::trace::trace_enabled("kp-cap") {
            eprintln!(
                "[exact-boundary-binonic] feasible-start ring decline: {}",
                crossing.diagnostic_label()
            );
        }
        return false;
    }
    true
}

/// Accept feasibility: per-patch fold gate plus the full ring-level
/// crossing semantics (PK-283 — a class the per-patch fold gates and
/// the per-step `fold_safe` cannot see). Run once per acceptance, not
/// per LM step.
fn acceptable_ring(binonic: &[RectangularGregoryBiNonic]) -> bool {
    !any_self_intersects(binonic) && ring_crossing_decline_label(binonic).is_none()
}

/// PK-283 accept-gate shared arm: the full ring-crossing scan, with the
/// codebase-standard the `kp-cap` trace topic decline trace. Returns the
/// crossing's diagnostic label so optimizer callers can record it as a
/// structured rung outcome.
pub(super) fn ring_crossing_decline_label(binonic: &[RectangularGregoryBiNonic]) -> Option<String> {
    let crossing = binonic_ring_crossing(binonic, RING_CROSSING_TOLERANCE)?;
    let label = crossing.diagnostic_label();
    if popkern_core::trace::trace_enabled("kp-cap") {
        eprintln!("[kp-ring-gate] accept decline: {label}");
    }
    Some(label)
}

fn any_self_intersects(binonic: &[RectangularGregoryBiNonic]) -> bool {
    binonic.iter().any(|patch| {
        GregoryPatch::RectangularBiNonic(Box::new(patch.clone()))
            .self_intersects_with_tolerance(1.0e-6)
    })
}

fn residual_vector(patches: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let n = patches.len();
    let mut out = Vec::with_capacity(n * 50);
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let cross = xtan_a.cross(&xtan_b);
            let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
            let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
            let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
            let ns = (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
            let denom = na * nb * ns;
            out.push(if denom > 1.0e-300 { numer / denom } else { 0.0 });
        }
    }
    out
}

fn measure_ring_residual(patches: &[RectangularGregoryBiNonic]) -> f64 {
    residual_vector(patches).iter().map(|v| v * v).sum()
}

#[cfg(test)]
#[path = "kp_exact_boundary_ring_gate_tests.rs"]
mod kp_exact_boundary_ring_gate_tests;
