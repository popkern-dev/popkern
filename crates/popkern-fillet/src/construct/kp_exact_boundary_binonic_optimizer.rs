//! Direct exact-boundary K-P binonic spoke-row optimizer.

use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};

pub(super) struct ExactBoundaryBinonicResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    pub final_residual: f64,
}

pub(super) fn optimize_exact_boundary_binonic(
    initial: &[RectangularGregoryBiNonic],
) -> Option<ExactBoundaryBinonicResult> {
    if initial.is_empty() || any_self_intersects(initial) {
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
    if any_self_intersects(&binonic) {
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
