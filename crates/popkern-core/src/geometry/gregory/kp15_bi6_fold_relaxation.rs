//! Fold relaxation for the KP15 shared-net cap over its two
//! boundary-neutral degrees of freedom.
//!
//! A boundary-faithful ring fit can fold the cap's deep interior
//! (PK-087 stage-6 measurement). Two families of perturbations leave
//! the fitted boundary rows **exactly** unchanged while re-placing
//! that interior:
//!
//! 1. the **kernel of the boundary design matrix** — ring shifts with
//!    `A·k = 0` (the boundary map has row-rank ~4n over 6n+1 ring
//!    unknowns, so this family is ~(2n+1)-dimensional per coordinate);
//! 2. the **paper free set** `P` (`p00`, tangent/quadratic central
//!    data, `p220`, `p210`, per-sector `p23`/`p33`) — free by the KP15
//!    construction, and the eq (2)/(3)/`sal` closure is independent of
//!    its values (oracle symbolic identity).
//!
//! This module searches both families with deterministic multistart
//! coordinate descent, scoring candidates by the production Jacobian
//! fold report and accepting only configurations that pass the full
//! production self-intersection detector (`self_intersects_with_
//! tolerance`, Jacobian tier + NURBS subdivision tier). On the
//! hexagon6_variable apex this turns the boundary-faithful fit from
//! folded (min normal cos −0.99) to detector-clean while the
//! whole-edge boundary residual stays bit-identical (kernel-probe
//! measurement, 2026-07-15).
//!
//! Invariants: boundary rows are invariant under both families (kernel
//! by `A·k = 0`; free set by construction order — the border rows are
//! assigned before `P` is consulted). Cross-sector eq (2)/(3) and
//! `sal` closure are structural for any shared net, so the only gate
//! this search can move is the fold verdict. The caller re-measures
//! fold and boundary independently before committing (all-or-nothing
//! acceptance unchanged).

use super::super::SharedControlRing;
use super::{ring_get, write_channel};
use crate::geometry::Point3;
use crate::geometry::gregory::{
    Bi6DeepInteriorFree, Bi6Reparams, GregoryPatch, RectangularGregoryBiNonic,
    assemble_bi6_polar_cap_from_net,
};

/// Fold tolerance for the full production detector — mirrors the
/// fillet promote path's `SELF_INT_TOL`.
const SELF_INT_TOL: f64 = 1e-6;

/// `sal` (deferred eq (8) closure) tolerance during candidate
/// assembly — structural machine-zero for a shared net; mirrors the
/// fillet wire's gate.
const SAL_TOLERANCE: f64 = 1e-6;

/// Deterministic multistart count. Start 0 descends from the fitted
/// configuration; further starts jitter the kernel coordinates with a
/// fixed-seed LCG.
const STARTS: u64 = 3;

/// Coordinate-descent rounds per start, with per-round step ladders.
const ROUNDS: usize = 3;

/// Relax the fold of the cap assembled from `(params, ring, free)` by
/// searching the boundary-kernel and free-set families. Returns
/// `Some((ring, free))` for the first configuration that passes the
/// full production self-intersection detector on every sector patch
/// (possibly the input itself), or `None` when the budgeted search
/// finds no detector-clean configuration.
pub(in super::super) fn relax_net_fold(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    n: usize,
    nu: usize,
) -> Option<(SharedControlRing, Bi6DeepInteriorFree)> {
    let m = kernel.len();
    let kdims = 3 * m;
    let fdims = free_dim(free);

    let kzero = vec![0.0; kdims];
    let fzero = vec![0.0; fdims];
    // Fast path: the fitted configuration may already be clean.
    if full_detector_clean(params, ring, free, kernel, &kzero, &fzero, n, nu) == Some(true) {
        return Some((ring.clone(), free.clone()));
    }

    for start in 0..STARTS {
        let mut kc = kzero.clone();
        let mut fc = fzero.clone();
        if start > 0 {
            let mut state = 0x9E37_79B9_7F4A_7C15_u64.wrapping_mul(start + 1);
            for c in kc.iter_mut() {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                *c = ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 0.02;
            }
        }
        let Some(mut best) = fast_score(params, ring, free, kernel, &kc, &fc, n, nu) else {
            continue;
        };
        for round in 0..ROUNDS {
            let steps: [f64; 2] = match round {
                0 => [0.02, 0.006],
                1 => [0.006, 0.002],
                _ => [0.002, 0.0006],
            };
            for &step in &steps {
                descend_kernel(
                    params, ring, free, kernel, &mut kc, &fc, step, &mut best, n, nu,
                );
                descend_free(
                    params, ring, free, kernel, &kc, &mut fc, step, &mut best, n, nu,
                );
            }
        }
        if full_detector_clean(params, ring, free, kernel, &kc, &fc, n, nu) == Some(true) {
            return Some(apply_shifts(ring, free, kernel, &kc, &fc, n, nu));
        }
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn descend_kernel(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &mut [f64],
    fc: &[f64],
    step: f64,
    best: &mut f64,
    n: usize,
    nu: usize,
) {
    for d in 0..kc.len() {
        for sign in [1.0, -1.0] {
            let saved = kc[d];
            kc[d] = saved + sign * step;
            match fast_score(params, ring, free, kernel, kc, fc, n, nu) {
                Some(score) if score > *best => *best = score,
                _ => kc[d] = saved,
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn descend_free(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &[f64],
    fc: &mut [f64],
    step: f64,
    best: &mut f64,
    n: usize,
    nu: usize,
) {
    for d in 0..fc.len() {
        for sign in [1.0, -1.0] {
            let saved = fc[d];
            fc[d] = saved + sign * step;
            match fast_score(params, ring, free, kernel, kc, fc, n, nu) {
                Some(score) if score > *best => *best = score,
                _ => fc[d] = saved,
            }
        }
    }
}

/// Fast fold objective: the minimum over sector patches of the
/// production Jacobian report's `min_jacobian_normal_cos` (larger is
/// less folded). `None` when the candidate fails to assemble.
#[allow(clippy::too_many_arguments)]
fn fast_score(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &[f64],
    fc: &[f64],
    n: usize,
    nu: usize,
) -> Option<f64> {
    let patches = build_patches(params, ring, free, kernel, kc, fc, n, nu)?;
    let mut score = f64::INFINITY;
    for gp in &patches {
        score = score.min(
            gp.jacobian_self_intersection_report()
                .min_jacobian_normal_cos,
        );
    }
    Some(score)
}

/// Full production verdict: `Some(true)` when every sector patch
/// passes `self_intersects_with_tolerance` (both tiers).
#[allow(clippy::too_many_arguments)]
fn full_detector_clean(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &[f64],
    fc: &[f64],
    n: usize,
    nu: usize,
) -> Option<bool> {
    let patches = build_patches(params, ring, free, kernel, kc, fc, n, nu)?;
    Some(
        !patches
            .iter()
            .any(|gp| gp.self_intersects_with_tolerance(SELF_INT_TOL)),
    )
}

/// Assemble the candidate's sector patches in the two-stripe slot
/// orientation (the frame every production gate measures).
#[allow(clippy::too_many_arguments)]
fn build_patches(
    params: &Bi6Reparams,
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &[f64],
    fc: &[f64],
    n: usize,
    nu: usize,
) -> Option<Vec<GregoryPatch>> {
    let (shifted_ring, shifted_free) = apply_shifts(ring, free, kernel, kc, fc, n, nu);
    let net: Vec<[[Point3; 4]; 4]> = (0..n)
        .map(|s| super::super::shared_net_layout(&shifted_ring, s, n))
        .collect();
    let cap = assemble_bi6_polar_cap_from_net(params, &net, &shifted_free, SAL_TOLERANCE).ok()?;
    if cap.patches_bi9.len() != n {
        return None;
    }
    Some(
        cap.patches_bi9
            .iter()
            .map(|grid| {
                let oriented: [[Point3; 10]; 10] = core::array::from_fn(|i| grid[9 - i]);
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                    grid: oriented,
                }))
            })
            .collect(),
    )
}

/// Apply kernel coordinates to the ring and free-set coordinates to
/// the free set, returning shifted clones.
fn apply_shifts(
    ring: &SharedControlRing,
    free: &Bi6DeepInteriorFree,
    kernel: &[Vec<f64>],
    kc: &[f64],
    fc: &[f64],
    n: usize,
    nu: usize,
) -> (SharedControlRing, Bi6DeepInteriorFree) {
    let m = kernel.len();
    let mut shifted = ring.clone();
    for coord in 0..3 {
        let mut channel: Vec<f64> = (0..nu).map(|u| ring_get(&shifted, u, n, coord)).collect();
        for (k, basis) in kernel.iter().enumerate() {
            let c = kc[coord * m + k];
            if c != 0.0 {
                for (u, b) in basis.iter().enumerate() {
                    channel[u] += c * b;
                }
            }
        }
        write_channel(&mut shifted, &channel, n, coord);
    }
    (shifted, shift_free(free, fc))
}

/// Number of scalar degrees of freedom in the free set.
fn free_dim(free: &Bi6DeepInteriorFree) -> usize {
    3 * (1 + 2 + 3 + 1 + 1 + free.p23.len() + free.p33.len())
}

/// Apply a flat scalar shift vector to a clone of the free set.
fn shift_free(free: &Bi6DeepInteriorFree, fc: &[f64]) -> Bi6DeepInteriorFree {
    let mut out = free.clone();
    let mut idx = 0;
    let mut shift = |p: &mut Point3| {
        p.x += fc[idx];
        p.y += fc[idx + 1];
        p.z += fc[idx + 2];
        idx += 3;
    };
    shift(&mut out.p00);
    for t in out.tan.iter_mut() {
        shift(t);
    }
    for t in out.two.iter_mut() {
        shift(t);
    }
    shift(&mut out.p220);
    shift(&mut out.p210);
    for p in out.p23.iter_mut() {
        shift(p);
    }
    for p in out.p33.iter_mut() {
        shift(p);
    }
    out
}
