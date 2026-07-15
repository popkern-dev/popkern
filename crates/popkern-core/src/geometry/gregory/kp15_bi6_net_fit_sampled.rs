//! Sampled-metric boundary fit for the KP15 shared control ring.
//!
//! The production cap-outer gate measures the *sampled* whole-edge
//! distance between the assembled cap's outer edges and the real
//! fillet boundary curves (11 points per edge), not the BB-CP
//! distance. Fitting the ring against the sampled points directly
//! therefore optimises the metric the acceptance actually enforces:
//! on the hexagon6_variable apex the CP-space fit floors at
//! ~5.09e-3 whole-edge while the sampled fit reaches ~2.62e-3 —
//! inside the 5e-3 gate (PK-087 kernel-probe measurement,
//! 2026-07-15).
//!
//! Also provides the orthonormal kernel basis of the boundary design
//! matrix: ring perturbations in the kernel leave every fitted
//! boundary row bit-identical while re-placing the net interior —
//! the search space of the fold relaxation
//! (`kp15_bi6_fold_relaxation`).

use super::super::SharedControlRing;
use super::{SectorOuterEdgeTarget, ring_get, write_channel};

/// Cubic Bernstein basis value `B³_i(t)`.
fn b3(i: usize, t: f64) -> f64 {
    let u = 1.0 - t;
    match i {
        0 => u * u * u,
        1 => 3.0 * u * u * t,
        2 => 3.0 * u * t * t,
        _ => t * t * t,
    }
}

/// Least-squares fit of the shared ring against the *sampled*
/// whole-edge points of both outer-edge targets per sector, with
/// Tikhonov ridge `lambda` toward the coarse seed.
///
/// `a_cp` is the CP-space boundary design matrix
/// (`boundary_design_matrix`); the sampled design rows are its
/// Bernstein combinations at 11 uniform parameters per edge. Returns
/// `None` when the ridge-regularised normal-equation solve is
/// singular.
pub(in super::super) fn sampled_boundary_fit(
    coarse: &SharedControlRing,
    target_edges: &[SectorOuterEdgeTarget],
    n: usize,
    nu: usize,
    a_cp: &[Vec<f64>],
    lambda: f64,
) -> Option<SharedControlRing> {
    let ts: Vec<f64> = (0..=10).map(|k| k as f64 / 10.0).collect();
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut targets: Vec<[f64; 3]> = Vec::new();
    for (s, sector_targets) in target_edges.iter().enumerate() {
        for (edge, tgt) in sector_targets.iter().enumerate() {
            let base = s * 8 + edge * 4;
            for &t in &ts {
                let mut row = vec![0.0; nu];
                for i in 0..4 {
                    let w = b3(i, t);
                    for (u, cell) in row.iter_mut().enumerate() {
                        *cell += w * a_cp[base + i][u];
                    }
                }
                rows.push(row);
                targets.push([
                    (0..4).map(|i| b3(i, t) * tgt[i].x).sum::<f64>(),
                    (0..4).map(|i| b3(i, t) * tgt[i].y).sum::<f64>(),
                    (0..4).map(|i| b3(i, t) * tgt[i].z).sum::<f64>(),
                ]);
            }
        }
    }
    let mut ring = coarse.clone();
    #[allow(clippy::needless_range_loop)]
    // `coord` indexes both `targets` cells and the ring channel
    for coord in 0..3 {
        let coarse_ch: Vec<f64> = (0..nu).map(|u| ring_get(coarse, u, n, coord)).collect();
        let mut ata = vec![vec![0.0; nu]; nu];
        let mut atr = vec![0.0; nu];
        for (r_idx, row) in rows.iter().enumerate() {
            let coarse_val: f64 = row.iter().zip(&coarse_ch).map(|(x, y)| x * y).sum();
            let resid = targets[r_idx][coord] - coarse_val;
            for i in 0..nu {
                atr[i] += row[i] * resid;
                for (j, cell) in ata[i].iter_mut().enumerate() {
                    *cell += row[i] * row[j];
                }
            }
        }
        for (i, row) in ata.iter_mut().enumerate() {
            row[i] += lambda * lambda;
        }
        let delta = crate::geometry::gregory::direct_linalg::solve_dense(ata, atr)?;
        let solved: Vec<f64> = (0..nu).map(|u| coarse_ch[u] + delta[u]).collect();
        write_channel(&mut ring, &solved, n, coord);
    }
    Some(ring)
}

/// Orthonormal basis of the kernel of the boundary design matrix
/// (rows = CP constraints) via Gram–Schmidt: orthonormalize the row
/// space, then complete standard basis vectors against it. Returns
/// `(row_rank, kernel_basis)`; every returned vector `k` satisfies
/// `|A·k| ≲ 1e-14`, so shifting the ring along `k` leaves the fitted
/// boundary rows unchanged.
pub(in super::super) fn boundary_kernel_basis(a: &[Vec<f64>], nu: usize) -> (usize, Vec<Vec<f64>>) {
    let mut q: Vec<Vec<f64>> = Vec::new();
    for row in a {
        let mut v = row.clone();
        for basis in &q {
            let dot: f64 = basis.iter().zip(&v).map(|(x, y)| x * y).sum();
            for (vi, bi) in v.iter_mut().zip(basis) {
                *vi -= dot * bi;
            }
        }
        let nrm = norm(&v);
        if nrm > 1e-10 {
            for vi in v.iter_mut() {
                *vi /= nrm;
            }
            q.push(v);
        }
    }
    let rank = q.len();
    let mut kernel: Vec<Vec<f64>> = Vec::new();
    for u in 0..nu {
        let mut v = vec![0.0; nu];
        v[u] = 1.0;
        for basis in q.iter().chain(kernel.iter()) {
            let dot: f64 = basis.iter().zip(&v).map(|(x, y)| x * y).sum();
            for (vi, bi) in v.iter_mut().zip(basis) {
                *vi -= dot * bi;
            }
        }
        let nrm = norm(&v);
        if nrm > 1e-10 {
            for vi in v.iter_mut() {
                *vi /= nrm;
            }
            kernel.push(v);
        }
    }
    (rank, kernel)
}

/// Euclidean norm of a scalar vector.
pub(in super::super) fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
