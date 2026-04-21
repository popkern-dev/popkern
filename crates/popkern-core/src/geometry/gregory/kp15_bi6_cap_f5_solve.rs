//! Phase I.10b — F_5 minimisation over the bi-6 cap's free
//! parameter set.
//!
//! Given a `Bi6CapInput` with **placeholder** free parameters
//! (the scaffold sets them to caller values or zeros), find the
//! free-parameter values that minimise F_5 across all N bi-6
//! patches of the cap.
//!
//! # Free parameter set
//!
//! Per [KP15 §3.6 step 4], the free parameter set `P` contains:
//! - `p⁰_{ij}` for `0 ≤ i + j ≤ 2` (excluding `(0, 0)` which is
//!   Catmull-Clark-fixed)
//! - `p⁰_{22}`
//! - `p^k_{23}`, `p^k_{33}` for k = 0..n − 1
//! - `p⁰_{21}` at `n ∈ {3, 6}`
//!
//! This scaffold-compatible solver covers the **interior subset**:
//! `p⁰_{22}`, `p^k_{23}`, `p^k_{33}`. The near-star `p⁰_{ij}` for
//! `i + j ≤ 2` are treated as caller-supplied `sectors[k].p_01`,
//! `p_02`, `p_11` for each sector (no (S1) propagation in the
//! scaffold). Full paper-compliant (S1) propagation is deferred
//! to I.10c when we wire into the fillet crate.
//!
//! Solver cardinality `|P|_scaffold = 2n + 1` (per coordinate;
//! 3 independent SPD systems of this size).
//!
//! # Algorithm
//!
//! 1. Assemble the cap with all free params = 0 → bias `b`.
//! 2. For each free param index `i ∈ 0..|P|`, assemble with the
//!    i-th probe direction (1 at index i, 0 elsewhere) and
//!    subtract the bias → column `M[:, i]` of the linear map.
//! 3. Compute `H = Σ_k M_k^T · Q · M_k` (same for all coords).
//! 4. Per coord `d ∈ {x, y, z}`, compute
//!    `g_d = Σ_k M_k^T · Q · (b_k)_d` and solve `H · p_d* = −g_d`.
//! 5. Assemble final cap with optimal `P`.

#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use super::kp15_bi6_cap::{Bi6CapInput, Bi6CapOutput, assemble_bi6_polar_cap};
use super::kp15_bi6_f5::f5_quadratic_form;
use crate::geometry::Point3;

/// Scaffold-compatible free parameter vector. Covers the interior
/// subset of the KP15 §3.6 step 4 set (see module docs).
#[derive(Debug, Clone)]
pub struct Bi6CapFreeParams {
    /// `p⁰_{22}` — center-diagonal CP at sector 0.
    pub p0_22: Point3,
    /// `p^k_{23}` for k = 0..n − 1 — row-2 column-3 interior CP.
    pub p_23: Vec<Point3>,
    /// `p^k_{33}` for k = 0..n − 1 — row-3 column-3 interior CP.
    pub p_33: Vec<Point3>,
}

impl Bi6CapFreeParams {
    /// All zero, for use as bias in linear-map extraction.
    pub fn zero(n: usize) -> Self {
        Self {
            p0_22: Point3::new(0.0, 0.0, 0.0),
            p_23: vec![Point3::new(0.0, 0.0, 0.0); n],
            p_33: vec![Point3::new(0.0, 0.0, 0.0); n],
        }
    }

    /// Length of the flat vector representation (per coordinate).
    pub fn flat_len(n: usize) -> usize {
        1 + 2 * n
    }

    /// Probe direction: flat-index `i` is `(1, 0, 0)`, others
    /// `(0, 0, 0)`. The linear map is the same for all three
    /// coordinates by F_5's per-coord independence, so we only
    /// probe with the x-axis unit vector and reuse the resulting
    /// `M` for y and z.
    ///
    /// Layout: `[p0_22, p_23[0], p_33[0], p_23[1], p_33[1], ...,
    /// p_23[n-1], p_33[n-1]]`.
    fn probe(n: usize, i: usize) -> Self {
        let mut out = Self::zero(n);
        let unit = Point3::new(1.0, 0.0, 0.0);
        if i == 0 {
            out.p0_22 = unit;
        } else {
            let k = (i - 1) / 2;
            if (i - 1).is_multiple_of(2) {
                out.p_23[k] = unit;
            } else {
                out.p_33[k] = unit;
            }
        }
        out
    }
}

/// Apply a free-parameter vector to a base input: overwrites
/// `p0_22` and builds `interior_filler` so `p^k_{23} =
/// free.p_23[k]` and `p^k_{33} = free.p_33[k]`.
fn apply_free_params(base: &Bi6CapInput, free: &Bi6CapFreeParams) -> Bi6CapInput {
    let mut out = base.clone();
    out.p0_22 = free.p0_22;
    let mut filler: Vec<[[Option<Point3>; 7]; 7]> = vec![[[None; 7]; 7]; out.n];
    for k in 0..out.n {
        filler[k][2][3] = Some(free.p_23[k]);
        filler[k][3][3] = Some(free.p_33[k]);
    }
    out.interior_filler = Some(filler);
    out
}

/// Minimise F_5 over the scaffold's free parameter subset. Returns
/// the optimal free-param values and the corresponding fully-
/// assembled cap.
pub fn f5_solve_optimal_bi6_cap(
    base: &Bi6CapInput,
) -> Result<(Bi6CapFreeParams, Bi6CapOutput), String> {
    let n = base.n;
    let p_len = Bi6CapFreeParams::flat_len(n);
    let cp_per_patch = 49;
    let total_cps = cp_per_patch * n;

    // Step 1: bias. Assemble with all free params = 0.
    let zero = Bi6CapFreeParams::zero(n);
    let bias_input = apply_free_params(base, &zero);
    let bias_out = assemble_bi6_polar_cap(&bias_input)?;

    // Flatten bias patches into per-coord vectors of length total_cps.
    let mut bias_x = vec![0.0_f64; total_cps];
    let mut bias_y = vec![0.0_f64; total_cps];
    let mut bias_z = vec![0.0_f64; total_cps];
    for k in 0..n {
        for i in 0..7 {
            for j in 0..7 {
                let idx = k * cp_per_patch + i * 7 + j;
                let p = bias_out.patches_bi6[k][i][j];
                bias_x[idx] = p.x;
                bias_y[idx] = p.y;
                bias_z[idx] = p.z;
            }
        }
    }

    // Step 2: build linear map M (total_cps × p_len). Use x-coord
    // probes; map is same for y/z by linearity.
    let mut m: Vec<Vec<f64>> = vec![vec![0.0_f64; p_len]; total_cps];
    for p_idx in 0..p_len {
        let probe = Bi6CapFreeParams::probe(n, p_idx);
        let probe_input = apply_free_params(base, &probe);
        let probe_out = assemble_bi6_polar_cap(&probe_input)?;
        for k in 0..n {
            for i in 0..7 {
                for j in 0..7 {
                    let idx = k * cp_per_patch + i * 7 + j;
                    m[idx][p_idx] = probe_out.patches_bi6[k][i][j].x - bias_x[idx];
                }
            }
        }
    }

    // Step 3: build H = Σ_k M_k^T · Q · M_k. Iterate per sector
    // because Q is per-patch.
    let q = f5_quadratic_form();
    let mut h = vec![vec![0.0_f64; p_len]; p_len];
    // Per sector, extract the 49 × p_len sub-block of M and
    // accumulate M_k^T · Q · M_k into H.
    for k in 0..n {
        // Extract M_k (49 × p_len).
        let mut m_k = vec![vec![0.0_f64; p_len]; cp_per_patch];
        for cp_idx in 0..cp_per_patch {
            let global_idx = k * cp_per_patch + cp_idx;
            for p_idx in 0..p_len {
                m_k[cp_idx][p_idx] = m[global_idx][p_idx];
            }
        }
        // QM = Q · M_k (49 × p_len).
        let mut qm = vec![vec![0.0_f64; p_len]; cp_per_patch];
        for r in 0..cp_per_patch {
            for c in 0..p_len {
                let mut s = 0.0;
                for q_idx in 0..cp_per_patch {
                    s += q[r][q_idx] * m_k[q_idx][c];
                }
                qm[r][c] = s;
            }
        }
        // H += M_k^T · QM (p_len × p_len).
        for r in 0..p_len {
            for c in 0..p_len {
                let mut s = 0.0;
                for m_row in 0..cp_per_patch {
                    s += m_k[m_row][r] * qm[m_row][c];
                }
                h[r][c] += s;
            }
        }
    }

    // Step 4: per coord, compute g_d = Σ_k M_k^T · Q · b_k_d and
    // solve H · p* = -g_d.
    let mut solutions: [Vec<f64>; 3] = [
        vec![0.0_f64; p_len],
        vec![0.0_f64; p_len],
        vec![0.0_f64; p_len],
    ];
    for (d_idx, bias) in [&bias_x, &bias_y, &bias_z].iter().enumerate() {
        let mut g = vec![0.0_f64; p_len];
        for k in 0..n {
            // Extract b_k (length 49).
            let mut b_k = vec![0.0_f64; cp_per_patch];
            for cp_idx in 0..cp_per_patch {
                b_k[cp_idx] = bias[k * cp_per_patch + cp_idx];
            }
            // Q · b_k.
            let mut qb = vec![0.0_f64; cp_per_patch];
            for r in 0..cp_per_patch {
                let mut s = 0.0;
                for c in 0..cp_per_patch {
                    s += q[r][c] * b_k[c];
                }
                qb[r] = s;
            }
            // g += M_k^T · qb.
            for p_idx in 0..p_len {
                let mut s = 0.0;
                for m_row in 0..cp_per_patch {
                    s += m[k * cp_per_patch + m_row][p_idx] * qb[m_row];
                }
                g[p_idx] += s;
            }
        }
        // Solve H · p* = -g.
        let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
        let p_star = solve_dense(h.clone(), neg_g)
            .ok_or_else(|| format!("F_5 SPD solve singular at coord {d_idx}"))?;
        solutions[d_idx] = p_star;
    }

    // Step 5: assemble optimal free params from per-coord solutions.
    let optimal = assemble_free_params(n, &solutions[0], &solutions[1], &solutions[2]);

    // Step 6: re-assemble cap with optimal free params.
    let final_input = apply_free_params(base, &optimal);
    let final_out = assemble_bi6_polar_cap(&final_input)?;

    Ok((optimal, final_out))
}

/// Pack three per-coord flat solutions into a `Bi6CapFreeParams`.
fn assemble_free_params(n: usize, sx: &[f64], sy: &[f64], sz: &[f64]) -> Bi6CapFreeParams {
    let p0_22 = Point3::new(sx[0], sy[0], sz[0]);
    let mut p_23 = Vec::with_capacity(n);
    let mut p_33 = Vec::with_capacity(n);
    for k in 0..n {
        p_23.push(Point3::new(sx[1 + 2 * k], sy[1 + 2 * k], sz[1 + 2 * k]));
        p_33.push(Point3::new(
            sx[1 + 2 * k + 1],
            sy[1 + 2 * k + 1],
            sz[1 + 2 * k + 1],
        ));
    }
    Bi6CapFreeParams { p0_22, p_23, p_33 }
}

#[cfg(test)]
#[path = "kp15_bi6_cap_f5_solve_tests.rs"]
mod tests;
