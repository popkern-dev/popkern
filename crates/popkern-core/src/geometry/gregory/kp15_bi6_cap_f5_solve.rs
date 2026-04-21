//! Phase I.10b + II.A — F_5 minimisation over the bi-6 cap's
//! free parameter set.
//!
//! Given a `Bi6CapInput` with **placeholder** free parameters
//! (the scaffold sets them to caller values or zeros), find the
//! free-parameter values that minimise F_5 across all N bi-6
//! patches of the cap.
//!
//! # Free parameter set
//!
//! Per [KP15 §3.6 step 4], the paper's free parameter set `P`
//! contains `p⁰_{22}`, `p^k_{23}`, `p^k_{33}` for k = 0..n − 1,
//! plus `p⁰_{21}` at `n ∈ {3, 6}`. The near-star CPs `p⁰_{01}`,
//! `p⁰_{02}`, `p⁰_{11}` and their (S1) propagates to k > 0 are
//! *derived* in step 3 from tangent-plane projection of Catmull-
//! Clark offsets — not free parameters in the paper.
//!
//! **This implementation extends the paper's free set** by
//! promoting the sector-0 near-star CPs to F_5-optimised free
//! parameters:
//! - Paper's set: `p⁰_{22}`, `p^k_{23}`, `p^k_{33}`
//! - **This code's set**: paper's set + `p⁰_{01}`, `p⁰_{02}`,
//!   `p⁰_{11}`
//!
//! Rationale (f5-tuning-plan.md option A): the step-3 / (S1)
//! propagation isn't implemented in this scaffold — callers
//! supply scaffold values via `sectors[k].p_01 / p_02 / p_11`
//! (typically from a tensor-border offset heuristic). Those
//! scaffold values can push the cap into a self-intersecting
//! configuration on sharp apex geometries. Promoting the
//! sector-0 near-star CPs to free lets the F_5 solver pick
//! values that avoid folds, at the cost of potentially breaking
//! paper's G² guarantee at the sector 0 / sector 1 boundary
//! (since sectors k > 0 still use scaffold values). When Phase
//! II option (B) lands, this extension can be revisited.
//!
//! Solver cardinality: `|P| = 2n + 4` per coordinate (3
//! independent SPD systems). Layout in the flat vector:
//!
//! ```text
//! [p0_01, p0_02, p0_11, p0_22, p_23[0], p_33[0], ..., p_23[n-1], p_33[n-1]]
//!    0       1       2      3     4        5               4+2(n-1)  5+2(n-1)
//! ```
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
use super::kp15_bi6_f5::{f_kappa_quadratic_form, f_kappa_quadratic_form_weighted};
use crate::geometry::Point3;

/// Free parameter vector. Covers the paper's interior-free set
/// (`p⁰_{22}`, `p^k_{23}`, `p^k_{33}`) plus the Phase II option-A
/// extension (sector-0 near-star CPs). See module docs for the
/// rationale.
#[derive(Debug, Clone)]
pub struct Bi6CapFreeParams {
    /// `p⁰_{01}` — sector 0 near-star CP, paper-derived in
    /// step 3. Promoted to free in option A to give F_5 the
    /// flexibility to avoid folds when scaffold tensor-border
    /// values push the cap into self-intersection.
    pub p0_01: Point3,
    /// `p⁰_{02}` — sector 0 near-star CP. See `p0_01`.
    pub p0_02: Point3,
    /// `p⁰_{11}` — sector 0 near-star CP. See `p0_01`.
    pub p0_11: Point3,
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
        let z = Point3::new(0.0, 0.0, 0.0);
        Self {
            p0_01: z,
            p0_02: z,
            p0_11: z,
            p0_22: z,
            p_23: vec![z; n],
            p_33: vec![z; n],
        }
    }

    /// Length of the flat vector representation (per coordinate).
    /// Layout: `[p0_01, p0_02, p0_11, p0_22, p_23[0], p_33[0],
    /// p_23[1], p_33[1], ..., p_23[n-1], p_33[n-1]]`.
    pub fn flat_len(n: usize) -> usize {
        4 + 2 * n
    }

    /// Probe direction: flat-index `i` is `(1, 0, 0)`, others
    /// `(0, 0, 0)`. The linear map is the same for all three
    /// coordinates by F_5's per-coord independence, so we only
    /// probe with the x-axis unit vector and reuse the resulting
    /// `M` for y and z.
    fn probe(n: usize, i: usize) -> Self {
        let mut out = Self::zero(n);
        let unit = Point3::new(1.0, 0.0, 0.0);
        match i {
            0 => out.p0_01 = unit,
            1 => out.p0_02 = unit,
            2 => out.p0_11 = unit,
            3 => out.p0_22 = unit,
            _ => {
                let offset = i - 4;
                let k = offset / 2;
                if offset.is_multiple_of(2) {
                    out.p_23[k] = unit;
                } else {
                    out.p_33[k] = unit;
                }
            }
        }
        out
    }
}

/// Apply a free-parameter vector to a base input: overwrites
/// `p0_22` + sector-0 near-star CPs, propagates sector 0's values
/// forward through sectors 1..n-1 via rule (S1), and builds
/// `interior_filler` for per-sector interior CPs at `(2, 3)` and
/// `(3, 3)`.
///
/// **Option B — (S1) propagation** (KP15 eq (6)): sector k > 0's
/// `(p_01, p_11, p_02)` are derived from sector (k-1)'s values
/// via `rule_s1_star_propagation`. This replaces the previous
/// behaviour where each sector's near-star CPs came from its own
/// tensor-border scaffold and were independently computed — a
/// pattern that broke the paper's cross-sector consistency
/// because the formulas for (S2) et al. use `p̀_{0j}` (previous
/// sector) to compute `ṕ_{ij}` (current sector).
///
/// The shared-edge CPs `(p_00, p_10, p_20, p_40, p_50, p_60)` are
/// left at their tensor-border values (not paper-strictly shared
/// across sectors because our scaffold has them per-sector; the
/// paper assumes all sectors' `p_i0` are identical). Closing that
/// gap is option B' — a separate refinement.
fn apply_free_params(base: &Bi6CapInput, free: &Bi6CapFreeParams) -> Bi6CapInput {
    let mut out = base.clone();
    out.p0_22 = free.p0_22;
    // Sector-0 near-star CPs (option A): overwrite with F_5 free
    // values.
    if !out.sectors.is_empty() {
        out.sectors[0].p_01 = free.p0_01;
        out.sectors[0].p_02 = free.p0_02;
        out.sectors[0].p_11 = free.p0_11;
    }
    // (S1) forward propagation for k = 1..n (option B).
    let n = out.sectors.len();
    for k in 1..n {
        let prev = out.sectors[k - 1];
        let s1 = super::kp15_bi6_rules::rule_s1_star_propagation(
            &out.params,
            prev.p_00,
            prev.p_01,
            prev.p_02,
            prev.p_10,
            prev.p_11,
            prev.p_20,
        );
        out.sectors[k].p_01 = s1.p01;
        out.sectors[k].p_11 = s1.p11;
        out.sectors[k].p_02 = s1.p02;
    }
    let mut filler: Vec<[[Option<Point3>; 7]; 7]> = vec![[[None; 7]; 7]; out.n];
    for k in 0..out.n {
        filler[k][2][3] = Some(free.p_23[k]);
        filler[k][3][3] = Some(free.p_33[k]);
    }
    out.interior_filler = Some(filler);
    out
}

/// Minimise F_5 over the free parameter set (see module docs).
/// Thin wrapper over [`f_kappa_solve_optimal_bi6_cap`]`(base, 5)`
/// kept for API compatibility.
pub fn f5_solve_optimal_bi6_cap(
    base: &Bi6CapInput,
) -> Result<(Bi6CapFreeParams, Bi6CapOutput), String> {
    f_kappa_solve_optimal_bi6_cap(base, 5)
}

/// Minimise F_κ with caller-supplied per-term weight multipliers
/// (option E of the F_5 tuning plan — per-valence tuning). Lets
/// the caller re-weight individual `(i, j)` terms away from the
/// paper's `κ!/(i!·j!)` defaults. See
/// [`f_kappa_quadratic_form_weighted`] for the weight layout.
pub fn f_kappa_solve_with_weights_optimal_bi6_cap(
    base: &Bi6CapInput,
    kappa: usize,
    weight_multipliers: &[f64],
) -> Result<(Bi6CapFreeParams, Bi6CapOutput), String> {
    let q = f_kappa_quadratic_form_weighted(kappa, weight_multipliers);
    f_kappa_solve_with_q(base, q)
}

/// Minimise F_κ for κ ∈ {1, 3, 4, 5} over the free parameter set.
/// Same as [`f5_solve_optimal_bi6_cap`] but parameterised on the
/// functional order.
///
/// Used by the F_κ fallback ladder (option C of the F_5 tuning
/// plan): when κ = 5 output self-intersects, retry with κ ∈
/// {4, 3, 1} to accept localised oscillation in exchange for
/// less aggressive smoothness penalisation. Lower κ → less
/// smoothness constraint → more flexibility for the solver to
/// avoid folds.
pub fn f_kappa_solve_optimal_bi6_cap(
    base: &Bi6CapInput,
    kappa: usize,
) -> Result<(Bi6CapFreeParams, Bi6CapOutput), String> {
    let q = f_kappa_quadratic_form(kappa);
    f_kappa_solve_with_q(base, q)
}

/// Internal: F_κ SPD solve given a precomputed 49×49 quadratic-
/// form matrix `q`. Shared by the default and weighted variants.
fn f_kappa_solve_with_q(
    base: &Bi6CapInput,
    q: [[f64; 49]; 49],
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
/// Layout mirrors `Bi6CapFreeParams::probe`.
fn assemble_free_params(n: usize, sx: &[f64], sy: &[f64], sz: &[f64]) -> Bi6CapFreeParams {
    let pt = |i: usize| Point3::new(sx[i], sy[i], sz[i]);
    let p0_01 = pt(0);
    let p0_02 = pt(1);
    let p0_11 = pt(2);
    let p0_22 = pt(3);
    let mut p_23 = Vec::with_capacity(n);
    let mut p_33 = Vec::with_capacity(n);
    for k in 0..n {
        p_23.push(pt(4 + 2 * k));
        p_33.push(pt(4 + 2 * k + 1));
    }
    Bi6CapFreeParams {
        p0_01,
        p0_02,
        p0_11,
        p0_22,
        p_23,
        p_33,
    }
}

#[cfg(test)]
#[path = "kp15_bi6_cap_f5_solve_tests.rs"]
mod tests;
