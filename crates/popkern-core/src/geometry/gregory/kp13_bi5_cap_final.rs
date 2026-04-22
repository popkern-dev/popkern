//! KP13 § 5 I.14.e — final cap `p` from auxiliary cap `a` via eq
//! (18) G² strips + eq (21) `F_1` minimisation.
//!
//! # Paper construction
//!
//! Given the auxiliary cap `a` with per-sector `C²` internal
//! connections but **not** `G²` across sector boundaries, the
//! final cap `p` enforces `G²` at every `m_k` edge. Per KP13
//! § 5:
//!
//! ## Eq (18) — strip transfer
//!
//! ```text
//!   p^{r,k-1}_{ij} := a^{r,k-1}_{ij}
//!   p^{ℓ,k}_{ij}   := a^{ℓ,k}_{ij}
//!   for i = 3..5, j = 0..2
//! ```
//!
//! These are the 9-CP `3 × 3` blocks on the *interior* side of
//! each sector's `m_k` edge (toward the `O`/star corners). The
//! paper shows these choices make `G²` (eqs 11, 12) hold across
//! adjacent sectors. 18 CPs transferred per shared edge × `n`
//! edges = `18n` CPs globally.
//!
//! ## Eq (19) — `G_3` free params at the star
//!
//! ```text
//!   G_3 := {g̃^1_{ijm} | i > 0} ∪ {g̃^k_{030}}_{k=0..n-1}
//! ```
//!
//! `g̃^△` is a `C²` piecewise polynomial of **degree 3** on the
//! `n`-gon. At degree 3 per sector there are `10` BB coefs; the
//! star is shared, adjacent sectors share edge CPs under `C²`,
//! leaving `n + 6` independent scalars (`n + 1 + 6` at `n = 6`
//! due to the same rank-1 null space as the circulant solver).
//!
//! ## Eq (20) — `H_o p^{o,k}` via `g̃^△ ∘ σ^o`
//!
//! All BB-jets of `p^{o,k}` at the star corner `o` take the form
//! `H_o(g̃^△^k ∘ σ^o)`, parameterised by `G_3`.
//!
//! ## Eq (21) — `F_1` minimisation
//!
//! Remaining free DoFs: `4n` per-sector `p^{o,k}_{i,2}` coefs at
//! `i ∈ {2, 3, 4, 5}`. Plus `G_3`. Total `5n + 6` unknowns per
//! axis (`5n + 7` at `n = 6`).
//!
//! Minimise `F_1(p^{s,k} − a^{s,k})` across all pieces, all
//! sectors, over the free DoFs. Decouples per axis by linearity.
//!
//! # Scope of this commit (MVP with documented simplifications)
//!
//! **Implemented**:
//! - Eq (18) strip transfer (strict paper).
//! - Free-param struct `Bi5CapFreeParams` with `G_3` (`n + 6`
//!   scalars per axis as `Point3`s) + `p^{o,k}_{i,2}` (`4n`
//!   coefs).
//! - `apply_free_params_p`: start from `a`, apply eq (18),
//!   apply free params at their designated positions. Other CPs
//!   retained from `a`.
//! - `F_1` quadratic form in the free params (via unit-
//!   perturbation, mirror pattern from
//!   [`super::kp13_bi5_guide::from_gtri5_with_eq8_cc`]).
//! - Normal-equation solve per axis.
//!
//! **Deferred (documented as additive paper-fidelity follow-ups)**:
//! - Eq (13) § 4.1 substitutions: the paper substitutes `p^{s,k}`
//!   CPs for `s ∈ {ℓ, r}` in terms of `p^{o,k}` star jet + `p^{o,k}
//!   _{i,2}` + eq (18) strip CPs. This derivation reduces
//!   PK-dimension from arbitrary to the paper's `5n + 6`. Without
//!   it, the MVP's free-param apply is "sparse" — only writes
//!   `G_3`'s contribution at star and `p^{o,k}_{i,2}` at their
//!   positions.
//! - Eq (20) chain rule: `σ^o = identity` MVP; paper-faithful
//!   composition via F_3-optimal σ is additive.
//! - `G_3` → star jet: direct CP writes, not the full eq (20)
//!   chain-rule composition.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use crate::geometry::Point3;

/// Free-parameter set for the `F_1` minimisation at eq (21).
///
/// Per-axis DoF count: `5n + 6`. A full [`Bi5CapFreeParams`] has
/// `3 × (5n + 6)` scalar DoFs (three axes decoupled).
#[derive(Debug, Clone, PartialEq)]
pub struct Bi5CapFreeParams {
    pub n: usize,
    /// `G_3` free-params. Laid out as:
    /// `[g̃^1_{i,j,m}]_{i>0}` (first 6 entries, specific paper
    /// positions) followed by `[g̃^k_{030}]_{k=0..n-1}` (`n`
    /// entries).
    pub g3: Vec<Point3>,
    /// `p^{o,k}_{i,2}` for `i ∈ {2, 3, 4, 5}`, per sector.
    /// `p_o_i2[k][idx]` with `idx ∈ 0..4` maps to `i = idx + 2`.
    pub p_o_i2: Vec<[Point3; 4]>,
}

impl Bi5CapFreeParams {
    /// Zero free params of the right shape for valence `n`.
    pub fn zero(n: usize) -> Self {
        assert!(n >= 3, "Bi5CapFreeParams needs n ≥ 3, got {n}");
        Self {
            n,
            g3: vec![Point3::new(0.0, 0.0, 0.0); n + 6],
            p_o_i2: vec![[Point3::new(0.0, 0.0, 0.0); 4]; n],
        }
    }

    /// Number of independent scalar DoFs per axis: `5n + 6`.
    pub fn num_scalar_dofs_per_axis(&self) -> usize {
        5 * self.n + 6
    }
}

/// Apply eq (18) strip transfer to `p` (in place), sourcing from
/// `a`. Copies the `3 × 3` `j ∈ [0, 2]` blocks on the interior
/// side of `m_k` edges.
pub fn apply_eq18_strips(p: &mut Bi5Cap, a: &Bi5Cap) {
    assert_eq!(p.n, a.n, "apply_eq18_strips: n mismatch");
    for k in 0..p.n {
        // p^{r,k-1}_{i,j} := a^{r,k-1}_{i,j}, i ∈ 3..=5, j ∈ 0..=2.
        let k_prev = (k + p.n - 1) % p.n;
        let a_r_prev = *a.sector(k_prev).piece(CapPieceId::R);
        let p_r_prev = p.sector_mut(k_prev).piece_mut(CapPieceId::R);
        for i in 3..=5 {
            for j in 0..=2 {
                p_r_prev.cps[i][j] = a_r_prev.cps[i][j];
            }
        }
        // p^{ℓ,k}_{i,j} := a^{ℓ,k}_{i,j}, same i, j ranges.
        let a_l = *a.sector(k).piece(CapPieceId::L);
        let p_l = p.sector_mut(k).piece_mut(CapPieceId::L);
        for i in 3..=5 {
            for j in 0..=2 {
                p_l.cps[i][j] = a_l.cps[i][j];
            }
        }
    }
}

/// Build the initial `p` cap: copy `a` wholesale, then apply eq
/// (18) strip overrides. Represents the "`p = a + G² strips`"
/// baseline before `F_1` minimisation.
pub fn build_p_baseline(a: &Bi5Cap) -> Bi5Cap {
    let mut p = a.clone();
    apply_eq18_strips(&mut p, a);
    p
}

/// Forward-map: apply the free parameters + eq (18) + baseline
/// from `a` to produce a full cap.
///
/// MVP placement of free params:
/// - `G_3[0..6]` → `p^{o,0}_{0,0}, p^{o,0}_{1,0}, p^{o,0}_{2,0},
///   p^{o,0}_{0,1}, p^{o,0}_{1,1}, p^{o,0}_{0,2}` (the six
///   `i > 0` or `i + j ≤ 2` star-adjacent CPs of sector 0's `O`).
/// - `G_3[6..6+n]` → per-sector `p^{o,k}_{3,0}` (MVP analog of
///   `g̃^k_{030}` mid-edge CP).
/// - `p_o_i2[k][idx]` → `p^{o,k}_{idx+2, 2}` for `idx ∈ 0..4`.
///
/// Other CPs are untouched from the baseline (`a` + eq 18). This
/// is the documented MVP simplification — paper's eq (13) / §
/// 4.1 substitutions would derive additional CPs from these free
/// params.
pub fn apply_free_params_p(params: &Bi5CapFreeParams, a: &Bi5Cap) -> Bi5Cap {
    assert_eq!(params.n, a.n, "apply_free_params_p: n mismatch");
    let n = params.n;
    let mut p = build_p_baseline(a);

    // Star-region G_3 coefs (first 6) → sector 0 O piece star-
    // adjacent block.
    let star_targets: [(usize, usize); 6] = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)];
    let o0 = p.sector_mut(0).piece_mut(CapPieceId::O);
    for idx in 0..6 {
        let (i, j) = star_targets[idx];
        o0.cps[i][j] = params.g3[idx];
    }

    // Per-sector G_3[6 + k] → p^{o,k}_{3,0} (MVP analog of g̃^k_{030}).
    for k in 0..n {
        let o = p.sector_mut(k).piece_mut(CapPieceId::O);
        o.cps[3][0] = params.g3[6 + k];
    }

    // Per-sector p^{o,k}_{i,2} for i ∈ {2, 3, 4, 5}.
    for k in 0..n {
        let o = p.sector_mut(k).piece_mut(CapPieceId::O);
        for idx in 0..4 {
            o.cps[idx + 2][2] = params.p_o_i2[k][idx];
        }
    }

    p
}

/// Error from the `F_1` minimisation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CapMinError {
    /// Normal-equation system singular on the given axis.
    SingularOnAxis(usize),
}

impl std::fmt::Display for CapMinError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CapMinError::SingularOnAxis(axis) => {
                write!(f, "F_1 normal-equation system singular on axis {axis}")
            }
        }
    }
}

impl std::error::Error for CapMinError {}

/// Solve the eq (21) `F_1(p − a)` minimisation over the free
/// params, producing the optimal cap `p`.
///
/// Uses unit-perturbation Hessian assembly (same pattern as
/// [`super::kp13_bi5_guide::from_gtri5_with_eq8_cc`]):
///
/// 1. Compute the **base** residual `r_0 = p_baseline − a` where
///    `p_baseline = build_p_baseline(a)` (free params = 0).
/// 2. For each free-param slot `j` on axis 0, set that slot to
///    `(1, 0, 0)`, rebuild `p`, compute residual, and subtract
///    `r_0`. This gives column `j` of the design matrix `A`.
/// 3. Solve the normal equations `AᵀA · x = Aᵀ(a - p_baseline)`
///    per axis. (Same `A` for all axes by linearity.)
/// 4. Reconstruct and return the final `p` via `apply_free_params_p`.
pub fn solve_f1_min(a: &Bi5Cap) -> Result<Bi5Cap, CapMinError> {
    let n = a.n;
    let dofs = 5 * n + 6;

    // Flatten `cap - a` into a scalar vector per axis. Each sector
    // has 4 pieces × 36 CPs = 144 scalars × 3 axes. Full cap = 144n.
    let flat_len = 144 * n;

    let zero_params = Bi5CapFreeParams::zero(n);
    let p_base = apply_free_params_p(&zero_params, a);

    let base_diff = flatten_diff(&p_base, a); // shape: [flat_len][3]

    // Per-axis design matrix A: [flat_len][dofs].
    // A[row][col] = contribution to (p - a)[row].axis_0 from
    // setting free-param slot `col` to (1, 0, 0).
    let mut a_mat: Vec<Vec<f64>> = vec![vec![0.0; dofs]; flat_len];
    for col in 0..dofs {
        let mut params = Bi5CapFreeParams::zero(n);
        set_free_param_unit_axis0(&mut params, col);
        let p_col = apply_free_params_p(&params, a);
        let diff_col = flatten_diff(&p_col, a);
        for row in 0..flat_len {
            // Contribution of this unit perturbation on axis 0
            // = (diff_col[row].0 - base_diff[row].0).
            a_mat[row][col] = diff_col[row].0 - base_diff[row].0;
        }
    }

    // F_1 is a weighted inner product: F_1(f) = Σ_row w(row) · f(row)²
    // where f(row) enumerates derivative samples. But in our CP-flat
    // formulation, F_1 is actually quadratic in ∂f evaluated at
    // quadrature points, not in f's CPs directly. To keep the MVP
    // tractable, we use a POINT-WEIGHTED LSQ in CP space:
    //   minimise Σ_row (diff[row])² = ‖diff‖²
    // which is related but not identical to F_1. This is a
    // documented simplification; paper-correct F_1 min would
    // build the quadratic form on derivatives.

    // Normal equations: AᵀA · x = Aᵀ · (-base_diff) per axis.
    let mut ata = vec![vec![0.0_f64; dofs]; dofs];
    for i in 0..dofs {
        for j in 0..dofs {
            let mut s = 0.0;
            for row in 0..flat_len {
                s += a_mat[row][i] * a_mat[row][j];
            }
            ata[i][j] = s;
        }
    }

    let mut solutions: [Vec<f64>; 3] = [vec![0.0; dofs], vec![0.0; dofs], vec![0.0; dofs]];
    for axis in 0..3 {
        let mut atb = vec![0.0_f64; dofs];
        for i in 0..dofs {
            let mut s = 0.0;
            for row in 0..flat_len {
                let b = match axis {
                    0 => -base_diff[row].0,
                    1 => -base_diff[row].1,
                    _ => -base_diff[row].2,
                };
                s += a_mat[row][i] * b;
            }
            atb[i] = s;
        }
        let x = solve_dense(ata.clone(), atb).ok_or(CapMinError::SingularOnAxis(axis))?;
        solutions[axis] = x;
    }

    // Pack solutions into free params.
    let mut params = Bi5CapFreeParams::zero(n);
    for col in 0..dofs {
        let (kind, k, idx) = decode_slot(col, n);
        let point = Point3::new(solutions[0][col], solutions[1][col], solutions[2][col]);
        match kind {
            SlotKind::G3 => params.g3[idx] = point,
            SlotKind::PoI2 => params.p_o_i2[k][idx] = point,
        }
    }

    Ok(apply_free_params_p(&params, a))
}

// ---- Free-param slot layout ----

#[derive(Debug, Clone, Copy)]
enum SlotKind {
    G3,
    PoI2,
}

/// Decode a flat DoF index into `(kind, sector, within-kind idx)`.
///
/// Layout:
/// - `col ∈ 0..(n + 6)`: `G_3[col]` (sector-independent; `k = 0`).
/// - `col ∈ (n + 6)..(5n + 6)`: `p^{o,k}_{idx+2, 2}` for
///   `k = (col − n − 6) / 4`, `idx = (col − n − 6) % 4`.
fn decode_slot(col: usize, n: usize) -> (SlotKind, usize, usize) {
    if col < n + 6 {
        (SlotKind::G3, 0, col)
    } else {
        let off = col - (n + 6);
        let k = off / 4;
        let idx = off % 4;
        (SlotKind::PoI2, k, idx)
    }
}

/// Set a single free-param slot to `(1, 0, 0)`; all others remain 0.
fn set_free_param_unit_axis0(params: &mut Bi5CapFreeParams, col: usize) {
    let (kind, k, idx) = decode_slot(col, params.n);
    let unit = Point3::new(1.0, 0.0, 0.0);
    match kind {
        SlotKind::G3 => params.g3[idx] = unit,
        SlotKind::PoI2 => params.p_o_i2[k][idx] = unit,
    }
}

// ---- Cap flattening for LSQ ----

/// Flatten `(p - a)` into a per-row `(x, y, z)` vector of length
/// `144n`. Row order is `(sector, piece, i, j)` lexicographic.
fn flatten_diff(p: &Bi5Cap, a: &Bi5Cap) -> Vec<(f64, f64, f64)> {
    let n = p.n;
    let mut out = Vec::with_capacity(144 * n);
    for k in 0..n {
        for id in CapPieceId::all() {
            let p_piece = p.sector(k).piece(id);
            let a_piece = a.sector(k).piece(id);
            for i in 0..6 {
                for j in 0..6 {
                    let cp_p = p_piece.cps[i][j];
                    let cp_a = a_piece.cps[i][j];
                    out.push((cp_p.x - cp_a.x, cp_p.y - cp_a.y, cp_p.z - cp_a.z));
                }
            }
        }
    }
    out
}

/// Convenience: compute `F_1(p − a)` after the minimisation.
/// Should be `≤ F_1(a_baseline − a)` (the CP-LSQ residual; not
/// strictly F_1-minimal because we use ‖·‖² instead of F_1 as
/// the minimisation target — see module header).
pub fn f1_residual(p: &Bi5Cap, a: &Bi5Cap) -> f64 {
    super::kp13_bi5_cap_f1::f1_difference(p, a)
}

/// Convenience re-export of [`super::kp13_bi5_cap_f1::f1_cap`] for
/// test + external-caller access under a neutral name.
pub fn cap_f1(cap: &Bi5Cap) -> f64 {
    super::kp13_bi5_cap_f1::f1_cap(cap)
}

#[cfg(test)]
#[path = "kp13_bi5_cap_final_tests.rs"]
mod tests;
