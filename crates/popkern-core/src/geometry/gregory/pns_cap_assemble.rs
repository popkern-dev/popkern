//! Phase IV PnS — End-to-end bi-4 cap assembly.
//!
//! Composes the PnS pipeline primitives into a single function
//! [`assemble_pns_bi4_cap`] that turns a boundary ring + pole into
//! `N` bi-4 sector patches:
//!
//! ```text
//! (Bi2TensorBorderRing, β, σ, pole) → N × [[Point3; 5]; 5]
//! ```
//!
//! # Architecture
//!
//! Storage: `patches[k][u_c][v_c]` uses **canonical** indexing
//! (sector's own parameterization). Each sector's 25 BB coefs are
//! classified by source:
//!
//! | Source                          | Canonical `(u_c, v_c)`             | Count |
//! |---------------------------------|------------------------------------|-------|
//! | Pole (constant)                 | `(0, 0)`                           | 1     |
//! | Central `a∘σ^k`                 | `u_c+v_c ≤ 2` (excl. pole)         | 5     |
//! | Boundary `b∘β` (v_c=3)          | `(0..4, 3)`                        | 5     |
//! | Boundary `b∘β` (v_c=4)          | `(0..4, 4)`                        | 5     |
//! | Prev-sharing (cur-prev seam)    | `(1..4, 0)`                        | 4     |
//! | Eq (2) closed-forms             | `(2, 1), (3, 1), (4, 1)`           | 3     |
//! | Free (F_κ min)                  | `(1, 2), (2, 2), (3, 2), (4, 2)`   | 4     |
//!
//! Sector-sharing relations (in canonical coords):
//! - cur-nxt seam: `cur[0][v_c] = nxt[v_c][0]` for `v_c ∈ 0..4`
//! - cur-prev seam: `cur[u_c][0] = prev[0][u_c]` for `u_c ∈ 0..4`
//!
//! Paper vs canonical axis swap (for cur side only):
//! - paper `p_{ij}` (cur) = canonical `p^k_{j, i}` (swapped)
//! - paper `´p_{ij}` (nxt) = canonical `p^{k+1}_{i, j}` (no swap)
//!
//! Eq (2) closed-forms (derived from paper's 5 outputs):
//! - paper cur_{10}, cur_{20}: automatic via central chain rule +
//!   paper's σ_{20,x} = (1 − 4/(3c))·d_1 + (4·cos(α/2)/(3c))·d_2
//!   (transcription fix 2026-04-24: was (1 − 4c/3), should be
//!   (1 − 4/(3c)))
//! - paper cur_{40} → canonical (4, 1) for sector k:
//!   `p^k_{4, 1} = 2·p^{k-1}_{0, 4} − p^{k-1}_{1, 4}`
//! - paper cur_{30} → canonical (3, 1) for sector k:
//!   `p^k_{3, 1} = [2·(4-c)·p^{k-1}_{0, 3} + c·(p^{k-1}_{1, 4} + p^k_{4, 1}) − 4·p^{k-1}_{1, 3}] / 4`
//! - paper nxt_{21} → canonical (2, 1) for sector k:
//!   `p^k_{2, 1} = −p^{k-1}_{1, 2} + (2-c)·p^{k-1}_{0, 2} + c·p^{k-1}_{0, 3}`
//!
//! # DOFs and the F_κ solve
//!
//! The truly free variables, determined jointly by the F_κ
//! minimisation (paper eq 6), are:
//!
//! - `a`'s 5 Point3 components `(a_x, a_y, a_xx, a_xy, a_yy)` — 15
//!   scalars (shared across sectors; `a_0` pinned to pole)
//! - `4 × N` free interior CPs per sector — `12·N` scalars
//!
//! Total: `15 + 12·N` scalars. Problem is separable per coordinate
//! (x, y, z), so we solve three `(5 + 4·N)`-dimensional symmetric
//! PSD systems.
//!
//! # Valence support
//!
//! Supported: `n ∈ [3, 12] \ {4}`. The bi-4 construction is
//! numerically singular at `n=4` (σ_{20,x} has a `1/c` factor
//! with `c=cos(π/2)=0`); for that case callers should use the
//! bi-3 variant or a different cap construction.
//!
//! Even `n ∈ {6, 8, 10, 12}` is supported through the asymmetric
//! closed-form treatment: instead of symmetrically solving
//! p_{2,1} via a circulant pair-sum (rank `n−1` for even n), we
//! treat canonical `(1, 2)` as a free DOF and derive canonical
//! `(2, 1)` one-way from the previous sector's data. No rank
//! deficiency to resolve.
//!
//! # Auto-compute pole
//!
//! [`pole_from_boundary`] derives the central pole from the
//! boundary ring via paper's eq (5) (see [`super::pns_center`]).
//! Useful for callers who don't want to compute the pole
//! manually. Assumes bi-4 e_def weight.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §4 (eqs 2, 4, 5, 6), Fig. 8c. The
//! paper-vs-canonical axis swap was clarified by clean-room
//! feedback from the KK-Peters reference implementation author
//! (2026-04-24); see memory note `project_pns_sigma_bug.md`.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use super::pns_beta::BetaReparam;
use super::pns_boundary::{Bi2TensorBorderRing, compose_bi4_outer_layers};
use super::pns_center::{central_point, e_def_bi4};
use super::pns_central_quad::{QuadMapJet, central_bb_from_a_and_sigma, rotate_sigma_jet};
use super::pns_gram::{bi4_flat_index, gram_bi4_f_kappa};
use super::pns_sigma::SigmaCenterJet;
use crate::geometry::Point3;

/// Inputs for the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleInput {
    /// Per-sector bi-2 tensor-border strips.
    pub boundary: Bi2TensorBorderRing,
    /// β reparameterisation for `n` sectors.
    pub beta: BetaReparam,
    /// Canonical σ centre jet.
    pub sigma: SigmaCenterJet,
    /// Pole (central point, eq 5).
    pub pole: Point3,
    /// F_κ order, `κ ∈ 1..=4`.
    pub kappa: usize,
}

/// Result of the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleResult {
    /// Fitted quadratic map `a` (with `a.a_0 = pole`).
    pub a: QuadMapJet,
    /// `N` bi-4 sector patches.
    pub patches: Vec<[[Point3; 5]; 5]>,
}

/// Cap assembly error cases.
#[derive(Debug, thiserror::Error)]
pub enum CapAssembleError {
    /// Valence is out of the supported range `[3, 12]`.
    #[error("PnS cap assemble: n must be in [3, 12], got {0}")]
    InvalidValence(usize),
    /// Valence `n = 4` is not supported — σ has a `1/c` singularity
    /// with `c = cos(π/2) = 0`.
    #[error(
        "PnS cap assemble: n=4 is not supported by the bi-4 construction (σ_{{20,x}} has a 1/c singularity with c=cos(π/2)=0). Use the bi-3 variant or a different cap."
    )]
    N4NotSupported,
    /// `κ` is out of the supported range `1..=4`.
    #[error("PnS cap assemble: kappa must be 1..=4, got {0}")]
    InvalidKappa(usize),
    /// The F_κ linear system was singular (shouldn't happen in
    /// practice for well-posed input).
    #[error("PnS cap assemble: F_κ linear system is singular")]
    SingularSystem,
}

/// The 4 "free" CPs per sector, in fixed order for DOF indexing.
/// **Canonical** `(u_c, v_c)` indexing. These are the deep-interior
/// CPs not pinned by:
/// - pole `(0, 0)`
/// - central `u_c + v_c ≤ 2` (via a∘σ^k)
/// - prev-sharing on `v_c = 0` (cur's seam with prev)
/// - boundary `v_c ∈ {3, 4}` (via b∘β)
/// - eq (2) closed-forms for canonical `(2, 1), (3, 1), (4, 1)` —
///   pinned by prev sector's `(1, 2)` + prev's boundary (see
///   [`build_coef_map`] Pass 3).
const FREE_CPS: [(usize, usize); 4] = [(1, 2), (2, 2), (3, 2), (4, 2)];

/// Index within [`FREE_CPS`] for a given `(i, j)`. Returns `None`
/// if the CP is not free.
fn free_cp_index(i: usize, j: usize) -> Option<usize> {
    FREE_CPS.iter().position(|&(a, b)| a == i && b == j)
}

/// Total scalar DOFs per coordinate: 5 (a's 5 Point3 coefs) + 4·n
/// (free CPs per sector, per [`FREE_CPS`]).
fn dof_count(n: usize) -> usize {
    5 + 4 * n
}

/// DOF index for a's `a_x` (0), `a_y` (1), `a_xx` (2), `a_xy` (3),
/// `a_yy` (4). 0..=4.
fn dof_a(a_idx: usize) -> usize {
    a_idx
}

/// DOF index for sector `k`'s free CP #`f` (0..4 per `FREE_CPS`).
fn dof_free(k: usize, f: usize) -> usize {
    5 + 4 * k + f
}

/// Compute the central pole via paper's eq (5) from the boundary
/// ring alone. The DS-net points `d^k_1, d^k_2, d^k_4` are
/// interpreted from the bi-2 tensor-border's inner row:
///
/// ```text
///   b^k_{i, 1} = strip.pos[i] + strip.deriv[i] / 2    (standard bi-2 derivative form)
///
///   d^k_1 = b^k_{1, 1}    (inner-row middle)
///   d^k_2 = b^k_{2, 1}    (inner-row right)
///   d^k_4 = b^k_{0, 1}    (inner-row left)
/// ```
///
/// Returns `pole = e·avg_k(d^k_1) + (1−e)·avg_k((d^k_2+d^k_4)/2)`
/// with `e = e_def_bi4(n)`.
pub fn pole_from_boundary(boundary: &Bi2TensorBorderRing) -> Point3 {
    let n = boundary.n();
    let e = e_def_bi4(n);
    let mut d1 = Vec::with_capacity(n);
    let mut d24_half = Vec::with_capacity(n);
    for strip in &boundary.sectors {
        // b^k_{i, 1} = pos[i] + deriv[i] / 2 (bi-2 derivative: row 1 = row 0 + deriv/2).
        let b01 = Point3::new(
            strip.pos[0].x + 0.5 * strip.deriv[0].x,
            strip.pos[0].y + 0.5 * strip.deriv[0].y,
            strip.pos[0].z + 0.5 * strip.deriv[0].z,
        );
        let b11 = Point3::new(
            strip.pos[1].x + 0.5 * strip.deriv[1].x,
            strip.pos[1].y + 0.5 * strip.deriv[1].y,
            strip.pos[1].z + 0.5 * strip.deriv[1].z,
        );
        let b21 = Point3::new(
            strip.pos[2].x + 0.5 * strip.deriv[2].x,
            strip.pos[2].y + 0.5 * strip.deriv[2].y,
            strip.pos[2].z + 0.5 * strip.deriv[2].z,
        );
        d1.push(b11);
        d24_half.push(Point3::new(
            0.5 * (b01.x + b21.x),
            0.5 * (b01.y + b21.y),
            0.5 * (b01.z + b21.z),
        ));
    }
    central_point(&d1, &d24_half, e)
}

/// Assemble a PnS bi-4 cap from boundary + pole + σ jet.
pub fn assemble_pns_bi4_cap(
    input: &CapAssembleInput,
) -> Result<CapAssembleResult, CapAssembleError> {
    let n = input.boundary.n();
    if !(3..=12).contains(&n) {
        return Err(CapAssembleError::InvalidValence(n));
    }
    // n=4 is numerically singular: c = cos(2π/4) = 0 appears in
    // σ_{20,x}'s `1/(3c)` factor. The paper's bi-4 construction
    // doesn't apply to n=4; callers should use the bi-3 variant
    // (`pns_bi3_twists`) or a different cap entirely.
    if n == 4 {
        return Err(CapAssembleError::N4NotSupported);
    }
    // Other even n (6, 8, 10, 12) work fine — my asymmetric
    // closed-form treatment (treating canonical (1, 2) as the
    // free DOF and (2, 1), (3, 1), (4, 1) as one-way derivations
    // from prev's data) avoids the rank-deficient symmetric
    // pair-sum issue that would otherwise plague even n.
    if !(1..=4).contains(&input.kappa) {
        return Err(CapAssembleError::InvalidKappa(input.kappa));
    }

    // Precompute per-sector σ^k and boundary (row0, row1).
    //
    // σ convention (refined 2026-04-24 per follow-up clarification):
    // σ^k = R_{kα} · σ^0 coefficient-by-coefficient (image-space
    // rotation). Each sector's BB-coefs genuinely differ — σ¹_{ij}
    // = R_α·σ⁰_{ij}. The "common recipe" phrase referred to the
    // FORMULAS in d₁, d₂, d₃ being identical per sector, not the
    // BB-coef values in ambient ℝ².
    let mut sigma_k = Vec::with_capacity(n);
    let mut row0 = Vec::with_capacity(n);
    let mut row1 = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        sigma_k.push(rotate_sigma_jet(&input.sigma, theta));
        let (r0, r1) = compose_bi4_outer_layers(&input.boundary.sectors[k], &input.beta);
        row0.push(r0);
        row1.push(r1);
    }

    // Build the coefficient map: for each sector k and each CP
    // (i, j), express it as a linear combination of DOFs plus a
    // constant Point3.
    let coef_map = build_coef_map(input, &sigma_k, &row0, &row1);

    // Assemble the DOF-space Gram matrix and RHS, solve per coord,
    // reconstruct patches.
    let dofs = solve_dof_system(&coef_map, input.kappa, n)?;
    let patches = reconstruct_patches(&coef_map, &dofs, n);
    let a = reconstruct_a(&dofs, input.pole);

    Ok(CapAssembleResult { a, patches })
}

/// Linear representation of a single BB-coef: `cp = Σ_d coef[d]·dof[d] + constant`.
#[derive(Debug, Clone)]
struct CpRepr {
    /// Scalar coefficients (same across x/y/z by separability).
    /// Length `dof_count(n)`.
    coef: Vec<f64>,
    /// Constant Point3 (contribution from pole + boundary).
    constant: Point3,
}

impl CpRepr {
    fn zeros(n: usize) -> Self {
        Self {
            coef: vec![0.0; dof_count(n)],
            constant: Point3::new(0.0, 0.0, 0.0),
        }
    }
}

/// Full coefficient map: `sector_cps[k][5·i + j]` gives the
/// representation for `P[k][i, j]`.
#[derive(Debug, Clone)]
struct CoefMap {
    n: usize,
    sector_cps: Vec<Vec<CpRepr>>,
}

fn build_coef_map(
    input: &CapAssembleInput,
    sigma_k: &[SigmaCenterJet],
    row0: &[[Point3; 5]],
    row1: &[[Point3; 5]],
) -> CoefMap {
    let n = input.boundary.n();
    let mut sector_cps: Vec<Vec<CpRepr>> = (0..n).map(|_| vec![CpRepr::zeros(n); 25]).collect();

    // Pass 1: fill in pole, central (a∘σ^k), boundary, free CPs.
    for k in 0..n {
        for i in 0..5 {
            for j in 0..5 {
                let repr = compute_non_shared_cp(input, sigma_k, row0, row1, k, i, j);
                if let Some(r) = repr {
                    sector_cps[k][bi4_flat_index(i, j)] = r;
                }
            }
        }
    }

    // Pass 2: prev-sharing on canonical u_c=4 edge. Cur's canonical
    // (4, v_c) = prev's canonical (0, v_c) via the cur-prev seam.
    // Derivation: for the pair (prev, cur) = (k-1, k), the shared
    // edge sits at prev's canonical u_c=0 and cur's canonical v_c=0
    // — BUT actually we want the OTHER seam. Let me re-derive.
    //
    // For pair (cur=k, nxt=k+1): seam sits at cur's canonical u_c=0.
    // Parametric identification along this seam: cur's (0, v) = nxt's (v, 0)
    // for v=0..4 (pole at v=0, outer end at v=4).
    //
    // Rewritten: patches[cur][0][v] = patches[nxt][v][0].
    // Equivalently (renaming): patches[k][u][0] = patches[prev][0][u].
    //
    // We pin cur's canonical v_c=0 edge from prev's u_c=0 edge.
    // (Except (0, 0) which is pole, already pinned.)
    for k in 0..n {
        let prev = (k + n - 1) % n;
        for u in 1..5 {
            // patches[k][u][0] = patches[prev][0][u]
            let src = sector_cps[prev][bi4_flat_index(0, u)].clone();
            sector_cps[k][bi4_flat_index(u, 0)] = src;
        }
    }

    // Pass 3: eq (2) closed-form for canonical (2, 1) per sector.
    // From the paper's pair-sum for pair (prev=k-1, cur=k):
    //   P[k][2, 1]_paper + P[k-1][2, 1]_paper
    //     = (2-c)·P[k-1][2, 0]_paper + c·P[k-1][3, 0]_paper
    // In canonical (applying cur-side swap for the "cur" side = k-1):
    //   canonical p^k_{2, 1} + canonical p^{k-1}_{1, 2}
    //     = (2-c)·p^{k-1}_{0, 2} + c·p^{k-1}_{0, 3}
    // Solve for canonical p^k_{2, 1}:
    //   p^k_{2, 1} = -p^{k-1}_{1, 2} + (2-c)·p^{k-1}_{0, 2} + c·p^{k-1}_{0, 3}
    // Additional closed-forms for canonical (3, 1) and (4, 1):
    //   cur_{40} (paper) = (cur_{41} + nxt_{41}) / 2
    //     → p^k_{4, 1} = 2·p^{k-1}_{0, 4} - p^{k-1}_{1, 4}
    //   cur_{30} (paper) = [4·(cur_{31} + nxt_{31})
    //                       - c·(cur_{41} + nxt_{41})] / (2·(4-c))
    //     → 2·(4-c)·p^{k-1}_{0, 3}
    //         = 4·(p^{k-1}_{1, 3} + p^k_{3, 1})
    //           - c·(p^{k-1}_{1, 4} + p^k_{4, 1})
    //     → p^k_{3, 1} = [2·(4-c)·p^{k-1}_{0, 3}
    //                    + c·(p^{k-1}_{1, 4} + p^k_{4, 1})
    //                    - 4·p^{k-1}_{1, 3}] / 4
    let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
    let two_minus_c = 2.0 - c;
    let four_minus_c = 4.0 - c;
    for k in 0..n {
        let prev = (k + n - 1) % n;
        // (4, 1): from prev's (0, 4), (1, 4) — pure boundary.
        let prev_04 = sector_cps[prev][bi4_flat_index(0, 4)].clone();
        let prev_14 = sector_cps[prev][bi4_flat_index(1, 4)].clone();
        let mut p_41 = CpRepr::zeros(n);
        for d in 0..p_41.coef.len() {
            p_41.coef[d] = 2.0 * prev_04.coef[d] - prev_14.coef[d];
        }
        p_41.constant = Point3::new(
            2.0 * prev_04.constant.x - prev_14.constant.x,
            2.0 * prev_04.constant.y - prev_14.constant.y,
            2.0 * prev_04.constant.z - prev_14.constant.z,
        );
        sector_cps[k][bi4_flat_index(4, 1)] = p_41.clone();

        // (3, 1): from prev's (0, 3), (1, 3), (1, 4) + just-
        // computed (4, 1).
        let prev_03 = sector_cps[prev][bi4_flat_index(0, 3)].clone();
        let prev_13 = sector_cps[prev][bi4_flat_index(1, 3)].clone();
        let mut p_31 = CpRepr::zeros(n);
        for d in 0..p_31.coef.len() {
            p_31.coef[d] = (2.0 * four_minus_c * prev_03.coef[d]
                + c * (prev_14.coef[d] + p_41.coef[d])
                - 4.0 * prev_13.coef[d])
                / 4.0;
        }
        p_31.constant = Point3::new(
            (2.0 * four_minus_c * prev_03.constant.x + c * (prev_14.constant.x + p_41.constant.x)
                - 4.0 * prev_13.constant.x)
                / 4.0,
            (2.0 * four_minus_c * prev_03.constant.y + c * (prev_14.constant.y + p_41.constant.y)
                - 4.0 * prev_13.constant.y)
                / 4.0,
            (2.0 * four_minus_c * prev_03.constant.z + c * (prev_14.constant.z + p_41.constant.z)
                - 4.0 * prev_13.constant.z)
                / 4.0,
        );
        sector_cps[k][bi4_flat_index(3, 1)] = p_31;

        // (2, 1): from prev's (1, 2) [free DOF], (0, 2) [central],
        // (0, 3) [boundary].
        let prev_12 = &sector_cps[prev][bi4_flat_index(1, 2)];
        let prev_02 = &sector_cps[prev][bi4_flat_index(0, 2)];
        let mut p_21 = CpRepr::zeros(n);
        for d in 0..p_21.coef.len() {
            p_21.coef[d] = -prev_12.coef[d] + two_minus_c * prev_02.coef[d] + c * prev_03.coef[d];
        }
        p_21.constant = Point3::new(
            -prev_12.constant.x + two_minus_c * prev_02.constant.x + c * prev_03.constant.x,
            -prev_12.constant.y + two_minus_c * prev_02.constant.y + c * prev_03.constant.y,
            -prev_12.constant.z + two_minus_c * prev_02.constant.z + c * prev_03.constant.z,
        );
        sector_cps[k][bi4_flat_index(2, 1)] = p_21;
    }

    CoefMap { n, sector_cps }
}

/// Compute the CANONICAL `(u_c, v_c)` representation for
/// `p^k_{u, v}` EXCLUDING prev-shared u_c=4 entries and eq (2)
/// closed-form canonical `(2, 1)` (both handled by later passes).
fn compute_non_shared_cp(
    input: &CapAssembleInput,
    sigma_k: &[SigmaCenterJet],
    row0: &[[Point3; 5]],
    row1: &[[Point3; 5]],
    k: usize,
    u: usize,
    v: usize,
) -> Option<CpRepr> {
    let n = input.boundary.n();

    // Deferred: prev-sharing on v_c=0 (pass 2, except pole at (0,0)),
    // eq (2) closed-forms for (2, 1), (3, 1), (4, 1) (pass 3).
    if v == 0 && u >= 1 {
        return None;
    }
    if matches!((u, v), (2, 1) | (3, 1) | (4, 1)) {
        return None;
    }

    let mut repr = CpRepr::zeros(n);

    // Pole at canonical (0, 0).
    if (u, v) == (0, 0) {
        repr.constant = input.pole;
        return Some(repr);
    }

    // Boundary from pns_boundary (paper i=4 row → canonical v_c=4;
    // paper i=3 row → canonical v_c=3). row0[u] = p^k_{u, 4},
    // row1[u] = 4·(p^k_{u, 3} - p^k_{u, 4}) so
    // p^k_{u, 3} = row0[u] + row1[u]/4.
    if v == 4 {
        repr.constant = row0[k][u];
        return Some(repr);
    }
    if v == 3 {
        repr.constant = Point3::new(
            row0[k][u].x + 0.25 * row1[k][u].x,
            row0[k][u].y + 0.25 * row1[k][u].y,
            row0[k][u].z + 0.25 * row1[k][u].z,
        );
        return Some(repr);
    }

    // Central (a∘σ^k) at canonical (u, v) with u+v ≤ 2 and (u, v) ≠ (0, 0).
    if u + v <= 2 {
        central_cp_repr(&mut repr, &sigma_k[k], input.pole, u, v);
        return Some(repr);
    }

    // Free CP.
    if let Some(f) = free_cp_index(u, v) {
        repr.coef[dof_free(k, f)] = 1.0;
        return Some(repr);
    }

    unreachable!(
        "build_coef_map: unexpected canonical CP ({u}, {v}) in sector {k} (u+v={}, free? {})",
        u + v,
        free_cp_index(u, v).is_some()
    );
}

/// Fill `repr` for a central CP `(i, j)` with `i+j ≤ 2` (excluding
/// the pole). Coefficients depend on `sigma_k` via the chain rule
/// in [`central_bb_from_a_and_sigma`]; constants = pole.
fn central_cp_repr(repr: &mut CpRepr, sigma: &SigmaCenterJet, pole: Point3, i: usize, j: usize) {
    // Use unit-basis evaluation: for each a-DOF, set that
    // component of `a` to (1, 0, 0), others zero, evaluate
    // central_bb_from_a_and_sigma; the resulting Point3.x = coef.
    //
    // (By separability the coef on a's x-component equals the coef
    // on y- and z-components, so we need only the x-component.)
    let zero = Point3::new(0.0, 0.0, 0.0);
    let basis = Point3::new(1.0, 0.0, 0.0);

    // Baseline (a = zero, a_0 = pole): central CP = pole.
    // Already handled as constant.
    repr.constant = pole;

    // a_x
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: basis,
            a_y: zero,
            a_xx: zero,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(0)] = cp.x - pole.x;
    }
    // a_y
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: basis,
            a_xx: zero,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(1)] = cp.x - pole.x;
    }
    // a_xx
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: basis,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(2)] = cp.x - pole.x;
    }
    // a_xy
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: zero,
            a_xy: basis,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(3)] = cp.x - pole.x;
    }
    // a_yy
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: zero,
            a_xy: zero,
            a_yy: basis,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(4)] = cp.x - pole.x;
    }
}

/// Extract central BB coef `(i, j)` from a `CentralSixBB`.
fn cp_at_ij(cbb: super::pns_central_quad::CentralSixBB, i: usize, j: usize) -> Point3 {
    match (i, j) {
        (0, 0) => cbb.p00,
        (1, 0) => cbb.p10,
        (0, 1) => cbb.p01,
        (2, 0) => cbb.p20,
        (1, 1) => cbb.p11,
        (0, 2) => cbb.p02,
        _ => panic!("cp_at_ij: not a central CP ({i}, {j})"),
    }
}

/// Solve the joint F_κ system in DOF space, per coordinate.
fn solve_dof_system(
    coef_map: &CoefMap,
    kappa: usize,
    n: usize,
) -> Result<Vec<[f64; 3]>, CapAssembleError> {
    let gram = gram_bi4_f_kappa(kappa);
    let n_dofs = dof_count(n);

    // Build G_dof = Σ_k M_k^T · G · M_k  (size n_dofs × n_dofs).
    // M_k is the (25) × n_dofs matrix of scalar coefficients for
    // sector k (using coef_map.sector_cps[k][idx].coef).
    //
    // Also build b_dof = Σ_k M_k^T · G · c_k (three components:
    // b_dof_x, b_dof_y, b_dof_z) where c_k is sector k's constants.
    let mut g_dof = vec![vec![0.0_f64; n_dofs]; n_dofs];
    let mut b_dof_x = vec![0.0_f64; n_dofs];
    let mut b_dof_y = vec![0.0_f64; n_dofs];
    let mut b_dof_z = vec![0.0_f64; n_dofs];

    for k in 0..coef_map.n {
        let cps = &coef_map.sector_cps[k];
        for a in 0..25 {
            for b in 0..25 {
                let gab = gram[a][b];
                if gab == 0.0 {
                    continue;
                }
                let ma = &cps[a];
                let mb = &cps[b];
                // M_k^T · G · M_k contribution:
                //   for each (d1, d2): G_dof[d1][d2] += ma.coef[d1] · gab · mb.coef[d2]
                for d1 in 0..n_dofs {
                    let ca = ma.coef[d1];
                    if ca == 0.0 {
                        continue;
                    }
                    for d2 in 0..n_dofs {
                        let cb = mb.coef[d2];
                        if cb == 0.0 {
                            continue;
                        }
                        g_dof[d1][d2] += ca * gab * cb;
                    }
                }
                // M_k^T · G · c_k contribution:
                //   b_dof[d] += ma.coef[d] · gab · mb.constant.?
                // Accumulate for each coord.
                for d in 0..n_dofs {
                    let ca = ma.coef[d];
                    if ca == 0.0 {
                        continue;
                    }
                    b_dof_x[d] += ca * gab * mb.constant.x;
                    b_dof_y[d] += ca * gab * mb.constant.y;
                    b_dof_z[d] += ca * gab * mb.constant.z;
                }
            }
        }
    }

    // Solve G_dof · x = -b_dof per coord. Gram is symmetric PSD;
    // if the free-DOF subspace has no kernel, solve_dense works.
    let mut dofs = vec![[0.0_f64; 3]; n_dofs];

    let xs =
        solve_dense(g_dof.clone(), negate(&b_dof_x)).ok_or(CapAssembleError::SingularSystem)?;
    let ys =
        solve_dense(g_dof.clone(), negate(&b_dof_y)).ok_or(CapAssembleError::SingularSystem)?;
    let zs = solve_dense(g_dof, negate(&b_dof_z)).ok_or(CapAssembleError::SingularSystem)?;
    for d in 0..n_dofs {
        dofs[d] = [xs[d], ys[d], zs[d]];
    }
    Ok(dofs)
}

fn negate(v: &[f64]) -> Vec<f64> {
    v.iter().map(|x| -x).collect()
}

/// Reconstruct the cap patches from solved DOF values.
fn reconstruct_patches(coef_map: &CoefMap, dofs: &[[f64; 3]], n: usize) -> Vec<[[Point3; 5]; 5]> {
    let mut patches = vec![[[Point3::new(0.0, 0.0, 0.0); 5]; 5]; n];
    for k in 0..n {
        for i in 0..5 {
            for j in 0..5 {
                let repr = &coef_map.sector_cps[k][bi4_flat_index(i, j)];
                let mut p = repr.constant;
                for d in 0..repr.coef.len() {
                    let c = repr.coef[d];
                    if c == 0.0 {
                        continue;
                    }
                    p.x += c * dofs[d][0];
                    p.y += c * dofs[d][1];
                    p.z += c * dofs[d][2];
                }
                patches[k][i][j] = p;
            }
        }
    }
    patches
}

/// Reconstruct `a` from the DOF solution (first 5 DOFs, per coord).
fn reconstruct_a(dofs: &[[f64; 3]], pole: Point3) -> QuadMapJet {
    let dof_to_p3 = |d: usize| Point3::new(dofs[d][0], dofs[d][1], dofs[d][2]);
    QuadMapJet {
        a_0: pole,
        a_x: dof_to_p3(dof_a(0)),
        a_y: dof_to_p3(dof_a(1)),
        a_xx: dof_to_p3(dof_a(2)),
        a_xy: dof_to_p3(dof_a(3)),
        a_yy: dof_to_p3(dof_a(4)),
    }
}

#[cfg(test)]
#[path = "pns_cap_assemble_tests.rs"]
mod tests;
