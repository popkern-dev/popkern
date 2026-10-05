//! Reference-grade Catmull-Clark characteristic ring.
//!
//! Constructs, for a valence-`n` extraordinary point, the exact
//! frequency-one characteristic map data of Catmull-Clark subdivision:
//! the extended control net, one subdivided net, and the three regular
//! 4×4 bicubic B-spline control grids of one characteristic-ring
//! sector. This is the source substrate required by KP20 §4.3 (the
//! reparameterization τ is G²-connected to "the tensor-border of the
//! characteristic ring of Catmull-Clark subdivision") and by KP18's
//! guide fit against corner jets of `g ∘ χ̃`.
//!
//! Provenance and verification. The construction follows Stam 1998
//! (`docs/papers/sig98.pdf`, extended subdivision matrices and picking
//! vectors) with the ring seed in closed form (PR08 §6.1 frequency-one
//! eigenstructure). It is pinned against an independent Python rebuild
//! from mesh connectivity (`tools/reference/cc-characteristic/`) that
//! also cross-checks PR08's characteristic-ring patch grids to ≤4e-15.
//! Two unpublished misprints in PR08 (p.112 DFT row k′=4; p.113
//! eigenvector table) make transcribing that book unsafe; this module
//! never transcribes the misprinted objects. See
//! `docs/architecture/kp20-paper-faithful-recovery-plan.md`, WP1/RC5.
//!
//! Conventions (Stam Fig. 1/3/5):
//! - Control-vertex labels are 1-based; label `1` is the extraordinary
//!   vertex. The one-ring is interleaved: even labels `2j` are
//!   edge-adjacent (spoke) vertices, odd labels `2j+1` are
//!   face-diagonal vertices, `j = 1..n`.
//! - The evaluated sector face is `{1, 4, 5, 6}`; `u` runs from
//!   vertex 1 toward 6, `v` from 1 toward 4.
//! - Picked grids are stored `grid[r][c]` with `r` the v-row and `c`
//!   the u-column (Stam Fig. 1 order: entry m = 1 + c + 4r).
//! - Net values are complex (`[re, im]`); sector `j` of the ring is
//!   the sector-0 data rotated by `e^{2πij/n}`.

use super::direct_linalg::solve_dense;
use popkern_math::PortableMath;

/// Complex scalar stored as `[re, im]`.
pub type CcComplex = [f64; 2];

/// Exact characteristic-ring data for one extraordinary valence.
#[derive(Debug, Clone, PartialEq)]
pub struct CcCharacteristicRing {
    /// Valence (`n >= 3`; `n = 4` is the regular sanity case).
    pub n: usize,
    /// Subdominant eigenvalue `λ = (c + 5 + √((c+9)(c+1)))/16`.
    pub lambda: f64,
    /// Extended level-0 net `C` (`2n+8` values; Stam label `k` at
    /// index `k−1`). Satisfies `A·C = λ·C` exactly.
    pub extended_net: Vec<CcComplex>,
    /// Once-subdivided net `C̄ = Ā·C` (`2n+17` values).
    pub refined_net: Vec<CcComplex>,
    /// The three picked 4×4 B-spline grids of the level-1 ring
    /// (Stam subpatches Ω₁, Ω₂, Ω₃ = PR08 ring patches C, B, A).
    pub grids: [[[CcComplex; 4]; 4]; 3],
    /// PR08 normalization value `f(1,1,0)`: the ring-surface value at
    /// the outer corner of the Ω₂ grid.
    pub f110: CcComplex,
}

/// Catmull-Clark subdominant eigenvalue for valence `n`.
pub fn cc_subdominant_lambda(n: usize) -> Result<f64, String> {
    if n < 3 {
        return Err(format!("CC characteristic ring requires n >= 3; got {n}"));
    }
    let c = (core::f64::consts::TAU / n as f64).pcos();
    Ok((c + 5.0 + ((c + 9.0) * (c + 1.0)).sqrt()) / 16.0)
}

/// Number of control vertices in the level-0 extended net.
pub fn cc_extended_net_len(n: usize) -> usize {
    2 * n + 8
}

/// Number of control vertices after one subdivision.
pub fn cc_refined_net_len(n: usize) -> usize {
    2 * n + 17
}

/// Build the exact characteristic ring for valence `n`.
pub fn cc_characteristic_ring(n: usize) -> Result<CcCharacteristicRing, String> {
    let lambda = cc_subdominant_lambda(n)?;
    let ring = cc_ring_seed(n, lambda);
    let outer = solve_outer_extension(n, lambda, &ring)?;

    let mut extended_net = ring;
    extended_net.extend_from_slice(&outer);

    let a_bar = cc_extended_subdivision_rows(n);
    let refined_net: Vec<CcComplex> = a_bar
        .iter()
        .map(|row| {
            let mut acc = [0.0_f64; 2];
            for (weight, value) in row.iter().zip(&extended_net) {
                acc[0] += weight * value[0];
                acc[1] += weight * value[1];
            }
            acc
        })
        .collect();

    let picks = cc_picking_vectors(n);
    let mut grids = [[[[0.0_f64; 2]; 4]; 4]; 3];
    for (grid, pick) in grids.iter_mut().zip(&picks) {
        for (m, &label) in pick.iter().enumerate() {
            let r = m / 4;
            let c = m % 4;
            grid[r][c] = refined_net[label - 1];
        }
    }

    // f(1,1,0): bicubic B-spline value at the outer corner of the Ω₂
    // grid — tensor product of the uniform-cubic endpoint mask
    // [1/6, 2/3, 1/6] over the 3×3 control points nearest grid corner
    // (3,3).
    let mask = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
    let mut f110 = [0.0_f64; 2];
    for (r, wr) in mask.iter().enumerate() {
        for (c, wc) in mask.iter().enumerate() {
            let value = grids[1][r + 1][c + 1];
            f110[0] += wr * wc * value[0];
            f110[1] += wr * wc * value[1];
        }
    }

    Ok(CcCharacteristicRing {
        n,
        lambda,
        extended_net,
        refined_net,
        grids,
        f110,
    })
}

/// Closed-form frequency-one ring seed (2n+1 values): the
/// extraordinary vertex maps to 0, edge vertex `2j` to
/// `(4λ−1)·e^{2πi(j−1)/n}`, face vertex `2j+1` to
/// `(1+e^{2πi/n})·e^{2πi(j−1)/n}`. No eigensolver is involved, so the
/// result is deterministic and phase-anchored at spoke vertex 2.
fn cc_ring_seed(n: usize, lambda: f64) -> Vec<CcComplex> {
    let alpha = core::f64::consts::TAU / n as f64;
    let edge_mag = 4.0 * lambda - 1.0;
    let face = [1.0 + alpha.pcos(), alpha.psin()];
    let mut seed = vec![[0.0_f64; 2]; 2 * n + 1];
    for j in 1..=n {
        let phase = (j as f64 - 1.0) * alpha;
        let rot = [phase.pcos(), phase.psin()];
        seed[2 * j - 1] = [edge_mag * rot[0], edge_mag * rot[1]];
        seed[2 * j] = complex_mul(face, rot);
    }
    seed
}

fn complex_mul(a: CcComplex, b: CcComplex) -> CcComplex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

/// Extend the ring seed to the seven outer control vertices by solving
/// `(S₁₂ − λI)·x = −S₁₁·ring` (the unique extension making the whole
/// extended net a λ-eigenvector of the extended subdivision matrix;
/// `S₁₂ − λI` is invertible because `λ > 1/8 ≥ spec(S₁₂)`).
fn solve_outer_extension(
    n: usize,
    lambda: f64,
    ring: &[CcComplex],
) -> Result<[CcComplex; 7], String> {
    let s11 = cc_s11_rows(n);
    let s12 = cc_s12_rows();

    let mut matrix = vec![vec![0.0_f64; 7]; 7];
    for (r, row) in matrix.iter_mut().enumerate() {
        for (c, entry) in row.iter_mut().enumerate() {
            *entry = s12[r][c] - if r == c { lambda } else { 0.0 };
        }
    }

    let mut solution = [[0.0_f64; 2]; 7];
    for axis in 0..2 {
        let rhs: Vec<f64> = s11
            .iter()
            .map(|row| {
                -row.iter()
                    .zip(ring)
                    .map(|(weight, value)| weight * value[axis])
                    .sum::<f64>()
            })
            .collect();
        let axis_solution = solve_dense(matrix.clone(), rhs).ok_or_else(|| {
            format!("CC outer-extension system is singular for n = {n} (unexpected: λ > 1/8)")
        })?;
        for (slot, value) in solution.iter_mut().zip(axis_solution) {
            slot[axis] = value;
        }
    }
    Ok(solution)
}

const D: f64 = 3.0 / 8.0;
const E: f64 = 1.0 / 16.0;
const F: f64 = 1.0 / 4.0;
const VA: f64 = 9.0 / 16.0;
const VB: f64 = 3.0 / 32.0;
const VC: f64 = 1.0 / 64.0;

/// Stam `S` (the (2n+1)×(2n+1) Catmull-Clark matrix around the
/// extraordinary vertex), built from the subdivision stencils applied
/// to the labeled star — not transcribed from a printed matrix.
fn cc_star_subdivision_rows(n: usize) -> Vec<Vec<f64>> {
    let k1 = 2 * n + 1;
    // Ring labels 2..=2n+1 are cyclic in j.
    let idx = |label: isize| -> usize {
        let folded = (label - 2).rem_euclid(2 * n as isize) + 2;
        folded as usize - 1
    };
    let mut s = vec![vec![0.0_f64; k1]; k1];

    // Extraordinary vertex row: a_n on itself, b_n on each edge
    // neighbor, c_n on each face neighbor (Stam Appendix A; the 1/n²
    // absorbs the β/n, γ/n of the vertex stencil).
    let an = 1.0 - 7.0 / (4.0 * n as f64);
    let bn = 3.0 / (2.0 * (n as f64) * (n as f64));
    let cn = 1.0 / (4.0 * (n as f64) * (n as f64));
    s[0][0] = an;
    for j in 1..=n as isize {
        s[0][idx(2 * j)] += bn;
        s[0][idx(2 * j + 1)] += cn;
    }

    for j in 1..=n as isize {
        // New edge point of spoke {1, 2j}.
        let row = idx(2 * j);
        s[row][0] += D;
        s[row][idx(2 * j)] += D;
        for v in [2 * j - 2, 2 * j - 1, 2 * j + 1, 2 * j + 2] {
            s[row][idx(v)] += E;
        }
        // New face point of sector face {1, 2j, 2j+1, 2j+2}.
        let row = idx(2 * j + 1);
        s[row][0] += F;
        for v in [2 * j, 2 * j + 1, 2 * j + 2] {
            s[row][idx(v)] += F;
        }
    }
    s
}

/// Stam `S₁₁` (7 × (2n+1)): the new one-ring-exterior vertices
/// `2n+2..2n+8` as stencils on the old star. Verified against mesh
/// connectivity by the reference oracle; `n = 3` folds label 8 onto
/// label 2 (`c₈ ≡ c₂`).
fn cc_s11_rows(n: usize) -> Vec<Vec<f64>> {
    let k1 = 2 * n + 1;
    let rows: [[f64; 8]; 7] = [
        [VC, 0.0, 0.0, VB, VA, VB, 0.0, 0.0],
        [E, 0.0, 0.0, E, D, D, 0.0, 0.0],
        [VB, 0.0, 0.0, VC, VB, VA, VB, VC],
        [E, 0.0, 0.0, 0.0, 0.0, D, D, E],
        [E, 0.0, 0.0, D, D, E, 0.0, 0.0],
        [VB, VC, VB, VA, VB, VC, 0.0, 0.0],
        [E, E, D, D, 0.0, 0.0, 0.0, 0.0],
    ];
    let mut s11 = vec![vec![0.0_f64; k1]; 7];
    for (out, row) in s11.iter_mut().zip(rows) {
        for (label0, weight) in row.iter().enumerate() {
            if *weight == 0.0 {
                continue;
            }
            let mut column = label0; // label = label0 + 1
            if n == 3 && label0 == 7 {
                column = 1; // c₈ ≡ c₂ at valence 3
            }
            out[column] += weight;
        }
    }
    s11
}

/// Stam `S₁₂` (7 × 7): the same new vertices on the old outer seven.
fn cc_s12_rows() -> [[f64; 7]; 7] {
    [
        [VC, VB, VC, 0.0, VB, VC, 0.0],
        [0.0, E, E, 0.0, 0.0, 0.0, 0.0],
        [0.0, VC, VB, VC, 0.0, 0.0, 0.0],
        [0.0, 0.0, E, E, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, E, E, 0.0],
        [0.0, 0.0, 0.0, 0.0, VC, VB, VC],
        [0.0, 0.0, 0.0, 0.0, 0.0, E, E],
    ]
}

/// Stam `S₂₁` (9 × (2n+1)) and `S₂₂` (9 × 7): the nine additional
/// level-1 vertices `2n+9..2n+17` (plain regular midpoint stencils;
/// note the printed S₂₁ has only seven explicit columns).
fn cc_s21_s22_rows(n: usize) -> (Vec<Vec<f64>>, [[f64; 7]; 9]) {
    let k1 = 2 * n + 1;
    let s21_cols: [[f64; 7]; 9] = [
        [0.0, 0.0, 0.0, 0.0, F, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, D, E, 0.0],
        [0.0, 0.0, 0.0, 0.0, F, F, 0.0],
        [0.0, 0.0, 0.0, 0.0, E, D, E],
        [0.0, 0.0, 0.0, 0.0, 0.0, F, F],
        [0.0, 0.0, 0.0, E, D, 0.0, 0.0],
        [0.0, 0.0, 0.0, F, F, 0.0, 0.0],
        [0.0, 0.0, E, D, E, 0.0, 0.0],
        [0.0, 0.0, F, F, 0.0, 0.0, 0.0],
    ];
    let mut s21 = vec![vec![0.0_f64; k1]; 9];
    for (out, row) in s21.iter_mut().zip(s21_cols) {
        for (label0, weight) in row.iter().enumerate() {
            out[label0] = *weight;
        }
    }
    let s22: [[f64; 7]; 9] = [
        [F, F, 0.0, 0.0, F, 0.0, 0.0],
        [E, D, E, 0.0, E, 0.0, 0.0],
        [0.0, F, F, 0.0, 0.0, 0.0, 0.0],
        [0.0, E, D, E, 0.0, 0.0, 0.0],
        [0.0, 0.0, F, F, 0.0, 0.0, 0.0],
        [E, E, 0.0, 0.0, D, E, 0.0],
        [0.0, 0.0, 0.0, 0.0, F, F, 0.0],
        [0.0, 0.0, 0.0, 0.0, E, D, E],
        [0.0, 0.0, 0.0, 0.0, 0.0, F, F],
    ];
    (s21, s22)
}

/// The extended level-0 subdivision matrix `A = [[S, 0], [S₁₁, S₁₂]]`
/// ((2n+8) × (2n+8)).
pub fn cc_level0_subdivision_rows(n: usize) -> Vec<Vec<f64>> {
    let k = cc_extended_net_len(n);
    let k1 = 2 * n + 1;
    let mut a = vec![vec![0.0_f64; k]; k];
    for (out, row) in a.iter_mut().zip(cc_star_subdivision_rows(n)) {
        out[..k1].copy_from_slice(&row);
    }
    let s11 = cc_s11_rows(n);
    let s12 = cc_s12_rows();
    for r in 0..7 {
        a[k1 + r][..k1].copy_from_slice(&s11[r]);
        a[k1 + r][k1..].copy_from_slice(&s12[r]);
    }
    a
}

/// The full one-step matrix `Ā = [[S, 0], [S₁₁, S₁₂], [S₂₁, S₂₂]]`
/// ((2n+17) × (2n+8)).
pub fn cc_extended_subdivision_rows(n: usize) -> Vec<Vec<f64>> {
    let k = cc_extended_net_len(n);
    let k1 = 2 * n + 1;
    let mut a_bar = cc_level0_subdivision_rows(n);
    let (s21, s22) = cc_s21_s22_rows(n);
    for r in 0..9 {
        let mut row = vec![0.0_f64; k];
        row[..k1].copy_from_slice(&s21[r]);
        row[k1..].copy_from_slice(&s22[r]);
        a_bar.push(row);
    }
    a_bar
}

/// Stam's three picking vectors (1-based labels into the refined net,
/// Fig. 1 row-major order `m = 1 + c + 4r`). `n = 3` folds the first
/// entry of `q¹` onto label 2.
pub fn cc_picking_vectors(n: usize) -> [[usize; 16]; 3] {
    let b = 2 * n; // labels 2n+k written as b + k
    let q1_first = if n == 3 { 2 } else { 8 };
    [
        [
            q1_first,
            7,
            b + 5,
            b + 13,
            1,
            6,
            b + 4,
            b + 12,
            4,
            5,
            b + 3,
            b + 11,
            b + 7,
            b + 6,
            b + 2,
            b + 10,
        ],
        [
            1,
            6,
            b + 4,
            b + 12,
            4,
            5,
            b + 3,
            b + 11,
            b + 7,
            b + 6,
            b + 2,
            b + 10,
            b + 16,
            b + 15,
            b + 14,
            b + 9,
        ],
        [
            2,
            1,
            6,
            b + 4,
            3,
            4,
            5,
            b + 3,
            b + 8,
            b + 7,
            b + 6,
            b + 2,
            b + 17,
            b + 16,
            b + 15,
            b + 14,
        ],
    ]
}

#[cfg(test)]
#[path = "cc_characteristic_ring_tests.rs"]
mod tests;
