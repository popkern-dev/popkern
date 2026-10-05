//! Exact KP20 §4.4 step (iii) sampling: `h := g ∘ τ` converted to
//! bi-6 pieces by full Fig-9 corner Hermite blocks.
//!
//! Per sector and quadrant piece, the (3,3)-jet of `h` is computed at
//! each of the four piece corners EXACTLY: the planar τ piece's
//! corner Taylor (from BB differences, all sixteen `∂^p_u ∂^q_v`,
//! p, q ≤ 3) is composed with the total-degree-6 guide through
//! truncated bivariate Taylor-jet arithmetic (closed under (3,3)
//! truncation), then converted with the verified
//! `hermite_to_bi6_corner_block` and merged with averaged
//! coefficients at overlapping locations (Fig 9c). Because a bi-6
//! polynomial is exactly reconstructed from its four (3,3) corner
//! blocks, planar inputs whose guide is the identity reproduce τ's
//! pieces exactly — the module gate.
//!
//! Charts: τ pieces in the WP2 quadrant convention (piece indices
//! globally aligned, P1 owns `o`); the guide in its own chart
//! (guide-u along spoke `k`, guide-v along spoke `k+1` — the
//! transpose of the τ-local convention, see `kp20_paper_guide`).

use super::kp20_bi6_cap::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use super::kp20_bi6_math::hermite_to_bi6_corner_block;
use crate::geometry::Point3;
use popkern_math::PortableMath;

/// Truncated bivariate Taylor jet through order (3, 3): entry
/// `[p][q]` is the Taylor coefficient of `s^p t^q`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Jet33(pub(crate) [[f64; 4]; 4]);

impl Jet33 {
    pub(crate) fn constant(value: f64) -> Self {
        let mut out = [[0.0_f64; 4]; 4];
        out[0][0] = value;
        Jet33(out)
    }

    pub(crate) fn mul(&self, other: &Jet33) -> Jet33 {
        let mut out = [[0.0_f64; 4]; 4];
        for p1 in 0..4 {
            for q1 in 0..4 {
                if self.0[p1][q1] == 0.0 {
                    continue;
                }
                for p2 in 0..(4 - p1) {
                    for q2 in 0..(4 - q1) {
                        out[p1 + p2][q1 + q2] += self.0[p1][q1] * other.0[p2][q2];
                    }
                }
            }
        }
        Jet33(out)
    }

    pub(crate) fn add_scaled(&mut self, other: &Jet33, scale: f64) {
        for p in 0..4 {
            for q in 0..4 {
                self.0[p][q] += scale * other.0[p][q];
            }
        }
    }

    /// Derivatives `h[p][q] = ∂^p_s ∂^q_t` from the coefficients.
    pub(crate) fn derivatives(&self) -> [[f64; 4]; 4] {
        let factorial = [1.0, 1.0, 2.0, 6.0];
        let mut out = [[0.0_f64; 4]; 4];
        for p in 0..4 {
            for q in 0..4 {
                out[p][q] = self.0[p][q] * factorial[p] * factorial[q];
            }
        }
        out
    }
}

const FACTORIAL: [f64; 4] = [1.0, 1.0, 2.0, 6.0];

/// All sixteen corner partials `∂^p_u ∂^q_v P` (p, q ≤ 3) of one
/// bi-6 piece component at the corner selected by `(ci, cj) ∈ {0,6}²`,
/// with both parameters pointing INTO the piece.
fn bi6_corner_jet33(piece: &Kp20Bi6Patch, ci: usize, cj: usize, axis: usize) -> Jet33 {
    let falling = |order: usize| -> f64 { (0..order).map(|k| (6 - k) as f64).product() };
    let binom = |order: usize, k: usize| -> f64 {
        let mut value = 1.0;
        for step in 0..k {
            value = value * (order - step) as f64 / (step + 1) as f64;
        }
        value
    };
    let read = |i: usize, j: usize| -> f64 {
        let si = if ci == 0 { i } else { 6 - i };
        let sj = if cj == 0 { j } else { 6 - j };
        let cp = piece.cps[si][sj];
        match axis {
            0 => cp.x,
            1 => cp.y,
            _ => cp.z,
        }
    };
    let mut out = [[0.0_f64; 4]; 4];
    for (a, row) in out.iter_mut().enumerate() {
        for (b, slot) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for p in 0..=a {
                for q in 0..=b {
                    let sign = if (a - p + b - q) % 2 == 0 { 1.0 } else { -1.0 };
                    acc += sign * binom(a, p) * binom(b, q) * read(p, q);
                }
            }
            // Taylor coefficient = derivative / (a! b!).
            *slot = acc * falling(a) * falling(b) / (FACTORIAL[a] * FACTORIAL[b]);
        }
    }
    Jet33(out)
}

/// Guide-chart conversion for one sector: global plane → (guide-u,
/// guide-v) with guide-u along spoke `k`, guide-v along spoke `k+1`.
fn guide_local(n: usize, sector: usize, global: [f64; 2]) -> [f64; 2] {
    let alpha = core::f64::consts::TAU / n as f64;
    let a0 = sector as f64 * alpha;
    let a1 = (sector + 1) as f64 * alpha;
    let spoke0 = [a0.pcos(), a0.psin()];
    let spoke1 = [a1.pcos(), a1.psin()];
    let det = spoke1[0] * spoke0[1] - spoke1[1] * spoke0[0];
    [
        (spoke1[0] * global[1] - spoke1[1] * global[0]) / det,
        (global[0] * spoke0[1] - global[1] * spoke0[0]) / det,
    ]
}

/// Compose the guide with one τ-piece corner jet: returns the
/// composed (3,3) derivative block per axis, exact.
fn composed_corner_block(
    n: usize,
    sector: usize,
    guide: &Kp20Tri6GuideSector,
    tau_piece: &Kp20Bi6Patch,
    ci: usize,
    cj: usize,
) -> [[Point3; 4]; 4] {
    // τ's planar corner Taylor in the guide chart.
    let tau_x = bi6_corner_jet33(tau_piece, ci, cj, 0);
    let tau_y = bi6_corner_jet33(tau_piece, ci, cj, 1);
    let mut u_jet = Jet33([[0.0; 4]; 4]);
    let mut v_jet = Jet33([[0.0; 4]; 4]);
    for p in 0..4 {
        for q in 0..4 {
            let local = guide_local(n, sector, [tau_x.0[p][q], tau_y.0[p][q]]);
            u_jet.0[p][q] = local[0];
            v_jet.0[p][q] = local[1];
        }
    }
    let mut w_jet = Jet33::constant(1.0);
    w_jet.add_scaled(&u_jet, -1.0);
    w_jet.add_scaled(&v_jet, -1.0);

    let powers = |base: &Jet33| -> [Jet33; 7] {
        let mut out = [Jet33::constant(1.0); 7];
        for k in 1..7 {
            out[k] = out[k - 1].mul(base);
        }
        out
    };
    let u_pow = powers(&u_jet);
    let v_pow = powers(&v_jet);
    let w_pow = powers(&w_jet);

    let factorial = |k: usize| -> f64 { (1..=k).map(|x| x as f64).product::<f64>().max(1.0) };
    let mut h = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    let mut accum = [
        Jet33::constant(0.0),
        Jet33::constant(0.0),
        Jet33::constant(0.0),
    ];
    for q_idx in 0..=6_usize {
        for r_idx in 0..=(6 - q_idx) {
            let cp = super::tri6_flat_index(6, q_idx, r_idx);
            let coeff =
                factorial(6) / (factorial(q_idx) * factorial(r_idx) * factorial(6 - q_idx - r_idx));
            let basis = u_pow[q_idx]
                .mul(&v_pow[r_idx])
                .mul(&w_pow[6 - q_idx - r_idx]);
            let value = guide.cps[cp];
            accum[0].add_scaled(&basis, coeff * value.x);
            accum[1].add_scaled(&basis, coeff * value.y);
            accum[2].add_scaled(&basis, coeff * value.z);
        }
    }
    let dx = accum[0].derivatives();
    let dy = accum[1].derivatives();
    let dz = accum[2].derivatives();
    for p in 0..4 {
        for q in 0..4 {
            h[p][q] = Point3::new(dx[p][q], dy[p][q], dz[p][q]);
        }
    }
    h
}

/// Sample one sector's 2×2 bi-6 macro patch from `g ∘ τ` (KP20 §4.4
/// step (iii), Fig 9): per piece, four exact corner Hermite blocks
/// merged with averaged coefficients at overlaps.
pub fn kp20_sample_paper_cap_sector(
    n: usize,
    sector: usize,
    guide: &Kp20Tri6GuideSector,
    tau_sector: &Kp20MacroPatch,
) -> Kp20MacroPatch {
    let mut out = Kp20MacroPatch::zero();
    for piece_id in Kp20PieceId::all() {
        let tau_piece = tau_sector.piece(piece_id);
        let mut sums = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
        let mut counts = [[0.0_f64; 7]; 7];
        for &(ci, cj) in &[(0_usize, 0_usize), (0, 6), (6, 0), (6, 6)] {
            let h = composed_corner_block(n, sector, guide, tau_piece, ci, cj);
            let mut block = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
            for axis in 0..3_usize {
                let mut scalar = [[0.0_f64; 4]; 4];
                for p in 0..4 {
                    for q in 0..4 {
                        scalar[p][q] = match axis {
                            0 => h[p][q].x,
                            1 => h[p][q].y,
                            _ => h[p][q].z,
                        };
                    }
                }
                let converted = hermite_to_bi6_corner_block(scalar);
                for (a, row) in converted.iter().enumerate() {
                    for (b, value) in row.iter().enumerate() {
                        match axis {
                            0 => block[a][b].x = *value,
                            1 => block[a][b].y = *value,
                            _ => block[a][b].z = *value,
                        }
                    }
                }
            }
            for (a, block_row) in block.iter().enumerate() {
                for (b, value) in block_row.iter().enumerate() {
                    let i = if ci == 0 { a } else { 6 - a };
                    let j = if cj == 0 { b } else { 6 - b };
                    sums[i][j] = Point3::new(
                        sums[i][j].x + value.x,
                        sums[i][j].y + value.y,
                        sums[i][j].z + value.z,
                    );
                    counts[i][j] += 1.0;
                }
            }
        }
        let piece = out.piece_mut(piece_id);
        for i in 0..7 {
            for j in 0..7 {
                let count = counts[i][j].max(1.0);
                piece.cps[i][j] = Point3::new(
                    sums[i][j].x / count,
                    sums[i][j].y / count,
                    sums[i][j].z / count,
                );
            }
        }
    }
    out
}

#[cfg(test)]
#[path = "kp20_paper_cap_sample_tests.rs"]
mod tests;
