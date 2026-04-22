//! KP13 § 5 F_1 functional for the bi-5 cap (I.14.e math core).
//!
//! Paper eq (16) defines the smoothness functional
//!
//! ```text
//!   F_κ f := ∫₀¹ ∫₀¹ Σ_{i+j=κ} (κ! / (i! j!)) · (∂^i_u ∂^j_v f)² du dv
//! ```
//!
//! For `κ = 1` this reduces to the classical Dirichlet / membrane
//! energy
//!
//! ```text
//!   F_1 f = ∫₀¹ ∫₀¹ (∂_u f)² + (∂_v f)² du dv
//! ```
//!
//! Paper eq (21) minimises `F_1(p − a)` over the free DoFs of the
//! final cap `p`, with `a` the auxiliary cap from I.14.a–d.
//!
//! # Evaluation strategy
//!
//! We use tensor-product Gauss–Legendre quadrature with 7 nodes
//! per axis. `∂_u p` and `∂_v p` for a bi-5 patch are polynomials
//! of bi-degree `(4, 5)` and `(5, 4)` respectively; squaring gives
//! bi-degree `(8, 10)` and `(10, 8)`. 7-point Gauss-Legendre is
//! exact for polynomials of degree `≤ 2 · 7 − 1 = 13` in one
//! direction — comfortably above our needs.
//!
//! # F_1 on a Bi5Cap
//!
//! A cap has `n` sectors, each with four `1/2 × 1/2` bi-5 sub-
//! patches (`O, L, R, M`). For each sub-patch with CP grid `P`,
//!
//! ```text
//!   F_1(P, scale_u, scale_v) = ∫₀¹ ∫₀¹ scale_u² (∂_u P)² + scale_v² (∂_v P)² du' dv'
//! ```
//!
//! where `(u', v')` are the piece-local coords and `scale_u =
//! scale_v = 2` (piece-local parameter advances 2× as fast as
//! sector-local). The total `F_1(cap) = Σ_{k, s} F_1(piece^{s,k})
//! · 1/(2·2)` — with the `1/4` Jacobian accounting for each
//! piece's `1/4` of the sector's parameter area.
//!
//! Net: `F_1(cap) = Σ F_1_piece(P)` where `F_1_piece` uses the
//! piece-local `(u', v')` parameterisation directly (the `2·2 ·
//! 1/4 = 1` scalings cancel).

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::bernstein_bi5::{bernstein_bi5_all, bernstein_bi5_deriv_all};
use super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, CapPieceId};

/// Gauss–Legendre nodes on `[0, 1]` (7-point).
/// Derived by mapping the standard `[-1, 1]` nodes via `t = (x +
/// 1) / 2`.
const GL7_NODES: [f64; 7] = [
    0.025_446_043_828_620_74,
    0.129_234_407_200_302_8,
    0.297_077_424_311_301_4,
    0.5,
    0.702_922_575_688_698_6,
    0.870_765_592_799_697_2,
    0.974_553_956_171_379_3,
];

/// Gauss–Legendre weights on `[0, 1]` (7-point). Sum to 1.
const GL7_WEIGHTS: [f64; 7] = [
    0.064_742_483_084_434_85,
    0.139_852_695_744_638_3,
    0.190_915_025_252_559_5,
    0.208_979_591_836_734_7,
    0.190_915_025_252_559_5,
    0.139_852_695_744_638_3,
    0.064_742_483_084_434_85,
];

/// `F_1` of a single bi-5 piece, over its piece-local `[0, 1]²`.
/// Returns the **sum of squared partials** for all three axes:
/// `F_1 = ∫∫ (∂u p)² + (∂v p)² = Σ_axis (axis-only F_1)`.
pub fn f1_piece(piece: &Bi5CapPiece) -> f64 {
    let mut total = 0.0_f64;
    for (ui, &u) in GL7_NODES.iter().enumerate() {
        let bu = bernstein_bi5_all(u);
        let bu_d = bernstein_bi5_deriv_all(u);
        for (vi, &v) in GL7_NODES.iter().enumerate() {
            let bv = bernstein_bi5_all(v);
            let bv_d = bernstein_bi5_deriv_all(v);

            // Evaluate ∂u p and ∂v p at (u, v).
            let mut du_x = 0.0_f64;
            let mut du_y = 0.0_f64;
            let mut du_z = 0.0_f64;
            let mut dv_x = 0.0_f64;
            let mut dv_y = 0.0_f64;
            let mut dv_z = 0.0_f64;
            for p in 0..6 {
                for q in 0..6 {
                    let cp = piece.cps[p][q];
                    let du_w = bu_d[p] * bv[q];
                    let dv_w = bu[p] * bv_d[q];
                    du_x += du_w * cp.x;
                    du_y += du_w * cp.y;
                    du_z += du_w * cp.z;
                    dv_x += dv_w * cp.x;
                    dv_y += dv_w * cp.y;
                    dv_z += dv_w * cp.z;
                }
            }
            let sq_du = du_x * du_x + du_y * du_y + du_z * du_z;
            let sq_dv = dv_x * dv_x + dv_y * dv_y + dv_z * dv_z;
            let w = GL7_WEIGHTS[ui] * GL7_WEIGHTS[vi];
            total += w * (sq_du + sq_dv);
        }
    }
    total
}

/// `F_1` of the whole cap: sum of `f1_piece` across every sector
/// and every sub-patch.
pub fn f1_cap(cap: &Bi5Cap) -> f64 {
    let mut total = 0.0_f64;
    for k in 0..cap.n {
        for id in CapPieceId::all() {
            total += f1_piece(cap.sector(k).piece(id));
        }
    }
    total
}

/// `F_1(p − a)`: the paper's minimisation target. Computes the
/// CP-wise difference `p − a` as a `Bi5Cap` and then applies
/// [`f1_cap`].
pub fn f1_difference(p: &Bi5Cap, a: &Bi5Cap) -> f64 {
    assert_eq!(p.n, a.n, "f1_difference: n mismatch");
    let mut diff = Bi5Cap::zero(p.n);
    for k in 0..p.n {
        for id in CapPieceId::all() {
            let p_piece = p.sector(k).piece(id);
            let a_piece = a.sector(k).piece(id);
            let d_piece = diff.sector_mut(k).piece_mut(id);
            for i in 0..6 {
                for j in 0..6 {
                    let cp_p = p_piece.cps[i][j];
                    let cp_a = a_piece.cps[i][j];
                    d_piece.cps[i][j].x = cp_p.x - cp_a.x;
                    d_piece.cps[i][j].y = cp_p.y - cp_a.y;
                    d_piece.cps[i][j].z = cp_p.z - cp_a.z;
                }
            }
        }
    }
    f1_cap(&diff)
}

#[cfg(test)]
#[path = "kp13_bi5_cap_f1_tests.rs"]
mod tests;
