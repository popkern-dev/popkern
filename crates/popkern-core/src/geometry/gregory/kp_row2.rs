//! Phase H.4.e — K-P row-2 CP derivation for the bidegree-9
//! patch (G2 continuity).
//!
//! Given a [`KpGuidePatch`] + 4 degree-9 boundary curves + the
//! H.4.d [`Row1Frame`], compute the **row-2 frame** — rows 2, 7
//! and cols 2, 7 of the bidegree-9 patch. These are the CPs that
//! make the bidegree-9 patch's second cross-derivative at each
//! boundary match the guide's second cross-derivative,
//! delivering **G2 continuity across shared spokes** by
//! construction.
//!
//! # Formula
//!
//! Per ADR
//! `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`
//! §2.3:
//!
//! ```text
//!   P[2][j] = P[1][j] + (1/9)  · [∂G/∂u(0, v)]_j^9
//!                     + (1/72) · [∂²G/∂u²(0, v)]_j^9
//!   P[7][j] = P[8][j] − (1/9)  · [∂G/∂u(1, v)]_j^9
//!                     + (1/72) · [∂²G/∂u²(1, v)]_j^9
//!   P[i][2] = P[i][1] + (1/9)  · [∂G/∂v(u, 0)]_i^9
//!                     + (1/72) · [∂²G/∂v²(u, 0)]_i^9
//!   P[i][7] = P[i][8] − (1/9)  · [∂G/∂v(u, 1)]_i^9
//!                     + (1/72) · [∂²G/∂v²(u, 1)]_i^9
//! ```
//!
//! where `[F(t)]_j^9` is the j-th Bernstein CP of `F(t)` at
//! degree 9 (obtained via
//! [`elevate_cubic_vec_to_nine`](super::elevate_cubic_vec_to_nine)).
//! Derivation: the degree-9 Bezier boundary identity
//! `B''(0) = 72·(P_2 − 2·P_1 + P_0)` combined with
//! `B'(0) = 9·(P_1 − P_0)` solved for `P_2` and rewritten
//! relative to `P_1` gives exactly `P_2 = P_1 + B'(0)/9 +
//! B''(0)/72`. Applied column-wise along each boundary.
//!
//! # Corner 3 × 3 blocks — automatic (validated in H.4.f)
//!
//! Early drafts of this module feared that the row-2 / col-2
//! formulas might disagree on `grid[2][2]` (and its mirrors at
//! the other 3 corners), requiring explicit tensor-product
//! mixed-partial formulas involving third- and fourth-order
//! guide partials (`G_uuv`, `G_uvv`, `G_uuvv`). **This turns out
//! to be unnecessary.**
//!
//! The polynomial-identity construction of H.4.d + H.4.e
//! preserves cross-derivative identities pointwise along each
//! boundary. Differentiating in the tangent direction then
//! propagates the identity to higher-order mixed partials: if
//! `∂S/∂u(0, v) ≡ G_u(0, v)` and `∂²S/∂u²(0, v) ≡ G_uu(0, v)`
//! as polynomial identities in `v`, then `∂²S/∂u∂v(0, v) ≡
//! G_uv(0, v)`, `∂³S/∂u²∂v(0, v) ≡ G_uuv(0, v)`, and
//! `∂⁴S/∂u²∂v²(0, v) ≡ G_uuvv(0, v)` follow by differentiation.
//! In particular they hold at `v = 0`, which is the "hard" corner
//! of the 3 × 3 block.
//!
//! `kp_row2_corner_tests.rs` validates this rigorously on an
//! arbitrary non-trivial bicubic guide: every mixed partial up
//! to fourth order at each of the 4 corners matches the guide
//! to `< 1e-8`. G2 at corners is **automatic, not deferred**.

use super::kp_degree_elevate::elevate_cubic_vec_to_nine;
use super::kp_guide_patch::KpGuidePatch;
use super::kp_row1::Row1Frame;
use crate::geometry::{Point3, Vec3};

/// The four G2-frame rows of a bidegree-9 patch, derived by the
/// K-P row-2 formula.
///
/// Indexing convention:
///
/// ```text
///   row2[j] = grid[2][j]   (two rows inside u = 0 boundary)
///   row7[j] = grid[7][j]   (two rows inside u = 1 boundary)
///   col2[i] = grid[i][2]   (two cols inside v = 0 boundary)
///   col7[i] = grid[i][7]   (two cols inside v = 1 boundary)
/// ```
///
/// The `col*` arrays agree with `row*` on overlapping indices by
/// construction (both use the same tensor-product guide data).
#[derive(Debug, Clone)]
pub struct Row2Frame {
    /// `grid[2][j]` for `j ∈ 0..=9`.
    pub row2: [Point3; 10],
    /// `grid[7][j]` for `j ∈ 0..=9`.
    pub row7: [Point3; 10],
    /// `grid[i][2]` for `i ∈ 0..=9`.
    pub col2: [Point3; 10],
    /// `grid[i][7]` for `i ∈ 0..=9`.
    pub col7: [Point3; 10],
}

/// Compute the row-2 frame from a guide + boundaries + H.4.d
/// row-1 frame.
///
/// See the module doc for the formulas and the corner 3 × 3
/// block limitation.
pub fn derive_row2_frame(
    guide: &KpGuidePatch,
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    row1_frame: &Row1Frame,
) -> Row2Frame {
    // Degree-9 representations of the guide's cross-derivative
    // cubics (first partial) and cross-second-derivative cubics
    // (second partial) on each boundary.
    let g_u_u0 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_u0());
    let g_u_u1 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_u1());
    let g_v_v0 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_v0());
    let g_v_v1 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_v1());
    let g_uu_u0 = elevate_cubic_vec_to_nine(guide.cross_second_deriv_cps_uu_u0());
    let g_uu_u1 = elevate_cubic_vec_to_nine(guide.cross_second_deriv_cps_uu_u1());
    let g_vv_v0 = elevate_cubic_vec_to_nine(guide.cross_second_deriv_cps_vv_v0());
    let g_vv_v1 = elevate_cubic_vec_to_nine(guide.cross_second_deriv_cps_vv_v1());

    // Row 2: P[2][j] = P[1][j] + (1/9) · G_u + (1/72) · G_uu
    //   j = 0, 9 → boundary (col 0, col 9).
    //   j ∈ 1..=8 → row-2 formula.
    let mut row2 = [Point3::new(0.0, 0.0, 0.0); 10];
    row2[0] = col0[2];
    for j in 1..=8 {
        row2[j] = row2_formula(row1_frame.row1[j], g_u_u0[j], g_uu_u0[j], 1.0);
    }
    row2[9] = col9[2];

    // Row 7: P[7][j] = P[8][j] − (1/9) · G_u + (1/72) · G_uu
    let mut row7 = [Point3::new(0.0, 0.0, 0.0); 10];
    row7[0] = col0[7];
    for j in 1..=8 {
        row7[j] = row2_formula(row1_frame.row8[j], g_u_u1[j], g_uu_u1[j], -1.0);
    }
    row7[9] = col9[7];

    // Col 2: skip overlaps with row 2 (at i = 2) and row 7 (at
    // i = 7). `col2[i]` for `i ∈ 3..=6` from the col-2 formula;
    // other indices carry values from row 2, row 7, row 1, row 8,
    // row 0, row 9 as appropriate.
    let mut col2 = [Point3::new(0.0, 0.0, 0.0); 10];
    col2[0] = row0[2];
    col2[1] = row1_frame.row1[2];
    col2[2] = row2[2];
    for i in 3..=6 {
        col2[i] = row2_formula(row1_frame.col1[i], g_v_v0[i], g_vv_v0[i], 1.0);
    }
    col2[7] = row7[2];
    col2[8] = row1_frame.row8[2];
    col2[9] = row9[2];

    // Col 7: mirror of col 2 at v = 1.
    let mut col7 = [Point3::new(0.0, 0.0, 0.0); 10];
    col7[0] = row0[7];
    col7[1] = row1_frame.row1[7];
    col7[2] = row2[7];
    for i in 3..=6 {
        col7[i] = row2_formula(row1_frame.col8[i], g_v_v1[i], g_vv_v1[i], -1.0);
    }
    col7[7] = row7[7];
    col7[8] = row1_frame.row8[7];
    col7[9] = row9[7];

    Row2Frame {
        row2,
        row7,
        col2,
        col7,
    }
}

/// Row-2 K-P formula. `sign = +1` for the `u = 0` / `v = 0`
/// boundaries (row 2, col 2); `sign = −1` for the `u = 1` /
/// `v = 1` boundaries (row 7, col 7). The second-derivative
/// term keeps its `+` coefficient either way — the sign flip
/// only affects the first-derivative term because the boundary
/// Bezier identity reverses direction.
#[inline]
fn row2_formula(p1: Point3, g_first: Vec3, g_second: Vec3, sign: f64) -> Point3 {
    let inv9 = sign / 9.0;
    let inv72 = 1.0 / 72.0;
    Point3::new(
        p1.x + g_first.x * inv9 + g_second.x * inv72,
        p1.y + g_first.y * inv9 + g_second.y * inv72,
        p1.z + g_first.z * inv9 + g_second.z * inv72,
    )
}

#[cfg(test)]
#[path = "kp_row2_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "kp_row2_corner_tests.rs"]
mod corner_tests;
