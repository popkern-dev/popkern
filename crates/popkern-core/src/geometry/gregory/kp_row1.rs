//! Phase H.4.d — K-P row-1 CP derivation for the bidegree-9
//! patch.
//!
//! Given a bicubic [`KpGuidePatch`] + four boundary curves at
//! degree 9 (`row 0`, `row 9`, `col 0`, `col 9`), compute the
//! four **frame rows** (`row 1`, `row 8`, `col 1`, `col 8`) of
//! the bidegree-9 patch. These are the row-1-equivalent CPs that
//! guarantee G1 across shared spokes **by construction** — see
//! ADR
//! `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`
//! §2.3.
//!
//! # Formula summary
//!
//! Per-edge mid-run (non-corner indices):
//!
//! ```text
//!   row1[j] = row0[j] + (1/9) · [∂G/∂u(0, v)]_j^9   for j ∈ 2..=7
//!   row8[j] = row9[j] − (1/9) · [∂G/∂u(1, v)]_j^9   for j ∈ 2..=7
//!   col1[i] = col0[i] + (1/9) · [∂G/∂v(u, 0)]_i^9   for i ∈ 2..=7
//!   col8[i] = col9[i] − (1/9) · [∂G/∂v(u, 1)]_i^9   for i ∈ 2..=7
//! ```
//!
//! where `[F(t)]_j^9` denotes the j-th Bernstein CP of `F(t)` at
//! degree 9 — i.e. `elevate_cubic_vec_to_nine(guide.cross_deriv_cps_*)`.
//!
//! Corner `2 × 2` blocks (`grid[1][1]`, `[1][8]`, `[8][1]`,
//! `[8][8]`) use the tensor-product mixed-partial identity
//!
//! ```text
//!   P[u_in][v_in] = P[u_c][v_c]
//!                   + u_sign · G_u / 9
//!                   + v_sign · G_v / 9
//!                   + u_sign · v_sign · G_uv / 81
//! ```
//!
//! where `u_sign = +1` at the `u = 0` corners, `-1` at `u = 1`;
//! analogously for `v_sign`. Derivatives are evaluated at the
//! corner.
//!
//! # Boundary consistency
//!
//! The caller-supplied boundary CPs must be **tangent-consistent
//! with the guide** at each corner. Specifically, the cubic-
//! Bezier tangent of each boundary curve at its endpoints must
//! equal the guide's corresponding directional partial:
//!
//! ```text
//!   9 · (col0[1] − col0[0]) = G_u(0, 0)
//!   9 · (row0[1] − row0[0]) = G_v(0, 0)
//!   ... (and 6 more at the other 3 corners)
//! ```
//!
//! If consistency holds, the mid-edge row-1 formula and the
//! `2 × 2` corner formula produce values compatible with the
//! boundary (i.e. `row1[0] = col0[1]` as predicted by the
//! formula = as given by the boundary). If consistency fails,
//! the frame has corner-row / corner-col conflicts — the
//! resulting patch will have small geometric artefacts at the
//! corners. Upstream code is responsible for constructing
//! boundaries that match the guide.

use super::kp_degree_elevate::elevate_cubic_vec_to_nine;
use super::kp_guide_patch::KpGuidePatch;
use crate::geometry::{Point3, Vec3};

/// The four frame rows of a bidegree-9 patch, derived by the K-P
/// row-1 formula + corner `2 × 2` blocks.
///
/// Indexing convention:
///
/// ```text
///   row1[j] = grid[1][j]   (adjacent to u = 0 boundary)
///   row8[j] = grid[8][j]   (adjacent to u = 1 boundary)
///   col1[i] = grid[i][1]   (adjacent to v = 0 boundary)
///   col8[i] = grid[i][8]   (adjacent to v = 1 boundary)
/// ```
///
/// Corner slots are consistent across `row1`/`col1` etc. by
/// construction (both the boundary-tangent assumption and the
/// `2 × 2` corner formula produce the same value there, modulo
/// input boundary consistency).
#[derive(Debug, Clone)]
pub struct Row1Frame {
    /// `grid[1][j]` for `j ∈ 0..=9` — the row adjacent to the
    /// `u = 0` boundary.
    pub row1: [Point3; 10],
    /// `grid[8][j]` for `j ∈ 0..=9` — the row adjacent to the
    /// `u = 1` boundary.
    pub row8: [Point3; 10],
    /// `grid[i][1]` for `i ∈ 0..=9` — the column adjacent to the
    /// `v = 0` boundary.
    pub col1: [Point3; 10],
    /// `grid[i][8]` for `i ∈ 0..=9` — the column adjacent to the
    /// `v = 1` boundary.
    pub col8: [Point3; 10],
}

/// Compute the row-1 frame from a guide + four degree-9 boundary
/// curves.
///
/// See the module doc for the formulas and boundary-consistency
/// requirement.
pub fn derive_row1_frame(
    guide: &KpGuidePatch,
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
) -> Row1Frame {
    // Degree-9 representations of the guide's cross-derivative
    // cubics on each boundary.
    let g_u_u0 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_u0());
    let g_u_u1 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_u1());
    let g_v_v0 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_v0());
    let g_v_v1 = elevate_cubic_vec_to_nine(guide.cross_deriv_cps_v1());

    // Per-corner data for the 2 × 2 blocks.
    let corner_00 = corner_2x2(
        row0[0],
        guide.partial_u(0.0, 0.0),
        guide.partial_v(0.0, 0.0),
        guide.second_partial_uv(0.0, 0.0),
        1.0,
        1.0,
    );
    let corner_10 = corner_2x2(
        row9[0],
        guide.partial_u(1.0, 0.0),
        guide.partial_v(1.0, 0.0),
        guide.second_partial_uv(1.0, 0.0),
        -1.0,
        1.0,
    );
    let corner_11 = corner_2x2(
        row9[9],
        guide.partial_u(1.0, 1.0),
        guide.partial_v(1.0, 1.0),
        guide.second_partial_uv(1.0, 1.0),
        -1.0,
        -1.0,
    );
    let corner_01 = corner_2x2(
        row0[9],
        guide.partial_u(0.0, 1.0),
        guide.partial_v(0.0, 1.0),
        guide.second_partial_uv(0.0, 1.0),
        1.0,
        -1.0,
    );

    // Row 1: j ∈ 2..=7 via row-1 mid-edge formula; j ∈ {0, 9}
    // take the boundary; j ∈ {1, 8} take the corner 2×2.
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
    row1[0] = col0[1];
    row1[1] = corner_00;
    for j in 2..=7 {
        row1[j] = add_scaled(row0[j], g_u_u0[j], 1.0 / 9.0);
    }
    row1[8] = corner_01;
    row1[9] = col9[1];

    // Row 8: mirror of row 1 at u = 1.
    let mut row8 = [Point3::new(0.0, 0.0, 0.0); 10];
    row8[0] = col0[8];
    row8[1] = corner_10;
    for j in 2..=7 {
        row8[j] = add_scaled(row9[j], g_u_u1[j], -1.0 / 9.0);
    }
    row8[8] = corner_11;
    row8[9] = col9[8];

    // Col 1: i ∈ 2..=7 via col-1 mid-edge formula.
    let mut col1 = [Point3::new(0.0, 0.0, 0.0); 10];
    col1[0] = row0[1];
    col1[1] = corner_00;
    for i in 2..=7 {
        col1[i] = add_scaled(col0[i], g_v_v0[i], 1.0 / 9.0);
    }
    col1[8] = corner_10;
    col1[9] = row9[1];

    // Col 8: i ∈ 2..=7 via col-8 mid-edge formula.
    let mut col8 = [Point3::new(0.0, 0.0, 0.0); 10];
    col8[0] = row0[8];
    col8[1] = corner_01;
    for i in 2..=7 {
        col8[i] = add_scaled(col9[i], g_v_v1[i], -1.0 / 9.0);
    }
    col8[8] = corner_11;
    col8[9] = row9[8];

    Row1Frame {
        row1,
        row8,
        col1,
        col8,
    }
}

/// Tensor-product mixed-partial formula for the corner 2 × 2
/// inner CP. `u_sign`, `v_sign` are `±1` depending on which
/// corner of the unit square the data is attached to (see module
/// doc).
#[inline]
fn corner_2x2(
    corner_pos: Point3,
    g_u: Vec3,
    g_v: Vec3,
    g_uv: Vec3,
    u_sign: f64,
    v_sign: f64,
) -> Point3 {
    let inv9 = 1.0 / 9.0;
    let inv81 = 1.0 / 81.0;
    Point3::new(
        corner_pos.x
            + u_sign * g_u.x * inv9
            + v_sign * g_v.x * inv9
            + u_sign * v_sign * g_uv.x * inv81,
        corner_pos.y
            + u_sign * g_u.y * inv9
            + v_sign * g_v.y * inv9
            + u_sign * v_sign * g_uv.y * inv81,
        corner_pos.z
            + u_sign * g_u.z * inv9
            + v_sign * g_v.z * inv9
            + u_sign * v_sign * g_uv.z * inv81,
    )
}

/// `p + v · s` as a [`Point3`]. The mid-edge row-1 formula is
/// exactly this shape.
#[inline]
fn add_scaled(p: Point3, v: Vec3, s: f64) -> Point3 {
    Point3::new(p.x + v.x * s, p.y + v.y * s, p.z + v.z * s)
}

#[cfg(test)]
#[path = "kp_row1_tests.rs"]
mod tests;
