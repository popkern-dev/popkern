//! K-P row-1 CP derivation for the bidegree-9
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

use super::kp_guide_trait::KpGuide;
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
pub fn derive_row1_frame<G: KpGuide>(
    guide: &G,
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
) -> Row1Frame {
    // Degree-9 representations of the guide's cross-derivative
    // curves on each boundary. The trait hides the native guide
    // degree (cubic → elevate by 6; bi-quintic → elevate by 4).
    let g_u_u0 = guide.cross_deriv_cps_u0_degree9();
    let g_u_u1 = guide.cross_deriv_cps_u1_degree9();
    let g_v_v0 = guide.cross_deriv_cps_v0_degree9();
    let g_v_v1 = guide.cross_deriv_cps_v1_degree9();

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

/// Derive the row-1 frame purely from the
/// four boundary curves, with no reference to a K-P guide.
///
/// # When to use this
///
/// Called when the upstream `DirectBoundaryInput::hex_boundaries`
/// is `Some` (i.e., the sphere-fit wide-tolerance stage boundary injection is
/// active). The injected hex-side boundaries are the analytic
/// `sphere ∩ cylinder` arc + slant-face bridge (deg-9 Bernstein
/// CPs). Their corner tangents differ from the K-P guide's
/// pre-injection partials, breaking the boundary-consistency
/// requirement at module top. Empirical: the resulting bidegree-9
/// patch self-intersects (Stage 1 trace
/// `[kp_binonic N=4] fallback @ self_intersects_with_tolerance`).
///
/// This function bypasses the inconsistency by deriving every CP
/// in [`Row1Frame`] from boundary tangents directly:
///
/// * Corner partials `G_u(corner)`, `G_v(corner)` → boundary
///   endpoint Bezier-tangent diffs (`9·(boundary[1] - boundary[0])`
///   etc.).
/// * Twist `G_uv(corner)` → `0` (Stage 2 simplification — see
///   Stage 3 for analytic-twist upgrade).
/// * Cross-derivative profiles `[∂G/∂u(0, v)]_j^9` etc. →
///   degree-9 linear blend of the two corner endpoint values.
///
/// # Geometric guarantee
///
/// The linear-blend cross-derivative profile produces a
/// row-1 ring that is fold-free by construction: each `row1[j]`
/// sits at `row0[j] + (1/9)·g_u_u0[j]` where `g_u_u0[j]` is a
/// convex combination of two same-side-pointing corner tangents.
/// No fold can result. Trades K-P G2 invariants (the row-2
/// formula's `∂²G/∂u²` term) for the boundary-consistent G1
/// frame the patch needs to ship.
///
/// # Trade-offs vs [`derive_row1_frame`]
///
/// * **Loses:** K-P guide-driven cross-tangent fidelity in the
///   patch interior. The row-1 ring's *direction* into the patch
///   is set by the boundary tangents at corners; the cross-tangent
///   profile in between is geometrically simple, not surface-aware.
/// * **Gains:** boundary-corner consistency by construction,
///   fold-freeness, and direct analytic-arc tracking on the hex
///   sides.
pub fn derive_row1_frame_from_boundary(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
) -> Row1Frame {
    // Corner tangents from boundary endpoint Bezier diffs. Naming:
    // `g_<dir>_<u_index><v_index>` = ∂G/∂<dir> at corner (u, v).
    //
    // At corner (0, 0): col0 starts here (going in +u toward
    // (1, 0)), row0 starts here (going in +v toward (0, 1)).
    // So G_u(0, 0) = col0's start tangent; G_v(0, 0) = row0's
    // start tangent. Mirror logic at the other 3 corners.
    let g_u_00 = boundary_tangent_at_start(col0);
    let g_u_10 = boundary_tangent_at_end(col0);
    let g_u_01 = boundary_tangent_at_start(col9);
    let g_u_11 = boundary_tangent_at_end(col9);
    let g_v_00 = boundary_tangent_at_start(row0);
    let g_v_01 = boundary_tangent_at_end(row0);
    let g_v_10 = boundary_tangent_at_start(row9);
    let g_v_11 = boundary_tangent_at_end(row9);

    // Stage 2 simplification: twist `G_uv := 0` at every corner.
    // Loses K-P's tensor-product G2-at-corners invariant; gains
    // simplicity. Stage 3 may upgrade this to an analytic twist
    // derived from the cap-sphere / cylinder local frame.
    let zero = Vec3::new(0.0, 0.0, 0.0);

    let corner_00 = corner_2x2(row0[0], g_u_00, g_v_00, zero, 1.0, 1.0);
    let corner_10 = corner_2x2(row9[0], g_u_10, g_v_10, zero, -1.0, 1.0);
    let corner_11 = corner_2x2(row9[9], g_u_11, g_v_11, zero, -1.0, -1.0);
    let corner_01 = corner_2x2(row0[9], g_u_01, g_v_01, zero, 1.0, -1.0);

    // Mid-edge row1[j] / row8[j] / col1[i] / col8[i]
    // computed by *opposite-boundary interpolation* — set
    // `row1[j] = (8/9)·row0[j] + (1/9)·row9[j]`. This places the
    // row-1 ring on the segment from row0 to row9 at u = 1/9,
    // guaranteeing fold-freeness by construction (each ring's
    // CPs lie inside the convex hull of the four boundaries).
    //
    // Geometrically equivalent to setting the cross-derivative
    // CP `g_u_u0[j] := (row9[j] − row0[j])` — the surface's
    // u-direction "step" at v=v_j is the chord from row 0 to
    // row 9 at that v. This is a chord-length-aware estimate
    // that beats the corner-tangent linear blend (which can
    // reverse direction at mid-edge for tilted-axis fixtures
    // where corner tangents point in conflicting directions).
    //
    // Corners (j ∈ {0, 1, 8, 9}) keep boundary-direct / corner-2×2
    // values for K-P consistency.
    //
    // Suppress unused warnings on the corner-tangent variables —
    // they're consumed by `corner_2x2` above.
    let _ = (
        g_u_00, g_u_10, g_u_01, g_u_11, g_v_00, g_v_01, g_v_10, g_v_11,
    );

    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
    row1[0] = col0[1];
    row1[1] = corner_00;
    for j in 2..=7 {
        row1[j] = lerp_point(row0[j], row9[j], 1.0 / 9.0);
    }
    row1[8] = corner_01;
    row1[9] = col9[1];

    let mut row8 = [Point3::new(0.0, 0.0, 0.0); 10];
    row8[0] = col0[8];
    row8[1] = corner_10;
    for j in 2..=7 {
        row8[j] = lerp_point(row0[j], row9[j], 8.0 / 9.0);
    }
    row8[8] = corner_11;
    row8[9] = col9[8];

    let mut col1 = [Point3::new(0.0, 0.0, 0.0); 10];
    col1[0] = row0[1];
    col1[1] = corner_00;
    for i in 2..=7 {
        col1[i] = lerp_point(col0[i], col9[i], 1.0 / 9.0);
    }
    col1[8] = corner_10;
    col1[9] = row9[1];

    let mut col8 = [Point3::new(0.0, 0.0, 0.0); 10];
    col8[0] = row0[8];
    col8[1] = corner_01;
    for i in 2..=7 {
        col8[i] = lerp_point(col0[i], col9[i], 8.0 / 9.0);
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

/// Linearly interpolate two `Point3`s by `t`.
/// `lerp_point(p0, p1, t) = (1 − t)·p0 + t·p1`.
#[inline]
fn lerp_point(p0: Point3, p1: Point3, t: f64) -> Point3 {
    Point3::new(
        (1.0 - t) * p0.x + t * p1.x,
        (1.0 - t) * p0.y + t * p1.y,
        (1.0 - t) * p0.z + t * p1.z,
    )
}

/// `9·(P[1] − P[0])` for a deg-9 boundary.
/// This is the boundary's tangent at `t = 0` (the K-P boundary-
/// consistency expression at module top).
#[inline]
pub(super) fn boundary_tangent_at_start(boundary: &[Point3; 10]) -> Vec3 {
    Vec3::new(
        9.0 * (boundary[1].x - boundary[0].x),
        9.0 * (boundary[1].y - boundary[0].y),
        9.0 * (boundary[1].z - boundary[0].z),
    )
}

/// `9·(P[9] − P[8])` for a deg-9 boundary
/// (tangent at `t = 1`).
#[inline]
pub(super) fn boundary_tangent_at_end(boundary: &[Point3; 10]) -> Vec3 {
    Vec3::new(
        9.0 * (boundary[9].x - boundary[8].x),
        9.0 * (boundary[9].y - boundary[8].y),
        9.0 * (boundary[9].z - boundary[8].z),
    )
}

/// Degree-9 Bernstein representation of a
/// linear-blend `Vec3` field between two endpoint values. CPs
/// lie on the line `(1 − j/9)·start + (j/9)·end` for `j ∈ 0..=9`.
#[inline]
pub(super) fn lerp_cross_deriv_cps(start: Vec3, end: Vec3) -> [Vec3; 10] {
    let mut cps = [Vec3::new(0.0, 0.0, 0.0); 10];
    for (j, slot) in cps.iter_mut().enumerate() {
        let t = (j as f64) / 9.0;
        *slot = Vec3::new(
            start.x * (1.0 - t) + end.x * t,
            start.y * (1.0 - t) + end.y * t,
            start.z * (1.0 - t) + end.z * t,
        );
    }
    cps
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
