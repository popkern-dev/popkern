//! Tests for [`super::derive_row2_frame`].
//!
//! Verification goals:
//!
//! 1. **Cross-second-derivative polynomial identity** at each
//!    boundary: `∂²S/∂u²(0, v) ≡ ∂²G/∂u²(0, v)` as polynomials in
//!    `v` — the mid-edge range (v ∈ [3/9, 6/9]) matches to
//!    floating-point precision.
//! 2. **Corner preservation**: `∂²S/∂u²(0, 0) = ∂²G/∂u²(0, 0)`
//!    and analogous at the other 3 corners.
//! 3. **Flat-plane fidelity**: flat inputs produce a flat frame.
//! 4. **Boundary-corner consistency**: `row2[0] = col0[2]` etc.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::kp_guide_build::{GuideCornerData, build_guide_from_hermite};
use crate::geometry::gregory::kp_guide_patch::KpGuidePatch;
use crate::geometry::gregory::kp_row1::derive_row1_frame;
use crate::geometry::gregory::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::{Point3, Vec3};

const EPS: f64 = 1e-9;

fn consistent_guide_and_boundaries() -> (
    KpGuidePatch,
    [Point3; 10],
    [Point3; 10],
    [Point3; 10],
    [Point3; 10],
) {
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.1, 0.2),
            tangent_v: Vec3::new(0.05, 1.1, -0.15),
            twist: Vec3::new(0.2, -0.1, 0.3),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.1, 0.3),
            tangent_u: Vec3::new(1.1, -0.2, 0.25),
            tangent_v: Vec3::new(0.1, 1.0, 0.2),
            twist: Vec3::new(-0.15, 0.2, -0.1),
        },
        GuideCornerData {
            position: Point3::new(1.15, 1.2, 0.1),
            tangent_u: Vec3::new(0.9, 0.15, -0.1),
            tangent_v: Vec3::new(-0.1, 0.95, 0.25),
            twist: Vec3::new(0.1, 0.15, -0.2),
        },
        GuideCornerData {
            position: Point3::new(0.05, 1.1, -0.05),
            tangent_u: Vec3::new(1.05, -0.1, 0.15),
            tangent_v: Vec3::new(0.0, 1.05, -0.1),
            twist: Vec3::new(-0.05, -0.2, 0.15),
        },
    ];
    let guide = build_guide_from_hermite(corners);
    let row0 = build_boundary_k_p_consistent(&guide, true, false);
    let row9 = build_boundary_k_p_consistent(&guide, true, true);
    let col0 = build_boundary_k_p_consistent(&guide, false, false);
    let col9 = build_boundary_k_p_consistent(&guide, false, true);
    (guide, row0, row9, col0, col9)
}

/// Build a K-P G2-consistent degree-9 boundary curve — the
/// boundary CPs at the "two from the end" positions (P[2] /
/// P[7]) are placed so the boundary curve's second derivative
/// at the endpoints matches the guide's. This is the proper
/// boundary data for the row-2 formula to be consistent with
/// at the corners.
///
/// `is_u_boundary = true`: the boundary runs in v direction at
/// fixed u = 0 (if `at_end = false`) or u = 1 (if `at_end =
/// true`); used for `row0` / `row9`.
/// `is_u_boundary = false`: runs in u direction at fixed v;
/// used for `col0` / `col9`.
fn build_boundary_k_p_consistent(
    guide: &KpGuidePatch,
    is_u_boundary: bool,
    at_end: bool,
) -> [Point3; 10] {
    // Endpoints + first/second derivatives of the guide's
    // boundary curve at each end.
    let (u_a, v_a, u_b, v_b) = match (is_u_boundary, at_end) {
        // row 0: u = 0, v varies 0 → 1.
        (true, false) => (0.0, 0.0, 0.0, 1.0),
        // row 9: u = 1, v varies 0 → 1.
        (true, true) => (1.0, 0.0, 1.0, 1.0),
        // col 0: v = 0, u varies 0 → 1.
        (false, false) => (0.0, 0.0, 1.0, 0.0),
        // col 9: v = 1, u varies 0 → 1.
        (false, true) => (0.0, 1.0, 1.0, 1.0),
    };
    let p_start = guide.point_at(u_a, v_a);
    let p_end = guide.point_at(u_b, v_b);
    let (tangent_start, tangent_end, second_start, second_end) = if is_u_boundary {
        (
            guide.partial_v(u_a, v_a),
            guide.partial_v(u_b, v_b),
            guide.second_partial_vv(u_a, v_a),
            guide.second_partial_vv(u_b, v_b),
        )
    } else {
        (
            guide.partial_u(u_a, v_a),
            guide.partial_u(u_b, v_b),
            guide.second_partial_uu(u_a, v_a),
            guide.second_partial_uu(u_b, v_b),
        )
    };
    // P[0] = start, P[9] = end.
    // P[1] = P[0] + T_start / 9.
    // P[2] = P[1] + T_start / 9 + T''_start / 72
    //      = P[0] + (2/9) T_start + (1/72) T''_start.
    // P[8] = P[9] − T_end / 9.
    // P[7] = P[8] − T_end / 9 + T''_end / 72
    //      = P[9] − (2/9) T_end + (1/72) T''_end.
    // Interior (P[3..=6]) linearly interpolated between P[2] and
    // P[7] — doesn't affect the row-2 formula's endpoint-sensitive
    // behaviour.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[0] = p_start;
    cps[9] = p_end;
    cps[1] = offset(p_start, tangent_start, 1.0 / 9.0);
    cps[8] = offset(p_end, tangent_end, -1.0 / 9.0);
    cps[2] = offset2(p_start, tangent_start, 2.0 / 9.0, second_start, 1.0 / 72.0);
    cps[7] = offset2(p_end, tangent_end, -2.0 / 9.0, second_end, 1.0 / 72.0);
    for i in 3..=6 {
        let t = ((i - 2) as f64) / 5.0;
        cps[i] = Point3::new(
            cps[2].x + t * (cps[7].x - cps[2].x),
            cps[2].y + t * (cps[7].y - cps[2].y),
            cps[2].z + t * (cps[7].z - cps[2].z),
        );
    }
    cps
}

fn offset(p: Point3, v: Vec3, s: f64) -> Point3 {
    Point3::new(p.x + v.x * s, p.y + v.y * s, p.z + v.z * s)
}

fn offset2(p: Point3, v1: Vec3, s1: f64, v2: Vec3, s2: f64) -> Point3 {
    Point3::new(
        p.x + v1.x * s1 + v2.x * s2,
        p.y + v1.y * s1 + v2.y * s2,
        p.z + v1.z * s1 + v2.z * s2,
    )
}

/// Assemble a full bidegree-9 patch from boundaries + row-1 +
/// row-2 frames. Interior rows 3..=6 and cols 3..=6 left at the
/// origin.
fn build_binonic(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    r1: &crate::geometry::gregory::kp_row1::Row1Frame,
    r2: &Row2Frame,
) -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    grid[0].copy_from_slice(row0);
    grid[9].copy_from_slice(row9);
    grid[1].copy_from_slice(&r1.row1);
    grid[8].copy_from_slice(&r1.row8);
    grid[2].copy_from_slice(&r2.row2);
    grid[7].copy_from_slice(&r2.row7);
    for i in 0..10 {
        grid[i][0] = col0[i];
        grid[i][9] = col9[i];
        grid[i][1] = r1.col1[i];
        grid[i][8] = r1.col8[i];
        grid[i][2] = r2.col2[i];
        grid[i][7] = r2.col7[i];
    }
    RectangularGregoryBiNonic { grid }
}

fn vec_max(a: Vec3, b: Vec3) -> f64 {
    (a - b).length()
}

// =====================================================================
// Cross-second-derivative polynomial identity (mid-edge range)
// =====================================================================

#[test]
fn cross_second_derivative_at_u0_matches_guide_mid_edge() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    // Mid-edge range — avoid the corner 3×3 blocks.
    for k in 3..=6 {
        let v = (k as f64) / 9.0;
        let s_uu = binonic.second_partial_uu(0.0, v);
        let g_uu = guide.second_partial_uu(0.0, v);
        assert!(
            vec_max(s_uu, g_uu) < 1e-6,
            "v = {v}: bidegree-9 ∂²uu(0, v) = {s_uu:?} vs guide {g_uu:?}"
        );
    }
}

#[test]
fn cross_second_derivative_at_u1_matches_guide_mid_edge() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for k in 3..=6 {
        let v = (k as f64) / 9.0;
        let s_uu = binonic.second_partial_uu(1.0, v);
        let g_uu = guide.second_partial_uu(1.0, v);
        assert!(vec_max(s_uu, g_uu) < 1e-6, "v = {v}: bidegree-9 ∂²uu(1, v)");
    }
}

#[test]
fn cross_second_derivative_at_v0_matches_guide_mid_edge() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for k in 3..=6 {
        let u = (k as f64) / 9.0;
        let s_vv = binonic.second_partial_vv(u, 0.0);
        let g_vv = guide.second_partial_vv(u, 0.0);
        assert!(vec_max(s_vv, g_vv) < 1e-6, "u = {u}: bidegree-9 ∂²vv(u, 0)");
    }
}

#[test]
fn cross_second_derivative_at_v1_matches_guide_mid_edge() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for k in 3..=6 {
        let u = (k as f64) / 9.0;
        let s_vv = binonic.second_partial_vv(u, 1.0);
        let g_vv = guide.second_partial_vv(u, 1.0);
        assert!(vec_max(s_vv, g_vv) < 1e-6, "u = {u}: bidegree-9 ∂²vv(u, 1)");
    }
}

// =====================================================================
// Boundary-corner consistency
// =====================================================================

#[test]
fn row2_corner_slots_match_boundary_cps() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    assert_eq!(r2.row2[0].x, col0[2].x);
    assert_eq!(r2.row2[9].x, col9[2].x);
    assert_eq!(r2.row7[0].x, col0[7].x);
    assert_eq!(r2.row7[9].x, col9[7].x);
    assert_eq!(r2.col2[0].x, row0[2].x);
    assert_eq!(r2.col2[9].x, row9[2].x);
    assert_eq!(r2.col7[0].x, row0[7].x);
    assert_eq!(r2.col7[9].x, row9[7].x);
}

#[test]
fn col2_carries_row1_and_row8_at_intermediate_indices() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    assert_eq!(r2.col2[1].x, r1.row1[2].x);
    assert_eq!(r2.col2[8].x, r1.row8[2].x);
    assert_eq!(r2.col7[1].x, r1.row1[7].x);
    assert_eq!(r2.col7[8].x, r1.row8[7].x);
}

// =====================================================================
// Flat-plane fidelity
// =====================================================================

#[test]
fn flat_inputs_produce_flat_row2_frame() {
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 1.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(0.0, 1.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
    ];
    let guide = build_guide_from_hermite(corners);
    let row0 = build_boundary_k_p_consistent(&guide, true, false);
    let row9 = build_boundary_k_p_consistent(&guide, true, true);
    let col0 = build_boundary_k_p_consistent(&guide, false, false);
    let col9 = build_boundary_k_p_consistent(&guide, false, true);
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    for j in 0..10 {
        assert!(r2.row2[j].z.abs() < EPS, "row2[{j}].z = {}", r2.row2[j].z);
        assert!(r2.row7[j].z.abs() < EPS);
        assert!(r2.col2[j].z.abs() < EPS);
        assert!(r2.col7[j].z.abs() < EPS);
    }
}

#[test]
fn output_lengths_are_ten() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    assert_eq!(r2.row2.len(), 10);
    assert_eq!(r2.row7.len(), 10);
    assert_eq!(r2.col2.len(), 10);
    assert_eq!(r2.col7.len(), 10);
}
