//! H.4.f validation — does the row-1 + row-2 construction
//! AUTOMATICALLY satisfy full G2 at the 4 corner 3 × 3 blocks,
//! or do we need explicit corner-block formulas?
//!
//! The theoretical argument: if H.4.d preserves the polynomial
//! identity `∂S/∂u(0, v) ≡ G_u(0, v)` in v, then differentiating
//! in v also preserves derived identities:
//!
//! ```text
//!   ∂²S/∂u∂v(0, v) ≡ G_uv(0, v)        (from H.4.d)
//!   ∂³S/∂u²∂v(0, v) ≡ G_uuv(0, v)      (from H.4.e, differentiated)
//!   ∂³S/∂u∂v²(u, 0) ≡ G_uvv(u, 0)      (from H.4.e, differentiated)
//!   ∂⁴S/∂u²∂v²(u, 0) ≡ G_uuvv(u, 0)    (from H.4.e differentiated twice)
//! ```
//!
//! These should hold at every point on the boundary — and in
//! particular at the corners. So the corner 3 × 3 block should
//! be satisfied by construction of H.4.d + H.4.e, **not**
//! requiring additional formulas.
//!
//! This test file validates that hypothesis. If every test
//! passes, the H.4.e documentation's note about "corner 3 × 3
//! deferred" is actually **unnecessary** — G1 and G2 at the
//! corners are inherited automatically. If tests fail, we need
//! explicit corner formulas.

#![allow(clippy::needless_range_loop)]

use crate::geometry::gregory::kp_guide_build::{GuideCornerData, build_guide_from_hermite};
use crate::geometry::gregory::kp_guide_patch::KpGuidePatch;
use crate::geometry::gregory::kp_row1::derive_row1_frame;
use crate::geometry::gregory::kp_row2::derive_row2_frame;
use crate::geometry::gregory::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::{Point3, Vec3};

const EPS: f64 = 1e-8;

fn vec_max(a: Vec3, b: Vec3) -> f64 {
    (a - b).length()
}

/// Build a guide with arbitrary, non-trivial Hermite data + the
/// K-P-consistent boundaries needed by H.4.d/e.
fn arbitrary_guide_with_boundaries() -> (
    KpGuidePatch,
    [Point3; 10],
    [Point3; 10],
    [Point3; 10],
    [Point3; 10],
) {
    let corners = [
        GuideCornerData {
            position: Point3::new(0.1, 0.2, 0.3),
            tangent_u: Vec3::new(1.2, 0.15, 0.25),
            tangent_v: Vec3::new(0.08, 1.15, -0.2),
            twist: Vec3::new(0.25, -0.15, 0.35),
        },
        GuideCornerData {
            position: Point3::new(1.3, 0.25, 0.45),
            tangent_u: Vec3::new(1.25, -0.25, 0.3),
            tangent_v: Vec3::new(0.12, 1.1, 0.25),
            twist: Vec3::new(-0.2, 0.25, -0.15),
        },
        GuideCornerData {
            position: Point3::new(1.5, 1.4, 0.2),
            tangent_u: Vec3::new(0.95, 0.2, -0.15),
            tangent_v: Vec3::new(-0.15, 1.0, 0.3),
            twist: Vec3::new(0.15, 0.2, -0.25),
        },
        GuideCornerData {
            position: Point3::new(0.2, 1.3, 0.05),
            tangent_u: Vec3::new(1.1, -0.15, 0.2),
            tangent_v: Vec3::new(-0.05, 1.1, -0.15),
            twist: Vec3::new(-0.1, -0.25, 0.2),
        },
    ];
    let guide = build_guide_from_hermite(corners);
    let row0 = k_p_boundary(&guide, true, false);
    let row9 = k_p_boundary(&guide, true, true);
    let col0 = k_p_boundary(&guide, false, false);
    let col9 = k_p_boundary(&guide, false, true);
    (guide, row0, row9, col0, col9)
}

fn k_p_boundary(guide: &KpGuidePatch, is_u: bool, at_end: bool) -> [Point3; 10] {
    let (ua, va, ub, vb) = match (is_u, at_end) {
        (true, false) => (0.0, 0.0, 0.0, 1.0),
        (true, true) => (1.0, 0.0, 1.0, 1.0),
        (false, false) => (0.0, 0.0, 1.0, 0.0),
        (false, true) => (0.0, 1.0, 1.0, 1.0),
    };
    let p_start = guide.point_at(ua, va);
    let p_end = guide.point_at(ub, vb);
    let (t_start, t_end, tt_start, tt_end) = if is_u {
        (
            guide.partial_v(ua, va),
            guide.partial_v(ub, vb),
            guide.second_partial_vv(ua, va),
            guide.second_partial_vv(ub, vb),
        )
    } else {
        (
            guide.partial_u(ua, va),
            guide.partial_u(ub, vb),
            guide.second_partial_uu(ua, va),
            guide.second_partial_uu(ub, vb),
        )
    };
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[0] = p_start;
    cps[9] = p_end;
    cps[1] = offset(p_start, t_start, 1.0 / 9.0);
    cps[8] = offset(p_end, t_end, -1.0 / 9.0);
    cps[2] = offset2(p_start, t_start, 2.0 / 9.0, tt_start, 1.0 / 72.0);
    cps[7] = offset2(p_end, t_end, -2.0 / 9.0, tt_end, 1.0 / 72.0);
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

fn build_binonic(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    r1: &crate::geometry::gregory::kp_row1::Row1Frame,
    r2: &crate::geometry::gregory::kp_row2::Row2Frame,
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

// =====================================================================
// Corner 3 × 3 preservation — the critical hypothesis tests
// =====================================================================

#[test]
fn third_partial_uuv_at_corners_matches_guide() {
    // ∂³S/∂u²∂v(corner) should equal G_uuv(corner) automatically,
    // since row-2 preserves ∂²S/∂u² polynomial identity in v.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for &(u, v) in &[(0.0_f64, 0.0_f64), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let s_uuv = binonic.third_partial_uuv(u, v);
        let g_uuv = guide.third_partial_uuv(u, v);
        assert!(
            vec_max(s_uuv, g_uuv) < EPS,
            "corner ({u}, {v}) ∂³uuv: S = {s_uuv:?}, G = {g_uuv:?}"
        );
    }
}

#[test]
fn third_partial_uvv_at_corners_matches_guide() {
    // Mirror of uuv — should hold because col-2 preserves
    // ∂²S/∂v² polynomial identity.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for &(u, v) in &[(0.0_f64, 0.0_f64), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let s_uvv = binonic.third_partial_uvv(u, v);
        let g_uvv = guide.third_partial_uvv(u, v);
        assert!(
            vec_max(s_uvv, g_uvv) < EPS,
            "corner ({u}, {v}) ∂³uvv: S = {s_uvv:?}, G = {g_uvv:?}"
        );
    }
}

#[test]
fn fourth_partial_uuvv_at_corners_matches_guide() {
    // The load-bearing test — ∂⁴S/∂u²∂v² at each corner. Involves
    // grid[2][2] (the corner-deep CP) which is set by BOTH
    // row-2 at j=2 AND col-2 at i=2; current H.4.e lets row-2
    // win. For this test to pass, the guide must be a valid
    // tensor-product bicubic (which it is) AND the polynomial
    // identities must compose tensor-product-wise.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for &(u, v) in &[(0.0_f64, 0.0_f64), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let s_uuvv = binonic.fourth_partial_uuvv(u, v);
        let g_uuvv = guide.fourth_partial_uuvv(u, v);
        assert!(
            vec_max(s_uuvv, g_uuvv) < EPS,
            "corner ({u}, {v}) ∂⁴uuvv: S = {s_uuvv:?}, G = {g_uuvv:?}"
        );
    }
}

#[test]
fn row2_at_j_equals_2_matches_col2_at_i_equals_2() {
    // Direct consistency check: does "row-2 wins at j=2" agree
    // with what col-2 would produce at i=2? If yes, the "row
    // wins" choice in H.4.e is without loss. If no, we have
    // ambiguity that requires explicit corner formulas.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);

    // Row-2 at j=2 produces r2.row2[2].
    let row_val = r2.row2[2];

    // Compute what col-2 formula at i=2 would produce.
    // From the col-2 formula: P[2][2] = col1[2] + [G_v(u,0)]_2^9 / 9 + [G_vv(u,0)]_2^9 / 72
    let col1_at_i2 = r1.col1[2]; // = P[2][1] from H.4.d's col-1 mid-edge formula
    let g_v_v0 = crate::geometry::gregory::kp_degree_elevate::elevate_cubic_vec_to_nine(
        guide.cross_deriv_cps_v0(),
    );
    let g_vv_v0 = crate::geometry::gregory::kp_degree_elevate::elevate_cubic_vec_to_nine(
        guide.cross_second_deriv_cps_vv_v0(),
    );
    let col_val = Point3::new(
        col1_at_i2.x + g_v_v0[2].x / 9.0 + g_vv_v0[2].x / 72.0,
        col1_at_i2.y + g_v_v0[2].y / 9.0 + g_vv_v0[2].y / 72.0,
        col1_at_i2.z + g_v_v0[2].z / 9.0 + g_vv_v0[2].z / 72.0,
    );

    let diff = ((row_val.x - col_val.x).powi(2)
        + (row_val.y - col_val.y).powi(2)
        + (row_val.z - col_val.z).powi(2))
    .sqrt();
    assert!(
        diff < EPS,
        "row-2 at j=2 disagrees with col-2 at i=2: row = {row_val:?}, col = {col_val:?}, diff = {diff:e}"
    );
}

// =====================================================================
// Higher-order partial FD validation (sanity)
// =====================================================================

#[test]
fn third_partial_uuv_polynomial_identity_along_u0_boundary() {
    // Stronger than corner-only: ∂³S/∂u²∂v(0, v) ≡ G_uuv(0, v)
    // as polynomials in v. If this holds for all v (not just
    // v ∈ {0, 1}), we have G2 at all boundary points including
    // the mid-edge range, automatically from H.4.e.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for k in 0..=10 {
        let v = (k as f64) / 10.0;
        let s = binonic.third_partial_uuv(0.0, v);
        let g = guide.third_partial_uuv(0.0, v);
        assert!(
            vec_max(s, g) < 1e-6,
            "v = {v}: ∂³uuv(0, v) S = {s:?}, G = {g:?}"
        );
    }
}

#[test]
fn fourth_partial_uuvv_matches_guide_at_boundary_sample_u0() {
    // ∂⁴S/∂u²∂v²(0, v) matches G_uuvv(0, v) — this is the G2
    // "cross-second-derivative second derivative" test.
    // Stronger than corner-only test; verifies the polynomial
    // identity at dense v samples along u=0.
    let (guide, row0, row9, col0, col9) = arbitrary_guide_with_boundaries();
    let r1 = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let r2 = derive_row2_frame(&guide, &row0, &row9, &col0, &col9, &r1);
    let binonic = build_binonic(&row0, &row9, &col0, &col9, &r1, &r2);
    for k in 0..=10 {
        let v = (k as f64) / 10.0;
        let s = binonic.fourth_partial_uuvv(0.0, v);
        let g = guide.fourth_partial_uuvv(0.0, v);
        assert!(vec_max(s, g) < 1e-6, "v = {v}: ∂⁴uuvv(0, v) S vs G");
    }
}

#[test]
fn guide_third_partial_uuv_matches_finite_difference() {
    let (guide, _, _, _, _) = arbitrary_guide_with_boundaries();
    let h = 1e-3;
    for &(u, v) in &[(0.3_f64, 0.4_f64), (0.6, 0.7), (0.5, 0.5)] {
        let analytic = guide.third_partial_uuv(u, v);
        // FD: ∂³uuv ≈ [S_uu(u, v+h) - S_uu(u, v-h)] / (2h)
        let uup = guide.second_partial_uu(u, v + h);
        let uum = guide.second_partial_uu(u, v - h);
        let fd = Vec3::new(
            (uup.x - uum.x) / (2.0 * h),
            (uup.y - uum.y) / (2.0 * h),
            (uup.z - uum.z) / (2.0 * h),
        );
        assert!(vec_max(analytic, fd) < 1e-3, "uuv FD at ({u}, {v})");
    }
}

#[test]
fn guide_fourth_partial_uuvv_matches_finite_difference() {
    let (guide, _, _, _, _) = arbitrary_guide_with_boundaries();
    let h = 1e-2;
    for &(u, v) in &[(0.3_f64, 0.4_f64), (0.5, 0.5)] {
        let analytic = guide.fourth_partial_uuvv(u, v);
        // FD: ∂⁴uuvv ≈ [S_uu(u, v+h) - 2 S_uu(u, v) + S_uu(u, v-h)] / h²
        let p = guide.second_partial_uu(u, v + h);
        let c = guide.second_partial_uu(u, v);
        let m = guide.second_partial_uu(u, v - h);
        let fd = Vec3::new(
            (p.x - 2.0 * c.x + m.x) / (h * h),
            (p.y - 2.0 * c.y + m.y) / (h * h),
            (p.z - 2.0 * c.z + m.z) / (h * h),
        );
        assert!(vec_max(analytic, fd) < 1e-1, "uuvv FD at ({u}, {v})");
    }
}
