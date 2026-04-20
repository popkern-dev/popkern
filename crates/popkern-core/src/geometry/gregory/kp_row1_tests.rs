//! Tests for [`super::derive_row1_frame`].
//!
//! Verification goals:
//!
//! 1. **Corner 2×2 preservation**: the four `grid[1][1]`,
//!    `[1][8]`, `[8][1]`, `[8][8]` CPs reproduce the guide's
//!    position + u/v tangents + mixed partial at each corner
//!    when the bidegree-9 patch is constructed from the frame.
//! 2. **Mid-edge polynomial identity**: the cross-tangent
//!    behaviour at the u = 0 / u = 1 boundaries exactly matches
//!    the guide's cross-tangent field in its degree-9
//!    representation.
//! 3. **Boundary-corner consistency**: with tangent-consistent
//!    boundary curves, `row1[0]` (= `col0[1]` from boundary)
//!    equals what the row-1 formula would predict.
//! 4. **Flat-plane fidelity**: flat guide + flat boundaries →
//!    flat frame (every CP in the z = 0 plane).

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::kp_guide_build::{GuideCornerData, build_guide_from_hermite};
use crate::geometry::gregory::kp_guide_patch::KpGuidePatch;
use crate::geometry::gregory::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::{Point3, Vec3};

const EPS: f64 = 1e-9;

/// Build a non-trivial guide with known Hermite data and the
/// matching bidegree-9 boundary CPs (tangent-consistent at all
/// four corners).
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

    // Build tangent-consistent boundaries by linearly
    // interpolating between the corner tangents at the endpoints:
    // row0 is a cubic Bezier at degree 9 that goes from c0.pos to
    // c3.pos with c0.tangent_v and c3.tangent_v at the endpoints.
    // To get a 10-CP curve with the right endpoint tangents, we
    // construct a degree-9 Bezier whose first two and last two CPs
    // encode the endpoint tangents via `P[1] = P[0] + T/9` and
    // `P[n-1] = P[n] − T/9`; interior CPs are taken along the
    // straight line between P[1] and P[n-1] (any choice works —
    // this keeps the test fixture simple).
    //
    // row0 runs along u = 0 boundary of guide, parameter v.
    let row0 = build_boundary(
        corners[0].position,
        corners[0].tangent_v,
        corners[3].tangent_v,
        corners[3].position,
    );
    // row9 runs along u = 1, parameter v.
    let row9 = build_boundary(
        corners[1].position,
        corners[1].tangent_v,
        corners[2].tangent_v,
        corners[2].position,
    );
    // col0 runs along v = 0, parameter u.
    let col0 = build_boundary(
        corners[0].position,
        corners[0].tangent_u,
        corners[1].tangent_u,
        corners[1].position,
    );
    // col9 runs along v = 1, parameter u.
    let col9 = build_boundary(
        corners[3].position,
        corners[3].tangent_u,
        corners[2].tangent_u,
        corners[2].position,
    );

    (guide, row0, row9, col0, col9)
}

/// Build a degree-9 boundary curve from endpoints + endpoint
/// tangents. `P[0] = start`, `P[9] = end`, `P[1] = start +
/// tangent_start / 9`, `P[8] = end − tangent_end / 9`. Interior
/// CPs (`P[2..=7]`) are linearly interpolated between `P[1]` and
/// `P[8]`.
fn build_boundary(
    start: Point3,
    tangent_start: Vec3,
    tangent_end: Vec3,
    end: Point3,
) -> [Point3; 10] {
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[0] = start;
    cps[9] = end;
    cps[1] = Point3::new(
        start.x + tangent_start.x / 9.0,
        start.y + tangent_start.y / 9.0,
        start.z + tangent_start.z / 9.0,
    );
    cps[8] = Point3::new(
        end.x - tangent_end.x / 9.0,
        end.y - tangent_end.y / 9.0,
        end.z - tangent_end.z / 9.0,
    );
    for i in 2..=7 {
        let t = ((i - 1) as f64) / 7.0;
        cps[i] = Point3::new(
            cps[1].x + t * (cps[8].x - cps[1].x),
            cps[1].y + t * (cps[8].y - cps[1].y),
            cps[1].z + t * (cps[8].z - cps[1].z),
        );
    }
    cps
}

/// Assemble a full [`RectangularGregoryBiNonic`] from boundaries
/// and frame, with interior rows 2..=7 / cols 2..=7 left at
/// origin (the test only inspects boundary-adjacent behaviour).
fn build_binonic_from_frame(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    frame: &Row1Frame,
) -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    grid[0].copy_from_slice(row0);
    grid[9].copy_from_slice(row9);
    grid[1].copy_from_slice(&frame.row1);
    grid[8].copy_from_slice(&frame.row8);
    for i in 0..10 {
        grid[i][0] = col0[i];
        grid[i][9] = col9[i];
        grid[i][1] = frame.col1[i];
        grid[i][8] = frame.col8[i];
    }
    RectangularGregoryBiNonic { grid }
}

fn vec_max(a: Vec3, b: Vec3) -> f64 {
    (a - b).length()
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

// =====================================================================
// Corner 2 × 2 tests — the bidegree-9 patch must reproduce the
// guide's position + partials + mixed partial at each corner.
// =====================================================================

#[test]
fn bidegree9_patch_reproduces_guide_corner_positions() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for &(u, v) in &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let bidegree9_pos = binonic.point_at(u, v);
        let guide_pos = guide.point_at(u, v);
        assert!(
            point_max(bidegree9_pos, guide_pos) < EPS,
            "corner ({u}, {v}): bidegree-9 {bidegree9_pos:?} vs guide {guide_pos:?}"
        );
    }
}

#[test]
fn bidegree9_patch_reproduces_guide_corner_partials() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for &(u, v) in &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let p_u = binonic.partial_u(u, v);
        let g_u = guide.partial_u(u, v);
        let p_v = binonic.partial_v(u, v);
        let g_v = guide.partial_v(u, v);
        assert!(
            vec_max(p_u, g_u) < EPS,
            "corner ({u}, {v}) ∂u: bidegree-9 {p_u:?} vs guide {g_u:?}"
        );
        assert!(
            vec_max(p_v, g_v) < EPS,
            "corner ({u}, {v}) ∂v: bidegree-9 {p_v:?} vs guide {g_v:?}"
        );
    }
}

#[test]
fn bidegree9_patch_reproduces_guide_corner_mixed_partials() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for &(u, v) in &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
        let p_uv = binonic.second_partial_uv(u, v);
        let g_uv = guide.second_partial_uv(u, v);
        assert!(
            vec_max(p_uv, g_uv) < EPS,
            "corner ({u}, {v}) ∂²uv: bidegree-9 {p_uv:?} vs guide {g_uv:?}"
        );
    }
}

// =====================================================================
// Mid-edge polynomial identity — the bidegree-9 patch's cross-
// tangent CURVE along the u = 0 boundary must equal the guide's
// cross-tangent curve ∂G/∂u(0, v) as a polynomial in v.
// =====================================================================

#[test]
fn cross_tangent_at_u0_matches_guide_polynomial_identity() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for k in 0..=10 {
        let v = (k as f64) / 10.0;
        let bidegree9_tan = binonic.partial_u(0.0, v);
        let guide_tan = guide.partial_u(0.0, v);
        assert!(
            vec_max(bidegree9_tan, guide_tan) < 1e-8,
            "v = {v}: bidegree-9 ∂u(0, v) = {bidegree9_tan:?} vs guide = {guide_tan:?}"
        );
    }
}

#[test]
fn cross_tangent_at_u1_matches_guide_polynomial_identity() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for k in 0..=10 {
        let v = (k as f64) / 10.0;
        let bidegree9_tan = binonic.partial_u(1.0, v);
        let guide_tan = guide.partial_u(1.0, v);
        assert!(
            vec_max(bidegree9_tan, guide_tan) < 1e-8,
            "v = {v}: bidegree-9 ∂u(1, v) vs guide"
        );
    }
}

#[test]
fn cross_tangent_at_v0_matches_guide_polynomial_identity() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for k in 0..=10 {
        let u = (k as f64) / 10.0;
        let bidegree9_tan = binonic.partial_v(u, 0.0);
        let guide_tan = guide.partial_v(u, 0.0);
        assert!(
            vec_max(bidegree9_tan, guide_tan) < 1e-8,
            "u = {u}: bidegree-9 ∂v(u, 0) vs guide"
        );
    }
}

#[test]
fn cross_tangent_at_v1_matches_guide_polynomial_identity() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    let binonic = build_binonic_from_frame(&row0, &row9, &col0, &col9, &frame);
    for k in 0..=10 {
        let u = (k as f64) / 10.0;
        let bidegree9_tan = binonic.partial_v(u, 1.0);
        let guide_tan = guide.partial_v(u, 1.0);
        assert!(
            vec_max(bidegree9_tan, guide_tan) < 1e-8,
            "u = {u}: bidegree-9 ∂v(u, 1) vs guide"
        );
    }
}

// =====================================================================
// Boundary-corner consistency — with tangent-consistent inputs,
// the row-1 formula's corner predictions match the boundary.
// =====================================================================

#[test]
fn row1_corner_slots_match_boundary_cps() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    // row1[0] comes from col0[1] (by construction in derive_row1_frame);
    // verify it matches our tangent-consistent boundary.
    assert_eq!(frame.row1[0].x, col0[1].x);
    assert_eq!(frame.row1[0].y, col0[1].y);
    assert_eq!(frame.row1[0].z, col0[1].z);
    // Similar for other corner slots.
    assert_eq!(frame.row1[9].x, col9[1].x);
    assert_eq!(frame.row8[0].x, col0[8].x);
    assert_eq!(frame.row8[9].x, col9[8].x);
    assert_eq!(frame.col1[0].x, row0[1].x);
    assert_eq!(frame.col1[9].x, row9[1].x);
    assert_eq!(frame.col8[0].x, row0[8].x);
    assert_eq!(frame.col8[9].x, row9[8].x);
}

#[test]
fn corner_2x2_consistency_between_row1_and_col1() {
    // grid[1][1] appears in both row1[1] and col1[1]. They must
    // be identical — both produced by the same corner_2x2
    // formula inside derive_row1_frame.
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    assert_eq!(frame.row1[1].x, frame.col1[1].x);
    assert_eq!(frame.row1[8].x, frame.col8[1].x);
    assert_eq!(frame.row8[1].x, frame.col1[8].x);
    assert_eq!(frame.row8[8].x, frame.col8[8].x);
}

// =====================================================================
// Flat-plane fidelity
// =====================================================================

#[test]
fn flat_guide_and_flat_boundaries_produce_flat_frame() {
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
    let row0 = build_boundary(
        corners[0].position,
        corners[0].tangent_v,
        corners[3].tangent_v,
        corners[3].position,
    );
    let row9 = build_boundary(
        corners[1].position,
        corners[1].tangent_v,
        corners[2].tangent_v,
        corners[2].position,
    );
    let col0 = build_boundary(
        corners[0].position,
        corners[0].tangent_u,
        corners[1].tangent_u,
        corners[1].position,
    );
    let col9 = build_boundary(
        corners[3].position,
        corners[3].tangent_u,
        corners[2].tangent_u,
        corners[2].position,
    );
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    // Every frame CP should have z ≈ 0.
    for j in 0..10 {
        assert!(
            frame.row1[j].z.abs() < EPS,
            "row1[{j}].z = {}",
            frame.row1[j].z
        );
        assert!(
            frame.row8[j].z.abs() < EPS,
            "row8[{j}].z = {}",
            frame.row8[j].z
        );
        assert!(
            frame.col1[j].z.abs() < EPS,
            "col1[{j}].z = {}",
            frame.col1[j].z
        );
        assert!(
            frame.col8[j].z.abs() < EPS,
            "col8[{j}].z = {}",
            frame.col8[j].z
        );
    }
}

#[test]
fn output_lengths_are_ten() {
    let (guide, row0, row9, col0, col9) = consistent_guide_and_boundaries();
    let frame = derive_row1_frame(&guide, &row0, &row9, &col0, &col9);
    assert_eq!(frame.row1.len(), 10);
    assert_eq!(frame.row8.len(), 10);
    assert_eq!(frame.col1.len(), 10);
    assert_eq!(frame.col8.len(), 10);
}
