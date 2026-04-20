//! Tests for [`super::extract_hex_corners`].
//!
//! Verification goals:
//!
//! 1. Output length = N.
//! 2. `c0.position == V_{2k}`, `c2.position == V_{2k+2}`,
//!    `c3.position == M_k` for every sub-patch.
//! 3. `M_k` equals the arithmetic midpoint of its two outer
//!    vertices (the phase-defined choice).
//! 4. Hex-edge tangent pattern for a **symmetric planar ring**:
//!    every hex-edge tangent lies in the `z = 0` plane; left-
//!    edge and top-edge tangents have exactly the expected
//!    magnitudes (`|M_k − V_{2k}|`, `|V_{2k+2} − M_k|`).
//! 5. **Shared-point invariant**: at each outer vertex `V_{2k+2}`
//!    (which is sub-patch `k`'s `c2` and sub-patch `(k+1)%N`'s
//!    `c0`), the two adjacent sub-patches' along-spoke tangents
//!    are exact negatives — the `extract_outer_corner_hermite`
//!    contract is preserved.
//! 6. `c3`'s tangents are raw (unprojected) — at `M_k` we have
//!    no input tangent plane, so the straight-line tangents go
//!    through unchanged.
//! 7. Twist is zero at every corner (phase invariant).

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

const EPS: f64 = 1e-10;

/// Re-use the `planar_ring` fixture shape from earlier tests.
fn planar_ring(n: usize, radius: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();

    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let v = Point3::new(radius * theta.cos(), radius * theta.sin(), 0.0);

        let mut cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for j in 0..8 {
            let t = (j as f64) / 7.0;
            cps[j] = Point3::new(
                c.x + t * (v.x - c.x),
                c.y + t * (v.y - c.y),
                c.z + t * (v.z - c.z),
            );
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps: cps,
            xtan_a_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_b_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_a_at_v: Vec3::new(0.0, 0.0, 1.0),
            xtan_b_at_v: Vec3::new(0.0, 0.0, 1.0),
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("e{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tps_xy,
        });
    }

    DirectBoundaryInput {
        central_vertex: c,
        spokes,
        outer_corners,
        central_vertex_tps: vec![tps_xy; n],
        seed_patches: None,
    }
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[test]
fn output_length_matches_n() {
    for &n in &[3_usize, 4, 5, 6, 8] {
        let input = planar_ring(n, 1.5);
        let corners = extract_hex_corners(&input);
        assert_eq!(corners.len(), n);
    }
}

#[test]
fn positions_match_expected_vertices() {
    // Since G5/G6: M_k is no longer the arithmetic midpoint —
    // it's a tangent-plane-consistent Hermite midpoint. We check
    // the outer-vertex positions only here; `m_k_lies_in_face_tangent_planes`
    // validates the refined midpoint construction.
    let input = planar_ring(4, 2.0);
    let corners = extract_hex_corners(&input);
    for k in 0..4 {
        let k_prev = (k + 4 - 1) % 4;
        let v_prev = input.outer_corners[k_prev].position;
        let v_next = input.outer_corners[k].position;
        assert!(
            point_max(corners[k].c0.position, v_prev) < EPS,
            "c0 pos k={k}"
        );
        assert!(
            point_max(corners[k].c2.position, v_next) < EPS,
            "c2 pos k={k}"
        );
    }
}

#[test]
fn m_k_lies_in_face_tangent_planes_after_refinement() {
    // G5 + G6: M_k is constructed so the cubic hex-edge Bezier's
    // endpoint tangents lie in the respective input face tangent
    // planes. For a planar ring with all face planes at z = 0, M_k
    // itself has z ≈ 0.
    let input = planar_ring(5, 1.5);
    let corners = extract_hex_corners(&input);
    for (k, c) in corners.iter().enumerate() {
        assert!(
            c.c3.position.z.abs() < EPS,
            "k={k}: M_k.z = {} (planar input expects ~0)",
            c.c3.position.z
        );
    }
}

#[test]
fn m_k_stays_between_adjacent_outer_vertices() {
    // Sanity: M_k should be in the convex hull of V_prev and V_next
    // (not wildly displaced). For a symmetric planar ring with radius
    // R, |M_k| ≤ R + small (bounded by the Hermite correction term).
    let radius = 2.0;
    let input = planar_ring(6, radius);
    let corners = extract_hex_corners(&input);
    for (k, c) in corners.iter().enumerate() {
        let d =
            (c.c3.position.x.powi(2) + c.c3.position.y.powi(2) + c.c3.position.z.powi(2)).sqrt();
        assert!(d <= radius + 0.1, "k={k}: |M_k| = {d} exceeds radius + 0.1");
    }
}

#[test]
fn hex_edge_tangents_lie_in_face_plane_at_outer_vertices() {
    // After projection onto the z=0 face plane, c0.tangent_v and
    // c2.tangent_u have z ≈ 0.
    let input = planar_ring(5, 1.5);
    let corners = extract_hex_corners(&input);
    for (k, c) in corners.iter().enumerate() {
        assert!(c.c0.tangent_v.z.abs() < EPS, "c0.tangent_v.z at k={k}");
        assert!(c.c2.tangent_u.z.abs() < EPS, "c2.tangent_u.z at k={k}");
    }
}

#[test]
fn c3_tangents_project_onto_blended_m_k_plane() {
    // G6: c3.tangent_u and c3.tangent_v are projected onto the
    // blended face tangent plane at M_k. For a planar ring the
    // blended plane is z = 0, so projection is a no-op and the
    // tangents equal the straight-line `V_next − M_k` and
    // `M_k − V_prev`. For a non-planar ring this test would
    // measure the projection effect (see the integration tests
    // in `kp_row1_tests.rs` / `kp_row2_tests.rs` for non-planar
    // cases).
    let input = planar_ring(4, 2.0);
    let corners = extract_hex_corners(&input);
    for k in 0..4 {
        let k_prev = (k + 4 - 1) % 4;
        let v_prev = input.outer_corners[k_prev].position;
        let v_next = input.outer_corners[k].position;
        let m = corners[k].hex_midpoint;
        let expected_tu = Vec3::new(v_next.x - m.x, v_next.y - m.y, v_next.z - m.z);
        let expected_tv = Vec3::new(m.x - v_prev.x, m.y - v_prev.y, m.z - v_prev.z);
        let diff_tu = (corners[k].c3.tangent_u - expected_tu).length();
        let diff_tv = (corners[k].c3.tangent_v - expected_tv).length();
        assert!(diff_tu < EPS, "c3 tangent_u k={k}: diff {diff_tu}");
        assert!(diff_tv < EPS, "c3 tangent_v k={k}: diff {diff_tv}");
    }
}

#[test]
fn c0_and_c2_along_spoke_tangents_agree_with_outer_extractor() {
    // c0.tangent_u should equal the along-spoke tangent at V_{2k}
    // (from outer-corner extraction at the previous index).
    // c2.tangent_v should equal the along-spoke tangent at V_{2k+2}.
    let input = planar_ring(4, 1.5);
    let outer = crate::geometry::gregory::extract_outer_corner_hermite(&input);
    let corners = extract_hex_corners(&input);
    for k in 0..4 {
        let k_prev = (k + 4 - 1) % 4;
        let c0_tu = corners[k].c0.tangent_u;
        let expected_c0_tu = outer[k_prev].patch_kp1_c0.tangent_u;
        let diff = (c0_tu - expected_c0_tu).length();
        assert!(diff < EPS, "c0.tangent_u k={k}: diff {diff}");

        let c2_tv = corners[k].c2.tangent_v;
        let expected_c2_tv = outer[k].patch_k_c2.tangent_v;
        let diff = (c2_tv - expected_c2_tv).length();
        assert!(diff < EPS, "c2.tangent_v k={k}: diff {diff}");
    }
}

#[test]
fn shared_vertex_along_spoke_tangents_remain_exact_negatives() {
    // At each outer vertex V_{2k+2}:
    //   sub-patch k's c2.tangent_v (along spoke) should equal
    //   -(sub-patch (k+1)%N's c0.tangent_u).
    let input = planar_ring(5, 1.8);
    let corners = extract_hex_corners(&input);
    let n = corners.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        let t_c2 = corners[k].c2.tangent_v;
        let t_c0 = corners[k_next].c0.tangent_u;
        assert!(
            (t_c2.x + t_c0.x).abs() < EPS,
            "x sum k={k}: {}",
            t_c2.x + t_c0.x
        );
        assert!((t_c2.y + t_c0.y).abs() < EPS, "y sum k={k}");
        assert!((t_c2.z + t_c0.z).abs() < EPS, "z sum k={k}");
    }
}

#[test]
fn twists_are_zero() {
    let input = planar_ring(4, 1.0);
    let corners = extract_hex_corners(&input);
    for c in &corners {
        assert!(c.c0.twist.length() < EPS);
        assert!(c.c2.twist.length() < EPS);
        assert!(c.c3.twist.length() < EPS);
    }
}

#[test]
fn planar_n6_hex_tangents_lie_in_xy_plane() {
    // For a planar ring (all face normals = +z), every hex-edge
    // tangent must lie in the z = 0 plane regardless of the
    // refined M_k position. (M_k itself lies in the z = 0 plane
    // by the G5 + G6 construction.)
    let input = planar_ring(6, 2.0);
    let corners = extract_hex_corners(&input);
    for k in 0..6 {
        assert!(corners[k].c0.tangent_v.z.abs() < EPS, "c0 k={k}");
        assert!(corners[k].c2.tangent_u.z.abs() < EPS, "c2 k={k}");
        assert!(corners[k].c3.tangent_u.z.abs() < EPS, "c3_u k={k}");
        assert!(corners[k].c3.tangent_v.z.abs() < EPS, "c3_v k={k}");
        assert!(corners[k].c3.position.z.abs() < EPS, "M_k.z k={k}");
    }
}
