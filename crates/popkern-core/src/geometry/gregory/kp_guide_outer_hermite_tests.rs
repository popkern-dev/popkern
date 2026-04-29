//! Tests for [`super::extract_outer_corner_hermite`].
//!
//! Verification goals:
//!
//! 1. Output length equals `N` (one entry per outer hex vertex).
//! 2. Every entry's `patch_k_c2.position` and
//!    `patch_kp1_c0.position` equal the outer vertex (trivial
//!    sanity).
//! 3. The along-spoke tangents at `c2` and `c0` are **negatives**
//!    of each other (both pre- and post-projection — a natural
//!    consequence of the parameter reversal on the shared spoke).
//! 4. Both along-spoke tangents lie in the input face tangent
//!    plane at `V` (projection zeros the out-of-plane component).
//! 5. Hex-edge tangents (`tangent_u` at `c2`, `tangent_v` at
//!    `c0`) are zero — placeholders for H.4.c.4.
//! 6. Twist is zero at every corner (phase invariant).

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

const EPS: f64 = 1e-10;

/// Build a symmetric N-sided planar input centred at the origin
/// on `z = 0`. Re-uses the same shape as the central-hermite
/// tests — straight spokes from `C` to each outer vertex, all
/// face tangent planes the `z = 0` plane.
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
        hex_boundaries: None,
    }
}

#[test]
fn output_length_matches_n() {
    for &n in &[3_usize, 4, 5, 6, 8] {
        let input = planar_ring(n, 1.5);
        let hermite = extract_outer_corner_hermite(&input);
        assert_eq!(hermite.len(), n, "N = {n}");
    }
}

#[test]
fn positions_match_outer_vertices() {
    let input = planar_ring(4, 2.0);
    let hermite = extract_outer_corner_hermite(&input);
    for (k, entry) in hermite.iter().enumerate() {
        let v_expected = input.outer_corners[k].position;
        assert!(
            (entry.outer_vertex.x - v_expected.x).abs() < EPS
                && (entry.outer_vertex.y - v_expected.y).abs() < EPS
                && (entry.outer_vertex.z - v_expected.z).abs() < EPS,
            "entry {k}: outer_vertex drift"
        );
        assert!(
            (entry.patch_k_c2.position.x - v_expected.x).abs() < EPS
                && (entry.patch_k_c2.position.y - v_expected.y).abs() < EPS,
            "entry {k}: c2 position drift"
        );
        assert!(
            (entry.patch_kp1_c0.position.x - v_expected.x).abs() < EPS
                && (entry.patch_kp1_c0.position.y - v_expected.y).abs() < EPS,
            "entry {k}: c0 position drift"
        );
        assert_eq!(entry.spoke_index, k);
    }
}

#[test]
fn along_spoke_tangents_are_negatives_of_each_other() {
    let input = planar_ring(5, 2.5);
    let hermite = extract_outer_corner_hermite(&input);
    for (k, entry) in hermite.iter().enumerate() {
        // patch_k_c2.tangent_v and patch_kp1_c0.tangent_u both
        // describe the same spoke at V — they must be negatives.
        let t_v = entry.patch_k_c2.tangent_v;
        let t_u = entry.patch_kp1_c0.tangent_u;
        assert!(
            (t_v.x + t_u.x).abs() < EPS,
            "entry {k}: x sum {}",
            t_v.x + t_u.x
        );
        assert!(
            (t_v.y + t_u.y).abs() < EPS,
            "entry {k}: y sum {}",
            t_v.y + t_u.y
        );
        assert!(
            (t_v.z + t_u.z).abs() < EPS,
            "entry {k}: z sum {}",
            t_v.z + t_u.z
        );
    }
}

#[test]
fn along_spoke_tangents_lie_in_face_tangent_plane() {
    let input = planar_ring(6, 1.5);
    let hermite = extract_outer_corner_hermite(&input);
    for (k, entry) in hermite.iter().enumerate() {
        // n_V = (0, 0, 1) for the planar ring; the z-components
        // must be zero after projection.
        assert!(
            entry.patch_k_c2.tangent_v.z.abs() < EPS,
            "entry {k}: c2 tangent_v.z = {}",
            entry.patch_k_c2.tangent_v.z
        );
        assert!(
            entry.patch_kp1_c0.tangent_u.z.abs() < EPS,
            "entry {k}: c0 tangent_u.z = {}",
            entry.patch_kp1_c0.tangent_u.z
        );
    }
}

#[test]
fn hex_edge_tangents_are_placeholders() {
    // tangent_u at c2 and tangent_v at c0 are zero in this
    // phase — filled in H.4.c.4.
    let input = planar_ring(4, 1.0);
    let hermite = extract_outer_corner_hermite(&input);
    for entry in &hermite {
        assert!(entry.patch_k_c2.tangent_u.length() < EPS);
        assert!(entry.patch_kp1_c0.tangent_v.length() < EPS);
    }
}

#[test]
fn twists_are_zero() {
    let input = planar_ring(4, 1.0);
    let hermite = extract_outer_corner_hermite(&input);
    for entry in &hermite {
        assert!(entry.patch_k_c2.twist.length() < EPS);
        assert!(entry.patch_kp1_c0.twist.length() < EPS);
    }
}

#[test]
fn planar_symmetric_magnitudes_match_hand_derivation() {
    // N = 4, radius 2. Each spoke k is a straight line from C
    // to V_k. Spoke k's cps[6] is at distance 6/7 along the
    // spoke from C, so:
    //   cps[6] = V_k · (6/7)   (for C at origin, radius = 2)
    //
    // along-spoke tangent at V for patch_k_c2:
    //   t_v = 7 · (V_k − cps[6]) = 7 · (V_k − (6/7) V_k) = V_k
    //
    // For k = 0: V_0 = (2, 0, 0), so t_v should be (2, 0, 0).
    let input = planar_ring(4, 2.0);
    let hermite = extract_outer_corner_hermite(&input);
    let entry0 = &hermite[0];
    assert!((entry0.patch_k_c2.tangent_v.x - 2.0).abs() < EPS);
    assert!(entry0.patch_k_c2.tangent_v.y.abs() < EPS);
    assert!(entry0.patch_k_c2.tangent_v.z.abs() < EPS);
    // Mirror: c0 of patch 1 has negated tangent_u.
    assert!((entry0.patch_kp1_c0.tangent_u.x + 2.0).abs() < EPS);
    assert!(entry0.patch_kp1_c0.tangent_u.y.abs() < EPS);
    assert!(entry0.patch_kp1_c0.tangent_u.z.abs() < EPS);
}

#[test]
fn projection_with_tilted_face_plane_removes_normal_component() {
    // Tilt one outer corner's face tangent plane so its normal
    // has a non-z component. Projection should still zero the
    // component along the tilted normal.
    let mut input = planar_ring(4, 2.0);
    // Tilt outer_corners[0]'s face plane: normal rotated 30° in
    // the x-z plane toward +x.
    let theta = 30.0_f64.to_radians();
    let tilted_normal = Vec3::new(theta.sin(), 0.0, theta.cos());
    input.outer_corners[0].tangent_plane = TangentPlane::from_normal(tilted_normal).unwrap();

    let hermite = extract_outer_corner_hermite(&input);
    let t_v = hermite[0].patch_k_c2.tangent_v;
    // Projection: dot product with the normal should be zero.
    let dot = t_v.x * tilted_normal.x + t_v.y * tilted_normal.y + t_v.z * tilted_normal.z;
    assert!(dot.abs() < EPS, "projection failed: dot = {dot}");
}
