//! Tests for [`super::extract_hex_corners_with_twist`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::gregory::kp_guide_build::build_guide_from_hermite;
use crate::geometry::gregory::kp_guide_central_twist::extract_central_vertex_hermite_with_twist;
use crate::geometry::gregory::kp_guide_hex_hermite::extract_hex_corners;

const EPS: f64 = 1e-9;

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

#[test]
fn position_and_tangents_unchanged_from_baseline() {
    let input = planar_ring(5, 2.0);
    let twisted = extract_hex_corners_with_twist(&input);
    let baseline = extract_hex_corners(&input);
    for k in 0..5 {
        // Positions and tangents should be bit-identical — only
        // twist changes.
        assert_eq!(twisted[k].c0.position.x, baseline[k].c0.position.x);
        assert_eq!(twisted[k].c2.position.x, baseline[k].c2.position.x);
        assert!((twisted[k].c0.tangent_u - baseline[k].c0.tangent_u).length() < EPS);
        assert!((twisted[k].c0.tangent_v - baseline[k].c0.tangent_v).length() < EPS);
        assert!((twisted[k].c2.tangent_u - baseline[k].c2.tangent_u).length() < EPS);
        assert!((twisted[k].c2.tangent_v - baseline[k].c2.tangent_v).length() < EPS);
        assert_eq!(twisted[k].c3.position.x, baseline[k].c3.position.x);
    }
}

#[test]
fn c3_twist_remains_zero() {
    // Hex midpoint twist is still zero — no pairwise G1
    // constraint applies there.
    let input = planar_ring(6, 2.5);
    let twisted = extract_hex_corners_with_twist(&input);
    for c in &twisted {
        assert!(c.c3.twist.length() < EPS);
    }
}

#[test]
fn twists_lie_in_face_tangent_plane_at_outer_vertex() {
    // Face tangent plane at V is z = 0 in the fixture, so
    // derived twists at c0 and c2 have z ≈ 0.
    let input = planar_ring(5, 1.8);
    let twisted = extract_hex_corners_with_twist(&input);
    for (k, c) in twisted.iter().enumerate() {
        assert!(
            c.c0.twist.z.abs() < EPS,
            "sub-patch {k}: c0.twist.z = {}",
            c.c0.twist.z
        );
        assert!(
            c.c2.twist.z.abs() < EPS,
            "sub-patch {k}: c2.twist.z = {}",
            c.c2.twist.z
        );
    }
}

#[test]
fn shared_vertex_adjacent_twists_satisfy_pairwise_constraint() {
    // For adjacent sub-patches k and (k+1)%N meeting at V:
    //   (Y_k − Y_{k+1}) · p_spoke = 0
    // in the face tangent plane at V. Recover Y_k, Y_{k+1}
    // from the twists and verify.
    let input = planar_ring(4, 2.0);
    let twisted = extract_hex_corners_with_twist(&input);
    let n = twisted.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        let v = input.outer_corners[k].position;
        let c2_k = twisted[k].c2;
        let c0_kn = twisted[k_next].c0;
        // grid[2][2] at c2: V − tu/3 − tv/3 + twist/9
        let y_k = Point3::new(
            v.x - c2_k.tangent_u.x / 3.0 - c2_k.tangent_v.x / 3.0 + c2_k.twist.x / 9.0,
            v.y - c2_k.tangent_u.y / 3.0 - c2_k.tangent_v.y / 3.0 + c2_k.twist.y / 9.0,
            v.z - c2_k.tangent_u.z / 3.0 - c2_k.tangent_v.z / 3.0 + c2_k.twist.z / 9.0,
        );
        // grid[1][1] at c0: V + tu/3 + tv/3 + twist/9
        let y_kn = Point3::new(
            v.x + c0_kn.tangent_u.x / 3.0 + c0_kn.tangent_v.x / 3.0 + c0_kn.twist.x / 9.0,
            v.y + c0_kn.tangent_u.y / 3.0 + c0_kn.tangent_v.y / 3.0 + c0_kn.twist.y / 9.0,
            v.z + c0_kn.tangent_u.z / 3.0 + c0_kn.tangent_v.z / 3.0 + c0_kn.twist.z / 9.0,
        );
        // Spoke direction at V (from c2_k.tangent_v, which is
        // the along-spoke at V pointing from C toward extension).
        // Perpendicular in 3D to the spoke AND in the face plane
        // at V = normal_V × spoke_direction (then normalised).
        let sigma = c2_k.tangent_v;
        let n_v = input.outer_corners[k].tangent_plane.normal();
        let p = n_v.cross(&sigma).normalized().unwrap();

        let diff = Vec3::new(y_k.x - y_kn.x, y_k.y - y_kn.y, y_k.z - y_kn.z);
        let dot = diff.dot(&p);
        assert!(
            dot.abs() < EPS,
            "k={k}: (Y_k − Y_{{k+1}}) · p_spoke = {dot} (should be 0)"
        );
    }
}

#[test]
fn twist_round_trips_through_build_guide_from_hermite() {
    // Build each guide using the twist-augmented Hermite data.
    // Evaluating second_partial_uv at c2 = (1, 1) and c0 =
    // (0, 0) should reproduce the input twists.
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let central = extract_central_vertex_hermite_with_twist(&input, &opts).unwrap();
    let hex = extract_hex_corners_with_twist(&input);
    for k in 0..4 {
        let corners = [hex[k].c0, central[k], hex[k].c2, hex[k].c3];
        let guide = build_guide_from_hermite(corners);
        let uv_c0 = guide.second_partial_uv(0.0, 0.0);
        let uv_c2 = guide.second_partial_uv(1.0, 1.0);
        assert!(
            (uv_c0 - hex[k].c0.twist).length() < EPS,
            "k={k} c0 twist round-trip"
        );
        assert!(
            (uv_c2 - hex[k].c2.twist).length() < EPS,
            "k={k} c2 twist round-trip"
        );
    }
}
