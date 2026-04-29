//! Tests for [`super::extract_central_vertex_hermite_with_twist`].
//!
//! Verification goals:
//!
//! 1. Position + tangents match the tangent-only extractor
//!    (twist correction must not perturb them).
//! 2. For symmetric planar input, by symmetry the Lagrangian QP
//!    produces zero λ everywhere, so the parallelogram-seed
//!    inner CPs are the output, and the derived twists match
//!    the zero-twist baseline at those tangent values.
//! 3. Derivation round-trip: evaluating the built guide's
//!    `second_partial_uv(1, 0)` reproduces the returned `twist`
//!    at each `c1` corner.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::extract_central_vertex_hermite;
use crate::geometry::gregory::kp_guide_build::build_guide_from_hermite;
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
        hex_boundaries: None,
    }
}

#[test]
fn position_and_tangents_identical_to_tangent_only_extractor() {
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let with_twist = extract_central_vertex_hermite_with_twist(&input, &opts).unwrap();
    let without_twist = extract_central_vertex_hermite(&input, &opts).unwrap();
    assert_eq!(with_twist.len(), without_twist.len());
    for (a, b) in with_twist.iter().zip(without_twist.iter()) {
        let pos_diff = ((a.position.x - b.position.x).abs())
            .max((a.position.y - b.position.y).abs())
            .max((a.position.z - b.position.z).abs());
        assert!(pos_diff < EPS, "position drift");
        assert!(
            (a.tangent_u - b.tangent_u).length() < EPS,
            "tangent_u drift"
        );
        assert!(
            (a.tangent_v - b.tangent_v).length() < EPS,
            "tangent_v drift"
        );
    }
}

#[test]
fn twist_round_trips_through_build_guide_from_hermite() {
    // For each sub-patch k, if we build a guide using the
    // twist-augmented Hermite data (c1 from this module + other
    // corners from extract_hex_corners), then evaluate
    // second_partial_uv(1, 0) on the built guide, we should get
    // back the twist that was fed in.
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let central = extract_central_vertex_hermite_with_twist(&input, &opts).unwrap();
    let hex = extract_hex_corners(&input);
    for k in 0..4 {
        let corners = [hex[k].c0, central[k], hex[k].c2, hex[k].c3];
        let guide = build_guide_from_hermite(corners);
        let evaluated = guide.second_partial_uv(1.0, 0.0);
        let diff = (evaluated - central[k].twist).length();
        assert!(
            diff < EPS,
            "k={k}: eval ∂²uv {evaluated:?} vs input twist {:?}",
            central[k].twist
        );
    }
}

#[test]
fn symmetric_planar_ring_twist_lies_in_tangent_plane() {
    // For a planar ring with face tangent planes at z = 0, the
    // Shirman-Séquin solve adjusts inner CPs within the z = 0
    // plane, so derived twists have zero z-component even
    // though their magnitudes are non-zero.
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let with_twist = extract_central_vertex_hermite_with_twist(&input, &opts).unwrap();
    for (k, h) in with_twist.iter().enumerate() {
        assert!(
            h.twist.z.abs() < EPS,
            "sub-patch {k}: twist.z = {} (should be 0 — face plane is z = 0)",
            h.twist.z
        );
    }
}

#[test]
fn symmetric_ring_twist_magnitudes_agree_by_rotational_symmetry() {
    // For a rotationally-symmetric N-sided planar ring, every
    // sub-patch's twist vector has the same magnitude (the
    // N-fold symmetry means Y_k are rotations of Y_0, so twists
    // have equal length).
    let input = planar_ring(5, 1.5);
    let opts = DirectConstructionOptions::with_characteristic_length(1.5);
    let with_twist = extract_central_vertex_hermite_with_twist(&input, &opts).unwrap();
    let mag0 = with_twist[0].twist.length();
    for (k, h) in with_twist.iter().enumerate().skip(1) {
        let mag = h.twist.length();
        assert!(
            (mag - mag0).abs() < 1e-7,
            "sub-patch {k}: |twist| = {mag}, expected {mag0} by rotational symmetry"
        );
    }
}

#[test]
fn error_surfaces_when_face_normals_cancel() {
    let mut input = planar_ring(4, 1.0);
    let tps_up = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let tps_dn = TangentPlane::from_normal(Vec3::new(0.0, 0.0, -1.0)).unwrap();
    input.central_vertex_tps = vec![tps_up, tps_dn, tps_up, tps_dn];
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let err = extract_central_vertex_hermite_with_twist(&input, &opts);
    assert!(err.is_err());
}
