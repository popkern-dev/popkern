//! Tests for [`super::extract_central_vertex_hermite`].
//!
//! Verification goals:
//!
//! 1. Every returned [`GuideCornerData`] has `position == C`
//!    (trivial sanity).
//! 2. After projection, `tangent_u` and `tangent_v` are both
//!    orthogonal to the common normal `n_C` (the tangent plane
//!    at `C` coincides across all N sub-patches).
//! 3. For a symmetric planar input, projection is a no-op and
//!    the tangent magnitudes reduce to the expected `7 · (cp1 −
//!    C)` forms.
//! 4. On a perturbed input (one spoke tilted out of the
//!    nominal plane), projection removes exactly the out-of-
//!    plane component and leaves the in-plane component
//!    untouched.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

const EPS: f64 = 1e-10;

/// Build a symmetric N-sided planar input centred at the origin
/// on `z = 0`. Each spoke is a straight line from `C = origin`
/// to `V_k = radius · (cos θ_k, sin θ_k, 0)` at
/// `θ_k = 2π·k / N`, expressed as 8 equidistant biseptic CPs.
/// All face tangent planes are the `z = 0` plane.
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
fn position_is_central_vertex_for_every_sub_patch() {
    let input = planar_ring(4, 1.0);
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");
    assert_eq!(hermite.len(), 4);
    for (k, h) in hermite.iter().enumerate() {
        assert!(
            (h.position.x - input.central_vertex.x).abs() < EPS
                && (h.position.y - input.central_vertex.y).abs() < EPS
                && (h.position.z - input.central_vertex.z).abs() < EPS,
            "sub-patch {k}: position drift {:?} vs C {:?}",
            h.position,
            input.central_vertex
        );
    }
}

#[test]
fn tangents_lie_in_common_plane_after_projection() {
    // For a planar input, n_C = ±(0, 0, 1). Every projected
    // tangent must have z ≈ 0 (the component along n_C).
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");
    for (k, h) in hermite.iter().enumerate() {
        // In this test n_C = (0, 0, 1) — check z-components.
        assert!(
            h.tangent_u.z.abs() < EPS,
            "sub-patch {k}: tangent_u.z = {} (should be 0 after projection)",
            h.tangent_u.z
        );
        assert!(
            h.tangent_v.z.abs() < EPS,
            "sub-patch {k}: tangent_v.z = {}",
            h.tangent_v.z
        );
    }
}

#[test]
fn twist_is_zero_in_this_phase() {
    // Twist is deferred to H.4.c.5 — initial implementation
    // sets it to zero.
    let input = planar_ring(4, 1.0);
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");
    for (k, h) in hermite.iter().enumerate() {
        assert!(h.twist.length() < EPS, "sub-patch {k}: twist {:?}", h.twist);
    }
}

#[test]
fn planar_symmetric_tangents_match_expected_magnitudes() {
    // For a planar N=4 ring with radius 2, spoke k's cps[1] is
    // at distance 2/7 along the spoke from C. At θ_k = 2π k / 4:
    //   spoke[0] direction = (+1, 0, 0) → cps[1] = (2/7, 0, 0)
    //   spoke[1] direction = (0, +1, 0) → cps[1] = (0, 2/7, 0)
    //   spoke[2] direction = (-1, 0, 0) → cps[1] = (-2/7, 0, 0)
    //   spoke[3] direction = (0, -1, 0) → cps[1] = (0, -2/7, 0)
    //
    // Guide 0's c1 Hermite (C = origin):
    //   tangent_v = 7 · (spoke[0].cps[1] − C) = (2, 0, 0)
    //   tangent_u = 7 · (C − spoke[3].cps[1]) = (0, +2, 0)
    //
    // (Guide 0 has spoke[-1 mod 4] = spoke[3] as its v=0 edge.
    // spoke[3] goes to V_3 = (0, -r, 0); reversed, its +u
    // direction at C points away from V_3 = +y direction,
    // magnitude = r = 2.)
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");

    let g0 = &hermite[0];
    assert!((g0.tangent_v.x - 2.0).abs() < EPS, "g0.tangent_v.x");
    assert!(g0.tangent_v.y.abs() < EPS);
    assert!(g0.tangent_v.z.abs() < EPS);
    assert!(g0.tangent_u.x.abs() < EPS);
    assert!((g0.tangent_u.y - 2.0).abs() < EPS, "g0.tangent_u.y");
    assert!(g0.tangent_u.z.abs() < EPS);
}

#[test]
fn projection_removes_only_out_of_plane_component() {
    // Start with the planar N=4 ring, then *tilt one spoke* out
    // of the xy-plane so its cps[1] has non-zero z. Projection
    // should kill the z component of tangents derived from this
    // spoke while leaving the in-plane part untouched.
    let mut input = planar_ring(4, 1.0);
    // Add +0.3 z to spoke 0's cps[1].
    input.spokes[0].boundary_cps[1].z += 0.3;

    // With tilted spoke 0, the tangent_v at guide 0's c1 is
    // 7 · (spoke[0].cps[1] - C) = 7 · ((1/7, 0, 0.3) - 0) =
    // (1, 0, 2.1) — RAW. After projection onto the common n_C
    // plane, z → 0. But note: the common n_C itself slightly
    // tilts because spoke 0's input-face tangent plane still
    // reports z=1 normal (we didn't alter central_vertex_tps),
    // so n_C stays ≈ (0, 0, 1) and projection is nearly
    // z-zeroing.
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");

    // Guide 1 uses spoke 0 as its v=0 boundary (reversed), so
    // tangent_u at c1 = 7·(C - spoke[0].cps[1]) = -(1, 0, 2.1).
    // Also guide 0 uses spoke 0 as its u=1 boundary → tangent_v.
    let g0_tv = hermite[0].tangent_v;
    let g1_tu = hermite[1].tangent_u;

    // The in-plane x component: g0's tangent_v should have
    // in-plane ~ (1, 0, 0) (length 1). g1's tangent_u should
    // have in-plane ~ (-1, 0, 0).
    assert!((g0_tv.x - 1.0).abs() < EPS, "g0_tv.x = {}", g0_tv.x);
    assert!(g0_tv.y.abs() < EPS);
    assert!(g0_tv.z.abs() < EPS, "g0_tv.z = {} (should be ~0)", g0_tv.z);

    assert!((g1_tu.x - (-1.0)).abs() < EPS, "g1_tu.x = {}", g1_tu.x);
    assert!(g1_tu.y.abs() < EPS);
    assert!(g1_tu.z.abs() < EPS, "g1_tu.z = {}", g1_tu.z);
}

#[test]
fn all_n_tangent_planes_coincide_after_projection() {
    // The load-bearing invariant: for every sub-patch k, the
    // span of (tangent_u, tangent_v) is the same plane (up to
    // orientation). Check this by verifying that the cross
    // product tangent_u × tangent_v is parallel to the same
    // normal for every k.
    let input = planar_ring(6, 3.5);
    let opts = DirectConstructionOptions::with_characteristic_length(3.5);
    let hermite = extract_central_vertex_hermite(&input, &opts).expect("should succeed");

    let first_cross = hermite[0].tangent_u.cross(&hermite[0].tangent_v);
    let first_normal = first_cross.normalized().expect("non-degenerate");

    for (k, h) in hermite.iter().enumerate().skip(1) {
        let cross = h.tangent_u.cross(&h.tangent_v);
        let normal = cross.normalized().expect("non-degenerate");
        // Up to ±1 sign (orientation). We're comparing unit
        // normals so `dot == ±1` iff planes coincide.
        let dot = normal.x * first_normal.x + normal.y * first_normal.y + normal.z * first_normal.z;
        assert!(
            (dot.abs() - 1.0).abs() < 1e-9,
            "sub-patch {k}: normal {normal:?} not parallel to first {first_normal:?} (dot={dot})"
        );
    }
}
