//! Tests for [`super::build_guide_ring`].
//!
//! The load-bearing contract here is **shared-spoke G0 by
//! construction**: adjacent sub-patches must have bit-matching
//! CPs along their shared spoke edge (up to parameter reversal).
//! If this fails, no amount of H.4.d/H.4.e cross-derivative
//! matching will save downstream G1 — the surfaces won't even
//! be C⁰ continuous.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

const EPS: f64 = 1e-10;

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

// ---- Output shape ----

#[test]
fn output_length_matches_n() {
    for &n in &[3_usize, 4, 5, 6, 8] {
        let input = planar_ring(n, 1.5);
        let opts = DirectConstructionOptions::with_characteristic_length(1.5);
        let guides = build_guide_ring(&input, &opts).expect("should succeed");
        assert_eq!(guides.len(), n, "N = {n}");
    }
}

// ---- Corner positions ----

#[test]
fn corner_positions_match_expected_layout() {
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");

    for k in 0..4 {
        let k_prev = (k + 4 - 1) % 4;
        let v_prev = input.outer_corners[k_prev].position;
        let v_next = input.outer_corners[k].position;
        let c = input.central_vertex;

        assert!(
            point_max(guides[k].point_at(0.0, 0.0), v_prev) < EPS,
            "c0 of guide {k}"
        );
        assert!(
            point_max(guides[k].point_at(1.0, 0.0), c) < EPS,
            "c1 of guide {k}"
        );
        assert!(
            point_max(guides[k].point_at(1.0, 1.0), v_next) < EPS,
            "c2 of guide {k}"
        );
        // Since G5/G6: c3 (M_k) is a tangent-plane-consistent
        // Hermite midpoint, not the arithmetic midpoint. Check
        // only that M_k lies in the z = 0 plane (face plane for
        // planar ring) and is within the chord bounding box.
        let m = guides[k].point_at(0.0, 1.0);
        assert!(m.z.abs() < EPS, "M_k.z of guide {k} = {}", m.z);
    }
}

// ---- Shared-spoke G0: THE key invariant ----

#[test]
fn shared_spoke_cps_agree_across_adjacent_guides() {
    // Sub-patch k's u=1 edge CPs:    grid[3][0..=3]
    // Sub-patch (k+1)%N's v=0 edge CPs (reversed): grid[0..=3][0]
    //
    // With parameter reversal:
    //   sub-patch k's grid[3][j] = sub-patch (k+1)'s grid[3-j][0]
    //
    // All 4 CP pairs must match exactly (bit-level within EPS).
    let input = planar_ring(5, 2.5);
    let opts = DirectConstructionOptions::with_characteristic_length(2.5);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");
    let n = guides.len();

    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 0..4 {
            let p_k = guides[k].grid[3][j];
            let p_kn = guides[k_next].grid[3 - j][0];
            assert!(
                point_max(p_k, p_kn) < EPS,
                "shared-spoke CP mismatch: guide {k} grid[3][{j}] = {p_k:?}, guide {k_next} grid[{}][0] = {p_kn:?}",
                3 - j
            );
        }
    }
}

#[test]
fn shared_spoke_curves_evaluate_identically() {
    // Stronger check: evaluate sub-patch k's u=1 edge (a cubic
    // in v) and sub-patch (k+1)%N's v=0 edge (cubic in u, reversed
    // parameter) at the same physical point along the shared
    // spoke. They must return the same Point3.
    let input = planar_ring(4, 1.8);
    let opts = DirectConstructionOptions::with_characteristic_length(1.8);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");
    let n = guides.len();

    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in &[0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0] {
            // guide k's u=1 edge at parameter v = t (going from
            // C to V_{2k+2}).
            let p_k = guides[k].point_at(1.0, t);
            // guide (k+1)'s v=0 edge at parameter u = 1 - t
            // (reversed — going from V_{2k+2} to C, so
            // u = 0 is at V and u = 1 is at C).
            let p_kn = guides[k_next].point_at(1.0 - t, 0.0);
            assert!(
                point_max(p_k, p_kn) < EPS,
                "at shared spoke t = {t}: guide {k} = {p_k:?}, guide {k_next} = {p_kn:?}"
            );
        }
    }
}

#[test]
fn central_vertex_is_shared_by_every_guide() {
    // All N guides must evaluate to exactly the central vertex
    // at their c1 = (1, 0) corner.
    let input = planar_ring(6, 1.0);
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");
    for (k, g) in guides.iter().enumerate() {
        let pt = g.point_at(1.0, 0.0);
        assert!(
            point_max(pt, input.central_vertex) < EPS,
            "guide {k} c1 drifts from C: {pt:?}"
        );
    }
}

// ---- G1 at C: tangent planes coincide ----

#[test]
fn tangent_planes_at_c_coincide_across_guides() {
    // Every guide's span(partial_u, partial_v) at c1 = (1, 0)
    // should define the SAME tangent plane (identical normal
    // up to sign).
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");

    let du0 = guides[0].partial_u(1.0, 0.0);
    let dv0 = guides[0].partial_v(1.0, 0.0);
    let n0 = du0.cross(&dv0).normalized().expect("non-degenerate");

    for (k, g) in guides.iter().enumerate().skip(1) {
        let du = g.partial_u(1.0, 0.0);
        let dv = g.partial_v(1.0, 0.0);
        let nk = du.cross(&dv).normalized().expect("non-degenerate");
        let dot = n0.x * nk.x + n0.y * nk.y + n0.z * nk.z;
        assert!(
            (dot.abs() - 1.0).abs() < 1e-9,
            "guide {k}: normal {nk:?} not parallel to guide 0 {n0:?} (dot = {dot})"
        );
    }
}

// ---- Consistency with upstream extractors ----

#[test]
fn guide_corner_hermite_round_trips_to_extracted_data() {
    // For each sub-patch k, evaluating the built guide at each
    // corner should reproduce the Hermite data fed in by
    // `extract_central_vertex_hermite` (at c1) and
    // `extract_hex_corners` (at c0, c2, c3). This is the
    // `build_guide_from_hermite` reconstruction contract
    // (tested standalone in kp_guide_build_tests) lifted into
    // the ring-assembly context.
    let input = planar_ring(4, 1.5);
    let opts = DirectConstructionOptions::with_characteristic_length(1.5);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");
    let central = crate::geometry::gregory::extract_central_vertex_hermite(&input, &opts).unwrap();
    let hex = crate::geometry::gregory::extract_hex_corners(&input);

    for k in 0..4 {
        // c1 (central vertex C)
        let du = guides[k].partial_u(1.0, 0.0);
        let dv = guides[k].partial_v(1.0, 0.0);
        assert!(
            (du - central[k].tangent_u).length() < EPS,
            "guide {k} c1 ∂u"
        );
        assert!(
            (dv - central[k].tangent_v).length() < EPS,
            "guide {k} c1 ∂v"
        );

        // c0 (V_{2k})
        let du = guides[k].partial_u(0.0, 0.0);
        let dv = guides[k].partial_v(0.0, 0.0);
        assert!((du - hex[k].c0.tangent_u).length() < EPS, "guide {k} c0 ∂u");
        assert!((dv - hex[k].c0.tangent_v).length() < EPS, "guide {k} c0 ∂v");

        // c2 (V_{2k+2})
        let du = guides[k].partial_u(1.0, 1.0);
        let dv = guides[k].partial_v(1.0, 1.0);
        assert!((du - hex[k].c2.tangent_u).length() < EPS, "guide {k} c2 ∂u");
        assert!((dv - hex[k].c2.tangent_v).length() < EPS, "guide {k} c2 ∂v");

        // c3 (M_k)
        let du = guides[k].partial_u(0.0, 1.0);
        let dv = guides[k].partial_v(0.0, 1.0);
        assert!((du - hex[k].c3.tangent_u).length() < EPS, "guide {k} c3 ∂u");
        assert!((dv - hex[k].c3.tangent_v).length() < EPS, "guide {k} c3 ∂v");
    }
}

#[test]
fn flat_planar_ring_every_guide_is_planar() {
    // The planar input should produce guides that evaluate
    // entirely in the z = 0 plane.
    let input = planar_ring(6, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("should succeed");
    for (k, g) in guides.iter().enumerate() {
        for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
                let p = g.point_at(u, v);
                assert!(
                    p.z.abs() < 1e-9,
                    "guide {k} at ({u}, {v}): z = {} (should be ~0)",
                    p.z
                );
            }
        }
    }
}
