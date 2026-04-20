//! Integration tests for [`super::build_binonic_ring`] — the
//! full H.4 pipeline end-to-end on a DirectBoundaryInput.
//!
//! Tests the composition of H.4.c + H.4.d + H.4.e + shared-spoke
//! CP averaging (G4), measuring exact G1/G2 across shared spokes
//! on the emitted bidegree-9 patches.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

const EPS: f64 = 1e-9;

fn planar_ring(n: usize, radius: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
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
            tangent_plane: tps,
        });
    }
    DirectBoundaryInput {
        central_vertex: c,
        spokes,
        outer_corners,
        central_vertex_tps: vec![tps; n],
        seed_patches: None,
    }
}

/// Non-planar consistent fixture: outer hex z-alternates, face
/// tangent planes at V tilted to contain the spoke direction.
fn consistent_3d_ring(n: usize, radius: f64, z_amp: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps_c = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let z_v = z_amp * if k % 2 == 0 { 1.0 } else { -1.0 };
        let v = Point3::new(radius * theta.cos(), radius * theta.sin(), z_v);
        let mut cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for j in 0..8 {
            let t = (j as f64) / 7.0;
            cps[j] = Point3::new(
                c.x + t * (v.x - c.x),
                c.y + t * (v.y - c.y),
                c.z + t * (v.z - c.z),
            );
        }
        let dir = Vec3::new(v.x - c.x, v.y - c.y, v.z - c.z)
            .normalized()
            .unwrap();
        let nv_raw = Vec3::new(
            -dir.z * dir.x,
            -dir.z * dir.y,
            dir.x * dir.x + dir.y * dir.y,
        );
        let nv = nv_raw.normalized().unwrap();
        let tps_v = TangentPlane::from_normal(nv).unwrap();
        spokes.push(SpokeBoundaryData {
            boundary_cps: cps,
            xtan_a_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_b_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_a_at_v: nv,
            xtan_b_at_v: nv,
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("e{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tps_v,
        });
    }
    DirectBoundaryInput {
        central_vertex: c,
        spokes,
        outer_corners,
        central_vertex_tps: vec![tps_c; n],
        seed_patches: None,
    }
}

// =====================================================================
// Shape + count + self-intersection
// =====================================================================

#[test]
fn output_length_matches_n_for_various_valences() {
    for &n in &[3_usize, 4, 5, 6, 8] {
        let input = planar_ring(n, 2.0);
        let opts = DirectConstructionOptions::with_characteristic_length(2.0);
        let patches = build_binonic_ring(&input, &opts).expect("build");
        assert_eq!(patches.len(), n, "N = {n}");
    }
}

// =====================================================================
// Shared-spoke G0 (position continuity)
// =====================================================================

#[test]
fn shared_spoke_g0_cps_match_across_adjacent_subpatches() {
    // Sub-patch k's u = 1 edge (grid[9][j]) and sub-patch (k+1)'s
    // v = 0 edge (grid[i][0]) represent the same physical curve.
    // They must match bit-level after reversal.
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    for k in 0..n {
        let kn = (k + 1) % n;
        for j in 0..10 {
            let a = patches[k].grid[9][j];
            let b = patches[kn].grid[9 - j][0];
            let d = ((a.x - b.x).abs())
                .max((a.y - b.y).abs())
                .max((a.z - b.z).abs());
            assert!(
                d < EPS,
                "k={k} j={j}: grid[9][{j}] ({a:?}) vs grid[{}][0] ({b:?})",
                9 - j
            );
        }
    }
}

// =====================================================================
// Shared-spoke G1 (cross-tangent coplanarity after averaging)
// =====================================================================

/// Scalar triple product `a · (b × c)`. Zero iff coplanar.
fn triple_product(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.dot(&b.cross(&c))
}

#[test]
fn shared_spoke_g1_coplanarity_on_planar() {
    // G1 at the shared spoke = spoke tangent, patch-k cross,
    // patch-(k+1) cross are coplanar. Measured by the triple
    // product; should be exactly zero for the planar ring (all
    // vectors in the xy plane → triple product always zero).
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    for k in 0..n {
        let kn = (k + 1) % n;
        for &t in &[0.05_f64, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95] {
            let spoke = patches[k].partial_v(1.0, t);
            let cross_k = patches[k].partial_u(1.0, t);
            let cross_kn = patches[kn].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_kn);
            assert!(
                tp.abs() < 1e-8,
                "k={k} t={t}: triple product = {tp:e} (should be ~0)"
            );
        }
    }
}

#[test]
fn shared_spoke_g1_coplanarity_on_consistent_3d_after_averaging() {
    // G4 averaging significantly reduces (but may not fully
    // close) the G1 gap on consistent 3D input. Records the
    // max triple product for regression tracking; acceptance
    // threshold is bounded by the commercial-grade G1 target
    // translated to triple product units.
    let input = consistent_3d_ring(4, 2.0, 0.3);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    let mut max_tp = 0.0_f64;
    for k in 0..n {
        let kn = (k + 1) % n;
        for &t in &[0.1_f64, 0.3, 0.5, 0.7, 0.9] {
            let spoke = patches[k].partial_v(1.0, t);
            let cross_k = patches[k].partial_u(1.0, t);
            let cross_kn = patches[kn].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_kn).abs();
            if tp > max_tp {
                max_tp = tp;
            }
        }
    }
    println!("consistent-3D G1 max triple product = {max_tp:e}");
    // Regression guard: bidegree-9 G1 residual on consistent-3D
    // with averaging should be well under 1.0 (at commit time
    // measured ~0.3-0.6). Shared-ribbon construction (#79)
    // would drive this to machine precision.
    assert!(max_tp < 1.0, "G1 regression on consistent-3D: {max_tp:e}");
    assert!(max_tp.is_finite());
}

// =====================================================================
// Shared-spoke G2 (cross-second-derivative matching after averaging)
// =====================================================================

#[test]
fn shared_spoke_g2_coplanarity_on_planar() {
    // G2 at the shared spoke = second-fundamental-form
    // Frobenius-norm agreement. For planar input, all second
    // partials are in the xy plane (z = 0 component for ∂²S/∂u²
    // and ∂²S/∂v²), so the z-components of both must be
    // approximately zero.
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    for k in 0..n {
        let kn = (k + 1) % n;
        for &t in &[0.05_f64, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95] {
            let s_uu_k = patches[k].second_partial_uu(1.0, t);
            let s_vv_kn = patches[kn].second_partial_vv(1.0 - t, 0.0);
            // Planar ring: z-component of second partials is 0.
            assert!(
                s_uu_k.z.abs() < 1e-8,
                "planar G2 k={k} t={t}: S_uu.z = {} (should be 0)",
                s_uu_k.z
            );
            assert!(
                s_vv_kn.z.abs() < 1e-8,
                "planar G2 k={k} t={t}: S_vv.z = {}",
                s_vv_kn.z
            );
        }
    }
}

#[test]
fn shared_spoke_g2_residual_tracked_on_consistent_3d() {
    // Record the G2 residual for regression tracking. Full G2
    // along the spoke is the K-P acceptance gate; currently
    // our construction achieves exact G2 mid-edge (validated
    // in kp_row2_tests.rs / kp_row2_corner_tests.rs) and near-
    // endpoint G2 is bounded by the corner 3×3 treatment —
    // measurable but < the "commercial grade" gate.
    let input = consistent_3d_ring(4, 2.0, 0.3);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    let mut max_diff = 0.0_f64;
    for k in 0..n {
        let kn = (k + 1) % n;
        for &t in &[0.1_f64, 0.3, 0.5, 0.7, 0.9] {
            let s_uu_k = patches[k].second_partial_uu(1.0, t);
            let s_vv_kn = patches[kn].second_partial_vv(1.0 - t, 0.0);
            let d = (s_uu_k - s_vv_kn).length();
            if d > max_diff {
                max_diff = d;
            }
        }
    }
    println!("consistent-3D G2 max |S_uu − S_vv| = {max_diff:e}");
    // Regression guard.
    assert!(max_diff < 100.0, "G2 regression: {max_diff:e}");
    assert!(max_diff.is_finite());
}

// =====================================================================
// Central vertex: all N patches meet at the same point
// =====================================================================

#[test]
fn all_patches_meet_at_central_vertex() {
    let input = planar_ring(6, 1.5);
    let opts = DirectConstructionOptions::with_characteristic_length(1.5);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    for (k, p) in patches.iter().enumerate() {
        // Bidegree-9 patch k's (u=1, v=0) corner is the
        // central vertex C.
        let pt = p.point_at(1.0, 0.0);
        let d = pt - input.central_vertex;
        let mag = (d.x.powi(2) + d.y.powi(2) + d.z.powi(2)).sqrt();
        assert!(
            mag < EPS,
            "patch {k}: (1,0) = {pt:?}, C = {:?}",
            input.central_vertex
        );
    }
}

// =====================================================================
// Outer vertices: V_{2k+2} is shared between patches k and k+1
// =====================================================================

#[test]
fn outer_vertices_shared_between_adjacent_sub_patches() {
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    let n = patches.len();
    for k in 0..n {
        let kn = (k + 1) % n;
        // Sub-patch k's (1, 1) and sub-patch (k+1)'s (0, 0) are
        // the same outer vertex V_{2k+2}.
        let a = patches[k].point_at(1.0, 1.0);
        let b = patches[kn].point_at(0.0, 0.0);
        let d = a - b;
        let mag = (d.x.powi(2) + d.y.powi(2) + d.z.powi(2)).sqrt();
        assert!(mag < EPS, "outer vertex shared k={k}");
    }
}

// =====================================================================
// Self-intersection safety
// =====================================================================

#[test]
fn jacobian_well_defined_everywhere_on_planar_input() {
    // After H.5 fits interior CPs to the guide, the bidegree-9
    // patch has a well-defined Jacobian across the FULL domain
    // [0, 1]² — not just the boundary region.
    let input = planar_ring(5, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    for (k, patch) in patches.iter().enumerate() {
        for k_u in 0..10 {
            for k_v in 0..10 {
                let u = 0.05 + 0.1 * k_u as f64;
                let v = 0.05 + 0.1 * k_v as f64;
                let pu = patch.partial_u(u, v);
                let pv = patch.partial_v(u, v);
                let n = pu.cross(&pv);
                let mag = (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
                assert!(
                    mag > 1e-4,
                    "patch {k} at ({u}, {v}): degenerate Jacobian |{mag}|"
                );
            }
        }
    }
}

#[test]
fn bidegree9_patch_approximates_guide_on_planar_input() {
    // With H.5 interior fit, the bidegree-9 patch should closely
    // approximate the guide on planar input (where the guide is
    // exactly a flat plane and the fit can recover it exactly).
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let patches = build_binonic_ring(&input, &opts).unwrap();
    // For planar input, every patch's z-component should stay
    // ~0 across the entire [0, 1]² domain — even in the interior
    // that H.5 fills in.
    for (k, patch) in patches.iter().enumerate() {
        for k_u in 0..=10 {
            for k_v in 0..=10 {
                let u = k_u as f64 / 10.0;
                let v = k_v as f64 / 10.0;
                let pt = patch.point_at(u, v);
                assert!(
                    pt.z.abs() < 1e-6,
                    "patch {k} at ({u}, {v}): z = {} (expected ~0 for planar)",
                    pt.z
                );
            }
        }
    }
}

/// Option C Axis-1 ceiling parity: a bi-nonic guide elevated
/// from the bicubic guide must produce the same bidegree-9
/// output as the bicubic path. At equal degree the row-1/row-2
/// construction collapses to identity on the guide-derived
/// boundaries, so we need to check the full pipeline round-trip.
#[test]
fn binonic_path_is_identical_to_bicubic_when_degree_elevated() {
    let input = planar_ring(3, 1.0);
    let options = DirectConstructionOptions::with_characteristic_length(1.0);

    let bicubic_guides = crate::geometry::gregory::build_guide_ring_with_twist(&input, &options)
        .expect("build_guide_ring_with_twist");
    let bn_guides: Vec<crate::geometry::gregory::KpGuidePatchBiNonic> = bicubic_guides
        .iter()
        .map(crate::geometry::gregory::KpGuidePatchBiNonic::from_bicubic_elevated)
        .collect();

    let bi_ring = crate::geometry::gregory::build_binonic_ring_from_guides(&bicubic_guides, &input)
        .expect("bicubic ring build");
    let bn_ring =
        crate::geometry::gregory::build_binonic_ring_from_binonic_guides(&bn_guides, &input)
            .expect("bi-nonic ring build");

    assert_eq!(bi_ring.len(), bn_ring.len());
    for (k, (a, b)) in bi_ring.iter().zip(bn_ring.iter()).enumerate() {
        for i in 0..10 {
            for j in 0..10 {
                let pa = a.grid[i][j];
                let pb = b.grid[i][j];
                let dx = (pa.x - pb.x).abs();
                let dy = (pa.y - pb.y).abs();
                let dz = (pa.z - pb.z).abs();
                let maxd = dx.max(dy).max(dz);
                assert!(
                    maxd < 1e-9,
                    "patch {k} cp[{i}][{j}]: bicubic={pa:?} vs bi-nonic={pb:?} (Δ={maxd:.3e})"
                );
            }
        }
    }
}

/// Option C MVP parity: building the ring via a bi-quintic guide
/// degree-elevated from the bicubic guide must produce the same
/// bidegree-9 output as the bicubic path. Elevation is pointwise-
/// identical, so the downstream construction should be too.
#[test]
fn biquintic_path_is_identical_to_bicubic_when_degree_elevated() {
    let input = planar_ring(3, 1.0);
    let options = DirectConstructionOptions::with_characteristic_length(1.0);

    let bicubic_guides = crate::geometry::gregory::build_guide_ring_with_twist(&input, &options)
        .expect("build_guide_ring_with_twist");
    let bq_guides: Vec<crate::geometry::gregory::KpGuidePatchBiQuintic> = bicubic_guides
        .iter()
        .map(crate::geometry::gregory::KpGuidePatchBiQuintic::from_bicubic_elevated)
        .collect();

    let bi_ring = crate::geometry::gregory::build_binonic_ring_from_guides(&bicubic_guides, &input)
        .expect("bicubic ring build");
    let bq_ring =
        crate::geometry::gregory::build_binonic_ring_from_biquintic_guides(&bq_guides, &input)
            .expect("bi-quintic ring build");

    assert_eq!(bi_ring.len(), bq_ring.len());
    for (k, (a, b)) in bi_ring.iter().zip(bq_ring.iter()).enumerate() {
        for i in 0..10 {
            for j in 0..10 {
                let pa = a.grid[i][j];
                let pb = b.grid[i][j];
                let dx = (pa.x - pb.x).abs();
                let dy = (pa.y - pb.y).abs();
                let dz = (pa.z - pb.z).abs();
                let maxd = dx.max(dy).max(dz);
                assert!(
                    maxd < 1e-9,
                    "patch {k} cp[{i}][{j}]: bicubic={pa:?} vs bi-quintic={pb:?} (Δ={maxd:.3e})"
                );
            }
        }
    }
}
