//! Gap-coverage tests for [`super::build_guide_ring`] (baseline,
//! twist = 0).
//!
//! Covers:
//! - Non-planar fixture + G0 on shared spokes, endpoint G1.
//! - Cross-derivative coplanarity along shared spokes (the
//!   triple-product test that QUANTIFIES the G1-along-spoke gap
//!   for the baseline construction).
//! - N = 3 and N = 8 edge cases.
//! - Error path: face normals at C cancel.
//! - H.4.d data-flow preview: cross_deriv_cps + degree-9 elevation.
//!
//! Twist-corrected equivalents live in `kp_guide_ring_twist_tests.rs`.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::elevate_cubic_vec_to_nine;

const EPS: f64 = 1e-9;

/// Planar N-sided fixture centred at origin on z = 0.
fn planar_ring(n: usize, radius: f64) -> DirectBoundaryInput {
    ring_with_spoke_bump(n, radius, 0.0)
}

/// Non-planar N-sided fixture with sin(πt) z-bump in the middle
/// of each spoke (zero at the endpoints — keeps C and V on z = 0).
fn non_planar_ring(n: usize, radius: f64, bump: f64) -> DirectBoundaryInput {
    ring_with_spoke_bump(n, radius, bump)
}

fn ring_with_spoke_bump(n: usize, radius: f64, bump: f64) -> DirectBoundaryInput {
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
            let z = bump * (std::f64::consts::PI * t).sin();
            cps[j] = Point3::new(c.x + t * (v.x - c.x), c.y + t * (v.y - c.y), z);
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

/// Truly-3D fixture — outer hex vertices have non-zero z
/// alternating in sign. Face tangent planes at V claimed to be
/// z = 0 (geometrically inconsistent — real input would tilt
/// them to match the spoke direction at V).
fn truly_3d_ring(n: usize, radius: f64, z_amp: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
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

fn triple_product(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.dot(&b.cross(&c))
}

// =====================================================================
// Non-planar G0 and endpoint G1
// =====================================================================

#[test]
fn non_planar_ring_shared_spoke_cps_agree() {
    let input = non_planar_ring(5, 2.0, 0.3);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    let n = guides.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 0..4 {
            assert!(
                point_max(guides[k].grid[3][j], guides[k_next].grid[3 - j][0]) < EPS,
                "non-planar shared-spoke CP at k={k} j={j}"
            );
        }
    }
}

#[test]
fn non_planar_ring_shared_spoke_curves_evaluate_identically() {
    let input = non_planar_ring(4, 1.5, 0.25);
    let opts = DirectConstructionOptions::with_characteristic_length(1.5);
    let guides = build_guide_ring(&input, &opts).expect("build");
    let n = guides.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in &[0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0] {
            let p_k = guides[k].point_at(1.0, t);
            let p_kn = guides[k_next].point_at(1.0 - t, 0.0);
            assert!(point_max(p_k, p_kn) < EPS, "at t={t} k={k}");
        }
    }
}

#[test]
fn non_planar_ring_tangent_planes_at_c_coincide() {
    let input = non_planar_ring(6, 2.0, 0.4);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    for (k, g) in guides.iter().enumerate() {
        let du = g.partial_u(1.0, 0.0);
        let dv = g.partial_v(1.0, 0.0);
        assert!(du.z.abs() < EPS, "guide {k} c1 ∂u.z");
        assert!(dv.z.abs() < EPS, "guide {k} c1 ∂v.z");
    }
}

#[test]
fn non_planar_ring_tangent_planes_at_v_coincide() {
    let input = non_planar_ring(5, 2.0, 0.35);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    let n = guides.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        for (name, v) in [
            ("g_k c2 ∂u", guides[k].partial_u(1.0, 1.0)),
            ("g_k c2 ∂v", guides[k].partial_v(1.0, 1.0)),
            ("g_{k+1} c0 ∂u", guides[k_next].partial_u(0.0, 0.0)),
            ("g_{k+1} c0 ∂v", guides[k_next].partial_v(0.0, 0.0)),
        ] {
            assert!(v.z.abs() < EPS, "{name} at k={k}: z={}", v.z);
        }
    }
}

// =====================================================================
// Cross-derivative coplanarity — baseline (twist = 0)
// =====================================================================

#[test]
fn cross_derivative_coplanarity_planar_ring() {
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    let n = guides.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let spoke = guides[k].partial_v(1.0, t);
            let cross_k = guides[k].partial_u(1.0, t);
            let cross_next = guides[k_next].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_next);
            assert!(tp.abs() < 1e-8, "planar k={k} t={t}: tp={tp:e}");
        }
    }
}

#[test]
fn cross_derivative_coplanarity_baseline_truly_3d_exposes_gap() {
    // Baseline (twist = 0) has significant triple product on
    // truly-3D input — measures ≈ 1.2 at commit time. The
    // twist-corrected variant in `kp_guide_ring_twist_tests.rs`
    // measurably reduces this.
    let input = truly_3d_ring(4, 2.0, 0.5);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    let n = guides.len();
    let mut max_abs_tp = 0.0_f64;
    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let spoke = guides[k].partial_v(1.0, t);
            let cross_k = guides[k].partial_u(1.0, t);
            let cross_next = guides[k_next].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_next);
            if tp.abs() > max_abs_tp {
                max_abs_tp = tp.abs();
            }
        }
    }
    println!("baseline truly-3D max |triple product| = {max_abs_tp:e}");
    // Endpoints must still be G1 (projection guarantees it).
    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in &[0.0, 1.0] {
            let spoke = guides[k].partial_v(1.0, t);
            let cross_k = guides[k].partial_u(1.0, t);
            let cross_next = guides[k_next].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_next);
            assert!(tp.abs() < 1e-8, "endpoint k={k} t={t}: tp={tp:e}");
        }
    }
}

// =====================================================================
// N = 3 and N = 8 edge cases
// =====================================================================

#[test]
fn n_equals_3_full_invariants() {
    let input = planar_ring(3, 1.0);
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let guides = build_guide_ring(&input, &opts).expect("N=3 build");
    assert_eq!(guides.len(), 3);
    for k in 0..3 {
        let k_next = (k + 1) % 3;
        for j in 0..4 {
            assert!(
                point_max(guides[k].grid[3][j], guides[k_next].grid[3 - j][0]) < EPS,
                "N=3 shared-spoke CP at k={k} j={j}"
            );
        }
    }
}

#[test]
fn n_equals_8_full_invariants() {
    let input = planar_ring(8, 3.0);
    let opts = DirectConstructionOptions::with_characteristic_length(3.0);
    let guides = build_guide_ring(&input, &opts).expect("N=8 build");
    assert_eq!(guides.len(), 8);
    for k in 0..8 {
        let k_next = (k + 1) % 8;
        for j in 0..4 {
            assert!(
                point_max(guides[k].grid[3][j], guides[k_next].grid[3 - j][0]) < EPS,
                "N=8 shared-spoke CP at k={k} j={j}"
            );
        }
    }
}

// =====================================================================
// Error paths
// =====================================================================

#[test]
fn cancelling_face_normals_at_c_returns_singular_error() {
    let mut input = planar_ring(4, 1.0);
    let tps_up = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let tps_dn = TangentPlane::from_normal(Vec3::new(0.0, 0.0, -1.0)).unwrap();
    input.central_vertex_tps = vec![tps_up, tps_dn, tps_up, tps_dn];
    let opts = DirectConstructionOptions::with_characteristic_length(1.0);
    let err = build_guide_ring(&input, &opts);
    assert!(
        matches!(
            err,
            Err(DirectConstructionError::CornerTangentPlaneSingular { .. })
        ),
        "expected CornerTangentPlaneSingular, got {err:?}"
    );
}

// =====================================================================
// H.4.d data-flow preview
// =====================================================================

#[test]
fn h4d_integration_preview_cross_deriv_u1_elevated_to_nine() {
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring(&input, &opts).expect("build");
    for g in &guides {
        let cubic_cps = g.cross_deriv_cps_u1();
        let elevated = elevate_cubic_vec_to_nine(cubic_cps);
        assert_eq!(elevated.len(), 10);
        assert!((elevated[0] - cubic_cps[0]).length() < EPS);
        assert!((elevated[9] - cubic_cps[3]).length() < EPS);
        // Pointwise polynomial identity.
        for k in 0..=10 {
            let t = (k as f64) / 10.0;
            let cubic_eval = eval_bezier_curve(&cubic_cps, t);
            let elev_eval = eval_bezier_curve(&elevated, t);
            assert!((cubic_eval - elev_eval).length() < 1e-10, "t={t}");
        }
    }
}

fn eval_bezier_curve(cps: &[Vec3], t: f64) -> Vec3 {
    let n = cps.len() - 1;
    let omt = 1.0 - t;
    let mut binom = vec![0.0_f64; n + 1];
    binom[0] = 1.0;
    for k in 1..=n {
        binom[k] = binom[k - 1] * ((n - k + 1) as f64) / (k as f64);
    }
    let mut out = Vec3::zero();
    for i in 0..=n {
        let w = binom[i] * t.powi(i as i32) * omt.powi((n - i) as i32);
        out.x += cps[i].x * w;
        out.y += cps[i].y * w;
        out.z += cps[i].z * w;
    }
    out
}
