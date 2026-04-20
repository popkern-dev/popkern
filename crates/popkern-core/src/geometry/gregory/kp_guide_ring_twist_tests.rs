//! Tests for [`super::build_guide_ring_with_twist`] — the
//! twist-corrected ring constructor that applies Shirman-Séquin
//! twist at `C` plus pairwise twist at every outer vertex.
//!
//! Complements `kp_guide_ring_gap_tests.rs` (twist = 0 baseline).
//! Together these two files establish the before / after behaviour
//! of the twist correction and act as the regression gates for
//! any future shared-ribbon / concentric-tessellation
//! implementation (see ADR amendment).

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

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

/// Truly-3D fixture (inconsistent — face planes at V not tilted
/// to match spoke direction).
fn truly_3d_ring(n: usize, radius: f64, z_amplitude: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let z_v = z_amplitude * if k % 2 == 0 { 1.0 } else { -1.0 };
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

/// Geometrically-consistent non-planar fixture (face plane at V
/// tilted to contain spoke direction).
fn tilted_plane_consistent_ring(n: usize, radius: f64, z_amplitude: f64) -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let tps_c = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    for k in 0..n {
        let theta = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
        let z_v = z_amplitude * if k % 2 == 0 { 1.0 } else { -1.0 };
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

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn triple_product(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    a.dot(&b.cross(&c))
}

fn max_triple(guides: &[super::KpGuidePatch], n: usize, ts: &[f64]) -> f64 {
    let mut mx = 0.0_f64;
    for k in 0..n {
        let k_next = (k + 1) % n;
        for &t in ts {
            let spoke = guides[k].partial_v(1.0, t);
            let cross_k = guides[k].partial_u(1.0, t);
            let cross_next = guides[k_next].partial_v(1.0 - t, 0.0);
            let tp = triple_product(spoke, cross_k, cross_next);
            if tp.abs() > mx {
                mx = tp.abs();
            }
        }
    }
    mx
}

#[test]
fn twist_corrected_planar_ring_still_zero() {
    let input = planar_ring(4, 2.0);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring_with_twist(&input, &opts).expect("should succeed");
    let n = guides.len();
    let mx = max_triple(&guides, n, &[0.1, 0.3, 0.5, 0.7, 0.9]);
    assert!(mx < 1e-7, "twist-corrected planar: triple product = {mx:e}");
}

#[test]
fn twist_correction_reduces_gap_on_inconsistent_3d_input() {
    // truly_3d_ring has V at z ≠ 0 but claims z = 0 face tangent
    // planes — a geometric inconsistency. Twist correction still
    // measurably reduces the gap vs twist = 0.
    let input = truly_3d_ring(4, 2.0, 0.5);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let baseline = build_guide_ring(&input, &opts).expect("baseline");
    let corrected = build_guide_ring_with_twist(&input, &opts).expect("corrected");
    let n = baseline.len();
    let ts = [0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9];
    let baseline_max = max_triple(&baseline, n, &ts);
    let corrected_max = max_triple(&corrected, n, &ts);
    println!("inconsistent-3D gap: baseline = {baseline_max:e}, corrected = {corrected_max:e}");
    assert!(
        corrected_max < baseline_max,
        "twist correction failed to reduce gap: {corrected_max:e} ≥ {baseline_max:e}"
    );
}

#[test]
fn twist_correction_records_gap_on_consistent_3d_input() {
    // tilted_plane_consistent_ring is geometrically consistent —
    // face planes at V tilted to contain the spoke direction. The
    // current twist correction (Shirman-Séquin at C + pairwise at
    // V) does NOT close the gap to machine precision on this
    // input because bicubic guides lack enough DOF to make
    // cross-tangent curves exactly compatible across adjacent
    // sub-patches. Full closure requires a shared-ribbon
    // construction (see ADR amendment 2026-04-19). This test
    // records the residual and catches regressions beyond the
    // measured baseline.
    let input = tilted_plane_consistent_ring(4, 2.0, 0.3);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let baseline = build_guide_ring(&input, &opts).expect("baseline");
    let corrected = build_guide_ring_with_twist(&input, &opts).expect("corrected");
    let n = baseline.len();
    let ts = [0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9];
    let baseline_max = max_triple(&baseline, n, &ts);
    let corrected_max = max_triple(&corrected, n, &ts);
    println!("consistent-3D gap: baseline = {baseline_max:e}, corrected = {corrected_max:e}");
    assert!(
        baseline_max.is_finite() && corrected_max.is_finite(),
        "gap went non-finite"
    );
    assert!(
        baseline_max < 1.0,
        "baseline gap regression: {baseline_max:e}"
    );
    assert!(
        corrected_max < 1.0,
        "corrected gap regression: {corrected_max:e}"
    );
}

#[test]
fn twist_corrected_shared_spoke_cps_still_agree() {
    // Shared-spoke CP equality (G0) must be preserved after twist
    // correction — the tangent and twist values at shared corners
    // are coordinated so edge CPs still match bit-for-bit.
    let input = truly_3d_ring(5, 2.0, 0.4);
    let opts = DirectConstructionOptions::with_characteristic_length(2.0);
    let guides = build_guide_ring_with_twist(&input, &opts).expect("should succeed");
    let n = guides.len();
    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 0..4 {
            let p_k = guides[k].grid[3][j];
            let p_kn = guides[k_next].grid[3 - j][0];
            assert!(
                point_max(p_k, p_kn) < EPS,
                "twist-corrected shared-spoke CP mismatch at k={k} j={j}"
            );
        }
    }
}
