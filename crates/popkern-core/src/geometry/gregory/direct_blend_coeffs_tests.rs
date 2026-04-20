//! Tests for [`super::compute_ring_blend_coeffs`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

fn demo_spoke_c_to_v(
    c: Point3,
    v: Point3,
    xtan_a_c: Vec3,
    xtan_b_c: Vec3,
    xtan_a_v: Vec3,
    xtan_b_v: Vec3,
) -> SpokeBoundaryData {
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 8];
    for j in 0..8 {
        let t = (j as f64) / 7.0;
        cps[j] = Point3::new(
            c.x + t * (v.x - c.x),
            c.y + t * (v.y - c.y),
            c.z + t * (v.z - c.z),
        );
    }
    SpokeBoundaryData {
        boundary_cps: cps,
        xtan_a_at_c: xtan_a_c,
        xtan_b_at_c: xtan_b_c,
        xtan_a_at_v: xtan_a_v,
        xtan_b_at_v: xtan_b_v,
        radius_at_c: 0.1,
        radius_at_v: 0.1,
        edge_id: "test".into(),
    }
}

fn tps_xy() -> TangentPlane {
    TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

#[test]
fn symmetric_flat_input_yields_alpha_minus_one_beta_zero() {
    // When xtan_A = xtan_B (= same geometric direction) and T
    // perpendicular to both, α = -1, β = 0 is the solution to
    // xtan_A + α·xtan_B + β·T = 0.
    let c = Point3::new(0.0, 0.0, 0.0);
    let v = Point3::new(1.0, 0.0, 0.0);
    // xtan_A = xtan_B = (0, 1, 0); T = (1, 0, 0).
    let spoke = demo_spoke_c_to_v(
        c,
        v,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    );
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![spoke.clone(), spoke.clone(), spoke],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            };
            3
        ],
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let coeffs = compute_ring_blend_coeffs(&input, &options).expect("solved");
    assert_eq!(coeffs.len(), 3);
    for c in &coeffs {
        assert!(
            (c.alpha_0 - (-1.0)).abs() < 1e-10,
            "alpha_0 = {} (want -1)",
            c.alpha_0
        );
        assert!(
            (c.alpha_1 - (-1.0)).abs() < 1e-10,
            "alpha_1 = {} (want -1)",
            c.alpha_1
        );
        assert!(c.beta_0.abs() < 1e-10, "beta_0 = {} (want 0)", c.beta_0);
        assert!(c.beta_1.abs() < 1e-10, "beta_1 = {} (want 0)", c.beta_1);
    }
}

#[test]
fn asymmetric_xtans_yield_finite_coeffs() {
    // xtan_A and xtan_B at angles to each other; as long as they
    // and T are not colinear, the solve produces finite α, β.
    let c = Point3::new(0.0, 0.0, 0.0);
    let v = Point3::new(1.0, 0.0, 0.0);
    let spoke = demo_spoke_c_to_v(
        c,
        v,
        Vec3::new(0.2, 1.0, 0.0).normalized().unwrap(),
        Vec3::new(-0.2, 1.0, 0.0).normalized().unwrap(),
        Vec3::new(0.1, 1.0, 0.0).normalized().unwrap(),
        Vec3::new(-0.1, 1.0, 0.0).normalized().unwrap(),
    );
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![spoke.clone(), spoke.clone(), spoke],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            };
            3
        ],
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let coeffs = compute_ring_blend_coeffs(&input, &options).expect("solved");
    for c in &coeffs {
        assert!(c.alpha_0.is_finite());
        assert!(c.alpha_1.is_finite());
        assert!(c.beta_0.is_finite());
        assert!(c.beta_1.is_finite());
        assert!(c.condition_number.is_finite());
    }
}

#[test]
fn degenerate_xtan_b_colinear_with_t_rejects() {
    // xtan_B parallel to T makes the 2×2 system singular —
    // should return CornerTangentPlaneSingular.
    let c = Point3::new(0.0, 0.0, 0.0);
    let v = Point3::new(1.0, 0.0, 0.0);
    let spoke = demo_spoke_c_to_v(
        c,
        v,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0), // parallel to T
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    );
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![spoke.clone(), spoke.clone(), spoke],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            };
            3
        ],
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    assert!(matches!(
        compute_ring_blend_coeffs(&input, &options),
        Err(crate::geometry::gregory::direct_errors::DirectConstructionError::CornerTangentPlaneSingular { .. })
    ));
}

#[test]
fn condition_number_reflects_geometry() {
    // Well-conditioned (perpendicular xtan_B and T) should yield
    // a larger condition number than degenerate (small angle).
    let c = Point3::new(0.0, 0.0, 0.0);
    let v = Point3::new(1.0, 0.0, 0.0);
    let good = demo_spoke_c_to_v(
        c,
        v,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    );
    // Near-degenerate: xtan_B almost parallel to T.
    let near_deg = demo_spoke_c_to_v(
        c,
        v,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.99, 0.14, 0.0).normalized().unwrap(),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.99, 0.14, 0.0).normalized().unwrap(),
    );
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let input_good = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![good.clone(), good.clone(), good],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            };
            3
        ],
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let input_bad = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![near_deg.clone(), near_deg.clone(), near_deg],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            };
            3
        ],
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let cg = compute_ring_blend_coeffs(&input_good, &options).expect("good solved");
    let cb = compute_ring_blend_coeffs(&input_bad, &options).expect("bad solved");
    // Good configuration's conditioning number should be at
    // least an order of magnitude better than the near-degenerate.
    assert!(
        cg[0].condition_number > cb[0].condition_number * 2.0,
        "good={} bad={}",
        cg[0].condition_number,
        cb[0].condition_number
    );
}
