//! Tests for `build()` and `can_build()` orchestration.
//!
//! End-to-end G1-exactness on realistic three-quad geometry
//! requires the popkern-fillet extraction helper (next Phase G
//! step) — tests here are structural: verify the pipeline runs,
//! the result struct is fully populated, and the error paths
//! fire as expected.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

fn flat_patch(
    origin: (f64, f64),
    u_dir: (f64, f64),
    v_dir: (f64, f64),
) -> RectangularGregoryTwoStripe {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            let s = (i as f64) / 7.0;
            let t = (j as f64) / 7.0;
            grid[i][j] = Point3::new(
                origin.0 + s * u_dir.0 + t * v_dir.0,
                origin.1 + s * u_dir.1 + t * v_dir.1,
                0.0,
            );
        }
    }
    RectangularGregoryTwoStripe { grid }
}

/// A 3-patch synthetic ring — not a physically realistic
/// three-quad but exercises the orchestration plumbing.
fn synthetic_3_patch_input() -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.5, 0.0);
    let vs = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 0.5, 0.0),
    ];
    let patches = vec![
        flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        flat_patch((-1.0, 1.0), (1.0, 0.0), (0.0, -1.0)),
    ];
    let mut spokes = Vec::with_capacity(3);
    for (k, &v) in vs.iter().enumerate() {
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
            // xtans pointed out-of-xy so the blend-coefficient
            // system is well-conditioned regardless of spoke
            // direction. See direct_row1_ring_tests for the
            // caveat about physical realism.
            xtan_a_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_b_at_c: Vec3::new(0.0, 0.0, 1.0),
            xtan_a_at_v: Vec3::new(0.0, 0.0, 1.0),
            xtan_b_at_v: Vec3::new(0.0, 0.0, 1.0),
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("e{k}"),
        });
    }
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    DirectBoundaryInput {
        central_vertex: c,
        spokes,
        outer_corners: vs
            .iter()
            .map(|&v| OuterCornerData {
                position: v,
                tangent_plane: tps_xy,
            })
            .collect(),
        central_vertex_tps: vec![tps_xy; 3],
        seed_patches: Some(patches),
    }
}

#[test]
fn can_build_returns_report_on_valid_input() {
    let input = synthetic_3_patch_input();
    let options = DirectConstructionOptions::with_characteristic_length(2.0);
    let report = can_build(&input, &options).expect("can_build");
    assert_eq!(report.achievable_continuity, ContinuityLevel::G1);
    assert!(report.near_degeneracy_score >= 0.0);
    assert!(report.near_degeneracy_score <= 1.0);
    assert!(report.estimated_construction_us > 0);
}

#[test]
fn build_returns_valid_result_struct_on_zero_validation() {
    let input = synthetic_3_patch_input();
    let options = DirectConstructionOptions::with_characteristic_length(2.0)
        .with_validation(ValidationLevel::None);
    let result = build(&input, &options).expect("build");
    assert_eq!(result.sub_patches.len(), 3);
    assert_eq!(result.achieved_continuity, ContinuityLevel::G1);
    assert_eq!(result.corner_diagnostics.len(), 3 * 2); // central + outer per patch
    assert_eq!(result.spoke_diagnostics.len(), 3);
    // At None validation, dense-check fields are NaN.
    assert!(result.dense_min_cos_g1.is_nan());
    assert!(result.max_position_gap_g0.is_nan());
}

#[test]
fn build_requires_seed_patches() {
    let mut input = synthetic_3_patch_input();
    input.seed_patches = None;
    let options = DirectConstructionOptions::with_characteristic_length(2.0);
    assert!(build(&input, &options).is_err());
}

#[test]
fn build_with_g2_target_still_delivers_g1_in_mvp() {
    // Phase G MVP is G1-only; can_build reports G1 even when
    // asked for G2. Documenting the current scope.
    let input = synthetic_3_patch_input();
    let options = DirectConstructionOptions::with_characteristic_length(2.0)
        .with_continuity(ContinuityLevel::G2);
    let report = can_build(&input, &options).expect("can_build");
    assert_eq!(report.achievable_continuity, ContinuityLevel::G1);
}

#[test]
fn diagnostics_populated_in_result() {
    let input = synthetic_3_patch_input();
    let options = DirectConstructionOptions::with_characteristic_length(2.0)
        .with_validation(ValidationLevel::None);
    let result = build(&input, &options).expect("build");
    // Per-corner diagnostics: formula family + near-degeneracy
    // score should be set.
    for diag in &result.corner_diagnostics {
        assert!(diag.near_degeneracy_score >= 0.0);
        assert!(diag.near_degeneracy_score <= 1.0);
        assert!(diag.tolerance_used > 0.0);
    }
    // Per-spoke diagnostics: residual + polynomial degrees.
    for diag in &result.spoke_diagnostics {
        assert!(diag.g1_residual.is_finite());
        assert_eq!(diag.alpha_polynomial_degree, 1);
        assert_eq!(diag.beta_polynomial_degree, 1);
    }
}
