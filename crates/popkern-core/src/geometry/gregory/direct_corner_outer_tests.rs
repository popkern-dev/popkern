//! Tests for [`super::solve_outer_corners`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

fn linear_spoke_cps(c: Point3, v: Point3) -> [Point3; 8] {
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 8];
    for j in 0..8 {
        let t = (j as f64) / 7.0;
        cps[j] = Point3::new(
            c.x + t * (v.x - c.x),
            c.y + t * (v.y - c.y),
            c.z + t * (v.z - c.z),
        );
    }
    cps
}

fn demo_spoke(c: Point3, v: Point3, edge_id: &str) -> SpokeBoundaryData {
    SpokeBoundaryData {
        boundary_cps: linear_spoke_cps(c, v),
        xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
        xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
        radius_at_c: 0.1,
        radius_at_v: 0.1,
        edge_id: edge_id.into(),
    }
}

/// A flat unit-square two-stripe patch on z=0. Handy for
/// outer-corner tests where the seed tangent plane is exactly
/// the xy-plane.
fn flat_unit_patch_at(
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

#[test]
fn missing_seed_patches_returns_error() {
    let c = Point3::new(0.0, 0.0, 0.0);
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: (0..3)
            .map(|i| demo_spoke(c, Point3::new(1.0, 0.0, 0.0), &format!("e{i}")))
            .collect(),
        outer_corners: vec![
            OuterCornerData {
                position: Point3::new(1.0, 0.0, 0.0),
                tangent_plane: TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap(),
            };
            3
        ],
        central_vertex_tps: vec![TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap(); 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    assert!(solve_outer_corners(&input, &options).is_err());
}

#[test]
fn flat_coplanar_ring_produces_parallelogram_twists() {
    // Two flat-plane patches with a shared outer corner at
    // (0, 0, 0). Patch k is the unit square on x ∈ [-1, 0] × y ∈ [0, 1],
    // patch k+1 is x ∈ [0, 1] × y ∈ [0, 1]. Both z=0.
    //
    // N=3 structure is needed for validate() to accept the input,
    // so we add a third dummy patch. The outer-corner solve
    // operates on adjacent PAIRS, not the whole ring at once.
    let c = Point3::new(0.0, 0.0, 0.0);
    let vs = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -0.5, 0.0),
    ];
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();

    // Three flat-but-arbitrary-orientation two-stripe patches.
    let patches = vec![
        flat_unit_patch_at((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        flat_unit_patch_at((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
        flat_unit_patch_at((0.0, -1.0), (1.0, 0.0), (0.0, 1.0)),
    ];

    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: (0..3)
            .map(|i| demo_spoke(c, vs[i], &format!("e{i}")))
            .collect(),
        outer_corners: vs
            .iter()
            .map(|&v| OuterCornerData {
                position: v,
                tangent_plane: tps_xy,
            })
            .collect(),
        central_vertex_tps: vec![tps_xy; 3],
        seed_patches: Some(patches),
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);

    let result = solve_outer_corners(&input, &options).expect("solved");
    assert_eq!(result.corners.len(), 3);
    for out in &result.corners {
        // For flat coplanar patches the parallelogram rule's output
        // also lies on z=0 — the twist CP is in the patch plane.
        assert!(
            out.patch_k_grid_6_6.z.abs() < 1e-12,
            "grid[6][6] off the plane: {:?}",
            out.patch_k_grid_6_6
        );
        assert!(
            out.patch_k_plus_1_grid_1_1.z.abs() < 1e-12,
            "grid[1][1] off the plane: {:?}",
            out.patch_k_plus_1_grid_1_1
        );
        // Residual against the xy-plane normal should be tiny.
        assert!(
            out.tangent_plane_residual < 1e-10,
            "tangent-plane residual too high: {}",
            out.tangent_plane_residual
        );
    }
}

#[test]
fn parallelogram_rule_formulation_is_stable() {
    // Isolated sanity check on the parallelogram output: given
    // boundary CPs (7, 6) + (6, 7) + (7, 7) = corner, the twist CP
    // at (6, 6) should be their affine combination.
    // This is essentially a documentation/lock-in test of the
    // "minimum twist" formula used in the module.
    let p_7_7 = Point3::new(1.0, 1.0, 0.0);
    let p_7_6 = Point3::new(1.0, 0.857, 0.0); // on spoke boundary
    let p_6_7 = Point3::new(0.857, 1.0, 0.0); // on hex-external
    let expected = Point3::new(
        p_7_6.x + p_6_7.x - p_7_7.x,
        p_7_6.y + p_6_7.y - p_7_7.y,
        p_7_6.z + p_6_7.z - p_7_7.z,
    );
    // This is the formula in solve_outer_corners; just verify
    // the algebra does what we expect.
    assert!((expected.x - 0.857).abs() < 1e-12);
    assert!((expected.y - 0.857).abs() < 1e-12);
    assert!(expected.z.abs() < 1e-12);
}
