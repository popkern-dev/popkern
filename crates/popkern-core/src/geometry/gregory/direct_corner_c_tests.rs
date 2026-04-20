//! Tests for [`super::solve_central_vertex_corner`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_linalg::orthonormal_tangent_basis;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

/// Build a degree-7 Bezier line from `c` to `v` as 8 equally-
/// spaced CPs. Handy for tests where the spoke curve is a
/// straight segment.
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

fn tps_xy() -> TangentPlane {
    TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

/// An N=3 symmetric planar configuration in the xy-plane:
/// central vertex at origin, three outer vertices at 120°
/// intervals on a circle of radius 1. All face tangent planes
/// are the xy-plane (normal +z).
fn symmetric_n3_input() -> DirectBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let angles = [
        0.0_f64,
        2.0 * std::f64::consts::PI / 3.0,
        4.0 * std::f64::consts::PI / 3.0,
    ];
    let vs: Vec<Point3> = angles
        .iter()
        .map(|&a| Point3::new(a.cos(), a.sin(), 0.0))
        .collect();
    DirectBoundaryInput {
        central_vertex: c,
        spokes: vs
            .iter()
            .enumerate()
            .map(|(i, &v)| demo_spoke(c, v, &format!("e{i}")))
            .collect(),
        outer_corners: vs
            .iter()
            .map(|&v| OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            })
            .collect(),
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    }
}

#[test]
fn symmetric_n3_returns_three_cps_on_xy_plane() {
    let input = symmetric_n3_input();
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let result = solve_central_vertex_corner(&input, &options).expect("solved");

    assert_eq!(result.inner_cps.len(), 3);
    // All three inner CPs must lie on the xy plane (z ≈ 0).
    for (k, cp) in result.inner_cps.iter().enumerate() {
        assert!(
            cp.z.abs() < 1e-10,
            "inner CP {k} not on xy plane: z = {}",
            cp.z
        );
    }
    // Common normal is +z (or -z equivalently).
    assert!(
        (result.common_normal.z.abs() - 1.0).abs() < 1e-10,
        "common normal should be ±z, got {:?}",
        result.common_normal
    );
}

#[test]
fn symmetric_n3_inner_cps_have_threefold_symmetry() {
    let input = symmetric_n3_input();
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let result = solve_central_vertex_corner(&input, &options).expect("solved");

    // Rotating each inner CP by 120° should map it to the next.
    let two_pi_over_3 = 2.0 * std::f64::consts::PI / 3.0;
    let (cos120, sin120) = (two_pi_over_3.cos(), two_pi_over_3.sin());
    for k in 0..3 {
        let this = result.inner_cps[k];
        let next = result.inner_cps[(k + 1) % 3];
        let rotated_x = cos120 * this.x - sin120 * this.y;
        let rotated_y = sin120 * this.x + cos120 * this.y;
        assert!(
            (rotated_x - next.x).abs() < 1e-9
                && (rotated_y - next.y).abs() < 1e-9
                && (this.z - next.z).abs() < 1e-10,
            "inner CP {k} rotation -> {next:?} expected, got ({rotated_x}, {rotated_y})"
        );
    }
}

#[test]
fn rejects_zero_normal_sum() {
    // Two N=2 faces with opposite normals — their sum is zero.
    // N=2 would already fail validate() but the rejection path
    // can still be exercised via a synthetic 3-patch input
    // where the normals sum to near-zero.
    let c = Point3::new(0.0, 0.0, 0.0);
    let angles = [
        0.0_f64,
        2.0 * std::f64::consts::PI / 3.0,
        4.0 * std::f64::consts::PI / 3.0,
    ];
    let vs: Vec<Point3> = angles
        .iter()
        .map(|&a| Point3::new(a.cos(), a.sin(), 0.0))
        .collect();
    // Three normals that sum to zero: pick +z, basis-rotated +z
    // by 120° about y, and the residual that cancels them.
    let tp_a = TangentPlane::from_normal(Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let tp_b = TangentPlane::from_normal(Vec3::new(-1.0, 0.0, 0.0)).unwrap();
    // Third plane's normal cancels the residual (0, 0, 0) exactly
    // when summed — use (0, 0, 0)-ish by picking a plane with
    // normal that zeroes the sum. A pair ±x sums to zero already;
    // add a plane with normal anti-parallel to what we have so
    // the remaining sum is zero. Here tp_a + tp_b = 0; we add
    // any plane and the result is that plane's normal only.
    //
    // To actually force a zero sum, use three coplanar-normal
    // planes whose normals sum to zero: {e^{2πi/3}} at k=0,1,2
    // in the xy plane. Each has normal in xy; sum = 0.
    let tp_c = TangentPlane::from_normal(Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let tps = [
        TangentPlane::from_normal(Vec3::new(1.0_f64.cos(), 1.0_f64.sin(), 0.0)).unwrap(),
        TangentPlane::from_normal(Vec3::new(
            (1.0_f64 + 2.0 * std::f64::consts::PI / 3.0).cos(),
            (1.0_f64 + 2.0 * std::f64::consts::PI / 3.0).sin(),
            0.0,
        ))
        .unwrap(),
        TangentPlane::from_normal(Vec3::new(
            (1.0_f64 + 4.0 * std::f64::consts::PI / 3.0).cos(),
            (1.0_f64 + 4.0 * std::f64::consts::PI / 3.0).sin(),
            0.0,
        ))
        .unwrap(),
    ];
    // Sanity: sum ≈ 0.
    let mut s = Vec3::new(0.0, 0.0, 0.0);
    for t in &tps {
        let nrm = t.normal();
        s = Vec3::new(s.x + nrm.x, s.y + nrm.y, s.z + nrm.z);
    }
    assert!(s.length() < 1e-9, "test setup: sum not zero");
    let _ = (tp_a, tp_b, tp_c);
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vs
            .iter()
            .enumerate()
            .map(|(i, &v)| demo_spoke(c, v, &format!("e{i}")))
            .collect(),
        outer_corners: vs
            .iter()
            .map(|&v| OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            })
            .collect(),
        central_vertex_tps: tps.to_vec(),
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let result = solve_central_vertex_corner(&input, &options);
    assert!(matches!(
        result,
        Err(crate::geometry::gregory::direct_errors::DirectConstructionError::CornerTangentPlaneSingular { .. })
    ));
}

#[test]
fn asymmetric_n3_returns_finite_cps() {
    // Three outer vertices at non-symmetric positions.
    let c = Point3::new(0.0, 0.0, 0.0);
    let vs = [
        Point3::new(1.0, 0.1, 0.0),
        Point3::new(-0.3, 0.9, 0.0),
        Point3::new(-0.7, -0.4, 0.0),
    ];
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vs
            .iter()
            .enumerate()
            .map(|(i, &v)| demo_spoke(c, v, &format!("e{i}")))
            .collect(),
        outer_corners: vs
            .iter()
            .map(|&v| OuterCornerData {
                position: v,
                tangent_plane: tps_xy(),
            })
            .collect(),
        central_vertex_tps: vec![tps_xy(); 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let result = solve_central_vertex_corner(&input, &options).expect("solved");
    for cp in &result.inner_cps {
        assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
    }
}

#[test]
fn orthonormal_tangent_basis_is_orthonormal() {
    let n = Vec3::new(0.3, -0.4, 0.5).normalized().unwrap();
    let (t1, t2) = orthonormal_tangent_basis(n);
    assert!((t1.length() - 1.0).abs() < 1e-10);
    assert!((t2.length() - 1.0).abs() < 1e-10);
    assert!(t1.dot(&t2).abs() < 1e-10);
    assert!(t1.dot(&n).abs() < 1e-10);
    assert!(t2.dot(&n).abs() < 1e-10);
}

#[test]
fn central_plane_from_input_averages_normals() {
    let input = symmetric_n3_input();
    let plane = central_plane_from_input(&input).expect("plane");
    // All three input normals are +z, average is +z.
    let n = plane.normal();
    assert!((n.z.abs() - 1.0).abs() < 1e-10);
}

#[test]
fn project_to_tangent_plane_zeros_normal_component() {
    let c = Point3::new(1.0, 2.0, 3.0);
    let n = Vec3::new(0.0, 0.0, 1.0);
    let p = Point3::new(1.5, 2.5, 5.0);
    let proj = project_to_tangent_plane(p, c, n);
    // Projected point has z = 3.0 (c.z), x and y unchanged.
    assert!((proj.x - p.x).abs() < 1e-12);
    assert!((proj.y - p.y).abs() < 1e-12);
    assert!((proj.z - c.z).abs() < 1e-12);
}
