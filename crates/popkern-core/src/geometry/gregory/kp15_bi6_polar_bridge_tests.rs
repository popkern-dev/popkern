//! Tests for the polar → bi-6 cap converter.
//!
//! Strategy:
//! 1. Build a deterministic pentagon `PolarBoundaryInput` +
//!    matching synthetic bi-3 `TensorBorder`.
//! 2. Convert via [`convert_polar_to_bi6_cap_input`].
//! 3. Verify the result passes `Bi6CapInput` structural checks
//!    and round-trips correctly through `assemble_bi6_polar_cap`.
//! 4. End-to-end: run `f5_solve_optimal_bi6_cap` on the
//!    converted input and verify finite, optimal output.

use super::*;
use crate::geometry::Vec3;
use crate::geometry::gregory::kp_polar_boundary_input::{
    PolarBoundaryInput, TensorBorder, TensorBorderBi3,
};
use crate::geometry::gregory::kp15_bi6_cap::assemble_bi6_polar_cap;
use crate::geometry::gregory::kp15_bi6_cap_f5_solve::f5_solve_optimal_bi6_cap;
use crate::geometry::gregory::{OuterCornerData, SpokeBoundaryData, TangentPlane};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

/// Build a pentagon `PolarBoundaryInput` with degree-7 straight-
/// line spokes from origin outward.
fn pentagon_polar_input(n: usize) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v = pt(theta.cos(), theta.sin(), 0.0);
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = pt(t * v.x, t * v.y, t * v.z);
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: 0.05,
            radius_at_v: 0.05,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_xy());
        per_spoke_radii.push(0.05);
    }
    PolarBoundaryInput {
        central_vertex: center,
        central_tangent_plane: tp_xy(),
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    }
}

/// Build a synthetic bi-3 `TensorBorder` with deterministic CPs.
/// Values are chosen to be finite and well-conditioned (no
/// near-degenerate configurations).
fn synthetic_bi3_tensor_border(n: usize) -> TensorBorder {
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        // Row 0: position curve (the outer arc) — 4 CPs along
        // the "outer face" of sector k.
        let row0: [Point3; 4] = [
            pt(dx, dy, 0.1),
            pt(dx * 1.1, dy * 1.1, 0.1),
            pt(dx * 1.2, dy * 1.2, 0.1),
            pt(dx * 1.3, dy * 1.3, 0.1),
        ];
        // Row 1: depth-1 cross-derivative — offset in z.
        let row1: [Point3; 4] = [
            pt(dx, dy, 0.2),
            pt(dx * 1.1, dy * 1.1, 0.2),
            pt(dx * 1.2, dy * 1.2, 0.2),
            pt(dx * 1.3, dy * 1.3, 0.2),
        ];
        // Row 2: depth-2 cross-derivative — further z offset.
        let row2: [Point3; 4] = [
            pt(dx, dy, 0.3),
            pt(dx * 1.1, dy * 1.1, 0.3),
            pt(dx * 1.2, dy * 1.2, 0.3),
            pt(dx * 1.3, dy * 1.3, 0.3),
        ];
        sectors.push([row0, row1, row2]);
    }
    TensorBorder::Bi3(TensorBorderBi3 { sectors })
}

// ---- Converter correctness ----

#[test]
fn pentagon_converts_to_bi6_cap_input() {
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).expect("pentagon converts cleanly");
    assert_eq!(cap_input.n, 5);
    assert_eq!(cap_input.sectors.len(), 5);
    assert_eq!(cap_input.c5s.len(), 5);
    assert_eq!(cap_input.c6s.len(), 5);
    cap_input.params.validate().expect("params valid");
}

#[test]
fn c7_equals_central_vertex() {
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    assert!((cap.c7.x - polar.central_vertex.x).abs() < 1e-12);
    assert!((cap.c7.y - polar.central_vertex.y).abs() < 1e-12);
    assert!((cap.c7.z - polar.central_vertex.z).abs() < 1e-12);
}

#[test]
fn rejects_bi5_tensor_border() {
    let polar = pentagon_polar_input(5);
    let mut sectors = Vec::with_capacity(5);
    for _ in 0..5 {
        let row: [Point3; 6] = core::array::from_fn(|_| pt(0.0, 0.0, 0.0));
        sectors.push([row, row, row]);
    }
    let tb = TensorBorder::Bi5(super::super::kp_polar_boundary_input::TensorBorderBi5 { sectors });
    let err = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap_err();
    assert!(err.contains("bi-3"));
}

#[test]
fn n_mismatch_is_detected() {
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(6); // wrong length
    let err = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap_err();
    assert!(err.contains("sector"));
}

#[test]
fn hexagon_pins_p0_21() {
    let polar = pentagon_polar_input(6);
    let tb = synthetic_bi3_tensor_border(6);
    let cap = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    assert!(cap.p0_21_pinned.is_some());
}

#[test]
fn pentagon_does_not_pin_p0_21() {
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    assert!(cap.p0_21_pinned.is_none());
}

#[test]
fn row0_bi6_elevation_preserves_endpoints() {
    // Degree elevation of a bi-3 curve preserves the 1st and
    // last CP. So after conversion, sector k's p_00 should equal
    // tensor_border row-0[0] and p_60 should equal row-0[3].
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    let bi3 = match &tb {
        TensorBorder::Bi3(b) => b,
        _ => unreachable!(),
    };
    for k in 0..5 {
        let row0 = bi3.sectors[k][0];
        assert!(points_close(cap.sectors[k].p_00, row0[0]));
        assert!(points_close(cap.sectors[k].p_60, row0[3]));
    }
}

// ---- End-to-end ----

#[test]
fn pentagon_pipeline_end_to_end_assembles() {
    // Full flow: polar input → bi-6 cap input → assembly.
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    let out = assemble_bi6_polar_cap(&cap_input).expect("pipeline assembles");
    assert_eq!(out.patches_bi6.len(), 5);
    assert_eq!(out.patches_bi9.len(), 5);
    for patch in &out.patches_bi9 {
        for row in patch {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn pentagon_pipeline_end_to_end_f5_optimises() {
    // Full flow including F_5 minimisation.
    let polar = pentagon_polar_input(5);
    let tb = synthetic_bi3_tensor_border(5);
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    let (_optimal, out) = f5_solve_optimal_bi6_cap(&cap_input).expect("F_5 converges");
    assert_eq!(out.patches_bi6.len(), 5);
    for patch in &out.patches_bi9 {
        for row in patch {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn hexagon_pipeline_end_to_end() {
    let polar = pentagon_polar_input(6);
    let tb = synthetic_bi3_tensor_border(6);
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    let (_optimal, out) = f5_solve_optimal_bi6_cap(&cap_input).expect("hexagon via F_5");
    assert_eq!(out.patches_bi6.len(), 6);
}

#[test]
fn higher_valences_pipeline() {
    for n in [7usize, 8, 10, 12] {
        let polar = pentagon_polar_input(n);
        let tb = synthetic_bi3_tensor_border(n);
        let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
        let (_opt, out) =
            f5_solve_optimal_bi6_cap(&cap_input).unwrap_or_else(|e| panic!("n={n} failed: {e}"));
        assert_eq!(out.patches_bi6.len(), n);
    }
}

// ---- helpers ----

fn points_close(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10 && (a.z - b.z).abs() < 1e-10
}
