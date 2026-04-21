//! Tests for the tensor-border extractor.

use super::*;
use crate::geometry::gregory::kp15_bi6_cap::assemble_bi6_polar_cap;
use crate::geometry::gregory::kp15_bi6_polar_bridge::convert_polar_to_bi6_cap_input;
use crate::geometry::gregory::{OuterCornerData, SpokeBoundaryData, TangentPlane};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

/// Per-sector tangent plane where v is radial-outward and u is
/// perpendicular (in the xy plane). This gives a geometrically
/// meaningful "into the hole" direction per sector for testing.
fn tp_radial(theta: f64) -> TangentPlane {
    let u = Vec3::new(-theta.sin(), theta.cos(), 0.0);
    let v = Vec3::new(theta.cos(), theta.sin(), 0.0);
    TangentPlane::from_basis(u, v).unwrap()
}

/// Build a regular N-gon polar input with unit radii and radial
/// spokes from origin.
fn regular_polar(n: usize) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v_outer = pt(theta.cos(), theta.sin(), 0.0);
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = pt(t * v_outer.x, t * v_outer.y, 0.0);
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v_outer,
            tangent_plane: tp_xy(),
        });
        // Sector k's face plane: radial-outward orientation at
        // the spoke direction (between spokes k and k+1, but we
        // simplify here).
        face_tangent_planes.push(tp_radial(theta));
        per_spoke_radii.push(0.1);
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

#[test]
fn pentagon_extracts_valid_tensor_border() {
    let polar = regular_polar(5);
    let tb =
        extract_bi3_tensor_border_from_polar_input(&polar).expect("pentagon extraction succeeds");
    tb.validate(5).expect("result is structurally valid");
    assert_eq!(tb.n_sectors(), 5);
    assert_eq!(
        tb.kind(),
        super::super::kp_polar_boundary_input::TensorBorderKind::Bi3
    );
}

#[test]
fn row_0_endpoints_match_spoke_endpoints() {
    // Row 0 of sector k should start at C = spoke[0] and end at
    // V_k = spoke[7].
    let polar = regular_polar(5);
    let tb = extract_bi3_tensor_border_from_polar_input(&polar).unwrap();
    for k in 0..5 {
        let row0 = tb.position_row(k);
        assert_eq!(row0.len(), 4);
        // Start point matches C.
        let c = polar.central_vertex;
        assert!(points_close(row0[0], c));
        // End point matches V_k.
        let v_k = polar.outer_corners[k].position;
        assert!(points_close(row0[3], v_k));
    }
}

#[test]
fn row_1_offset_matches_formula() {
    // For each CP, row1 should be row0 + (1/3) · r_at_cp · v_axis.
    let polar = regular_polar(5);
    let tb = extract_bi3_tensor_border_from_polar_input(&polar).unwrap();
    for k in 0..5 {
        let row0 = tb.position_row(k);
        let row1 = tb.first_cross_row(k);
        let n_hat = polar.face_tangent_planes[k].v;
        // Radii: constant 0.1 in this test.
        let r = 0.1;
        let offset_mag = r / 3.0;
        let expected_offset = Vec3::new(
            offset_mag * n_hat.x,
            offset_mag * n_hat.y,
            offset_mag * n_hat.z,
        );
        for i in 0..4 {
            let actual = Vec3::new(
                row1[i].x - row0[i].x,
                row1[i].y - row0[i].y,
                row1[i].z - row0[i].z,
            );
            assert!(vec_close(actual, expected_offset));
        }
    }
}

#[test]
fn row_2_is_twice_row_1_offset() {
    // row2 − row0 = 2 × (row1 − row0) in the linear-offset scheme.
    let polar = regular_polar(5);
    let tb = extract_bi3_tensor_border_from_polar_input(&polar).unwrap();
    for k in 0..5 {
        let row0 = tb.position_row(k);
        let row1 = tb.first_cross_row(k);
        let row2 = tb.second_cross_row(k);
        for i in 0..4 {
            let d1_x = row1[i].x - row0[i].x;
            let d2_x = row2[i].x - row0[i].x;
            assert!(
                (d2_x - 2.0 * d1_x).abs() < 1e-12,
                "sector {k} cp {i}: row2.x offset = {d2_x}, expected 2·row1 = {}",
                2.0 * d1_x
            );
        }
    }
}

#[test]
fn variable_radius_uses_linear_interp() {
    // radius_at_c = 0.05, radius_at_v = 0.1. CP i=0 should use
    // r = 0.05; CP i=3 should use r = 0.1; CPs i=1, 2 should
    // interpolate linearly.
    let mut polar = regular_polar(5);
    for spoke in polar.spokes.iter_mut() {
        spoke.radius_at_c = 0.05;
        spoke.radius_at_v = 0.1;
    }
    let tb = extract_bi3_tensor_border_from_polar_input(&polar).unwrap();
    for k in 0..5 {
        let row0 = tb.position_row(k);
        let row1 = tb.first_cross_row(k);
        let n_hat = polar.face_tangent_planes[k].v;

        let expected_r_at_cp = [
            0.05,
            (2.0 * 0.05 + 0.1) / 3.0,
            (0.05 + 2.0 * 0.1) / 3.0,
            0.1,
        ];
        for i in 0..4 {
            let expected_offset = n_hat * (expected_r_at_cp[i] / 3.0);
            let actual_offset = Vec3::new(
                row1[i].x - row0[i].x,
                row1[i].y - row0[i].y,
                row1[i].z - row0[i].z,
            );
            assert!(vec_close(actual_offset, expected_offset));
        }
    }
}

#[test]
fn n3_and_n6_extract_successfully() {
    for n in [3usize, 6] {
        let polar = regular_polar(n);
        let tb = extract_bi3_tensor_border_from_polar_input(&polar)
            .unwrap_or_else(|e| panic!("n={n}: {e}"));
        tb.validate(n).unwrap();
    }
}

#[test]
fn end_to_end_extract_convert_assemble() {
    // Full pipeline: polar input → extract tensor-border →
    // convert to Bi6CapInput → assemble. No explicit tensor-
    // border from the caller.
    let polar = regular_polar(5);
    let tb = extract_bi3_tensor_border_from_polar_input(&polar).unwrap();
    let cap_input = convert_polar_to_bi6_cap_input(&polar, &tb).unwrap();
    let out = assemble_bi6_polar_cap(&cap_input).expect("pipeline");
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
fn extraction_rejects_invalid_polar_input() {
    let mut polar = regular_polar(5);
    polar.per_spoke_radii[2] = 0.0; // invalid
    assert!(extract_bi3_tensor_border_from_polar_input(&polar).is_err());
}

// ---- helpers ----

fn points_close(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10 && (a.z - b.z).abs() < 1e-10
}

fn vec_close(a: Vec3, b: Vec3) -> bool {
    (a.x - b.x).abs() < 1e-12 && (a.y - b.y).abs() < 1e-12 && (a.z - b.z).abs() < 1e-12
}
