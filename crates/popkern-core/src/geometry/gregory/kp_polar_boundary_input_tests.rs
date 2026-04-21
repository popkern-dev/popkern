//! Tests for [`PolarBoundaryInput`] and [`TensorBorder`].
//!
//! Covers:
//! - `validate` accepts well-formed N=5 / N=6 / N=8 inputs.
//! - `validate` rejects mismatched vector lengths and non-finite
//!   radii.
//! - `is_constant_radius` distinguishes uniform from non-uniform.
//! - `TensorBorder` encoding round-trips bi-3 / bi-5 correctly
//!   and accessor methods return the right rows.

use super::*;
use crate::geometry::Vec3;
use crate::geometry::gregory::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn make_spoke(c: Point3, v: Point3, edge_id: &str) -> SpokeBoundaryData {
    // Straight-line degree-7 Bezier from C to V.
    let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
    for (i, cp) in boundary_cps.iter_mut().enumerate() {
        let t = i as f64 / 7.0;
        *cp = Point3::new(
            c.x + t * (v.x - c.x),
            c.y + t * (v.y - c.y),
            c.z + t * (v.z - c.z),
        );
    }
    SpokeBoundaryData {
        boundary_cps,
        xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
        xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
        radius_at_c: 0.05,
        radius_at_v: 0.05,
        edge_id: edge_id.to_string(),
    }
}

/// Build a regular N-sided polar input at the origin, radius 1.0
/// for the outer-corner circle, uniform fillet radius 0.05.
fn make_regular_polar(n: usize) -> PolarBoundaryInput {
    let c = Point3::new(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v = Point3::new(theta.cos(), theta.sin(), 0.0);
        spokes.push(make_spoke(c, v, &format!("edge_{k}")));
        outer_corners.push(OuterCornerData {
            position: v,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_xy());
        per_spoke_radii.push(0.05);
    }
    PolarBoundaryInput {
        central_vertex: c,
        central_tangent_plane: tp_xy(),
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    }
}

#[test]
fn regular_n5_validates() {
    let input = make_regular_polar(5);
    assert_eq!(input.n(), 5);
    input.validate().expect("regular N=5 should validate");
    assert!(input.is_constant_radius());
}

#[test]
fn regular_n6_validates() {
    let input = make_regular_polar(6);
    assert_eq!(input.n(), 6);
    input.validate().expect("regular N=6 should validate");
}

#[test]
fn regular_n8_validates() {
    let input = make_regular_polar(8);
    assert_eq!(input.n(), 8);
    input.validate().expect("regular N=8 should validate");
}

#[test]
fn rejects_below_n3() {
    let input = make_regular_polar(3);
    // N=3 is at the lower bound but accepted for the K-P 2013 §6
    // n=3 special case.
    input.validate().expect("N=3 should validate");

    // Mutate to N=2: rebuild a shorter input.
    let mut short = make_regular_polar(3);
    short.spokes.pop();
    short.outer_corners.pop();
    short.face_tangent_planes.pop();
    short.per_spoke_radii.pop();
    assert!(short.validate().is_err());
}

#[test]
fn rejects_length_mismatch() {
    let mut input = make_regular_polar(5);
    input.outer_corners.pop();
    let err = input.validate().expect_err("length mismatch must fail");
    assert!(err.contains("outer_corners"));

    let mut input = make_regular_polar(5);
    input.face_tangent_planes.pop();
    let err = input.validate().expect_err("face TP mismatch must fail");
    assert!(err.contains("face_tangent_planes"));

    let mut input = make_regular_polar(5);
    input.per_spoke_radii.pop();
    let err = input.validate().expect_err("radii mismatch must fail");
    assert!(err.contains("per_spoke_radii"));
}

#[test]
fn rejects_non_positive_radius() {
    let mut input = make_regular_polar(5);
    input.per_spoke_radii[2] = 0.0;
    assert!(input.validate().is_err());

    let mut input = make_regular_polar(5);
    input.per_spoke_radii[0] = -0.05;
    assert!(input.validate().is_err());

    let mut input = make_regular_polar(5);
    input.per_spoke_radii[3] = f64::NAN;
    assert!(input.validate().is_err());
}

#[test]
fn is_constant_radius_variable_case() {
    let mut input = make_regular_polar(5);
    input.per_spoke_radii = vec![0.05, 0.06, 0.05, 0.07, 0.05];
    assert!(!input.is_constant_radius());
    input.validate().expect("variable radii should validate");
}

#[test]
fn tensor_border_validate_matches_expected_n() {
    let mut sectors = Vec::with_capacity(5);
    for _ in 0..5 {
        let row: [Point3; 4] = core::array::from_fn(|_| Point3::new(0.0, 0.0, 0.0));
        sectors.push([row, row, row]);
    }
    let tb = TensorBorder::Bi3(TensorBorderBi3 { sectors });
    tb.validate(5)
        .expect("5-sector tensor-border should validate against n=5");
    assert!(tb.validate(6).is_err());
    assert!(tb.validate(4).is_err());
}

#[test]
fn tensor_border_validate_rejects_non_finite_cps() {
    let mut sectors = Vec::with_capacity(5);
    for _ in 0..5 {
        let row: [Point3; 4] = core::array::from_fn(|_| Point3::new(0.0, 0.0, 0.0));
        sectors.push([row, row, row]);
    }
    sectors[2][1][2] = Point3::new(f64::NAN, 0.0, 0.0);
    let tb = TensorBorder::Bi3(TensorBorderBi3 { sectors });
    let err = tb.validate(5).expect_err("NaN CP should fail validation");
    assert!(err.contains("sector 2"));
    assert!(err.contains("row 1"));
}

#[test]
fn from_direct_boundary_input_copies_fields() {
    let polar = make_regular_polar(5);
    // Build a DirectBoundaryInput that mirrors the polar one
    // (same spokes + outer corners + per-sector tangent planes).
    let direct = DirectBoundaryInput {
        central_vertex: polar.central_vertex,
        spokes: polar.spokes.clone(),
        outer_corners: polar.outer_corners.clone(),
        central_vertex_tps: polar.face_tangent_planes.clone(),
        seed_patches: None,
    };
    let polar_from_direct = PolarBoundaryInput::from_direct_boundary_input(&direct)
        .expect("converter succeeds for well-formed input");
    assert_eq!(polar_from_direct.n(), 5);
    assert_eq!(polar_from_direct.spokes.len(), 5);
    assert_eq!(polar_from_direct.outer_corners.len(), 5);
    assert_eq!(polar_from_direct.face_tangent_planes.len(), 5);
    assert_eq!(polar_from_direct.per_spoke_radii.len(), 5);
    // per_spoke_radii[k] = midpoint of spoke's c/v radii.
    for (k, r) in polar_from_direct.per_spoke_radii.iter().enumerate() {
        let expected = 0.5 * (polar.spokes[k].radius_at_c + polar.spokes[k].radius_at_v);
        assert!((r - expected).abs() < 1e-15);
    }
    polar_from_direct
        .validate()
        .expect("converted input validates");
}

#[test]
fn from_direct_rejects_mismatched_central_vertex_tps() {
    use crate::geometry::gregory::DirectBoundaryInput;
    let polar = make_regular_polar(5);
    let mut tps_short = polar.face_tangent_planes.clone();
    tps_short.pop();
    let direct = DirectBoundaryInput {
        central_vertex: polar.central_vertex,
        spokes: polar.spokes.clone(),
        outer_corners: polar.outer_corners.clone(),
        central_vertex_tps: tps_short,
        seed_patches: None,
    };
    let err = PolarBoundaryInput::from_direct_boundary_input(&direct).unwrap_err();
    assert!(err.contains("central_vertex_tps"));
}

#[test]
fn tensor_border_bi3_accessors() {
    let mut sectors = Vec::with_capacity(5);
    for k in 0..5 {
        let row0: [Point3; 4] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 0.0));
        let row1: [Point3; 4] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 1.0));
        let row2: [Point3; 4] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 2.0));
        sectors.push([row0, row1, row2]);
    }
    let tb = TensorBorder::Bi3(TensorBorderBi3 { sectors });
    assert_eq!(tb.n_sectors(), 5);
    assert_eq!(tb.kind(), TensorBorderKind::Bi3);
    assert_eq!(tb.position_row(0).len(), 4);
    assert_eq!(tb.first_cross_row(0).len(), 4);
    assert_eq!(tb.second_cross_row(0).len(), 4);
    // Spot-check: sector 3 row 0 CP 2 == Point3(3, 2, 0).
    let p = tb.position_row(3)[2];
    assert!((p.x - 3.0).abs() < 1e-12 && (p.y - 2.0).abs() < 1e-12 && p.z.abs() < 1e-12);
    // Sector 2 row 1 CP 1 should have z = 1.0.
    assert!((tb.first_cross_row(2)[1].z - 1.0).abs() < 1e-12);
    assert!((tb.second_cross_row(2)[1].z - 2.0).abs() < 1e-12);
}

#[test]
fn tensor_border_bi5_accessors() {
    let mut sectors = Vec::with_capacity(6);
    for k in 0..6 {
        let row0: [Point3; 6] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 0.0));
        let row1: [Point3; 6] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 1.0));
        let row2: [Point3; 6] = core::array::from_fn(|i| Point3::new(k as f64, i as f64, 2.0));
        sectors.push([row0, row1, row2]);
    }
    let tb = TensorBorder::Bi5(TensorBorderBi5 { sectors });
    assert_eq!(tb.n_sectors(), 6);
    assert_eq!(tb.kind(), TensorBorderKind::Bi5);
    assert_eq!(tb.position_row(0).len(), 6);
    assert_eq!(tb.first_cross_row(0).len(), 6);
    assert_eq!(tb.second_cross_row(0).len(), 6);
    let p = tb.position_row(4)[5];
    assert!((p.x - 4.0).abs() < 1e-12 && (p.y - 5.0).abs() < 1e-12);
}
