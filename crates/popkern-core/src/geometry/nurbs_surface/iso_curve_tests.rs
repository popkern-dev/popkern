//! Exact iso-parameter curve extraction tests (see `iso_curve.rs`).

use super::*;
use crate::geometry::{CurveEval, Point3, SurfaceEval};

/// Flat unit bilinear patch (duplicated from `tests.rs` — the test
/// modules are siblings and share no helper crate).
fn bilinear_patch() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    }
}

/// A rational quarter-cylinder strip (degree-2 arc section in v,
/// linear extrusion in u): radius 1, axis +X, the section arc weights
/// [1, cos45°, 1]. Exercises the weighted de Boor combination the
/// conic/circular fillet strips use.
fn rational_quarter_cylinder() -> NurbsSurface {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    NurbsSurface {
        degree_u: 1,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 1.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
        ],
        num_u: 2,
        num_v: 3,
        weights: vec![1.0, w, 1.0, 1.0, w, 1.0],
    }
}

#[test]
fn iso_curve_at_u_matches_surface_to_roundoff() {
    let surface = rational_quarter_cylinder();
    for &u in &[0.0, 0.17, 0.5, 0.83, 1.0] {
        let iso = surface.iso_curve_at_u(u);
        assert_eq!(iso.degree, surface.degree_v);
        assert_eq!(iso.knots, surface.knots_v);
        for i in 0..=32 {
            let v = i as f64 / 32.0;
            let d = iso.point_at(v).distance_to(&surface.point_at(u, v));
            assert!(d < 1e-13, "iso_u({u}) deviates {d:.3e} at v={v}");
        }
    }
}

#[test]
fn iso_curve_at_v_matches_surface_to_roundoff() {
    let surface = rational_quarter_cylinder();
    for &v in &[0.0, 0.25, 0.61, 1.0] {
        let iso = surface.iso_curve_at_v(v);
        assert_eq!(iso.degree, surface.degree_u);
        assert_eq!(iso.knots, surface.knots_u);
        for i in 0..=32 {
            let u = i as f64 / 32.0;
            let d = iso.point_at(u).distance_to(&surface.point_at(u, v));
            assert!(d < 1e-13, "iso_v({v}) deviates {d:.3e} at u={u}");
        }
    }
}

#[test]
fn iso_curve_interior_station_is_the_exact_circular_section() {
    // The u-iso at any interior station of the extruded arc is the
    // unit circle section: every sample at radius 1 from the axis.
    let surface = rational_quarter_cylinder();
    let iso = surface.iso_curve_at_u(0.37);
    for i in 0..=64 {
        let v = i as f64 / 64.0;
        let p = iso.point_at(v);
        let r = (p.y * p.y + p.z * p.z).sqrt();
        assert!((r - 1.0).abs() < 1e-13, "radius {r} at v={v}");
    }
}

#[test]
fn iso_curve_on_nonrational_patch() {
    let surface = bilinear_patch();
    let iso = surface.iso_curve_at_v(0.5);
    for i in 0..=8 {
        let u = i as f64 / 8.0;
        let d = iso.point_at(u).distance_to(&surface.point_at(u, 0.5));
        assert!(d < 1e-14);
    }
}
