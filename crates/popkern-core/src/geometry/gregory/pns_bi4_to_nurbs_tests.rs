//! Tests for `pns_bi4_to_nurbs` — bi-4 patch → NurbsSurface.

#![allow(clippy::needless_range_loop)]

use super::pns_bi4_to_nurbs;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn bi4_bernstein(i: usize, u: f64) -> f64 {
    let v = 1.0 - u;
    match i {
        0 => v * v * v * v,
        1 => 4.0 * u * v * v * v,
        2 => 6.0 * u * u * v * v,
        3 => 4.0 * u * u * u * v,
        4 => u * u * u * u,
        _ => panic!("i out of range"),
    }
}

/// Evaluate a bi-4 patch directly via the Bernstein basis.
fn eval_bi4(patch: &[[Point3; 5]; 5], u: f64, v: f64) -> Point3 {
    let mut acc = Point3::new(0.0, 0.0, 0.0);
    for i in 0..5 {
        for j in 0..5 {
            let w = bi4_bernstein(i, u) * bi4_bernstein(j, v);
            let p = patch[i][j];
            acc.x += w * p.x;
            acc.y += w * p.y;
            acc.z += w * p.z;
        }
    }
    acc
}

fn make_test_patch() -> [[Point3; 5]; 5] {
    let mut patch = [[pt(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            patch[i][j] = pt(
                i as f64 * 0.25,
                j as f64 * 0.25,
                (i as f64 * j as f64) * 0.1,
            );
        }
    }
    patch
}

#[test]
fn nurbs_conversion_matches_bi4_evaluation() {
    let patch = make_test_patch();
    let nurbs = pns_bi4_to_nurbs(&patch);
    for u_n in 0..=5 {
        for v_n in 0..=5 {
            let u = (u_n as f64) / 5.0;
            let v = (v_n as f64) / 5.0;
            let expected = eval_bi4(&patch, u, v);
            let actual = nurbs.point_at(u, v);
            let dx = (actual.x - expected.x).abs();
            let dy = (actual.y - expected.y).abs();
            let dz = (actual.z - expected.z).abs();
            assert!(
                dx < 1e-10 && dy < 1e-10 && dz < 1e-10,
                "at (u={u}, v={v}) bi-4 = {expected:?}, nurbs = {actual:?}"
            );
        }
    }
}

#[test]
fn nurbs_corners_match_patch_corners() {
    let patch = make_test_patch();
    let nurbs = pns_bi4_to_nurbs(&patch);
    assert_eq!(nurbs.point_at(0.0, 0.0), patch[0][0]);
    assert_eq!(nurbs.point_at(0.0, 1.0), patch[0][4]);
    assert_eq!(nurbs.point_at(1.0, 0.0), patch[4][0]);
    assert_eq!(nurbs.point_at(1.0, 1.0), patch[4][4]);
}

#[test]
fn nurbs_parameter_range_is_0_1() {
    let patch = make_test_patch();
    let nurbs = pns_bi4_to_nurbs(&patch);
    let ((u0, u1), (v0, v1)) = nurbs.parameter_range();
    assert_eq!((u0, u1), (0.0, 1.0));
    assert_eq!((v0, v1), (0.0, 1.0));
}
