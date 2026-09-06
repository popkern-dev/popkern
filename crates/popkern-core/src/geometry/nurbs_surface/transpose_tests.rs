// This file is only compiled under #[cfg(test)] (see mod.rs:
// `#[cfg(test)] #[path = "transpose_tests.rs"] mod transpose_tests;`).
use super::*;
use crate::geometry::{Point3, SurfaceEval};

/// The bilinear test patch (a private copy of the one in `tests.rs` —
/// the modules are siblings and keep their helpers local).
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
        weights: vec![1.0; 4],
    }
}

/// PK-456: transposition is exact — S'(u, v) = S(v, u) across a
/// non-symmetric rational grid, and the normal reverses.
#[test]
fn transposed_surface_evaluates_swapped_and_reverses_the_normal() {
    let mut surface = bilinear_patch();
    // Break the symmetry: lift one corner and weight another.
    surface.control_points[1].z = 0.7;
    surface.weights[2] = 2.5;
    let t = surface.transposed();
    for i in 0..=4 {
        for j in 0..=4 {
            let (u, v) = (i as f64 / 4.0, j as f64 / 4.0);
            let a = surface.point_at(u, v);
            let b = t.point_at(v, u);
            assert!(a.distance_to(&b) < 1e-12, "S({u},{v}) vs S'({v},{u})");
        }
    }
    let n = surface.normal_at(0.3, 0.6);
    let nt = t.normal_at(0.6, 0.3);
    assert!(
        (n.normalized().unwrap() + nt.normalized().unwrap()).length() < 1e-12,
        "transposition reverses the chart normal"
    );
}
