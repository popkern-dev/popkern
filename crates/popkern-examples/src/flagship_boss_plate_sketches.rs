//! Profile sketches and the extrusion path of the flagship boss plate.

use super::*;

/// Closed degree-1 rectangle wire in the z=0 plane, centered at the
/// origin: the plate sketch.
pub(super) fn rectangle_sketch(width: f64, depth: f64) -> NurbsCurve {
    let (hw, hd) = (0.5 * width, 0.5 * depth);
    let pts = vec![
        Point3::new(-hw, -hd, 0.0),
        Point3::new(hw, -hd, 0.0),
        Point3::new(hw, hd, 0.0),
        Point3::new(-hw, hd, 0.0),
        Point3::new(-hw, -hd, 0.0),
    ];
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n],
    }
}

/// Closed rational quadratic circle wire (standard 9-CP form) in the
/// plane z = `z`, centered on the z-axis: the boss sketch.
pub(super) fn circle_sketch(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = vec![
        Point3::new(r, 0.0, z),
        Point3::new(r, r, z),
        Point3::new(0.0, r, z),
        Point3::new(-r, r, z),
        Point3::new(-r, 0.0, z),
        Point3::new(-r, -r, z),
        Point3::new(0.0, -r, z),
        Point3::new(r, -r, z),
        Point3::new(r, 0.0, z),
    ];
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts,
        weights,
    }
}

/// Straight extrusion path from `z0` to `z1` along +Z.
pub(super) fn extrusion_path(z0: f64, z1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    }
}
