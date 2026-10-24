use super::*;
use popkern_core::geometry::{Point3, Vec3};
use std::f64::consts::{PI, TAU};

fn cone() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.5_f64.atan(),
    }
}

/// The notched cone wall's loop as the folded convention writes it:
/// rim 2π → 0, seam down to the apex, the other seam copy back up at
/// u = 0 around a notch at u ∈ [−π/2, 0].
fn folded_notched_wall() -> Vec<(f64, f64)> {
    vec![
        (TAU, 2.0),
        (1.5 * PI, 2.0),
        (PI, 2.0),
        (PI / 2.0, 2.0),
        (0.0, 2.0),
        (0.0, 1.0),
        (0.0, 0.0),
        (0.0, 0.2),
        (-PI / 2.0, 0.2),
        (-PI / 2.0, 1.2),
        (0.0, 1.2),
        (0.0, 2.0),
    ]
}

#[test]
fn a_loop_folded_through_the_apex_walks_the_apex_row() {
    let s = cone();
    let mut uv = folded_notched_wall();
    let mut xyz: Vec<_> = uv.iter().map(|&(u, v)| s.point_at(u, v)).collect();
    assert!(lift_through_pole(&mut uv, &mut xyz, &s, TAU, true));
    assert_eq!(uv.len(), xyz.len());
    let (first, last) = (uv[0], uv[uv.len() - 1]);
    assert!((first.0 - last.0).abs() < 1e-12 && (first.1 - last.1).abs() < 1e-12);
    // The wall less its notch: 2π·2 − (π/2)·1.
    let area2: f64 = (0..uv.len())
        .map(|i| {
            let (a, b) = (uv[i], uv[(i + 1) % uv.len()]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum();
    assert!((area2.abs() / 2.0 - (2.0 * TAU - PI / 2.0)).abs() < 1e-12);
}

#[test]
fn a_loop_closing_without_a_period_is_left_alone() {
    let s = cone();
    let mut uv = vec![(-0.3, 1.0), (0.3, 1.0), (0.3, 1.5), (-0.3, 1.5)];
    let mut xyz: Vec<_> = uv.iter().map(|&(u, v)| s.point_at(u, v)).collect();
    let before = uv.clone();
    assert!(!lift_through_pole(&mut uv, &mut xyz, &s, TAU, true));
    assert_eq!(uv, before);
}

#[test]
fn a_full_period_loop_with_no_pole_is_left_alone() {
    let s = cone();
    // Rim to rim with no apex visit: nowhere to absorb the period.
    let mut uv = vec![
        (TAU, 2.0),
        (1.5 * PI, 2.0),
        (PI, 2.0),
        (PI / 2.0, 2.0),
        (0.0, 2.0),
        (0.0, 1.0),
    ];
    let mut xyz: Vec<_> = uv.iter().map(|&(u, v)| s.point_at(u, v)).collect();
    assert!(!lift_through_pole(&mut uv, &mut xyz, &s, TAU, true));
}
