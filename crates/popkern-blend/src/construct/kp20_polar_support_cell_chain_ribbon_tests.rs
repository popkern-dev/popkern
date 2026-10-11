// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};

use super::super::BoundarySide;
use super::super::chain_author::boundary_controls;
use super::super::row_fit::bernstein9;
use super::{project_row, second_row_controls};

/// Patch whose UMin row runs along x with the interior row-1
/// controls lifted out of the target tangent plane (z), so the
/// initial cross-derivative violates perpendicularity to Z. The
/// row-1 endpoint controls stay in-plane — they are frozen by
/// the projection (they belong to transverse chain curves).
fn tilted_patch() -> RectangularGregoryBiNonic {
    let mut patch = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            let lifted = i == 1 && (1..=8).contains(&j);
            patch.grid[i][j] = Point3::new(
                j as f64 / 9.0,
                i as f64 * 0.1,
                if lifted { 0.05 } else { 0.0 },
            );
        }
    }
    patch
}

fn cross_derivative_dot(patch: &RectangularGregoryBiNonic, t: f64, normal: Vec3) -> f64 {
    let row0 = boundary_controls(patch, BoundarySide::UMin);
    let row1 = second_row_controls(patch, BoundarySide::UMin);
    let basis = bernstein9(t);
    let mut dot = 0.0;
    for i in 0..10 {
        let delta = Vec3::new(
            row1[i].x - row0[i].x,
            row1[i].y - row0[i].y,
            row1[i].z - row0[i].z,
        );
        dot += basis[i] * delta.dot(&normal);
    }
    dot
}

#[test]
fn projection_drives_cross_derivative_into_plane() {
    let mut patch = tilted_patch();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let before = cross_derivative_dot(&patch, 0.5, normal).abs();
    assert!(before > 1.0e-3, "fixture must start out of plane");
    project_row(&mut patch, BoundarySide::UMin, |_| normal).expect("projection must solve");
    // Perpendicular at interior samples (frozen endpoints keep
    // their old values, so measure away from the corners).
    for sample in 1..=19 {
        let t = sample as f64 / 20.0;
        if !(0.08..=0.92).contains(&t) {
            continue;
        }
        let after = cross_derivative_dot(&patch, t, normal).abs();
        assert!(
            after < 1.0e-6,
            "cross-derivative should be in-plane at t={t}, got {after}"
        );
    }
}

#[test]
fn projection_freezes_row0_and_row1_endpoints() {
    let mut patch = tilted_patch();
    let row0_before = boundary_controls(&patch, BoundarySide::UMin);
    let row1_before = second_row_controls(&patch, BoundarySide::UMin);
    project_row(&mut patch, BoundarySide::UMin, |_| Vec3::new(0.0, 0.0, 1.0))
        .expect("projection must solve");
    let row0_after = boundary_controls(&patch, BoundarySide::UMin);
    let row1_after = second_row_controls(&patch, BoundarySide::UMin);
    for i in 0..10 {
        assert!(
            row0_before[i].distance_to(&row0_after[i]) < 1.0e-15,
            "row 0 must not move"
        );
    }
    assert!(
        row1_before[0].distance_to(&row1_after[0]) < 1.0e-15
            && row1_before[9].distance_to(&row1_after[9]) < 1.0e-15,
        "row-1 endpoint controls are frozen (they belong to transverse chain curves)"
    );
}
