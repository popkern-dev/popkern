use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};

use super::super::row_fit::bernstein9;
use super::{derive_row1, refill_interior_coons, solve_leg};

#[test]
fn solve_leg_is_exact_for_coplanar_legs() {
    let u = Vec3::new(1.0, 0.0, 0.0);
    let v = Vec3::new(0.3, 1.0, 0.0);
    let leg = u * 2.0 + v * -0.7;
    let (alpha, beta) = solve_leg(leg, u, v).expect("solve must succeed");
    assert!((alpha - 2.0).abs() < 1.0e-12 && (beta + 0.7).abs() < 1.0e-12);
}

#[test]
fn solve_leg_declines_parallel_basis() {
    let u = Vec3::new(1.0, 0.0, 0.0);
    let v = Vec3::new(2.0, 0.0, 0.0);
    assert!(solve_leg(Vec3::new(1.0, 1.0, 0.0), u, v).is_none());
}

/// Planar configuration: row 0 along x in the z = 0 plane, end
/// legs and ribbon field in-plane. The derived row-1 band must
/// stay exactly in the plane (the cross-derivative lies in
/// span(C', R) pointwise), even though the OLD row-1 interior
/// controls are far out of plane.
#[test]
fn derived_band_stays_in_the_shared_plane() {
    let mut row0 = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in row0.iter_mut().enumerate() {
        *point = Point3::new(index as f64 / 9.0, 0.0, 0.0);
    }
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in row1.iter_mut().enumerate() {
        // Old band: lifted far out of plane in z except the
        // frozen endpoints, which carry exact in-plane legs.
        let z = if index == 0 || index == 9 { 0.0 } else { 0.4 };
        *point = Point3::new(index as f64 / 9.0, 0.1, z);
    }
    let field = [Vec3::new(0.0, 1.0, 0.0); 9];
    let derived = derive_row1(&row0, &row1, &field).expect("derivation must solve");
    for (index, point) in derived.iter().enumerate() {
        assert!(
            point.z.abs() < 1.0e-12,
            "derived control {index} must lie in the shared plane, got z={}",
            point.z
        );
    }
    assert!(
        derived[0].distance_to(&row1[0]) < 1.0e-15 && derived[9].distance_to(&row1[9]) < 1.0e-15,
        "frozen endpoint controls must be reproduced"
    );
    // The derived cross-derivative must interpolate the end legs.
    let delta_at = |row: &[Point3; 10], t: f64| -> Vec3 {
        let basis = bernstein9(t);
        let mut delta = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            delta = delta + (row[i] - row0[i]) * basis[i];
        }
        delta * 9.0
    };
    let leg0 = (row1[0] - row0[0]) * 9.0;
    let derived_leg0 = delta_at(&derived, 0.0);
    assert!(
        (derived_leg0 - leg0).dot(&(derived_leg0 - leg0)).sqrt() < 1.0e-12,
        "derived band must reproduce the exact start leg"
    );
}

#[test]
fn coons_fill_reproduces_bilinear_interior() {
    let mut patch = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            // Bilinear surface: any Coons fill of its bands must
            // reproduce it exactly.
            patch.grid[i][j] = Point3::new(i as f64, j as f64, 0.5 * i as f64 * j as f64);
        }
    }
    let expected = patch.grid;
    // Scramble the interior, keep the bands.
    for i in 2..8 {
        for j in 2..8 {
            patch.grid[i][j] = Point3::new(-1.0, -1.0, -1.0);
        }
    }
    refill_interior_coons(&mut patch);
    for i in 2..8 {
        for j in 2..8 {
            assert!(
                patch.grid[i][j].distance_to(&expected[i][j]) < 1.0e-12,
                "interior ({i},{j}) must match the bilinear surface"
            );
        }
    }
}
