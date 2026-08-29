use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};

use super::super::band_derive::{derive_row1, solve_leg};
use super::super::row_fit::bernstein9;
use super::refill_interior_coons;

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
    let mut clamped = 0usize;
    let mut one_sided = 0usize;
    let derived = derive_row1(&row0, &row1, &field, &mut clamped, &mut one_sided)
        .expect("derivation must solve");
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

/// A degenerate (chain-collinear) start leg must trigger
/// one-sided derivation: the band copies the healthy end's
/// transverse magnitude instead of going tangential, and stays
/// in the shared plane.
#[test]
fn degenerate_corner_derives_one_sided() {
    let mut row0 = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in row0.iter_mut().enumerate() {
        *point = Point3::new(index as f64 / 9.0, 0.0, 0.0);
    }
    let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in row1.iter_mut().enumerate() {
        *point = Point3::new(index as f64 / 9.0, 0.1, 0.0);
    }
    // Start leg chain-collinear (the T-stem 180-degree case):
    // row1[0] sits BEHIND row0[0] along the chain, no transverse
    // part. End leg healthy.
    row1[0] = Point3::new(-0.05, 0.0, 0.0);
    let field = [Vec3::new(0.0, 1.0, 0.0); 9];
    let mut clamped = 0usize;
    let mut one_sided = 0usize;
    let derived = derive_row1(&row0, &row1, &field, &mut clamped, &mut one_sided)
        .expect("derivation must solve");
    assert_eq!(one_sided, 1, "degenerate start leg must derive one-sided");
    // Transverse component must stay POSITIVE (band does not
    // collapse), monotonically NON-DECREASING from the degenerate
    // corner toward the healthy end (the borrow is scaled by the
    // leg-magnitude ratio, so the band's transverse content at
    // the degenerate corner matches the small frozen leg's scale
    // and grows toward the healthy end), and reach near 0.1 at
    // the healthy end's neighbor.
    for point in derived.iter().take(9).skip(1) {
        assert!(point.y > 0.0, "band must not collapse: {}", point.y);
        assert!(point.z.abs() < 1.0e-12, "must stay in the shared plane");
    }
    for window in derived
        .iter()
        .take(9)
        .skip(1)
        .collect::<Vec<_>>()
        .windows(2)
    {
        assert!(
            window[1].y >= window[0].y - 1.0e-12,
            "transverse must grow toward the healthy end, got {} then {}",
            window[0].y,
            window[1].y
        );
    }
    assert!(
        (derived[8].y - 0.1).abs() < 2.0e-2,
        "transverse near the healthy end must approach its scale, got {}",
        derived[8].y
    );
}

#[test]
// The (i, j) indices address `patch.grid` and `expected` in lockstep and
// appear in the assert message; iterator form would obscure the grid walk.
#[allow(clippy::needless_range_loop)]
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
