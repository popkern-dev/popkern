//! Unit tests for
//! [`super::rectangular_biseptic::RectangularGregoryBiseptic`].

use super::rectangular_biseptic::*;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};

fn flat_unit_square() -> RectangularGregoryBiseptic {
    let pt = |i: usize, j: usize| Point3::new((i as f64) / 7.0, (j as f64) / 7.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    let corner_inner_cps = [
        (pt(1, 1), pt(1, 1)),
        (pt(1, 6), pt(1, 6)),
        (pt(6, 6), pt(6, 6)),
        (pt(6, 1), pt(6, 1)),
    ];
    RectangularGregoryBiseptic {
        grid,
        corner_inner_cps,
    }
}

#[test]
fn flat_patch_evaluates_in_plane() {
    let p = flat_unit_square();
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = p.point_at(u, v);
            assert!(
                pt.z.abs() < 1e-11,
                "flat patch at ({u}, {v}) must have z ≈ 0; got {}",
                pt.z
            );
            assert!(
                (pt.x - u).abs() < 1e-10,
                "flat patch x at u={u} must equal u; got {}",
                pt.x
            );
            assert!(
                (pt.y - v).abs() < 1e-10,
                "flat patch y at v={v} must equal v; got {}",
                pt.y
            );
        }
    }
}

#[test]
fn flat_patch_partial_u_matches_finite_difference() {
    let p = flat_unit_square();
    let h = 1e-5;
    for &u in &[0.2, 0.4, 0.6, 0.8] {
        for &v in &[0.2, 0.4, 0.6, 0.8] {
            let analytic = p.partial_u(u, v);
            let p_plus = p.point_at(u + h, v);
            let p_minus = p.point_at(u - h, v);
            let fd = Vec3::new(
                (p_plus.x - p_minus.x) / (2.0 * h),
                (p_plus.y - p_minus.y) / (2.0 * h),
                (p_plus.z - p_minus.z) / (2.0 * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(
                diff < 1e-6,
                "partial_u mismatch at ({u}, {v}): analytic {analytic:?} vs FD {fd:?}"
            );
        }
    }
}

#[test]
fn flat_patch_partial_v_matches_finite_difference() {
    let p = flat_unit_square();
    let h = 1e-5;
    for &u in &[0.2, 0.4, 0.6, 0.8] {
        for &v in &[0.2, 0.4, 0.6, 0.8] {
            let analytic = p.partial_v(u, v);
            let p_plus = p.point_at(u, v + h);
            let p_minus = p.point_at(u, v - h);
            let fd = Vec3::new(
                (p_plus.x - p_minus.x) / (2.0 * h),
                (p_plus.y - p_minus.y) / (2.0 * h),
                (p_plus.z - p_minus.z) / (2.0 * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(diff < 1e-6, "partial_v mismatch at ({u}, {v})");
        }
    }
}

/// Local septic Bernstein — hermetic test reference.
fn bs(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s.powi(7),
        1 => 7.0 * t * s.powi(6),
        2 => 21.0 * t * t * s.powi(5),
        3 => 35.0 * t.powi(3) * s.powi(4),
        4 => 35.0 * t.powi(4) * s.powi(3),
        5 => 21.0 * t.powi(5) * s * s,
        6 => 7.0 * t.powi(6) * s,
        7 => t.powi(7),
        _ => 0.0,
    }
}

/// At `u = 0` the patch reduces to the u=0 boundary Bezier
/// (degree-7 over `j`). No interior contribution at the boundary.
#[test]
fn boundary_at_u_zero_matches_boundary_bezier() {
    let p = flat_unit_square();
    for &v in &[0.0, 0.2, 0.5, 0.8, 1.0] {
        let patch_val = p.point_at(0.0, v);
        let mut bez = Point3::new(0.0, 0.0, 0.0);
        for j in 0..8 {
            let b = bs(j, v);
            bez.x += b * p.grid[0][j].x;
            bez.y += b * p.grid[0][j].y;
            bez.z += b * p.grid[0][j].z;
        }
        let diff = ((patch_val.x - bez.x).powi(2)
            + (patch_val.y - bez.y).powi(2)
            + (patch_val.z - bez.z).powi(2))
        .sqrt();
        assert!(diff < 1e-12, "u=0 boundary deviation at v={v}: {diff:.4e}");
    }
}

/// Free interior CPs at row 1 j=2..=5 drive the cross-tangent at
/// the u=0 boundary. Moving `grid[1][3].z` should lift the u=0
/// cross-tangent's z component.
#[test]
fn free_interior_cp_affects_cross_tangent() {
    let mut p = flat_unit_square();
    let xt_before = p.partial_u(0.0, 0.5);
    p.grid[1][3].z = 0.3;
    let xt_after = p.partial_u(0.0, 0.5);
    assert!(
        (xt_after.z - xt_before.z).abs() > 0.1,
        "free interior (1,3) should shift cross-tangent z; before={:?} after={:?}",
        xt_before,
        xt_after
    );
}

/// Corner-adjacent position (1,1) is a rational-blend placeholder;
/// setting `grid[1][1]` directly should have no effect.
#[test]
fn corner_adjacent_cp_placeholder_is_ignored() {
    let mut p = flat_unit_square();
    let orig = p.point_at(0.2, 0.2);
    p.grid[1][1] = Point3::new(100.0, 100.0, 100.0);
    let after = p.point_at(0.2, 0.2);
    let diff =
        ((orig.x - after.x).powi(2) + (orig.y - after.y).powi(2) + (orig.z - after.z).powi(2))
            .sqrt();
    assert!(
        diff < 1e-12,
        "grid[1][1] must be a placeholder overwritten by rational blend; diff = {diff:.4e}"
    );
}

/// Flat patch should have +z normal everywhere.
#[test]
fn flat_patch_normal_is_z() {
    let p = flat_unit_square();
    let n = p.partial_u(0.3, 0.7).cross(&p.partial_v(0.3, 0.7));
    let n_unit = n.normalized().unwrap();
    assert!((n_unit.x).abs() < 1e-10, "normal.x should be 0");
    assert!((n_unit.y).abs() < 1e-10, "normal.y should be 0");
    assert!(
        (n_unit.z.abs() - 1.0).abs() < 1e-10,
        "normal.z should be ±1"
    );
}

/// Second partials are finite at interior samples — the polynomial-identity solver work
/// (G2 solver) will use these; regression against numerical
/// blowup at corner-adjacent rational-blend positions.
#[test]
fn second_partials_finite_at_interior_samples() {
    let p = flat_unit_square();
    for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        for &v in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let suu = p.second_partial_uu(u, v);
            let suv = p.second_partial_uv(u, v);
            let svv = p.second_partial_vv(u, v);
            for (label, vec) in [("suu", suu), ("suv", suv), ("svv", svv)] {
                assert!(
                    vec.x.is_finite() && vec.y.is_finite() && vec.z.is_finite(),
                    "{label} at ({u}, {v}) non-finite: {vec:?}"
                );
            }
        }
    }
}

/// Build a curved biseptic patch with a non-trivial interior CP
/// for round-trip validation of the exact NURBS conversion.
fn curved_unit_square() -> RectangularGregoryBiseptic {
    let mut p = flat_unit_square();
    // Lift a handful of free interior CPs so the round-trip
    // exercises the corner-adjacent rational-blend pathway.
    p.grid[2][3].z = 0.15;
    p.grid[3][4].z = 0.1;
    p.corner_inner_cps[0].0.z = 0.03;
    p.corner_inner_cps[0].1.z = 0.04;
    p.corner_inner_cps[2].0.z = 0.02;
    p.corner_inner_cps[2].1.z = 0.05;
    p
}

/// Exact bidegree-(11, 11) NURBS conversion must reproduce
/// `point_at` to floating-point precision at interior samples.
/// Polynomial identity is exact on the open interior; corner
/// ε-regularisation kicks in at `(u, v) ∈ {0, 1}²`.
#[test]
fn to_nurbs_exact_reproduces_point_at_on_curved_patch() {
    let p = curved_unit_square();
    let nurbs = p.to_nurbs_exact();
    for &u in &[0.15, 0.3, 0.5, 0.7, 0.85] {
        for &v in &[0.15, 0.3, 0.5, 0.7, 0.85] {
            let direct = p.point_at(u, v);
            let via = nurbs.point_at(u, v);
            let diff = ((direct.x - via.x).powi(2)
                + (direct.y - via.y).powi(2)
                + (direct.z - via.z).powi(2))
            .sqrt();
            assert!(
                diff < 1e-9,
                "biseptic NURBS round-trip mismatch at ({u}, {v}): direct {direct:?} vs NURBS {via:?} diff {diff:.4e}"
            );
        }
    }
}

/// At-corner reproduction via the ε-regularisation trick: the
/// four `(0,0), (0,1), (1,0), (1,1)` evaluations must match the
/// analytic Gregory corner positions.
#[test]
fn to_nurbs_exact_corner_points_match() {
    let p = curved_unit_square();
    let nurbs = p.to_nurbs_exact();
    for (u, v, expected) in [
        (0.0, 0.0, p.grid[0][0]),
        (0.0, 1.0, p.grid[0][7]),
        (1.0, 0.0, p.grid[7][0]),
        (1.0, 1.0, p.grid[7][7]),
    ] {
        let via = nurbs.point_at(u, v);
        let diff = ((via.x - expected.x).powi(2)
            + (via.y - expected.y).powi(2)
            + (via.z - expected.z).powi(2))
        .sqrt();
        assert!(
            diff < 1e-10,
            "biseptic corner ({u}, {v}): expected {expected:?}, got {via:?}, diff={diff:.4e}"
        );
    }
}

/// `from_biquintic_full` must preserve patch corners exactly
/// (degree-elevation is lossless at the endpoints of each
/// boundary Bezier) and preserve the 4 corner inner-CP pairs.
#[test]
fn from_biquintic_full_preserves_corners_and_inner_cps() {
    use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
    let biq = {
        let pt = |i: usize, j: usize| Point3::new((i as f64) / 5.0, (j as f64) / 5.0, 0.0);
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = pt(i, j);
            }
        }
        let corner_inner_cps = [
            (Point3::new(0.2, 0.2, 0.05), Point3::new(0.2, 0.2, 0.04)),
            (Point3::new(0.2, 0.8, 0.03), Point3::new(0.2, 0.8, 0.06)),
            (Point3::new(0.8, 0.8, 0.07), Point3::new(0.8, 0.8, 0.05)),
            (Point3::new(0.8, 0.2, 0.02), Point3::new(0.8, 0.2, 0.08)),
        ];
        RectangularGregoryBiquinticFull {
            grid,
            corner_inner_cps,
        }
    };
    let bs = RectangularGregoryBiseptic::from_biquintic_full(&biq);
    for (expected_idx, bs_idx) in [
        ((0usize, 0usize), (0usize, 0usize)),
        ((0, 5), (0, 7)),
        ((5, 0), (7, 0)),
        ((5, 5), (7, 7)),
    ] {
        let ep = biq.grid[expected_idx.0][expected_idx.1];
        let gp = bs.grid[bs_idx.0][bs_idx.1];
        let d = ((ep.x - gp.x).powi(2) + (ep.y - gp.y).powi(2) + (ep.z - gp.z).powi(2)).sqrt();
        assert!(
            d < 1e-14,
            "biquintic corner {expected_idx:?} must equal biseptic corner {bs_idx:?}: diff={d:.4e}"
        );
    }
    for k in 0..4 {
        let (bu, bv) = biq.corner_inner_cps[k];
        let (su, sv) = bs.corner_inner_cps[k];
        let du = ((bu.x - su.x).powi(2) + (bu.y - su.y).powi(2) + (bu.z - su.z).powi(2)).sqrt();
        let dv = ((bv.x - sv.x).powi(2) + (bv.y - sv.y).powi(2) + (bv.z - sv.z).powi(2)).sqrt();
        assert!(du < 1e-14 && dv < 1e-14, "corner {k} CP-pair mismatch");
    }
}

/// Degree elevation preserves the boundary Bezier curve: every
/// sample along the boundary matches the biquintic's boundary
/// sample.
#[test]
fn from_biquintic_full_preserves_boundary_curve() {
    use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
    let biq = {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
        for (i, row) in grid.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                let s = (i as f64) / 5.0;
                let t = (j as f64) / 5.0;
                // Quadratic bump on the u=0 edge.
                let z = if i == 0 { 0.3 * t * (1.0 - t) } else { 0.0 };
                *cell = Point3::new(s, t, z);
            }
        }
        let corner_inner_cps = [
            (Point3::new(0.2, 0.2, 0.0), Point3::new(0.2, 0.2, 0.0)),
            (Point3::new(0.2, 0.8, 0.06), Point3::new(0.2, 0.8, 0.06)),
            (Point3::new(0.8, 0.8, 0.0), Point3::new(0.8, 0.8, 0.0)),
            (Point3::new(0.8, 0.2, 0.0), Point3::new(0.8, 0.2, 0.0)),
        ];
        RectangularGregoryBiquinticFull {
            grid,
            corner_inner_cps,
        }
    };
    let bs = RectangularGregoryBiseptic::from_biquintic_full(&biq);
    for t_sample in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let p_biq = biq.point_at(0.0, t_sample);
        let p_bs = bs.point_at(0.0, t_sample);
        let diff =
            ((p_biq.x - p_bs.x).powi(2) + (p_biq.y - p_bs.y).powi(2) + (p_biq.z - p_bs.z).powi(2))
                .sqrt();
        assert!(
            diff < 1e-10,
            "u=0 boundary disagreement at v={t_sample}: diff={diff:.4e}"
        );
    }
}

// Surface trait stub so tests compile identically to biquintic-full's pattern.
#[allow(dead_code)]
fn _surface_eval_compat<S: SurfaceEval>(_s: &S) {}
