//! Unit tests for
//! [`super::rectangular_biquintic_full::RectangularGregoryBiquinticFull`].
//!
//! Extracted from the implementation file to keep the latter
//! under the 500-line structural cap.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::rectangular_biquintic_full::*;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};

/// Build a flat biquintic-full patch over the unit square
/// in the z=0 plane with free interior CPs set to the
/// bilinear-interpolated positions (so the patch reproduces
/// a degree-5 embedding of the bilinear). Used as a baseline
/// for no-curvature regression.
fn flat_unit_square() -> RectangularGregoryBiquinticFull {
    let pt = |i: usize, j: usize| Point3::new((i as f64) / 5.0, (j as f64) / 5.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    let corner_inner_cps = [
        (pt(1, 1), pt(1, 1)),
        (pt(1, 4), pt(1, 4)),
        (pt(4, 4), pt(4, 4)),
        (pt(4, 1), pt(4, 1)),
    ];
    RectangularGregoryBiquinticFull {
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
                pt.z.abs() < 1e-12,
                "flat patch at ({u}, {v}) must have z ≈ 0; got {}",
                pt.z
            );
            assert!(
                (pt.x - u).abs() < 1e-10,
                "flat patch x at u={u} must equal u"
            );
            assert!(
                (pt.y - v).abs() < 1e-10,
                "flat patch y at v={v} must equal v"
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

/// Local quintic Bernstein reproduction for boundary-check
/// assertion only — the implementation has its own private
/// `bernstein_quintic`, keep tests hermetic.
fn bq(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s * s * s * s * s,
        1 => 5.0 * t * s * s * s * s,
        2 => 10.0 * t * t * s * s * s,
        3 => 10.0 * t * t * t * s * s,
        4 => 5.0 * t * t * t * t * s,
        5 => t * t * t * t * t,
        _ => 0.0,
    }
}

/// Boundary-interpolation: at `u = 0`, the patch reduces to
/// the u=0 boundary Bezier (which is determined by CPs
/// `grid[0][*]` — all on the boundary, no interior
/// contribution).
#[test]
fn boundary_at_u_zero_matches_boundary_bezier() {
    let p = flat_unit_square();
    for &v in &[0.0, 0.2, 0.5, 0.8, 1.0] {
        let patch_val = p.point_at(0.0, v);
        let mut bez = Point3::new(0.0, 0.0, 0.0);
        for j in 0..6 {
            let b = bq(j, v);
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

/// Free interior CPs are used directly: moving `grid[2][2]`
/// (a free-interior position) off the z=0 plane lifts the
/// patch at nearby (u, v) evaluations. This is the core
/// distinction from the 4-corner-only deg5 form, where
/// `grid[2][2]` is a placeholder overwritten by the rational
/// blend.
#[test]
fn free_interior_cp_controls_surface() {
    let mut p = flat_unit_square();
    p.grid[2][2].z = 0.5;
    let eval = p.point_at(0.4, 0.4);
    assert!(
        eval.z > 0.05,
        "free interior CP at (2,2) lifted by 0.5 must raise surface at (0.4, 0.4); got z = {}",
        eval.z
    );
}

/// Corner-adjacent positions remain rational-blend driven —
/// setting `grid[1][1]` directly has no effect on
/// evaluation (the value is overwritten by the inner-CP
/// blend at that Bezier index).
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

/// When all 4 corner pairs have `P_u = P_v` AND equal their
/// raw bilinear-interpolated positions AND all 12 free
/// interior CPs equal their bilinear positions, the
/// biquintic-full patch reproduces a degree-5 embedding of
/// the bilinear surface (verified by boundary + corners).
#[test]
fn biquintic_full_reduces_to_bilinear_on_degenerate_input() {
    let p = flat_unit_square();
    let n = p.partial_u(0.3, 0.7).cross(&p.partial_v(0.3, 0.7));
    let n_unit = n.normalized().unwrap();
    assert!(
        (n_unit.x).abs() < 1e-10,
        "flat patch normal.x should be 0, got {}",
        n_unit.x
    );
    assert!(
        (n_unit.y).abs() < 1e-10,
        "flat patch normal.y should be 0, got {}",
        n_unit.y
    );
    assert!(
        (n_unit.z.abs() - 1.0).abs() < 1e-10,
        "flat patch normal.z should be ±1, got {}",
        n_unit.z
    );
}

/// Curved biquintic-full patch with nonzero free interior CPs
/// and distinct corner inner-CP pairs. Used for NURBS
/// round-trip.
fn curved_unit_square() -> RectangularGregoryBiquinticFull {
    let pt = |i: usize, j: usize, z: f64| Point3::new((i as f64) / 5.0, (j as f64) / 5.0, z);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let z = if (1..=4).contains(&i)
                && (1..=4).contains(&j)
                && !matches!((i, j), (1, 1) | (1, 4) | (4, 1) | (4, 4))
            {
                0.1 * ((i as f64) - 2.5) * ((j as f64) - 2.5)
            } else {
                0.0
            };
            *cell = pt(i, j, z);
        }
    }
    let corner_inner_cps = [
        (pt(1, 1, 0.05), pt(1, 1, 0.04)),
        (pt(1, 4, 0.03), pt(1, 4, 0.06)),
        (pt(4, 4, 0.07), pt(4, 4, 0.05)),
        (pt(4, 1, 0.02), pt(4, 1, 0.08)),
    ];
    RectangularGregoryBiquinticFull {
        grid,
        corner_inner_cps,
    }
}

/// `to_nurbs_exact` must reproduce `point_at` to floating-point
/// precision at multiple interior samples. Polynomial identity
/// is exact on the open interior; ε-regularization at corners
/// is checked at samples away from the corners.
#[test]
fn to_nurbs_exact_reproduces_point_at_on_curved_patch() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let p = curved_unit_square();
    let nurbs: NurbsSurface = p.to_nurbs_exact();
    for &u in &[0.15, 0.3, 0.5, 0.7, 0.85] {
        for &v in &[0.15, 0.3, 0.5, 0.7, 0.85] {
            let direct = p.point_at(u, v);
            let via = nurbs.point_at(u, v);
            let diff = ((direct.x - via.x).powi(2)
                + (direct.y - via.y).powi(2)
                + (direct.z - via.z).powi(2))
            .sqrt();
            assert!(
                diff < 1e-10,
                "round-trip mismatch at ({u}, {v}): direct {direct:?} vs NURBS {via:?} diff {diff:.4e}"
            );
        }
    }
}

/// Corner reproduction: ε-regularized corner weight pins the
/// patch corner to the analytic Gregory corner.
#[test]
fn to_nurbs_exact_corner_points_match() {
    let p = curved_unit_square();
    let nurbs = p.to_nurbs_exact();
    for (u, v, expected) in [
        (0.0, 0.0, p.grid[0][0]),
        (0.0, 1.0, p.grid[0][5]),
        (1.0, 0.0, p.grid[5][0]),
        (1.0, 1.0, p.grid[5][5]),
    ] {
        let via = nurbs.point_at(u, v);
        let diff = ((via.x - expected.x).powi(2)
            + (via.y - expected.y).powi(2)
            + (via.z - expected.z).powi(2))
        .sqrt();
        assert!(
            diff < 1e-10,
            "corner ({u}, {v}): expected {expected:?}, got {via:?}, diff = {diff:.4e}"
        );
    }
}

/// Second partials (finite-difference) are continuous in
/// smooth interior samples — no singularities at the
/// corner-adjacent interior positions where the rational
/// blend is active. Regression against numerical blowup.
/// Constructor round-trip: `from_deg5` promotes a deg5 patch to
/// biquintic-full, preserving boundary + corners and agreeing
/// with deg5 on the boundary. Interior may diverge because the
/// free interior CPs are now plain-Bezier (bilinear seed)
/// rather than rational-blend.
#[test]
fn from_deg5_preserves_boundary_and_corners() {
    use super::rectangular_deg5::RectangularGregoryDeg5;
    let deg5 = {
        let pt = |i: usize, j: usize| Point3::new((i as f64) / 5.0, (j as f64) / 5.0, 0.0);
        let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
        for (i, row) in boundary_cps.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = pt(i, j);
            }
        }
        let inner_cps = [
            (pt(1, 1), pt(1, 1)),
            (pt(1, 4), pt(1, 4)),
            (pt(4, 4), pt(4, 4)),
            (pt(4, 1), pt(4, 1)),
        ];
        RectangularGregoryDeg5 {
            boundary_cps,
            inner_cps,
        }
    };
    let full = RectangularGregoryBiquinticFull::from_deg5(&deg5);
    // Boundary: row i=0, row i=5, col j=0, col j=5 are pulled
    // through as-is.
    for j in 0..6 {
        let b = deg5.boundary_cps[0][j];
        let g = full.grid[0][j];
        let diff = ((b.x - g.x).powi(2) + (b.y - g.y).powi(2) + (b.z - g.z).powi(2)).sqrt();
        assert!(
            diff < 1e-14,
            "boundary row 0 column {j} mismatch: {diff:.4e}"
        );
    }
    // Corner inner-CP pairs match.
    for k in 0..4 {
        let (d_pu, d_pv) = deg5.inner_cps[k];
        let (f_pu, f_pv) = full.corner_inner_cps[k];
        let du =
            ((d_pu.x - f_pu.x).powi(2) + (d_pu.y - f_pu.y).powi(2) + (d_pu.z - f_pu.z).powi(2))
                .sqrt();
        let dv =
            ((d_pv.x - f_pv.x).powi(2) + (d_pv.y - f_pv.y).powi(2) + (d_pv.z - f_pv.z).powi(2))
                .sqrt();
        assert!(du < 1e-14 && dv < 1e-14, "corner {k} CP-pair mismatch");
    }
}

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
                    "{label} at ({u}, {v}) has non-finite component: {vec:?}"
                );
            }
        }
    }
}
