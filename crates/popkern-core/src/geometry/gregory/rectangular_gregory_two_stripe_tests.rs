//! Unit tests for
//! [`super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe`]
//! — the two-stripe construction patch type.

use super::rectangular_gregory_two_stripe::*;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};

fn flat_unit_square() -> RectangularGregoryTwoStripe {
    let pt = |i: usize, j: usize| Point3::new((i as f64) / 7.0, (j as f64) / 7.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    RectangularGregoryTwoStripe { grid }
}

#[test]
fn flat_patch_evaluates_in_plane() {
    let p = flat_unit_square();
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = p.point_at(u, v);
            assert!(pt.z.abs() < 1e-12, "z = {}", pt.z);
            assert!((pt.x - u).abs() < 1e-10);
            assert!((pt.y - v).abs() < 1e-10);
        }
    }
}

#[test]
fn corner_adjacent_cp_now_affects_evaluation() {
    // In two-stripe, (1, 1) is a PLAIN Bezier CP — moving it
    // MUST change the evaluated surface at nearby (u, v). This
    // is the structural change vs biseptic where (1, 1) is
    // an ignored placeholder.
    let mut p = flat_unit_square();
    let before = p.point_at(0.15, 0.15);
    p.grid[1][1].z = 0.5;
    let after = p.point_at(0.15, 0.15);
    assert!(
        (after.z - before.z).abs() > 0.01,
        "grid[1][1] must be a plain Bezier CP (no rational blend placeholder); \
         before={before:?} after={after:?}"
    );
}

#[test]
fn partial_u_matches_finite_difference() {
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
                "({u}, {v}): analytic {analytic:?} vs FD {fd:?}"
            );
        }
    }
}

#[test]
fn partial_v_matches_finite_difference() {
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
            assert!(diff < 1e-6);
        }
    }
}

/// Phase-10b prerequisite: `second_partial_uu` returns the
/// analytic value (no finite differencing). Verify against FD
/// — this is the G2 solver's linear coefficient foundation.
#[test]
fn second_partial_uu_matches_finite_difference() {
    // Pick a patch with curvature so the second partial is
    // nonzero.
    let mut p = flat_unit_square();
    p.grid[2][2].z = 0.4;
    p.grid[3][3].z = 0.3;
    p.grid[4][4].z = 0.2;
    let h = 1e-4;
    for &u in &[0.3, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.7] {
            let analytic = p.second_partial_uu(u, v);
            let du_p = p.partial_u(u + h, v);
            let du_m = p.partial_u(u - h, v);
            let fd = Vec3::new(
                (du_p.x - du_m.x) / (2.0 * h),
                (du_p.y - du_m.y) / (2.0 * h),
                (du_p.z - du_m.z) / (2.0 * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(
                diff < 1e-4,
                "suu ({u}, {v}): analytic {analytic:?} vs FD {fd:?} diff={diff:.4e}"
            );
        }
    }
}

/// NURBS round-trip: pure bidegree-7 Bezier converts losslessly
/// (unit weights, no corner ε-regularisation needed).
#[test]
fn to_nurbs_exact_reproduces_point_at() {
    let mut p = flat_unit_square();
    p.grid[2][3].z = 0.2;
    p.grid[4][5].z = 0.15;
    p.grid[1][1].z = 0.1; // the newly-free corner position
    let nurbs = p.to_nurbs_exact();
    for &u in &[0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0] {
        for &v in &[0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0] {
            let direct = p.point_at(u, v);
            let via = nurbs.point_at(u, v);
            let diff = ((direct.x - via.x).powi(2)
                + (direct.y - via.y).powi(2)
                + (direct.z - via.z).powi(2))
            .sqrt();
            assert!(diff < 1e-12, "({u}, {v}): diff={diff:.4e}");
        }
    }
}

/// `from_biseptic` must preserve the 28 boundary CPs verbatim
/// (row 0, row 7, col 0, col 7) and replace the 4 corner-
/// adjacent positions with the biseptic's rational-blend value
/// evaluated at `(i/7, j/7)`.
#[test]
fn from_biseptic_preserves_boundary() {
    use super::rectangular_biseptic::RectangularGregoryBiseptic;
    let pt = |i: usize, j: usize| Point3::new((i as f64) / 7.0, (j as f64) / 7.0, 0.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = pt(i, j);
        }
    }
    let corner_inner_cps = [
        (Point3::new(0.14, 0.14, 0.05), Point3::new(0.14, 0.14, 0.04)),
        (Point3::new(0.14, 0.86, 0.03), Point3::new(0.14, 0.86, 0.06)),
        (Point3::new(0.86, 0.86, 0.07), Point3::new(0.86, 0.86, 0.05)),
        (Point3::new(0.86, 0.14, 0.02), Point3::new(0.86, 0.14, 0.08)),
    ];
    let bs = RectangularGregoryBiseptic {
        grid,
        corner_inner_cps,
    };
    let ts = RectangularGregoryTwoStripe::from_biseptic(&bs);
    // Boundary positions must match verbatim.
    for j in 0..8 {
        for i in [0usize, 7usize] {
            let b = bs.grid[i][j];
            let t = ts.grid[i][j];
            let d = ((b.x - t.x).powi(2) + (b.y - t.y).powi(2) + (b.z - t.z).powi(2)).sqrt();
            assert!(d < 1e-14, "boundary row {i} col {j} diff={d:.4e}");
        }
    }
    for i in 0..8 {
        for j in [0usize, 7usize] {
            let b = bs.grid[i][j];
            let t = ts.grid[i][j];
            let d = ((b.x - t.x).powi(2) + (b.y - t.y).powi(2) + (b.z - t.z).powi(2)).sqrt();
            assert!(d < 1e-14, "boundary col {j} row {i} diff={d:.4e}");
        }
    }
    // Corner-adjacent positions must NOT equal the biseptic's
    // placeholder (which was the bilinear seed) — they should
    // now be rational-blend evaluated values, which depend on
    // the corner_inner_cps (lifted in z).
    for (i, j) in [(1usize, 1usize), (1, 6), (6, 1), (6, 6)] {
        let t = ts.grid[i][j];
        assert!(
            t.z.abs() > 1e-6,
            "({i}, {j}) should have non-flat z from blend evaluation; got z={}",
            t.z
        );
    }
}

#[allow(dead_code)]
fn _surface_eval_compat<S: SurfaceEval>(_s: &S) {}
