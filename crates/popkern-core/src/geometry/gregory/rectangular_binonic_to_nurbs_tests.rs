//! Tests for [`RectangularGregoryBiNonic::to_nurbs_exact`].
//!
//! A plain bidegree-(9, 9) Bezier patch is *definitionally* a
//! single-interval clamped B-spline of degree 9 with clamped knots
//! `[0;10, 1;10]` and unit weights. The conversion is therefore
//! exact by polynomial identity everywhere — no ε-regularisation,
//! no rational-denominator clearing, no corner special-casing.
//!
//! These tests pin the identity at f64 precision (`< 1e-12`) across
//! a dense `(u, v)` grid, shape-metadata, corner/boundary behaviour,
//! and flat-plane preservation.

use crate::geometry::Point3;
use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::surface::SurfaceEval;

const EPS: f64 = 1e-12;

/// Flat unit-square bidegree-9 Bezier patch on `z = 0`. The grid
/// is the degree-9 linear interpolation of the `(0, 0) → (1, 1)`
/// square, so `point_at(u, v) = (u, v, 0)` exactly.
fn flat_unit_square_binonic() -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.0);
        }
    }
    RectangularGregoryBiNonic { grid }
}

/// Non-trivial "bumpy" bidegree-9 patch — flat base grid with a
/// smooth z-bump in the interior so the NURBS image has real
/// curvature to compare against.
fn bumpy_binonic() -> RectangularGregoryBiNonic {
    let mut p = flat_unit_square_binonic();
    for i in 1..9 {
        for j in 1..9 {
            let s = i as f64 / 9.0;
            let t = j as f64 / 9.0;
            // A mild asymmetric bump so every interior CP moves.
            p.grid[i][j].z =
                0.35 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.12 * (s - 0.5) * (t - 0.5);
        }
    }
    p
}

fn max_abs(a: Point3, b: Point3) -> f64 {
    (a.x - b.x)
        .abs()
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[test]
fn shape_metadata_is_bidegree_9() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    assert_eq!(nurbs.degree_u, 9);
    assert_eq!(nurbs.degree_v, 9);
    assert_eq!(nurbs.num_u, 10);
    assert_eq!(nurbs.num_v, 10);
    assert_eq!(nurbs.control_points.len(), 100);
    assert_eq!(nurbs.weights.len(), 100);
    assert_eq!(nurbs.knots_u.len(), 20);
    assert_eq!(nurbs.knots_v.len(), 20);
}

#[test]
fn all_weights_are_unity() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    for &w in &nurbs.weights {
        assert!((w - 1.0).abs() < EPS, "weight should be 1.0, got {w}");
    }
}

#[test]
fn knot_vector_is_clamped_zero_one() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    // Both knot vectors: ten 0s then ten 1s.
    for k in 0..10 {
        assert_eq!(nurbs.knots_u[k], 0.0, "knots_u[{k}] should be 0.0");
        assert_eq!(nurbs.knots_v[k], 0.0, "knots_v[{k}] should be 0.0");
    }
    for k in 10..20 {
        assert_eq!(nurbs.knots_u[k], 1.0, "knots_u[{k}] should be 1.0");
        assert_eq!(nurbs.knots_v[k], 1.0, "knots_v[{k}] should be 1.0");
    }
}

#[test]
fn control_points_equal_grid_in_row_major_order() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    for (i, row) in p.grid.iter().enumerate() {
        for (j, expected) in row.iter().enumerate() {
            let idx = i * 10 + j;
            let cp = nurbs.control_points[idx];
            assert!(
                max_abs(cp, *expected) < EPS,
                "CP[{i}][{j}] mismatch: expected {expected:?}, got {cp:?}"
            );
        }
    }
}

#[test]
fn corner_evaluation_matches_grid_corners() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    for &(u, v, i, j) in &[
        (0.0_f64, 0.0_f64, 0usize, 0usize),
        (0.0, 1.0, 0, 9),
        (1.0, 0.0, 9, 0),
        (1.0, 1.0, 9, 9),
    ] {
        let expected = p.grid[i][j];
        let got = nurbs.point_at(u, v);
        assert!(
            max_abs(got, expected) < EPS,
            "corner ({u}, {v}): expected {expected:?}, got {got:?}"
        );
    }
}

#[test]
fn nurbs_matches_patch_at_dense_grid() {
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    let params = [0.0_f64, 0.1, 0.25, 0.37, 0.5, 0.63, 0.75, 0.9, 1.0];
    let mut max_err = 0.0_f64;
    for &u in &params {
        for &v in &params {
            let p_patch = p.point_at(u, v);
            let p_nurbs = nurbs.point_at(u, v);
            let err = max_abs(p_patch, p_nurbs);
            if err > max_err {
                max_err = err;
            }
        }
    }
    assert!(
        max_err < EPS,
        "patch vs NURBS disagree: max error {max_err:e} (limit {EPS:e})"
    );
}

#[test]
fn flat_patch_stays_planar_through_nurbs() {
    let p = flat_unit_square_binonic();
    let nurbs = p.to_nurbs_exact();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let pt = nurbs.point_at(u, v);
            assert!(
                pt.z.abs() < EPS,
                "flat patch NURBS should be planar; got z = {} at ({u}, {v})",
                pt.z
            );
            assert!(
                (pt.x - u).abs() < EPS && (pt.y - v).abs() < EPS,
                "flat patch at ({u}, {v}): expected ({u}, {v}, 0), got {pt:?}"
            );
        }
    }
}

#[test]
fn zero_patch_maps_to_origin_everywhere() {
    let p = RectangularGregoryBiNonic::zeros();
    let nurbs = p.to_nurbs_exact();
    for &u in &[0.0, 0.2, 0.5, 0.8, 1.0] {
        for &v in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let pt = nurbs.point_at(u, v);
            assert!(
                pt.x.abs() < EPS && pt.y.abs() < EPS && pt.z.abs() < EPS,
                "zeros patch NURBS should map to origin; got {pt:?} at ({u}, {v})"
            );
        }
    }
}

#[test]
fn boundary_evaluation_matches_boundary_cps() {
    // On u = 0 the NURBS is the degree-9 Bezier in v of the first
    // grid row; its corner values are exactly grid[0][0] and
    // grid[0][9]. Symmetric claim on u = 1, v = 0, v = 1.
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    // u = 0, vary v, compare to patch.
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let pn = nurbs.point_at(0.0, v);
        let pp = p.point_at(0.0, v);
        assert!(
            max_abs(pn, pp) < EPS,
            "u=0 boundary mismatch at v={v}: {pn:?} vs {pp:?}"
        );
    }
    // u = 1, vary v.
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let pn = nurbs.point_at(1.0, v);
        let pp = p.point_at(1.0, v);
        assert!(
            max_abs(pn, pp) < EPS,
            "u=1 boundary mismatch at v={v}: {pn:?} vs {pp:?}"
        );
    }
    // v = 0, vary u.
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let pn = nurbs.point_at(u, 0.0);
        let pp = p.point_at(u, 0.0);
        assert!(
            max_abs(pn, pp) < EPS,
            "v=0 boundary mismatch at u={u}: {pn:?} vs {pp:?}"
        );
    }
    // v = 1, vary u.
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let pn = nurbs.point_at(u, 1.0);
        let pp = p.point_at(u, 1.0);
        assert!(
            max_abs(pn, pp) < EPS,
            "v=1 boundary mismatch at u={u}: {pn:?} vs {pp:?}"
        );
    }
}

#[test]
fn partial_derivatives_agree_with_patch() {
    // NURBS degree-9 derivative must equal the patch's analytic
    // partial at representative interior points. This catches any
    // knot-multiplicity off-by-one in the clamped B-spline.
    let p = bumpy_binonic();
    let nurbs = p.to_nurbs_exact();
    let params = [0.15_f64, 0.4, 0.6, 0.85];
    for &u in &params {
        for &v in &params {
            let du_p = p.partial_u(u, v);
            let du_n = nurbs.partial_u(u, v);
            let dv_p = p.partial_v(u, v);
            let dv_n = nurbs.partial_v(u, v);
            // NURBS `partial_u` in this crate is a finite-difference
            // estimator (see `nurbs_surface/mod.rs`), so tolerate a
            // modest ULP budget on the gradient comparison.
            let err_u = ((du_p.x - du_n.x).abs())
                .max((du_p.y - du_n.y).abs())
                .max((du_p.z - du_n.z).abs());
            let err_v = ((dv_p.x - dv_n.x).abs())
                .max((dv_p.y - dv_n.y).abs())
                .max((dv_p.z - dv_n.z).abs());
            assert!(
                err_u < 1e-4,
                "∂u mismatch at ({u}, {v}): patch {du_p:?}, nurbs {du_n:?}"
            );
            assert!(
                err_v < 1e-4,
                "∂v mismatch at ({u}, {v}): patch {dv_p:?}, nurbs {dv_n:?}"
            );
        }
    }
}
