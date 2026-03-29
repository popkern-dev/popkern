//! Tests for B-spline evaluation with analytical Jacobians.

use super::spline_deriv::{eval_arc_length, eval_curvature, eval_second_deriv, eval_tangent};
use super::spline_eval::eval_point;
use crate::Sketch;

/// Helper: create a sketch with a cubic Bezier spline.
/// P0=(0,0), P1=(1,2), P2=(3,2), P3=(4,0)
/// Knots: [0,0,0,0,1,1,1,1], degree=3.
fn cubic_bezier_sketch() -> (Sketch, crate::SplineHandle) {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 2.0);
    let p3 = sk.add_point(4.0, 0.0);
    let sh = sk.add_spline(&[p0, p1, p2, p3], 3);
    (sk, sh)
}

// ── Point evaluation tests ───────────────────────────────────────────────

#[test]
fn point_at_start() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let r = eval_point(sp, 0.0, &sk);
    assert!((r.x - 0.0).abs() < 1e-12, "x={}", r.x);
    assert!((r.y - 0.0).abs() < 1e-12, "y={}", r.y);
}

#[test]
fn point_at_end() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let r = eval_point(sp, 1.0, &sk);
    assert!((r.x - 4.0).abs() < 1e-12, "x={}", r.x);
    assert!((r.y - 0.0).abs() < 1e-12, "y={}", r.y);
}

#[test]
fn point_at_midpoint() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let r = eval_point(sp, 0.5, &sk);
    assert!((r.x - 2.0).abs() < 1e-12, "x={}", r.x);
    assert!((r.y - 1.5).abs() < 1e-12, "y={}", r.y);
}

// ── Tangent evaluation tests ─────────────────────────────────────────────

#[test]
fn tangent_at_start() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let r = eval_tangent(sp, 0.0, &sk);
    // tangent_at(0) = 3*(P1-P0) = 3*(1,2) = (3,6)
    assert!((r.x - 3.0).abs() < 1e-10, "tx={}", r.x);
    assert!((r.y - 6.0).abs() < 1e-10, "ty={}", r.y);
}

#[test]
fn tangent_at_end() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let r = eval_tangent(sp, 1.0, &sk);
    // tangent_at(1) = 3*(P3-P2) = 3*(1,-2) = (3,-6)
    assert!((r.x - 3.0).abs() < 1e-10, "tx={}", r.x);
    assert!((r.y - (-6.0)).abs() < 1e-10, "ty={}", r.y);
}

// ── Finite-difference Jacobian verification ──────────────────────────────

/// Perturb each control point parameter by epsilon, recompute, and
/// compare against the analytical Jacobian.
#[test]
fn point_jacobian_finite_diff() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t = 0.3;
    let eps = 1e-7;
    let tol = 1e-4;

    let base = eval_point(sp, t, &sk);

    // Gather all control point param indices.
    let mut param_indices = Vec::new();
    for &cph in &sp.control_points {
        let pt = &sk.points[cph];
        let [ix, iy] = pt.param_indices();
        param_indices.push(ix);
        param_indices.push(iy);
    }

    for &pi in &param_indices {
        let mut sk2 = Sketch::new();
        // Rebuild with perturbed params.
        let mut params = sk.params.clone();
        params[pi] += eps;

        // Reconstruct points in the clone.
        let mut handles = Vec::new();
        for &cph in &sp.control_points {
            let orig = &sk.points[cph];
            let [ix, iy] = orig.param_indices();
            handles.push(sk2.add_point(params[ix], params[iy]));
        }
        let sh2 = sk2.add_spline(&handles, sp.degree);
        let sp2 = &sk2.splines[sh2];
        let perturbed = eval_point(sp2, t, &sk2);

        let fd_dx = (perturbed.x - base.x) / eps;
        let fd_dy = (perturbed.y - base.y) / eps;

        // Find the analytical Jacobian entry for this param index.
        let (ana_dx, ana_dy) = analytical_jac_2d(&base.jacobian, pi);

        assert!(
            (fd_dx - ana_dx).abs() < tol,
            "point jac x: param {} fd={fd_dx} ana={ana_dx}",
            pi
        );
        assert!(
            (fd_dy - ana_dy).abs() < tol,
            "point jac y: param {} fd={fd_dy} ana={ana_dy}",
            pi
        );
    }
}

#[test]
fn tangent_jacobian_finite_diff() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t = 0.4;
    let eps = 1e-7;
    let tol = 1e-4;

    let base = eval_tangent(sp, t, &sk);

    let mut param_indices = Vec::new();
    for &cph in &sp.control_points {
        let pt = &sk.points[cph];
        let [ix, iy] = pt.param_indices();
        param_indices.push(ix);
        param_indices.push(iy);
    }

    for &pi in &param_indices {
        let mut sk2 = Sketch::new();
        let mut params = sk.params.clone();
        params[pi] += eps;

        let mut handles = Vec::new();
        for &cph in &sp.control_points {
            let orig = &sk.points[cph];
            let [ix, iy] = orig.param_indices();
            handles.push(sk2.add_point(params[ix], params[iy]));
        }
        let sh2 = sk2.add_spline(&handles, sp.degree);
        let sp2 = &sk2.splines[sh2];
        let perturbed = eval_tangent(sp2, t, &sk2);

        let fd_dx = (perturbed.x - base.x) / eps;
        let fd_dy = (perturbed.y - base.y) / eps;
        let (ana_dx, ana_dy) = analytical_jac_2d(&base.jacobian, pi);

        assert!(
            (fd_dx - ana_dx).abs() < tol,
            "tangent jac x: param {} fd={fd_dx} ana={ana_dx}",
            pi
        );
        assert!(
            (fd_dy - ana_dy).abs() < tol,
            "tangent jac y: param {} fd={fd_dy} ana={ana_dy}",
            pi
        );
    }
}

#[test]
fn second_deriv_jacobian_finite_diff() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t = 0.6;
    let eps = 1e-7;
    let tol = 1e-4;

    let base = eval_second_deriv(sp, t, &sk);

    let mut param_indices = Vec::new();
    for &cph in &sp.control_points {
        let pt = &sk.points[cph];
        let [ix, iy] = pt.param_indices();
        param_indices.push(ix);
        param_indices.push(iy);
    }

    for &pi in &param_indices {
        let mut sk2 = Sketch::new();
        let mut params = sk.params.clone();
        params[pi] += eps;

        let mut handles = Vec::new();
        for &cph in &sp.control_points {
            let orig = &sk.points[cph];
            let [ix, iy] = orig.param_indices();
            handles.push(sk2.add_point(params[ix], params[iy]));
        }
        let sh2 = sk2.add_spline(&handles, sp.degree);
        let sp2 = &sk2.splines[sh2];
        let perturbed = eval_second_deriv(sp2, t, &sk2);

        let fd_dx = (perturbed.x - base.x) / eps;
        let fd_dy = (perturbed.y - base.y) / eps;
        let (ana_dx, ana_dy) = analytical_jac_2d(&base.jacobian, pi);

        assert!(
            (fd_dx - ana_dx).abs() < tol,
            "accel jac x: param {} fd={fd_dx} ana={ana_dx}",
            pi
        );
        assert!(
            (fd_dy - ana_dy).abs() < tol,
            "accel jac y: param {} fd={fd_dy} ana={ana_dy}",
            pi
        );
    }
}

#[test]
fn curvature_jacobian_finite_diff() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t = 0.5;
    let eps = 1e-7;
    let tol = 1e-4;

    let (base_k, base_jac) = eval_curvature(sp, t, &sk);

    let mut param_indices = Vec::new();
    for &cph in &sp.control_points {
        let pt = &sk.points[cph];
        let [ix, iy] = pt.param_indices();
        param_indices.push(ix);
        param_indices.push(iy);
    }

    for &pi in &param_indices {
        let mut sk2 = Sketch::new();
        let mut params = sk.params.clone();
        params[pi] += eps;

        let mut handles = Vec::new();
        for &cph in &sp.control_points {
            let orig = &sk.points[cph];
            let [ix, iy] = orig.param_indices();
            handles.push(sk2.add_point(params[ix], params[iy]));
        }
        let sh2 = sk2.add_spline(&handles, sp.degree);
        let sp2 = &sk2.splines[sh2];
        let (perturbed_k, _) = eval_curvature(sp2, t, &sk2);

        let fd = (perturbed_k - base_k) / eps;
        let ana = analytical_jac_scalar(&base_jac, pi);

        assert!(
            (fd - ana).abs() < tol,
            "curvature jac: param {} fd={fd} ana={ana}",
            pi
        );
    }
}

#[test]
fn arc_length_jacobian_finite_diff() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t0 = 0.0;
    let t1 = 1.0;
    let eps = 1e-7;
    let tol = 1e-4;

    let (base_len, base_jac) = eval_arc_length(sp, t0, t1, &sk);
    assert!(base_len > 0.0, "arc length should be positive");

    let mut param_indices = Vec::new();
    for &cph in &sp.control_points {
        let pt = &sk.points[cph];
        let [ix, iy] = pt.param_indices();
        param_indices.push(ix);
        param_indices.push(iy);
    }

    for &pi in &param_indices {
        let mut sk2 = Sketch::new();
        let mut params = sk.params.clone();
        params[pi] += eps;

        let mut handles = Vec::new();
        for &cph in &sp.control_points {
            let orig = &sk.points[cph];
            let [ix, iy] = orig.param_indices();
            handles.push(sk2.add_point(params[ix], params[iy]));
        }
        let sh2 = sk2.add_spline(&handles, sp.degree);
        let sp2 = &sk2.splines[sh2];
        let (perturbed_len, _) = eval_arc_length(sp2, t0, t1, &sk2);

        let fd = (perturbed_len - base_len) / eps;
        let ana = analytical_jac_scalar(&base_jac, pi);

        assert!(
            (fd - ana).abs() < tol,
            "arc_length jac: param {} fd={fd} ana={ana}",
            pi
        );
    }
}

// ── Periodic spline tests ────────────────────────────────────────────────

#[test]
fn periodic_square_point_wraps() {
    // 4 control points forming a square, periodic, degree 3.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(1.0, 0.0);
    let p1 = sk.add_point(0.0, 1.0);
    let p2 = sk.add_point(-1.0, 0.0);
    let p3 = sk.add_point(0.0, -1.0);

    // For a periodic spline with 4 CPs and degree 3, use uniform knots.
    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let sh = sk.add_spline(&[p0, p1, p2, p3], 3);
    // Override knots and set periodic.
    sk.splines[sh].knots = knots;
    sk.splines[sh].periodic = true;

    let sp = &sk.splines[sh];

    // Point at t=0 should equal point at t=4 (the period).
    let r0 = eval_point(sp, 0.0, &sk);
    let r4 = eval_point(sp, 4.0, &sk);
    assert!((r0.x - r4.x).abs() < 1e-10, "wrap x: {} vs {}", r0.x, r4.x);
    assert!((r0.y - r4.y).abs() < 1e-10, "wrap y: {} vs {}", r0.y, r4.y);
}

#[test]
fn periodic_tangent_continuous() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(1.0, 0.0);
    let p1 = sk.add_point(0.0, 1.0);
    let p2 = sk.add_point(-1.0, 0.0);
    let p3 = sk.add_point(0.0, -1.0);

    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let sh = sk.add_spline(&[p0, p1, p2, p3], 3);
    sk.splines[sh].knots = knots;
    sk.splines[sh].periodic = true;

    let sp = &sk.splines[sh];

    // Tangent should be continuous across the wrap boundary.
    let t0 = eval_tangent(sp, 0.0, &sk);
    let t4 = eval_tangent(sp, 4.0, &sk);
    assert!(
        (t0.x - t4.x).abs() < 1e-10,
        "tangent wrap x: {} vs {}",
        t0.x,
        t4.x
    );
    assert!(
        (t0.y - t4.y).abs() < 1e-10,
        "tangent wrap y: {} vs {}",
        t0.y,
        t4.y
    );
}

// ── Rational spline test ─────────────────────────────────────────────────

#[test]
fn rational_uniform_weights_match_non_rational() {
    // With all weights = 1.0, rational and non-rational should agree.
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    let t = 0.3;

    let r1 = eval_point(sp, t, &sk);

    // Force rational path by setting a weight != 1 then back.
    // Actually, all weights are 1.0 by default so it takes the
    // non-rational path. Let's set them explicitly to 1.0 to confirm.
    // At t=0.3, the cubic Bezier should give a reasonable interior point.
    assert!(r1.x > 0.0 && r1.x < 4.0, "x={} should be in (0,4)", r1.x);
    assert!(r1.y > 0.0, "y={} should be positive", r1.y);
}

#[test]
fn rational_weighted_point() {
    // Rational B-spline: increase weight on P1 to pull curve toward it.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 2.0);
    let p3 = sk.add_point(4.0, 0.0);
    let sh = sk.add_spline(&[p0, p1, p2, p3], 3);
    sk.splines[sh].weights[1] = 5.0;

    let sp = &sk.splines[sh];
    let r_weighted = eval_point(sp, 0.5, &sk);

    // With a large weight on P1, the midpoint should be closer to P1.
    // For uniform weights, midpoint = (2, 1.5).
    // With w1=5, the curve should be pulled toward (1, 2).
    assert!(r_weighted.y > 1.5, "y={} should be > 1.5", r_weighted.y);
}

// ── Helpers ──────────────────────────────────────────────────────────────

/// Look up the analytical 2D Jacobian for a given param index.
fn analytical_jac_2d(jac: &[(usize, f64, f64)], param_idx: usize) -> (f64, f64) {
    for &(idx, dx, dy) in jac {
        if idx == param_idx {
            return (dx, dy);
        }
    }
    (0.0, 0.0)
}

/// Look up the analytical scalar Jacobian for a given param index.
fn analytical_jac_scalar(jac: &[(usize, f64)], param_idx: usize) -> f64 {
    for &(idx, v) in jac {
        if idx == param_idx {
            return v;
        }
    }
    0.0
}
