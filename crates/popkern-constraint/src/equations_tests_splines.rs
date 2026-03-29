//! Tests for spline endpoint constraint equations.

use crate::constraint::{Constraint, ConstraintEntry, ConstraintKind, ConstraintStatus};
use crate::equations::{equation_count, evaluate_constraint};
use crate::Sketch;

/// Build a ConstraintEntry with default metadata.
fn entry(constraint: Constraint, fixed_values: Vec<f64>) -> ConstraintEntry {
    ConstraintEntry {
        constraint,
        kind: ConstraintKind::Driving,
        status: ConstraintStatus::Unknown,
        fixed_values,
    }
}

/// Create a sketch with a cubic Bezier spline P0=(0,0), P1=(1,2), P2=(3,2), P3=(4,0).
fn cubic_bezier_sketch() -> (Sketch, crate::SplineHandle) {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 2.0);
    let p3 = sk.add_point(4.0, 0.0);
    let sh = sk.add_spline(&[p0, p1, p2, p3], 3);
    (sk, sh)
}

/// Finite-difference Jacobian check for constraints that use aux params.
/// `aux_indices` lists the parameter indices of auxiliary variables.
fn check_jacobian_fd_with_aux(
    sketch: &mut Sketch,
    handle: Option<crate::ConstraintHandle>,
    entry: &ConstraintEntry,
    row_offset: usize,
    eps: f64,
) {
    let base = evaluate_constraint(sketch, handle, entry, row_offset);
    let n_params = sketch.params.len();

    for param_idx in 0..n_params {
        let orig = sketch.params[param_idx];
        sketch.params[param_idx] = orig + eps;
        let perturbed = evaluate_constraint(sketch, handle, entry, row_offset);
        sketch.params[param_idx] = orig;

        for (eq_idx, (base_r, pert_r)) in
            base.residuals.iter().zip(&perturbed.residuals).enumerate()
        {
            let fd = (pert_r - base_r) / eps;
            let analytic = base
                .jacobian
                .iter()
                .filter(|(r, c, _)| *r == row_offset + eq_idx && *c == param_idx)
                .map(|(_, _, v)| v)
                .sum::<f64>();
            let tol = eps * 200.0;
            assert!(
                (fd - analytic).abs() < tol,
                "Jacobian mismatch for {:?} at eq={eq_idx} param={param_idx}: \
                 fd={fd:.8} analytic={analytic:.8} (diff={:.2e})",
                entry.constraint,
                (fd - analytic).abs()
            );
        }
    }
}

// ── 1. PointOnSpline tests ────────────────────────────────────────────────

#[test]
fn point_on_spline_satisfied() {
    let (mut sk, sh) = cubic_bezier_sketch();
    // Place a point exactly at the spline midpoint: S(0.5) = (2, 1.5)
    let pt = sk.add_point(2.0, 1.5);
    let ch = sk.add_constraint(Constraint::PointOnSpline(pt, sh));
    // Set aux param t = 0.5
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.5;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 2);
    assert!(
        output.residuals[0].abs() < 1e-10,
        "F1 = {}",
        output.residuals[0]
    );
    assert!(
        output.residuals[1].abs() < 1e-10,
        "F2 = {}",
        output.residuals[1]
    );
}

#[test]
fn point_on_spline_violated() {
    let (mut sk, sh) = cubic_bezier_sketch();
    // Place a point NOT on the spline
    let pt = sk.add_point(2.0, 3.0);
    let ch = sk.add_constraint(Constraint::PointOnSpline(pt, sh));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.5;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    // S(0.5) = (2, 1.5), point = (2, 3) => F2 = 1.5 - 3 = -1.5
    assert!(output.residuals[0].abs() < 1e-10); // x matches
    assert!((output.residuals[1] - (-1.5)).abs() < 1e-10);
}

#[test]
fn point_on_spline_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let pt = sk.add_point(2.5, 1.0);
    let ch = sk.add_constraint(Constraint::PointOnSpline(pt, sh));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.4;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

#[test]
fn point_on_spline_equation_count() {
    assert_eq!(
        equation_count(&Constraint::PointOnSpline(
            crate::PointHandle::default(),
            crate::SplineHandle::default()
        )),
        2
    );
}

// ── 2. SplineEndTangent tests ─────────────────────────────────────────────

#[test]
fn spline_end_tangent_satisfied() {
    let (sk, sh) = cubic_bezier_sketch();
    // Tangent at start: 3*(P1-P0) = 3*(1,2) = (3,6), direction (1,2).
    // Create a line with direction (1,2): from (0,0) to (1,2).
    let l_start = sk.points.keys().next().unwrap_or_default(); // p0
    let mut sk = sk;
    let la = sk.add_point(0.0, 0.0);
    let lb = sk.add_point(1.0, 2.0);
    let line = sk.add_line(la, lb);

    let c = Constraint::SplineEndTangent(sh, 0, line);
    let e = entry(c, vec![]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 1);
    // F = 3 * 2 - 6 * 1 = 0
    assert!(
        output.residuals[0].abs() < 1e-10,
        "residual = {}",
        output.residuals[0]
    );

    let _ = l_start; // suppress warning
}

#[test]
fn spline_end_tangent_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let la = sk.add_point(0.5, 0.3);
    let lb = sk.add_point(2.0, 1.7);
    let line = sk.add_line(la, lb);

    let c = Constraint::SplineEndTangent(sh, 0, line);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

#[test]
fn spline_end_tangent_at_end_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let la = sk.add_point(1.0, 0.0);
    let lb = sk.add_point(3.0, -1.5);
    let line = sk.add_line(la, lb);

    let c = Constraint::SplineEndTangent(sh, 1, line);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── 3. SplineEndTangentToArc tests ────────────────────────────────────────

#[test]
fn spline_end_tangent_to_arc_satisfied() {
    // Spline start tangent is (3, 6), direction (1, 2).
    // For tangent to arc: S'_x * (py - cy) - S'_y * (px - cx) = 0
    // where (px, py) = S(0) = (0, 0) and (cx, cy) is arc center.
    // 3 * (0 - cy) - 6 * (0 - cx) = -3*cy + 6*cx = 0 => cy = 2*cx.
    // Put arc centre at (1, 2).
    let (mut sk, sh) = cubic_bezier_sketch();
    let arc = sk.add_arc(1.0, 2.0, 1.0, 0.0, std::f64::consts::PI);

    let c = Constraint::SplineEndTangentToArc(sh, 0, arc);
    let e = entry(c, vec![]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 1);
    assert!(
        output.residuals[0].abs() < 1e-10,
        "residual = {}",
        output.residuals[0]
    );
}

#[test]
fn spline_end_tangent_to_arc_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let arc = sk.add_arc(0.5, 1.5, 2.0, 0.0, std::f64::consts::PI);

    let c = Constraint::SplineEndTangentToArc(sh, 0, arc);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

#[test]
fn spline_end_tangent_to_arc_end_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let arc = sk.add_arc(3.0, -1.0, 1.5, 0.0, std::f64::consts::PI);

    let c = Constraint::SplineEndTangentToArc(sh, 1, arc);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── 4. SplineFixedEndTangent tests ────────────────────────────────────────

#[test]
fn spline_fixed_end_tangent_satisfied() {
    let (sk, sh) = cubic_bezier_sketch();
    // Tangent at start = (3, 6).
    let c = Constraint::SplineFixedEndTangent(sh, 0, 3.0, 6.0);
    let e = entry(c, vec![3.0, 6.0]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 2);
    assert!(
        output.residuals[0].abs() < 1e-10,
        "F1 = {}",
        output.residuals[0]
    );
    assert!(
        output.residuals[1].abs() < 1e-10,
        "F2 = {}",
        output.residuals[1]
    );
}

#[test]
fn spline_fixed_end_tangent_violated() {
    let (sk, sh) = cubic_bezier_sketch();
    // Tangent at start = (3, 6), but we fix to (1, 0).
    let c = Constraint::SplineFixedEndTangent(sh, 0, 1.0, 0.0);
    let e = entry(c, vec![1.0, 0.0]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert!((output.residuals[0] - 2.0).abs() < 1e-10); // 3 - 1 = 2
    assert!((output.residuals[1] - 6.0).abs() < 1e-10); // 6 - 0 = 6
}

#[test]
fn spline_fixed_end_tangent_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let c = Constraint::SplineFixedEndTangent(sh, 0, 1.0, 2.0);
    let e = entry(c, vec![1.0, 2.0]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

#[test]
fn spline_fixed_end_tangent_end_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let c = Constraint::SplineFixedEndTangent(sh, 1, 2.0, -3.0);
    let e = entry(c, vec![2.0, -3.0]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── 5. SplineFixedEndCurvature tests ──────────────────────────────────────

#[test]
fn spline_fixed_end_curvature_satisfied() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    // Compute the actual curvature at start.
    let (kappa, _) = crate::equations::spline_deriv::eval_curvature(sp, 0.0, &sk);

    let c = Constraint::SplineFixedEndCurvature(sh, 0, kappa);
    let e = entry(c, vec![]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 1);
    assert!(
        output.residuals[0].abs() < 1e-10,
        "residual = {}",
        output.residuals[0]
    );
}

#[test]
fn spline_fixed_end_curvature_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let c = Constraint::SplineFixedEndCurvature(sh, 0, 0.5);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

#[test]
fn spline_fixed_end_curvature_end_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let c = Constraint::SplineFixedEndCurvature(sh, 1, -0.3);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── Through-point solve integration ───────────────────────────────────────

#[test]
fn through_point_spline_solve() {
    let mut sk = Sketch::new();
    let pts = vec![(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 0.0)];
    let (sh, tp_handles) = sk.add_spline_through_points(&pts).unwrap();

    // Should solve without error.
    let result = sk.solve().unwrap();
    let norm = result.residual_norm;
    assert!(norm < 1e-8, "residual norm {norm} too large after solve");

    // Verify through-points are still on the spline.
    for (i, &tph) in tp_handles.iter().enumerate() {
        let (px, py) = sk.point_position(tph).unwrap();
        assert!(
            (px - pts[i].0).abs() < 1e-6 && (py - pts[i].1).abs() < 1e-6,
            "through-point {} moved: ({px}, {py}) vs ({}, {})",
            i,
            pts[i].0,
            pts[i].1
        );
    }

    let _ = sh; // suppress warning
}

// ── 6. SplineG1 tests ────────────────────────────────────────────────────

#[test]
fn spline_g1_satisfied() {
    // Two Bezier splines that share an endpoint and have matching tangent.
    // S1: P0=(0,0), P1=(1,2), P2=(3,2), P3=(4,0).  End tangent = (3, -6).
    // S2: P0=(4,0), P1=(5,-2), P2=(7,-2), P3=(8,0). Start tangent = 3*(5-4,-2-0) = (3, -6).
    // They share endpoint (4,0) and tangent direction (3,-6).
    let mut sk = Sketch::new();
    let p = [
        sk.add_point(0.0, 0.0),
        sk.add_point(1.0, 2.0),
        sk.add_point(3.0, 2.0),
        sk.add_point(4.0, 0.0),
    ];
    let sh1 = sk.add_spline(&p, 3);

    let q = [
        sk.add_point(4.0, 0.0),
        sk.add_point(5.0, -2.0),
        sk.add_point(7.0, -2.0),
        sk.add_point(8.0, 0.0),
    ];
    let sh2 = sk.add_spline(&q, 3);

    let c = Constraint::SplineG1(sh1, 1, sh2, 0);
    let e = entry(c, vec![]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 3);
    for (i, r) in output.residuals.iter().enumerate() {
        assert!(r.abs() < 1e-10, "G1 F{} = {}", i + 1, r);
    }
}

#[test]
fn spline_g1_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = [
        sk.add_point(0.0, 0.0),
        sk.add_point(1.0, 2.0),
        sk.add_point(3.0, 2.5),
        sk.add_point(4.0, 0.5),
    ];
    let sh1 = sk.add_spline(&p, 3);
    let q = [
        sk.add_point(4.5, 0.0),
        sk.add_point(5.0, -1.5),
        sk.add_point(7.0, -2.0),
        sk.add_point(8.0, 0.0),
    ];
    let sh2 = sk.add_spline(&q, 3);

    let c = Constraint::SplineG1(sh1, 1, sh2, 0);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── 7. SplineG2 tests ────────────────────────────────────────────────────

#[test]
fn spline_g2_equation_count() {
    let s = crate::SplineHandle::default();
    assert_eq!(equation_count(&Constraint::SplineG2(s, 0, s, 1)), 4);
}

#[test]
fn spline_g2_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = [
        sk.add_point(0.0, 0.0),
        sk.add_point(1.0, 2.0),
        sk.add_point(3.0, 2.5),
        sk.add_point(4.0, 0.5),
    ];
    let sh1 = sk.add_spline(&p, 3);
    let q = [
        sk.add_point(4.5, 0.0),
        sk.add_point(5.0, -1.5),
        sk.add_point(7.0, -2.0),
        sk.add_point(8.0, 0.0),
    ];
    let sh2 = sk.add_spline(&q, 3);

    let c = Constraint::SplineG2(sh1, 1, sh2, 0);
    let e = entry(c, vec![]);
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-7);
}

// ── 8. SplineLength tests ────────────────────────────────────────────────

#[test]
fn spline_length_satisfied() {
    let (sk, sh) = cubic_bezier_sketch();
    let sp = &sk.splines[sh];
    // Compute the actual arc length.
    let (actual_len, _) = crate::equations::spline_deriv::eval_arc_length(sp, 0.0, 1.0, &sk);

    let c = Constraint::SplineLength(sh, actual_len);
    let e = entry(c, vec![]);
    let output = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(output.residuals.len(), 1);
    assert!(
        output.residuals[0].abs() < 1e-10,
        "residual = {}",
        output.residuals[0]
    );
}

#[test]
fn spline_length_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let c = Constraint::SplineLength(sh, 5.0);
    let e = entry(c, vec![]);
    // Use wider tolerance since arc length is approximate (quadrature).
    check_jacobian_fd_with_aux(&mut sk, None, &e, 0, 1e-6);
}

#[test]
fn spline_g1_g2_length_equation_counts() {
    let s = crate::SplineHandle::default();
    assert_eq!(equation_count(&Constraint::SplineG1(s, 0, s, 1)), 3);
    assert_eq!(equation_count(&Constraint::SplineG2(s, 0, s, 1)), 4);
    assert_eq!(equation_count(&Constraint::SplineLength(s, 5.0)), 1);
}
