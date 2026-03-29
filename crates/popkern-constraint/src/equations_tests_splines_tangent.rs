//! Tests for spline interior tangency constraint equations.
//!
//! Split from `equations_tests_splines.rs` to stay within the 500-line
//! file limit.

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
        group: None,
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

// ── SplineTangentToLine tests ────────────────────────────────────────────

#[test]
fn spline_tangent_to_line_satisfied() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let la = sk.add_point(-1.0, -2.0);
    let lb = sk.add_point(1.0, 2.0);
    let line = sk.add_line(la, lb);

    let ch = sk.add_constraint(Constraint::SplineTangentToLine(sh, line));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.0;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 2);
    for (i, r) in output.residuals.iter().enumerate() {
        assert!(r.abs() < 1e-10, "F{} = {}", i + 1, r);
    }
}

#[test]
fn spline_tangent_to_line_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let la = sk.add_point(0.5, 0.3);
    let lb = sk.add_point(3.0, 2.7);
    let line = sk.add_line(la, lb);

    let ch = sk.add_constraint(Constraint::SplineTangentToLine(sh, line));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.4;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

// ── SplineTangentToArc tests ─────────────────────────────────────────────

#[test]
fn spline_tangent_to_arc_satisfied() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let arc = sk.add_arc(2.0, 3.0, 1.5, 0.0, std::f64::consts::TAU);

    let ch = sk.add_constraint(Constraint::SplineTangentToArc(sh, arc));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.5;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 2);
    for (i, r) in output.residuals.iter().enumerate() {
        assert!(r.abs() < 1e-10, "F{} = {}", i + 1, r);
    }
}

#[test]
fn spline_tangent_to_arc_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let arc = sk.add_arc(1.5, 2.0, 1.0, 0.0, std::f64::consts::TAU);

    let ch = sk.add_constraint(Constraint::SplineTangentToArc(sh, arc));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.3;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

// ── SplineTangentToEllipse tests ─────────────────────────────────────────

#[test]
fn spline_tangent_to_ellipse_satisfied() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let ell = sk.add_ellipse(2.0, 0.0, 2.0, 1.5, 0.0);

    let ch = sk.add_constraint(Constraint::SplineTangentToEllipse(sh, ell));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.5;
        sk.params[aux_start + 1] = std::f64::consts::FRAC_PI_2;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 3);
    for (i, r) in output.residuals.iter().enumerate() {
        assert!(r.abs() < 1e-10, "F{} = {}", i + 1, r);
    }
}

#[test]
fn spline_tangent_to_ellipse_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let ell = sk.add_ellipse(1.5, 0.5, 2.0, 1.0, 0.3);

    let ch = sk.add_constraint(Constraint::SplineTangentToEllipse(sh, ell));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.35;
        sk.params[aux_start + 1] = 1.2;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

// ── SplineTangentToParabola tests ────────────────────────────────────────

#[test]
fn spline_tangent_to_parabola_satisfied() {
    // Rotated parabola (theta=pi/2) so tangent at t_p=0 is (-1, 0),
    // which is parallel to spline tangent (3, 0) at t_s=0.5.
    // Position: E(0) = (vx, vy) = (2, 1.5) = S(0.5).
    let (mut sk, sh) = cubic_bezier_sketch();
    let para = sk.add_parabola(2.0, 1.5, 1.0, std::f64::consts::FRAC_PI_2);

    let ch = sk.add_constraint(Constraint::SplineTangentToParabola(sh, para));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.5;
        sk.params[aux_start + 1] = 0.0;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 3);
    for (i, r) in output.residuals.iter().enumerate() {
        assert!(r.abs() < 1e-10, "F{} = {}", i + 1, r);
    }
}

#[test]
fn spline_tangent_to_parabola_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let para = sk.add_parabola(1.0, 0.5, 0.8, 0.5);

    let ch = sk.add_constraint(Constraint::SplineTangentToParabola(sh, para));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.4;
        sk.params[aux_start + 1] = 0.7;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

// ── SplineTangentToHyperbola tests ───────────────────────────────────────

#[test]
fn spline_tangent_to_hyperbola_runs() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let hyp = sk.add_hyperbola(1.0, 0.5, 2.0, 1.0, 0.3);

    let ch = sk.add_constraint(Constraint::SplineTangentToHyperbola(sh, hyp));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.3;
        sk.params[aux_start + 1] = 0.5;
    }

    let e = &sk.constraints[ch].clone();
    let output = evaluate_constraint(&sk, Some(ch), e, 0);
    assert_eq!(output.residuals.len(), 3);
}

#[test]
fn spline_tangent_to_hyperbola_jacobian_fd() {
    let (mut sk, sh) = cubic_bezier_sketch();
    let hyp = sk.add_hyperbola(1.0, 0.5, 2.0, 1.0, 0.3);

    let ch = sk.add_constraint(Constraint::SplineTangentToHyperbola(sh, hyp));
    if let Some(&aux_start) = sk.aux_param_map.get(&ch) {
        sk.params[aux_start] = 0.3;
        sk.params[aux_start + 1] = 0.5;
    }

    let e = sk.constraints[ch].clone();
    check_jacobian_fd_with_aux(&mut sk, Some(ch), &e, 0, 1e-7);
}

// ── Equation counts ──────────────────────────────────────────────────────

#[test]
fn tangent_constraint_equation_counts() {
    let s = crate::SplineHandle::default();
    let l = crate::LineHandle::default();
    let a = crate::ArcHandle::default();
    let e = crate::EllipseHandle::default();
    let p = crate::ParabolaHandle::default();
    let h = crate::HyperbolaHandle::default();
    assert_eq!(equation_count(&Constraint::SplineTangentToLine(s, l)), 2);
    assert_eq!(equation_count(&Constraint::SplineTangentToArc(s, a)), 2);
    assert_eq!(equation_count(&Constraint::SplineTangentToEllipse(s, e)), 3);
    assert_eq!(
        equation_count(&Constraint::SplineTangentToParabola(s, p)),
        3
    );
    assert_eq!(
        equation_count(&Constraint::SplineTangentToHyperbola(s, h)),
        3
    );
}
