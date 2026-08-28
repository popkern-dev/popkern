//! Tests for arc-to-conic tangency constraint residuals and Jacobians.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::equations::{equation_count, evaluate_constraint};

// ── Helper for arc-tangent tests that need aux params ────────────────────

/// Build a sketch with an arc-to-ellipse tangency constraint.
/// The aux parameter is initialised to `t_init`.
#[allow(clippy::too_many_arguments)]
fn build_arc_ellipse(
    acx: f64,
    acy: f64,
    r: f64,
    cx: f64,
    cy: f64,
    a: f64,
    b: f64,
    theta: f64,
    t_init: f64,
) -> Sketch {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(acx, acy, r, 0.0, std::f64::consts::TAU);
    let ell = sk.add_ellipse(cx, cy, a, b, theta);
    let ch = sk.add_constraint(Constraint::TangentArcEllipse(arc, ell));
    if let Some(&aux_idx) = sk.aux_param_map.get(&ch) {
        sk.params[aux_idx] = t_init;
    }
    sk
}

/// Build a sketch with an arc-to-parabola tangency constraint.
#[allow(clippy::too_many_arguments)]
fn build_arc_parabola(
    acx: f64,
    acy: f64,
    r: f64,
    vx: f64,
    vy: f64,
    f: f64,
    theta: f64,
    t_init: f64,
) -> Sketch {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(acx, acy, r, 0.0, std::f64::consts::TAU);
    let par = sk.add_parabola(vx, vy, f, theta);
    let ch = sk.add_constraint(Constraint::TangentArcParabola(arc, par));
    if let Some(&aux_idx) = sk.aux_param_map.get(&ch) {
        sk.params[aux_idx] = t_init;
    }
    sk
}

/// Build a sketch with an arc-to-hyperbola tangency constraint.
#[allow(clippy::too_many_arguments)]
fn build_arc_hyperbola(
    acx: f64,
    acy: f64,
    r: f64,
    cx: f64,
    cy: f64,
    a: f64,
    b: f64,
    theta: f64,
    t_init: f64,
) -> Sketch {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(acx, acy, r, 0.0, std::f64::consts::TAU);
    let hyp = sk.add_hyperbola(cx, cy, a, b, theta);
    let ch = sk.add_constraint(Constraint::TangentArcHyperbola(arc, hyp));
    if let Some(&aux_idx) = sk.aux_param_map.get(&ch) {
        sk.params[aux_idx] = t_init;
    }
    sk
}

/// Check Jacobian via finite differences for a sketch with one driving
/// constraint (the first one). Uses the constraint handle to look up
/// aux params properly.
fn check_jacobian_fd_with_handle(sk: &mut Sketch, eps: f64) {
    let (handle, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    let entry_clone = entry_ref.clone();
    let base = evaluate_constraint(sk, Some(handle), &entry_clone, 0);
    for param_idx in 0..sk.params.len() {
        let orig = sk.params[param_idx];
        sk.params[param_idx] = orig + eps;
        let perturbed = evaluate_constraint(sk, Some(handle), &entry_clone, 0);
        sk.params[param_idx] = orig;
        for (eq_idx, (base_r, pert_r)) in
            base.residuals.iter().zip(&perturbed.residuals).enumerate()
        {
            let fd = (pert_r - base_r) / eps;
            let analytic = base
                .jacobian
                .iter()
                .find(|(r, c, _)| *r == eq_idx && *c == param_idx)
                .map(|(_, _, v)| *v)
                .unwrap_or(0.0);
            assert!(
                (fd - analytic).abs() < eps * 100.0,
                "Jacobian mismatch for {:?} at eq={eq_idx} param={param_idx}: \
                 fd={fd:.8} analytic={analytic:.8}",
                entry_clone.constraint
            );
        }
    }
}

// ── TangentArcEllipse ───────────────────────────────────────────────────

#[test]
fn tangent_arc_ellipse_equation_count() {
    let sk = build_arc_ellipse(4.0, 0.0, 1.0, 0.0, 0.0, 3.0, 2.0, 0.0, 0.0);
    let (_, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    assert_eq!(equation_count(&entry_ref.constraint), 2);
}

#[test]
fn tangent_arc_ellipse_satisfied() {
    // Ellipse a=3, b=2 at origin. At t=0: E=(3, 0), T=(0, 2).
    // Arc at (4, 0) with r=1: |d|=|(1,0)|=1=r, d.T=1*0+0*2=0.
    let sk = build_arc_ellipse(4.0, 0.0, 1.0, 0.0, 0.0, 3.0, 2.0, 0.0, 0.0);
    let (handle, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    let out = evaluate_constraint(&sk, Some(handle), entry_ref, 0);
    assert!(out.residuals[0].abs() < 1e-12, "F1={}", out.residuals[0]);
    assert!(out.residuals[1].abs() < 1e-12, "F2={}", out.residuals[1]);
}

#[test]
fn tangent_arc_ellipse_jacobian_fd() {
    let mut sk = build_arc_ellipse(5.0, 1.0, 2.0, 0.3, -0.2, 2.5, 1.2, 0.4, 0.7);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

#[test]
fn tangent_arc_ellipse_jacobian_fd_rotated() {
    let mut sk = build_arc_ellipse(3.0, 4.0, 1.5, 1.0, 1.0, 3.0, 1.0, 1.2, -0.5);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

// ── TangentArcParabola ──────────────────────────────────────────────────

#[test]
fn tangent_arc_parabola_equation_count() {
    let sk = build_arc_parabola(3.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    let (_, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    assert_eq!(equation_count(&entry_ref.constraint), 2);
}

#[test]
fn tangent_arc_parabola_satisfied() {
    // Parabola f=1, theta=0. At t=0: E=(0,0), T=(0,1).
    // Arc at (1, 0) with r=1: |d|=|(1,0)|=1=r, d.T=1*0+0*1=0.
    let sk = build_arc_parabola(1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    let (handle, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    let out = evaluate_constraint(&sk, Some(handle), entry_ref, 0);
    assert!(out.residuals[0].abs() < 1e-12, "F1={}", out.residuals[0]);
    assert!(out.residuals[1].abs() < 1e-12, "F2={}", out.residuals[1]);
}

#[test]
fn tangent_arc_parabola_jacobian_fd() {
    let mut sk = build_arc_parabola(5.0, 3.0, 2.0, 0.5, -0.3, 1.5, 0.3, 1.2);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

#[test]
fn tangent_arc_parabola_jacobian_fd_rotated() {
    let mut sk = build_arc_parabola(2.0, 5.0, 1.5, 1.0, 1.0, 2.0, 0.8, -0.9);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

// ── TangentArcHyperbola ─────────────────────────────────────────────────

#[test]
fn tangent_arc_hyperbola_equation_count() {
    let sk = build_arc_hyperbola(4.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
    let (_, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    assert_eq!(equation_count(&entry_ref.constraint), 2);
}

#[test]
fn tangent_arc_hyperbola_satisfied() {
    // Hyperbola a=2, b=1 at origin. At t=0: E=(2,0), T=(0,1).
    // Arc at (3, 0) with r=1: |d|=|(1,0)|=1=r, d.T=1*0+0*1=0.
    let sk = build_arc_hyperbola(3.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0);
    let (handle, entry_ref) = sk.constraints.iter().next().expect("no constraints");
    let out = evaluate_constraint(&sk, Some(handle), entry_ref, 0);
    assert!(out.residuals[0].abs() < 1e-12, "F1={}", out.residuals[0]);
    assert!(out.residuals[1].abs() < 1e-12, "F2={}", out.residuals[1]);
}

#[test]
fn tangent_arc_hyperbola_jacobian_fd() {
    let mut sk = build_arc_hyperbola(5.0, 1.0, 2.0, 0.2, -0.1, 1.2, 0.8, 0.3, 0.5);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

#[test]
fn tangent_arc_hyperbola_jacobian_fd_rotated() {
    let mut sk = build_arc_hyperbola(3.0, 4.0, 1.5, 1.0, 1.0, 2.0, 1.5, 1.0, -0.8);
    check_jacobian_fd_with_handle(&mut sk, 1e-7);
}

// ── Solver convergence tests ────────────────────────────────────────────

#[test]
fn tangent_arc_ellipse_solver_converges() {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(4.5, 0.3, 1.2, 0.0, std::f64::consts::TAU);
    let ell = sk.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.0);
    sk.add_constraint(Constraint::FixedEllipse(ell));
    sk.add_constraint(Constraint::TangentArcEllipse(arc, ell));
    let result = sk.solve().expect("solve failed");
    assert!(
        result.residual_norm < 1e-8,
        "residual_norm={}",
        result.residual_norm
    );
}

#[test]
fn tangent_arc_parabola_solver_converges() {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(1.5, 0.5, 1.0, 0.0, std::f64::consts::TAU);
    let par = sk.add_parabola(0.0, 0.0, 1.0, 0.0);
    sk.add_constraint(Constraint::FixedParabola(par));
    sk.add_constraint(Constraint::TangentArcParabola(arc, par));
    let result = sk.solve().expect("solve failed");
    assert!(
        result.residual_norm < 1e-8,
        "residual_norm={}",
        result.residual_norm
    );
}

#[test]
fn tangent_arc_hyperbola_solver_converges() {
    let mut sk = Sketch::new();
    let arc = sk.add_arc(3.5, 0.3, 1.0, 0.0, std::f64::consts::TAU);
    let hyp = sk.add_hyperbola(0.0, 0.0, 2.0, 1.0, 0.0);
    sk.add_constraint(Constraint::FixedHyperbola(hyp));
    sk.add_constraint(Constraint::TangentArcHyperbola(arc, hyp));
    let result = sk.solve().expect("solve failed");
    assert!(
        result.residual_norm < 1e-8,
        "residual_norm={}",
        result.residual_norm
    );
}
