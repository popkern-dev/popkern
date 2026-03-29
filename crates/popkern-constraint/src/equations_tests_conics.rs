//! Tests for conic constraint residuals and Jacobian correctness.

use super::{check_jacobian_fd, entry};
use crate::constraint::Constraint;
use crate::equations::{equation_count, evaluate_constraint};
use crate::Sketch;

// ── PointOnEllipse ────────────────────────────────────────────────────────

#[test]
fn point_on_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.0, 0.0);
    let e = sk.add_ellipse(0.0, 0.0, 2.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnEllipse(p, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn point_on_ellipse_satisfied_rotated() {
    let mut sk = Sketch::new();
    // Ellipse with a=3, b=1, rotated 90 degrees. Point at (0, 3) is on it.
    let theta = std::f64::consts::FRAC_PI_2;
    let p = sk.add_point(0.0, 3.0);
    let e = sk.add_ellipse(0.0, 0.0, 3.0, 1.0, theta);
    let ce = entry(Constraint::PointOnEllipse(p, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn point_on_ellipse_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 0.0);
    let e = sk.add_ellipse(0.0, 0.0, 2.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnEllipse(p, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert!(out.residuals[0].abs() > 0.5);
}

#[test]
fn point_on_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.5, 0.7);
    let e = sk.add_ellipse(0.3, -0.2, 2.5, 1.2, 0.4);
    let ce = entry(Constraint::PointOnEllipse(p, e), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── PointOnParabola ───────────────────────────────────────────────────────

#[test]
fn point_on_parabola_satisfied() {
    let mut sk = Sketch::new();
    // v^2 = 4*f*u, f=1 => at u=1, v=2: point (1, 2) in local frame
    let p = sk.add_point(1.0, 2.0);
    let pa = sk.add_parabola(0.0, 0.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnParabola(p, pa), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn point_on_parabola_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 1.0);
    let pa = sk.add_parabola(0.0, 0.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnParabola(p, pa), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert!(out.residuals[0].abs() > 0.5);
}

#[test]
fn point_on_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.0, 3.0);
    let pa = sk.add_parabola(0.5, -0.3, 1.5, 0.3);
    let ce = entry(Constraint::PointOnParabola(p, pa), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── PointOnHyperbola ──────────────────────────────────────────────────────

#[test]
fn point_on_hyperbola_satisfied() {
    let mut sk = Sketch::new();
    // u^2/a^2 - v^2/b^2 = 1 with a=2, b=1, at u=2, v=0
    let p = sk.add_point(2.0, 0.0);
    let h = sk.add_hyperbola(0.0, 0.0, 2.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnHyperbola(p, h), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn point_on_hyperbola_satisfied_off_axis() {
    let mut sk = Sketch::new();
    // a=1, b=1: u^2 - v^2 = 1. At u=cosh(1), v=sinh(1)
    let u = 1.0_f64.cosh();
    let v = 1.0_f64.sinh();
    let p = sk.add_point(u, v);
    let h = sk.add_hyperbola(0.0, 0.0, 1.0, 1.0, 0.0);
    let ce = entry(Constraint::PointOnHyperbola(p, h), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn point_on_hyperbola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 1.5);
    let h = sk.add_hyperbola(0.2, -0.1, 2.0, 1.5, 0.5);
    let ce = entry(Constraint::PointOnHyperbola(p, h), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── FixedEllipse ──────────────────────────────────────────────────────────

#[test]
fn fixed_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let e = sk.add_ellipse(1.0, 2.0, 3.0, 1.5, 0.5);
    let ce = entry(Constraint::FixedEllipse(e), vec![1.0, 2.0, 3.0, 1.5, 0.5]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 5);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn fixed_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    let e = sk.add_ellipse(1.0, 2.0, 3.0, 1.5, 0.5);
    let ce = entry(Constraint::FixedEllipse(e), vec![0.0, 0.0, 1.0, 1.0, 0.0]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── FixedParabola ─────────────────────────────────────────────────────────

#[test]
fn fixed_parabola_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_parabola(1.0, 2.0, 3.0, 0.5);
    let ce = entry(Constraint::FixedParabola(p), vec![1.0, 2.0, 3.0, 0.5]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 4);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn fixed_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_parabola(1.0, 2.0, 3.0, 0.5);
    let ce = entry(Constraint::FixedParabola(p), vec![0.0, 0.0, 1.0, 0.0]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── FixedHyperbola ────────────────────────────────────────────────────────

#[test]
fn fixed_hyperbola_satisfied() {
    let mut sk = Sketch::new();
    let h = sk.add_hyperbola(1.0, 2.0, 3.0, 1.5, 0.5);
    let ce = entry(Constraint::FixedHyperbola(h), vec![1.0, 2.0, 3.0, 1.5, 0.5]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 5);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn fixed_hyperbola_jacobian_fd() {
    let mut sk = Sketch::new();
    let h = sk.add_hyperbola(1.0, 2.0, 3.0, 1.5, 0.5);
    let ce = entry(Constraint::FixedHyperbola(h), vec![0.0, 0.0, 1.0, 1.0, 0.0]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── ConcentricEllipse ─────────────────────────────────────────────────────

#[test]
fn concentric_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(3.0, 4.0, 5.0, 2.0, 0.0);
    let e2 = sk.add_ellipse(3.0, 4.0, 2.0, 1.0, 1.0);
    let ce = entry(Constraint::ConcentricEllipse(e1, e2), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 2);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn concentric_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(1.0, 2.0, 3.0, 1.5, 0.3);
    let e2 = sk.add_ellipse(4.0, 5.0, 2.0, 1.0, 0.7);
    let ce = entry(Constraint::ConcentricEllipse(e1, e2), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── ConcentricArcEllipse ──────────────────────────────────────────────────

#[test]
fn concentric_arc_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let a = sk.add_arc(3.0, 4.0, 2.0, 0.0, std::f64::consts::PI);
    let e = sk.add_ellipse(3.0, 4.0, 5.0, 2.0, 0.0);
    let ce = entry(Constraint::ConcentricArcEllipse(a, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 2);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn concentric_arc_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    let a = sk.add_arc(1.0, 2.0, 3.0, 0.0, std::f64::consts::PI);
    let e = sk.add_ellipse(4.0, 5.0, 2.0, 1.0, 0.7);
    let ce = entry(Constraint::ConcentricArcEllipse(a, e), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── EqualSemiAxis ─────────────────────────────────────────────────────────

#[test]
fn equal_semi_axis_satisfied() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(0.0, 0.0, 5.0, 2.0, 0.0);
    let e2 = sk.add_ellipse(10.0, 10.0, 5.0, 3.0, 1.0);
    let ce = entry(Constraint::EqualSemiAxis(e1, e2), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn equal_semi_axis_jacobian_fd() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(0.0, 0.0, 5.0, 2.0, 0.0);
    let e2 = sk.add_ellipse(10.0, 10.0, 3.0, 1.0, 1.0);
    let ce = entry(Constraint::EqualSemiAxis(e1, e2), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── FocalDistance ──────────────────────────────────────────────────────────

#[test]
fn focal_distance_satisfied() {
    let mut sk = Sketch::new();
    // a=5, b=3, focal distance c = sqrt(25-9) = 4
    let e = sk.add_ellipse(0.0, 0.0, 5.0, 3.0, 0.0);
    let ce = entry(Constraint::FocalDistance(e, 4.0), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn focal_distance_jacobian_fd() {
    let mut sk = Sketch::new();
    let e = sk.add_ellipse(1.0, 2.0, 5.0, 3.0, 0.5);
    let ce = entry(Constraint::FocalDistance(e, 3.0), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── EllipseEccentricity ──────────────────────────────────────────────────

#[test]
fn ellipse_eccentricity_satisfied() {
    let mut sk = Sketch::new();
    // a=5, b=3, ecc = sqrt(1 - 9/25) = sqrt(16/25) = 0.8
    let e = sk.add_ellipse(0.0, 0.0, 5.0, 3.0, 0.0);
    let ce = entry(Constraint::EllipseEccentricity(e, 0.8), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn ellipse_eccentricity_jacobian_fd() {
    let mut sk = Sketch::new();
    let e = sk.add_ellipse(1.0, 2.0, 5.0, 3.0, 0.5);
    let ce = entry(Constraint::EllipseEccentricity(e, 0.6), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── HyperbolaEccentricity ────────────────────────────────────────────────

#[test]
fn hyperbola_eccentricity_satisfied() {
    let mut sk = Sketch::new();
    // a=3, b=4, ecc = sqrt(1 + 16/9) = sqrt(25/9) = 5/3
    let h = sk.add_hyperbola(0.0, 0.0, 3.0, 4.0, 0.0);
    let ce = entry(Constraint::HyperbolaEccentricity(h, 5.0 / 3.0), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-12, "res={}", out.residuals[0]);
}

#[test]
fn hyperbola_eccentricity_jacobian_fd() {
    let mut sk = Sketch::new();
    let h = sk.add_hyperbola(1.0, 2.0, 3.0, 4.0, 0.5);
    let ce = entry(Constraint::HyperbolaEccentricity(h, 1.5), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── TangentLineEllipse ───────────────────────────────────────────────────

#[test]
fn tangent_line_ellipse_satisfied() {
    let mut sk = Sketch::new();
    // Horizontal line y=1 tangent to ellipse a=2, b=1 centered at origin
    let p1 = sk.add_point(-5.0, 1.0);
    let p2 = sk.add_point(5.0, 1.0);
    let l = sk.add_line(p1, p2);
    let e = sk.add_ellipse(0.0, 0.0, 2.0, 1.0, 0.0);
    let ce = entry(Constraint::TangentLineEllipse(l, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-10, "res={}", out.residuals[0]);
}

#[test]
fn tangent_line_ellipse_violated() {
    let mut sk = Sketch::new();
    // Line y=0.5 is NOT tangent to unit circle (a=1, b=1)
    let p1 = sk.add_point(-5.0, 0.5);
    let p2 = sk.add_point(5.0, 0.5);
    let l = sk.add_line(p1, p2);
    let e = sk.add_ellipse(0.0, 0.0, 1.0, 1.0, 0.0);
    let ce = entry(Constraint::TangentLineEllipse(l, e), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert!(out.residuals[0].abs() > 0.1);
}

#[test]
fn tangent_line_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    // Use moderate coordinates to keep derivative magnitudes within FD
    // tolerance (the dual-conic residual is polynomial, so large coords
    // produce large second derivatives).
    let p1 = sk.add_point(-1.0, 1.2);
    let p2 = sk.add_point(1.5, 1.5);
    let l = sk.add_line(p1, p2);
    let e = sk.add_ellipse(0.2, -0.1, 1.5, 0.8, 0.3);
    let ce = entry(Constraint::TangentLineEllipse(l, e), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── TangentLineParabola ──────────────────────────────────────────────────

#[test]
fn tangent_line_parabola_satisfied() {
    let mut sk = Sketch::new();
    // Parabola v^2 = 4u (f=1), tangent line at u=1, v=2: slope dv/du = 2/(2v) = 1/2
    // In parametric: v = 2t, u = t^2. At t=1: u=1, v=2.
    // Tangent direction: (du/dt, dv/dt) = (2t, 2) = (2, 2) at t=1.
    // Line through (1,2) with direction (2,2): point + t*(2,2)
    // So another point: (3, 4)
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 4.0);
    let l = sk.add_line(p1, p2);
    let pa = sk.add_parabola(0.0, 0.0, 1.0, 0.0);
    let ce = entry(Constraint::TangentLineParabola(l, pa), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-10, "res={}", out.residuals[0]);
}

#[test]
fn tangent_line_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(-2.0, 1.0);
    let p2 = sk.add_point(3.0, 4.0);
    let l = sk.add_line(p1, p2);
    let pa = sk.add_parabola(0.5, -0.3, 1.5, 0.3);
    let ce = entry(Constraint::TangentLineParabola(l, pa), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}

// ── TangentLineHyperbola ─────────────────────────────────────────────────

#[test]
fn tangent_line_hyperbola_satisfied() {
    let mut sk = Sketch::new();
    // Hyperbola u^2 - v^2 = 1 (a=1, b=1). At point (cosh(t), sinh(t)).
    // Take t=1: u=cosh(1), v=sinh(1).
    // Tangent direction: (sinh(1), cosh(1)).
    let u0 = 1.0_f64.cosh();
    let v0 = 1.0_f64.sinh();
    let du = 1.0_f64.sinh();
    let dv = 1.0_f64.cosh();
    let p1 = sk.add_point(u0, v0);
    let p2 = sk.add_point(u0 + du, v0 + dv);
    let l = sk.add_line(p1, p2);
    let h = sk.add_hyperbola(0.0, 0.0, 1.0, 1.0, 0.0);
    let ce = entry(Constraint::TangentLineHyperbola(l, h), vec![]);
    let out = evaluate_constraint(&sk, None, &ce, 0);
    assert_eq!(equation_count(&ce.constraint), 1);
    assert!(out.residuals[0].abs() < 1e-10, "res={}", out.residuals[0]);
}

#[test]
fn tangent_line_hyperbola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(-1.0, 0.8);
    let p2 = sk.add_point(1.5, 0.5);
    let l = sk.add_line(p1, p2);
    let h = sk.add_hyperbola(0.2, -0.1, 1.2, 0.8, 0.3);
    let ce = entry(Constraint::TangentLineHyperbola(l, h), vec![]);
    check_jacobian_fd(&mut sk, &ce, 0, 1e-7);
}
