//! Tests for arc-related and tier-2 constraint residuals and Jacobians.

use super::{check_jacobian_fd, entry};
use crate::constraint::Constraint;
use crate::equations::{equation_count, evaluate_constraint};
use crate::Sketch;

// ── PointOnArc ─────────────────────────────────────────────────────────────

#[test]
fn point_on_arc_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(5.0, 0.0);
    let a = sk.add_arc(0.0, 0.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PointOnArc(p, a), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn point_on_arc_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 0.0);
    let a = sk.add_arc(0.0, 0.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PointOnArc(p, a), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 1.0);
}

#[test]
fn point_on_arc_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 4.0);
    let a = sk.add_arc(1.0, 1.0, 3.5, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PointOnArc(p, a), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── PointLineDistance ──────────────────────────────────────────────────────

#[test]
fn point_line_distance_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 3.0);
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(10.0, 0.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::PointLineDistance(p, l, 3.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn point_line_distance_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.0, 5.0);
    let a = sk.add_point(1.0, 1.0);
    let b = sk.add_point(6.0, 3.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::PointLineDistance(p, l, 2.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Angle ──────────────────────────────────────────────────────────────────

#[test]
fn angle_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(1.0, 0.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(0.0, 1.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(
        Constraint::Angle(l1, l2, std::f64::consts::FRAC_PI_2),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-14);
}

#[test]
fn angle_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(1.0, 2.0);
    let b1 = sk.add_point(4.0, 3.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(-1.0, 2.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Angle(l1, l2, 1.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Tangent ────────────────────────────────────────────────────────────────

#[test]
fn tangent_satisfied() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 5.0);
    let b = sk.add_point(10.0, 5.0);
    let l = sk.add_line(a, b);
    let arc = sk.add_arc(3.0, 0.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Tangent(l, arc), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn tangent_jacobian_fd() {
    let mut sk = Sketch::new();
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(5.0, 4.0);
    let l = sk.add_line(a, b);
    let arc = sk.add_arc(3.0, 1.0, 2.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Tangent(l, arc), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── ArcTangent ─────────────────────────────────────────────────────────────

#[test]
fn arc_tangent_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 3.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(5.0, 0.0, 2.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::ArcTangent(a1, a2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn arc_tangent_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 3.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(7.0, 5.0, 2.5, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::ArcTangent(a1, a2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── EqualRadius ────────────────────────────────────────────────────────────

#[test]
fn equal_radius_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(10.0, 10.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::EqualRadius(a1, a2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn equal_radius_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(10.0, 10.0, 3.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::EqualRadius(a1, a2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Symmetric ──────────────────────────────────────────────────────────────

#[test]
fn symmetric_satisfied() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 3.0);
    let p2 = sk.add_point(1.0, -1.0);
    let la = sk.add_point(0.0, 1.0);
    let lb = sk.add_point(5.0, 1.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::Symmetric(p1, p2, l), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(equation_count(&e.constraint), 2);
    assert!(
        out.residuals[0].abs() < 1e-14,
        "midpoint-on-line residual: {}",
        out.residuals[0]
    );
    assert!(
        out.residuals[1].abs() < 1e-14,
        "perpendicular residual: {}",
        out.residuals[1]
    );
}

#[test]
fn symmetric_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(2.0, 5.0);
    let p2 = sk.add_point(4.0, -1.0);
    let la = sk.add_point(1.0, 1.0);
    let lb = sk.add_point(6.0, 3.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::Symmetric(p1, p2, l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Midpoint ───────────────────────────────────────────────────────────────

#[test]
fn midpoint_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.5, 3.5);
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(4.0, 5.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Midpoint(p, l), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-14));
}

#[test]
fn midpoint_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 4.0);
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(6.0, 7.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Midpoint(p, l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Concentric ─────────────────────────────────────────────────────────────

#[test]
fn concentric_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(3.0, 4.0, 2.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Concentric(a1, a2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn concentric_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(1.0, 2.0, 2.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Concentric(a1, a2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Collinear ──────────────────────────────────────────────────────────────

#[test]
fn collinear_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(2.0, 2.0);
    let a2 = sk.add_point(3.0, 3.0);
    let b2 = sk.add_point(5.0, 5.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Collinear(l1, l2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-14);
}

#[test]
fn collinear_violated() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(2.0, 2.0);
    let a2 = sk.add_point(3.0, 4.0);
    let b2 = sk.add_point(5.0, 5.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Collinear(l1, l2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 0.5);
}

#[test]
fn collinear_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(1.0, 2.0);
    let b1 = sk.add_point(4.0, 6.0);
    let a2 = sk.add_point(3.0, 5.0);
    let b2 = sk.add_point(7.0, 8.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Collinear(l1, l2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Coradial ───────────────────────────────────────────────────────────────

#[test]
fn coradial_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Coradial(a1, a2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert_eq!(equation_count(&e.constraint), 3);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn coradial_violated() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 5.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(3.0, 4.0, 3.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Coradial(a1, a2), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 0.5, "cx diff expected");
    assert!(out.residuals[1].abs() > 0.5, "cy diff expected");
    assert!(out.residuals[2].abs() > 0.5, "r diff expected");
}

#[test]
fn coradial_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 4.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(3.0, 5.0, 2.5, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::Coradial(a1, a2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── OffsetLine ─────────────────────────────────────────────────────────────

#[test]
fn offset_line_satisfied() {
    // L1 at y=3, L2 at y=0, horizontal: distance = 3.
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 3.0);
    let b1 = sk.add_point(4.0, 3.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(4.0, 0.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::OffsetLine(l1, l2, 3.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn offset_line_violated() {
    // L1 at y=5, L2 at y=0, distance is 5 but we ask for 3.
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 5.0);
    let b1 = sk.add_point(4.0, 5.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(4.0, 0.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::OffsetLine(l1, l2, 3.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 1.0);
}

#[test]
fn offset_line_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(1.0, 4.0);
    let b1 = sk.add_point(5.0, 4.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(6.0, 0.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::OffsetLine(l1, l2, 2.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── OffsetArc ──────────────────────────────────────────────────────────────

#[test]
fn offset_arc_satisfied() {
    // r1=7, r2=4 → diff=3.
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 7.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(0.0, 0.0, 4.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::OffsetArc(a1, a2, 3.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn offset_arc_violated() {
    // r1=7, r2=4 → diff=3, but we ask for d=1.
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 7.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(0.0, 0.0, 4.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::OffsetArc(a1, a2, 1.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 1.0);
}

#[test]
fn offset_arc_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 6.0, 0.0, std::f64::consts::PI);
    let a2 = sk.add_arc(3.0, 4.0, 3.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::OffsetArc(a1, a2, 2.5), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── OrdinateX ──────────────────────────────────────────────────────────────

#[test]
fn ordinate_x_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(4.0, 7.0);
    let e = entry(Constraint::OrdinateX(p, 4.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn ordinate_x_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(4.0, 7.0);
    let e = entry(Constraint::OrdinateX(p, 10.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!((out.residuals[0] + 6.0).abs() < 1e-15);
}

#[test]
fn ordinate_x_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 5.0);
    let e = entry(Constraint::OrdinateX(p, 7.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── OrdinateY ──────────────────────────────────────────────────────────────

#[test]
fn ordinate_y_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 8.0);
    let e = entry(Constraint::OrdinateY(p, 8.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn ordinate_y_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 8.0);
    let e = entry(Constraint::OrdinateY(p, 2.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!((out.residuals[0] - 6.0).abs() < 1e-15);
}

#[test]
fn ordinate_y_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.0, 9.0);
    let e = entry(Constraint::OrdinateY(p, 5.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}
