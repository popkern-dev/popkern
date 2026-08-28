//! Tests for pattern constraint residuals and Jacobians.

use super::{check_jacobian_fd, entry};
use crate::Sketch;
use crate::constraint::Constraint;
use crate::elements::ElementHandle;
use crate::equations::{equation_count, evaluate_constraint};

// ── RelativePosition ──────────────────────────────────────────────────────

#[test]
fn relative_position_satisfied() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(4.0, 6.0);
    let e = entry(Constraint::RelativePosition(p1, p2, 3.0, 4.0), vec![]);
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-15),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn relative_position_violated() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(5.0, 6.0);
    let e = entry(Constraint::RelativePosition(p1, p2, 3.0, 4.0), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!((out.residuals[0] - 1.0).abs() < 1e-15);
    assert!(out.residuals[1].abs() < 1e-15);
}

#[test]
fn relative_position_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.5, 2.7);
    let p2 = sk.add_point(4.3, 6.1);
    let e = entry(Constraint::RelativePosition(p1, p2, 2.0, 3.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── RotatedPosition ───────────────────────────────────────────────────────

#[test]
fn rotated_position_satisfied() {
    let mut sk = Sketch::new();
    // Rotate (1,0) about (0,0) by 90 degrees -> (0,1)
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(0.0, 1.0);
    let center = sk.add_point(0.0, 0.0);
    let e = entry(
        Constraint::RotatedPosition(p1, p2, center, std::f64::consts::FRAC_PI_2),
        vec![],
    );
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn rotated_position_satisfied_offset_center() {
    let mut sk = Sketch::new();
    // Rotate (3,2) about (2,2) by 90 degrees: offset=(1,0), rotated offset=(0,1)
    // -> (2, 3)
    let p1 = sk.add_point(3.0, 2.0);
    let p2 = sk.add_point(2.0, 3.0);
    let center = sk.add_point(2.0, 2.0);
    let e = entry(
        Constraint::RotatedPosition(p1, p2, center, std::f64::consts::FRAC_PI_2),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn rotated_position_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(2.5, 1.3);
    let p2 = sk.add_point(0.7, 3.2);
    let center = sk.add_point(1.0, 1.0);
    let e = entry(Constraint::RotatedPosition(p1, p2, center, 0.8), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── ReflectedPosition ─────────────────────────────────────────────────────

#[test]
fn reflected_position_satisfied_horizontal() {
    let mut sk = Sketch::new();
    // Reflect (1, 3) across y=1 -> (1, -1)
    let p1 = sk.add_point(1.0, 3.0);
    let p2 = sk.add_point(1.0, -1.0);
    let la = sk.add_point(0.0, 1.0);
    let lb = sk.add_point(5.0, 1.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::ReflectedPosition(p1, p2, l), vec![]);
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn reflected_position_satisfied_vertical() {
    let mut sk = Sketch::new();
    // Reflect (3, 2) across x=1 -> (-1, 2)
    let p1 = sk.add_point(3.0, 2.0);
    let p2 = sk.add_point(-1.0, 2.0);
    let la = sk.add_point(1.0, 0.0);
    let lb = sk.add_point(1.0, 5.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::ReflectedPosition(p1, p2, l), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn reflected_position_satisfied_diagonal() {
    let mut sk = Sketch::new();
    // Reflect (2, 0) across y=x -> (0, 2)
    let p1 = sk.add_point(2.0, 0.0);
    let p2 = sk.add_point(0.0, 2.0);
    let la = sk.add_point(0.0, 0.0);
    let lb = sk.add_point(1.0, 1.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::ReflectedPosition(p1, p2, l), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-13),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn reflected_position_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(2.5, 3.7);
    let p2 = sk.add_point(-0.5, 1.2);
    let la = sk.add_point(1.0, 0.5);
    let lb = sk.add_point(4.0, 3.5);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::ReflectedPosition(p1, p2, l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── PinToCenter ───────────────────────────────────────────────────────────

#[test]
fn pin_to_center_arc_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 4.0);
    let a = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PinToCenter(p, ElementHandle::Arc(a)), vec![]);
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn pin_to_center_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    let el = sk.add_ellipse(1.0, 2.0, 3.0, 2.0, 0.0);
    let e = entry(
        Constraint::PinToCenter(p, ElementHandle::Ellipse(el)),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn pin_to_center_parabola_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(5.0, 6.0);
    let pa = sk.add_parabola(5.0, 6.0, 2.0, 0.0);
    let e = entry(
        Constraint::PinToCenter(p, ElementHandle::Parabola(pa)),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn pin_to_center_hyperbola_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(7.0, 8.0);
    let h = sk.add_hyperbola(7.0, 8.0, 3.0, 2.0, 0.0);
    let e = entry(
        Constraint::PinToCenter(p, ElementHandle::Hyperbola(h)),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn pin_to_center_arc_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(0.0, 0.0);
    let a = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PinToCenter(p, ElementHandle::Arc(a)), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!((out.residuals[0] + 3.0).abs() < 1e-15);
    assert!((out.residuals[1] + 4.0).abs() < 1e-15);
}

#[test]
fn pin_to_center_arc_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    let a = sk.add_arc(3.0, 4.0, 5.0, 0.0, std::f64::consts::PI);
    let e = entry(Constraint::PinToCenter(p, ElementHandle::Arc(a)), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

#[test]
fn pin_to_center_ellipse_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    let el = sk.add_ellipse(3.0, 4.0, 5.0, 3.0, 0.5);
    let e = entry(
        Constraint::PinToCenter(p, ElementHandle::Ellipse(el)),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

#[test]
fn pin_to_center_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 2.0);
    let pa = sk.add_parabola(3.0, 4.0, 2.0, 0.5);
    let e = entry(
        Constraint::PinToCenter(p, ElementHandle::Parabola(pa)),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── AngleOffsetArc ────────────────────────────────────────────────────────

#[test]
fn angle_offset_arc_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.3, 1.5);
    let a2 = sk.add_arc(0.0, 0.0, 5.0, 0.8, 2.0);
    let e = entry(Constraint::AngleOffsetArc(a1, a2, 0.5), vec![]);
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn angle_offset_arc_violated() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.3, 1.5);
    let a2 = sk.add_arc(0.0, 0.0, 5.0, 1.3, 2.5);
    let e = entry(Constraint::AngleOffsetArc(a1, a2, 0.5), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() > 0.1);
}

#[test]
fn angle_offset_arc_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 3.0, 0.5, 1.8);
    let a2 = sk.add_arc(4.0, 5.0, 2.0, 1.2, 2.5);
    let e = entry(Constraint::AngleOffsetArc(a1, a2, 0.7), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── AngleOffsetConic ──────────────────────────────────────────────────────

#[test]
fn angle_offset_conic_ellipse_satisfied() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.5);
    let e2 = sk.add_ellipse(0.0, 0.0, 3.0, 2.0, 1.3);
    let e = entry(
        Constraint::AngleOffsetConic(ElementHandle::Ellipse(e1), ElementHandle::Ellipse(e2), 0.8),
        vec![],
    );
    assert_eq!(equation_count(&e.constraint), 1);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals[0].abs() < 1e-14,
        "residual: {}",
        out.residuals[0]
    );
}

#[test]
fn angle_offset_conic_parabola_satisfied() {
    let mut sk = Sketch::new();
    let p1 = sk.add_parabola(0.0, 0.0, 2.0, 0.3);
    let p2 = sk.add_parabola(0.0, 0.0, 2.0, 1.1);
    let e = entry(
        Constraint::AngleOffsetConic(
            ElementHandle::Parabola(p1),
            ElementHandle::Parabola(p2),
            0.8,
        ),
        vec![],
    );
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(out.residuals[0].abs() < 1e-14);
}

#[test]
fn angle_offset_conic_jacobian_fd() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(1.0, 2.0, 3.0, 2.0, 0.5);
    let e2 = sk.add_ellipse(4.0, 5.0, 2.0, 1.0, 1.3);
    let e = entry(
        Constraint::AngleOffsetConic(ElementHandle::Ellipse(e1), ElementHandle::Ellipse(e2), 0.8),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

#[test]
fn angle_offset_conic_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_parabola(1.0, 2.0, 3.0, 0.5);
    let p2 = sk.add_parabola(4.0, 5.0, 2.0, 1.3);
    let e = entry(
        Constraint::AngleOffsetConic(
            ElementHandle::Parabola(p1),
            ElementHandle::Parabola(p2),
            0.8,
        ),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── MirrorAnglesArc ───────────────────────────────────────────────────────

#[test]
fn mirror_angles_arc_satisfied() {
    let mut sk = Sketch::new();
    // Mirror line along x-axis (angle = 0).
    // Arc1 angles: start=0.3, end=1.0
    // Mirrored: arc2.start = 2*0 - 1.0 = -1.0, arc2.end = 2*0 - 0.3 = -0.3
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.3, 1.0);
    let a2 = sk.add_arc(0.0, 0.0, 5.0, -1.0, -0.3);
    let la = sk.add_point(0.0, 0.0);
    let lb = sk.add_point(1.0, 0.0);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::MirrorAnglesArc(a1, a2, l), vec![]);
    assert_eq!(equation_count(&e.constraint), 2);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}",
        out.residuals
    );
}

#[test]
fn mirror_angles_arc_satisfied_45deg() {
    let mut sk = Sketch::new();
    // Mirror line at 45 degrees (angle = pi/4).
    // Arc1: start=0.0, end=pi/2
    // Mirrored: arc2.start = 2*(pi/4) - pi/2 = 0, arc2.end = 2*(pi/4) - 0 = pi/2
    // So for a symmetric case, the constraint is self-consistent.
    let pi4 = std::f64::consts::FRAC_PI_4;
    let pi2 = std::f64::consts::FRAC_PI_2;
    let a1 = sk.add_arc(0.0, 0.0, 5.0, 0.0, pi2);
    let a2 = sk.add_arc(0.0, 0.0, 5.0, 0.0, pi2);
    let la = sk.add_point(0.0, 0.0);
    let lb = sk.add_point(1.0, 1.0); // line angle = pi/4
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::MirrorAnglesArc(a1, a2, l), vec![]);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals.iter().all(|r| r.abs() < 1e-14),
        "residuals: {:?}, pi4={pi4}",
        out.residuals
    );
}

#[test]
fn mirror_angles_arc_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(1.0, 2.0, 3.0, 0.5, 1.8);
    let a2 = sk.add_arc(4.0, 5.0, 2.0, -0.7, 0.3);
    let la = sk.add_point(1.0, 0.5);
    let lb = sk.add_point(4.0, 3.5);
    let l = sk.add_line(la, lb);
    let e = entry(Constraint::MirrorAnglesArc(a1, a2, l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── MirrorAnglesConic ─────────────────────────────────────────────────────

#[test]
fn mirror_angles_conic_ellipse_satisfied() {
    let mut sk = Sketch::new();
    // Mirror line along x-axis (angle=0). rot1=0.5, mirrored rot2 = 2*0 - 0.5 = -0.5
    let e1 = sk.add_ellipse(0.0, 0.0, 3.0, 2.0, 0.5);
    let e2 = sk.add_ellipse(0.0, 0.0, 3.0, 2.0, -0.5);
    let la = sk.add_point(0.0, 0.0);
    let lb = sk.add_point(1.0, 0.0);
    let l = sk.add_line(la, lb);
    let e = entry(
        Constraint::MirrorAnglesConic(ElementHandle::Ellipse(e1), ElementHandle::Ellipse(e2), l),
        vec![],
    );
    assert_eq!(equation_count(&e.constraint), 1);
    let out = evaluate_constraint(&sk, None, &e, 0);
    assert!(
        out.residuals[0].abs() < 1e-14,
        "residual: {}",
        out.residuals[0]
    );
}

#[test]
fn mirror_angles_conic_jacobian_fd() {
    let mut sk = Sketch::new();
    let e1 = sk.add_ellipse(1.0, 2.0, 3.0, 2.0, 0.5);
    let e2 = sk.add_ellipse(4.0, 5.0, 2.0, 1.0, -0.3);
    let la = sk.add_point(1.0, 0.5);
    let lb = sk.add_point(4.0, 3.5);
    let l = sk.add_line(la, lb);
    let e = entry(
        Constraint::MirrorAnglesConic(ElementHandle::Ellipse(e1), ElementHandle::Ellipse(e2), l),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

#[test]
fn mirror_angles_conic_parabola_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_parabola(1.0, 2.0, 3.0, 0.5);
    let p2 = sk.add_parabola(4.0, 5.0, 2.0, -0.3);
    let la = sk.add_point(1.0, 0.5);
    let lb = sk.add_point(4.0, 3.5);
    let l = sk.add_line(la, lb);
    let e = entry(
        Constraint::MirrorAnglesConic(ElementHandle::Parabola(p1), ElementHandle::Parabola(p2), l),
        vec![],
    );
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}
