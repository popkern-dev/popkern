//! Tests for tier-1 constraint residuals and Jacobian correctness.

#[path = "equations_tests_tier2.rs"]
mod tier2;

#[path = "equations_tests_conics.rs"]
mod conics;

use crate::constraint::{Constraint, ConstraintEntry, ConstraintKind, ConstraintStatus};
use crate::equations::{equation_count, evaluate_constraint};
use crate::Sketch;

/// Check that the analytic Jacobian matches a finite-difference approximation.
pub(super) fn check_jacobian_fd(
    sketch: &mut Sketch,
    entry: &ConstraintEntry,
    row_offset: usize,
    eps: f64,
) {
    let base = evaluate_constraint(sketch, entry, row_offset);
    for param_idx in 0..sketch.params.len() {
        let orig = sketch.params[param_idx];
        sketch.params[param_idx] = orig + eps;
        let perturbed = evaluate_constraint(sketch, entry, row_offset);
        sketch.params[param_idx] = orig;
        for (eq_idx, (base_r, pert_r)) in
            base.residuals.iter().zip(&perturbed.residuals).enumerate()
        {
            let fd = (pert_r - base_r) / eps;
            let analytic = base
                .jacobian
                .iter()
                .find(|(r, c, _)| *r == row_offset + eq_idx && *c == param_idx)
                .map(|(_, _, v)| *v)
                .unwrap_or(0.0);
            assert!(
                (fd - analytic).abs() < eps * 100.0,
                "Jacobian mismatch for {:?} at eq={eq_idx} param={param_idx}: \
                 fd={fd:.8} analytic={analytic:.8}",
                entry.constraint
            );
        }
    }
}

/// Build a ConstraintEntry with default metadata.
pub(super) fn entry(constraint: Constraint, fixed_values: Vec<f64>) -> ConstraintEntry {
    ConstraintEntry {
        constraint,
        kind: ConstraintKind::Driving,
        status: ConstraintStatus::Unknown,
        fixed_values,
    }
}

// ── Coincident ─────────────────────────────────────────────────────────────

#[test]
fn coincident_satisfied() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(3.0, 4.0);
    let p2 = sk.add_point(3.0, 4.0);
    let e = entry(Constraint::Coincident(p1, p2), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert_eq!(equation_count(&e.constraint), 2);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn coincident_violated() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 4.0);
    let e = entry(Constraint::Coincident(p1, p2), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() > 1e-10);
}

#[test]
fn coincident_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.5, 2.7);
    let p2 = sk.add_point(3.1, 4.8);
    let e = entry(Constraint::Coincident(p1, p2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Fixed ──────────────────────────────────────────────────────────────────

#[test]
fn fixed_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(5.0, 6.0);
    let e = entry(Constraint::Fixed(p), vec![5.0, 6.0]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn fixed_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(5.0, 6.0);
    let e = entry(Constraint::Fixed(p), vec![1.0, 2.0]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!((out.residuals[0] - 4.0).abs() < 1e-15);
}

#[test]
fn fixed_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(5.0, 6.0);
    let e = entry(Constraint::Fixed(p), vec![1.0, 2.0]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── FixedLine ──────────────────────────────────────────────────────────────

#[test]
fn fixed_line_satisfied() {
    let mut sk = Sketch::new();
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(3.0, 4.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::FixedLine(l), vec![1.0, 2.0, 3.0, 4.0]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert_eq!(equation_count(&e.constraint), 4);
    assert!(out.residuals.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn fixed_line_jacobian_fd() {
    let mut sk = Sketch::new();
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(3.0, 4.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::FixedLine(l), vec![0.0, 0.0, 0.0, 0.0]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Horizontal ─────────────────────────────────────────────────────────────

#[test]
fn horizontal_satisfied() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 5.0);
    let b = sk.add_point(10.0, 5.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Horizontal(l), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn horizontal_violated() {
    let mut sk = Sketch::new();
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(10.0, 5.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Horizontal(l), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!((out.residuals[0] - 5.0).abs() < 1e-15);
}

#[test]
fn horizontal_jacobian_fd() {
    let mut sk = Sketch::new();
    let a = sk.add_point(1.0, 2.0);
    let b = sk.add_point(3.0, 7.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Horizontal(l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Vertical ───────────────────────────────────────────────────────────────

#[test]
fn vertical_satisfied() {
    let mut sk = Sketch::new();
    let a = sk.add_point(5.0, 0.0);
    let b = sk.add_point(5.0, 10.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::Vertical(l), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

// ── Distance ───────────────────────────────────────────────────────────────

#[test]
fn distance_satisfied() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(3.0, 4.0);
    let e = entry(Constraint::Distance(p1, p2, 5.0), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn distance_violated() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(3.0, 4.0);
    let e = entry(Constraint::Distance(p1, p2, 10.0), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() > 1.0);
}

#[test]
fn distance_jacobian_fd() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(4.0, 6.0);
    let e = entry(Constraint::Distance(p1, p2, 5.0), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── PointOnLine ────────────────────────────────────────────────────────────

#[test]
fn point_on_line_satisfied() {
    let mut sk = Sketch::new();
    let p = sk.add_point(2.0, 2.0);
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(4.0, 4.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::PointOnLine(p, l), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn point_on_line_violated() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.0, 0.0);
    let a = sk.add_point(0.0, 0.0);
    let b = sk.add_point(0.0, 4.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::PointOnLine(p, l), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() > 0.5);
}

#[test]
fn point_on_line_jacobian_fd() {
    let mut sk = Sketch::new();
    let p = sk.add_point(1.5, 3.7);
    let a = sk.add_point(0.5, 1.2);
    let b = sk.add_point(4.0, 6.0);
    let l = sk.add_line(a, b);
    let e = entry(Constraint::PointOnLine(p, l), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Parallel ───────────────────────────────────────────────────────────────

#[test]
fn parallel_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(2.0, 1.0);
    let a2 = sk.add_point(1.0, 3.0);
    let b2 = sk.add_point(3.0, 4.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Parallel(l1, l2), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn parallel_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(2.3, 1.5);
    let a2 = sk.add_point(1.0, 3.0);
    let b2 = sk.add_point(3.7, 4.1);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Parallel(l1, l2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Perpendicular ──────────────────────────────────────────────────────────

#[test]
fn perpendicular_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(1.0, 0.0);
    let a2 = sk.add_point(0.0, 0.0);
    let b2 = sk.add_point(0.0, 1.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Perpendicular(l1, l2), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-15);
}

#[test]
fn perpendicular_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(1.0, 2.0);
    let b1 = sk.add_point(3.0, 5.0);
    let a2 = sk.add_point(0.0, 1.0);
    let b2 = sk.add_point(4.0, 2.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::Perpendicular(l1, l2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── EqualLength ────────────────────────────────────────────────────────────

#[test]
fn equal_length_satisfied() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(3.0, 4.0);
    let a2 = sk.add_point(10.0, 10.0);
    let b2 = sk.add_point(13.0, 14.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::EqualLength(l1, l2), vec![]);
    let out = evaluate_constraint(&sk, &e, 0);
    assert!(out.residuals[0].abs() < 1e-12);
}

#[test]
fn equal_length_jacobian_fd() {
    let mut sk = Sketch::new();
    let a1 = sk.add_point(0.0, 0.0);
    let b1 = sk.add_point(3.0, 4.0);
    let a2 = sk.add_point(1.0, 1.0);
    let b2 = sk.add_point(5.0, 7.0);
    let l1 = sk.add_line(a1, b1);
    let l2 = sk.add_line(a2, b2);
    let e = entry(Constraint::EqualLength(l1, l2), vec![]);
    check_jacobian_fd(&mut sk, &e, 0, 1e-7);
}

// ── Row offset ─────────────────────────────────────────────────────────────

#[test]
fn row_offset_applied_correctly() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(1.0, 2.0);
    let p2 = sk.add_point(3.0, 4.0);
    let e = entry(Constraint::Coincident(p1, p2), vec![]);
    let out = evaluate_constraint(&sk, &e, 10);
    for (r, _, _) in &out.jacobian {
        assert!(*r >= 10);
    }
}
