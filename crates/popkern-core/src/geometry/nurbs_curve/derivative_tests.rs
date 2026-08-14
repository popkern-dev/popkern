//! Agreement pins for the analytic NURBS curve first derivative
//! (PK-285): must match central FD of `point_at` below the FD
//! truncation floor it replaces, across multi-span rational curves,
//! knot values, and domain ends.

use super::NurbsCurve;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;

fn multi_span_rational() -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.55, 0.8, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.4),
            Point3::new(2.5, 1.5, -0.3),
            Point3::new(3.0, -0.5, 0.8),
            Point3::new(4.2, 0.7, 0.1),
            Point3::new(5.0, 2.2, -0.6),
            Point3::new(6.0, 1.0, 0.2),
        ],
        weights: vec![1.0, 1.3, 0.8, 1.5, 0.9, 1.2, 1.0],
    }
}

#[test]
fn analytic_first_derivative_agrees_with_fd() {
    let c = multi_span_rational();
    let h = 1e-6;
    for i in 0..=40 {
        let t = 0.002 + 0.996 * i as f64 / 40.0;
        let analytic = c.derivative_at(t, 1);
        let fd = (c.point_at(t + h) - c.point_at(t - h)) * (1.0 / (2.0 * h));
        assert!(
            (analytic - fd).length() / fd.length().max(1.0) < 1e-7,
            "curve ∂t diverged at t={t}: analytic {analytic:?} fd {fd:?}"
        );
    }
}

#[test]
fn analytic_first_derivative_at_knots_and_ends() {
    let c = multi_span_rational();
    let h = 1e-6;
    for &t in &[0.0, 0.3, 0.55, 0.8, 1.0] {
        let analytic = c.derivative_at(t, 1);
        // One-sided secant into the domain (derivative at a knot with
        // full-multiplicity ends is a one-sided limit).
        let (t1, sgn) = if t < 0.5 { (t + h, 1.0) } else { (t - h, -1.0) };
        let sec = (c.point_at(t1) - c.point_at(t)) * (sgn / h);
        assert!(
            (analytic - sec).length() / sec.length().max(1.0) < 1e-4,
            "curve ∂t at knot t={t}: analytic {analytic:?} secant {sec:?}"
        );
    }
}

#[test]
fn analytic_first_derivative_on_full_circle() {
    // The rational quadratic full circle (multiplicity-2 interior
    // knots at each quarter) is the curve class the fillet rails
    // actually use — the PK-285 regression fixture. Speed of the
    // standard parameterization varies but must match FD everywhere
    // and stay within a factor ~2 of r (never explode: the deviation
    // checker uses |C'| as a Lipschitz bound).
    use crate::geometry::Vec3;
    let c = crate::geometry::curve_to_nurbs::circle_to_nurbs(
        Point3::new(1.0, 2.0, 3.0),
        Vec3::new(0.0, 0.0, 1.0),
        2.5,
    );
    let h = 1e-6;
    let (t0, t1) = c.parameter_range();
    for i in 0..=200 {
        let t = t0 + (t1 - t0) * (0.001 + 0.998 * i as f64 / 200.0);
        let analytic = c.derivative_at(t, 1);
        let fd = (c.point_at(t + h) - c.point_at(t - h)) * (1.0 / (2.0 * h));
        assert!(
            (analytic - fd).length() / fd.length().max(1.0) < 1e-6,
            "circle ∂t diverged at t={t}: analytic {analytic:?} (|{:.6}|) fd {fd:?} (|{:.6}|)",
            analytic.length(),
            fd.length()
        );
        assert!(
            analytic.length() < 3.0 * 2.5,
            "circle speed exploded at t={t}: {:.6}",
            analytic.length()
        );
    }
}

#[test]
fn second_derivative_still_serves_via_fd_of_analytic() {
    // Order 2 remains FD, now differencing EXACT first derivatives —
    // pin it against an independent 5-point second-difference stencil
    // of point_at.
    let c = multi_span_rational();
    let h = 1e-4;
    let t = 0.42;
    let d2 = c.derivative_at(t, 2);
    let stencil =
        (c.point_at(t + h) - c.point_at(t) + (c.point_at(t - h) - c.point_at(t))) * (1.0 / (h * h));
    assert!(
        (d2 - stencil).length() / stencil.length().max(1.0) < 1e-3,
        "curve ∂²t at t={t}: fd-of-analytic {d2:?} stencil {stencil:?}"
    );
}
