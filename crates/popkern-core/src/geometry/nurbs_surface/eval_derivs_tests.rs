//! Agreement pins for the fused analytic point+partials evaluation
//! (PK-281): the point must match `point_at` (same surface, different
//! contraction — tolerance-level agreement), and the analytic
//! derivatives must match high-order central differences of
//! `point_at` well below the FD truncation floor they replace.

use super::super::NurbsSurface;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;

fn multi_span_rational() -> NurbsSurface {
    let num_u = 7;
    let num_v = 5;
    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);
    for i in 0..num_u {
        for j in 0..num_v {
            let (fi, fj) = (i as f64, j as f64);
            control_points.push(Point3::new(
                fi + 0.3 * (fj * 1.7).sin(),
                fj + 0.2 * (fi * 2.3).cos(),
                0.5 * ((fi - 3.0) * (fj - 2.0)) + 0.1 * fi,
            ));
            weights.push(1.0 + 0.4 * ((i * 5 + j * 3) % 7) as f64 / 7.0);
        }
    }
    NurbsSurface {
        degree_u: 3,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.4, 0.7, 1.0, 1.0, 1.0],
        control_points,
        num_u,
        num_v,
        weights,
    }
}

fn bezier_quadratic() -> NurbsSurface {
    NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: (0..3)
            .flat_map(|i| {
                (0..3).map(move |j| Point3::new(i as f64, j as f64, ((i * j) as f64 * 0.7).sin()))
            })
            .collect(),
        num_u: 3,
        num_v: 3,
        weights: vec![1.0, 1.2, 1.0, 0.9, 1.5, 0.9, 1.0, 1.2, 1.0],
    }
}

fn check_agreement(s: &NurbsSurface) {
    let n = 17;
    // Central FD of point_at with h chosen for ~1e-10 truncation.
    let h = 1e-6;
    for iu in 0..=n {
        for iv in 0..=n {
            let u = 0.001 + 0.998 * iu as f64 / n as f64;
            let v = 0.001 + 0.998 * iv as f64 / n as f64;
            let (p, su, sv) = s.point_and_first_partials(u, v);

            let p_ref = s.point_at(u, v);
            assert!(
                p.distance_to(&p_ref) < 1e-12,
                "fused point diverged at ({u}, {v}): {} vs {:?}",
                p.distance_to(&p_ref),
                p_ref
            );

            let fd_u = (s.point_at(u + h, v) - s.point_at(u - h, v)) * (1.0 / (2.0 * h));
            let fd_v = (s.point_at(u, v + h) - s.point_at(u, v - h)) * (1.0 / (2.0 * h));
            let scale_u = fd_u.length().max(1.0);
            let scale_v = fd_v.length().max(1.0);
            assert!(
                (su - fd_u).length() / scale_u < 1e-7,
                "analytic ∂u diverged from FD at ({u}, {v}): analytic {su:?} fd {fd_u:?}"
            );
            assert!(
                (sv - fd_v).length() / scale_v < 1e-7,
                "analytic ∂v diverged from FD at ({u}, {v}): analytic {sv:?} fd {fd_v:?}"
            );
        }
    }
}

#[test]
fn fused_eval_agrees_on_multi_span_rational() {
    check_agreement(&multi_span_rational());
}

#[test]
fn fused_eval_agrees_on_rational_bezier() {
    check_agreement(&bezier_quadratic());
}

#[test]
fn fused_eval_agrees_at_domain_ends_and_knots() {
    let s = multi_span_rational();
    // Domain ends and interior knots exercise the span-clamp branches;
    // derivatives there are one-sided limits, so compare against
    // one-sided FD just inside the domain.
    for &(u, v) in &[
        (0.0, 0.0),
        (1.0, 1.0),
        (0.25, 0.4),
        (0.5, 0.7),
        (0.75, 0.0),
        (0.0, 1.0),
    ] {
        let (p, su, sv) = s.point_and_first_partials(u, v);
        assert!(p.distance_to(&s.point_at(u, v)) < 1e-12);
        // Directional secant into the domain.
        let h = 1e-6;
        let (u1, sgn_u) = if u < 0.5 { (u + h, 1.0) } else { (u - h, -1.0) };
        let (v1, sgn_v) = if v < 0.5 { (v + h, 1.0) } else { (v - h, -1.0) };
        let sec_u = (s.point_at(u1, v) - s.point_at(u, v)) * (sgn_u / h);
        let sec_v = (s.point_at(u, v1) - s.point_at(u, v)) * (sgn_v / h);
        assert!(
            (su - sec_u).length() / sec_u.length().max(1.0) < 1e-4,
            "∂u at knot ({u}, {v}): analytic {su:?} secant {sec_u:?}"
        );
        assert!(
            (sv - sec_v).length() / sec_v.length().max(1.0) < 1e-4,
            "∂v at knot ({u}, {v}): analytic {sv:?} secant {sec_v:?}"
        );
    }
}

#[test]
fn fused_eval_degree_one_direction() {
    // Ruled strip: degree 1 in v — exercises the lower-row capture
    // special case (the degree-0 row is trivial).
    let s = NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.3),
            Point3::new(1.0, 0.0, 0.5),
            Point3::new(1.0, 1.0, 0.1),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.4),
        ],
        num_u: 3,
        num_v: 2,
        weights: vec![1.0, 1.1, 1.3, 1.0, 1.0, 0.9],
    };
    let h = 1e-6;
    for &(u, v) in &[(0.3, 0.5), (0.7, 0.2), (0.5, 0.9)] {
        let (p, su, sv) = s.point_and_first_partials(u, v);
        assert!(p.distance_to(&s.point_at(u, v)) < 1e-12);
        let fd_u = (s.point_at(u + h, v) - s.point_at(u - h, v)) * (1.0 / (2.0 * h));
        let fd_v = (s.point_at(u, v + h) - s.point_at(u, v - h)) * (1.0 / (2.0 * h));
        assert!((su - fd_u).length() < 1e-6);
        assert!((sv - fd_v).length() < 1e-6);
    }
}
