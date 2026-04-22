//! Tests for the Taylor jet + chain-rule composition.

use super::super::kp13_bi5_chi::{BilinearChi, ChiCorner};
use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- Basic Taylor jet ----

#[test]
fn zero_jet_has_zero_norm() {
    assert_eq!(TaylorJet2::zero().sq_norm(), 0.0);
}

#[test]
fn jet_sub_self_is_zero() {
    let j = TaylorJet2 {
        value: Point3::new(1.0, 2.0, 3.0),
        d_u: Point3::new(0.5, 0.0, 0.0),
        d_v: Point3::new(0.0, 0.5, 0.0),
        d_uu: Point3::new(1.0, 0.0, 0.0),
        d_uv: Point3::new(0.0, 0.0, 1.0),
        d_vv: Point3::new(0.0, 1.0, 0.0),
    };
    let diff = j.sub(&j);
    assert_eq!(diff.sq_norm(), 0.0);
}

// ---- Chain rule: identity inner map ----

#[test]
fn identity_chi_passes_outer_jet_through() {
    // χ(u, v) = (u, v): x = u, y = v. Derivatives: x_u = 1, y_v = 1,
    // everything else zero. Chain rule should give outer's partials
    // directly.
    let chi = ChiJet {
        pos: (0.25, 0.75),
        d_u: (1.0, 0.0),
        d_v: (0.0, 1.0),
        d_uu: (0.0, 0.0),
        d_uv: (0.0, 0.0),
        d_vv: (0.0, 0.0),
    };
    let outer = OuterHessian {
        value: Point3::new(1.0, 2.0, 3.0),
        d_x: Point3::new(4.0, 5.0, 6.0),
        d_y: Point3::new(7.0, 8.0, 9.0),
        d_xx: Point3::new(10.0, 11.0, 12.0),
        d_xy: Point3::new(13.0, 14.0, 15.0),
        d_yy: Point3::new(16.0, 17.0, 18.0),
    };
    let j = compose_chain_rule(&outer, &chi);
    assert!(approx_p(j.value, outer.value));
    assert!(approx_p(j.d_u, outer.d_x));
    assert!(approx_p(j.d_v, outer.d_y));
    assert!(approx_p(j.d_uu, outer.d_xx));
    assert!(approx_p(j.d_uv, outer.d_xy));
    assert!(approx_p(j.d_vv, outer.d_yy));
}

// ---- Chain rule: scaled inner map ----

#[test]
fn scaled_chi_scales_partials_correctly() {
    // χ(u, v) = (a·u, b·v): x_u = a, y_v = b, else 0.
    // Chain rule:
    //   ∂u = f_x · a
    //   ∂v = f_y · b
    //   ∂uu = f_xx · a²
    //   ∂uv = f_xy · a · b
    //   ∂vv = f_yy · b²
    let a = 2.0;
    let b = 3.0;
    let chi = ChiJet {
        pos: (0.0, 0.0),
        d_u: (a, 0.0),
        d_v: (0.0, b),
        d_uu: (0.0, 0.0),
        d_uv: (0.0, 0.0),
        d_vv: (0.0, 0.0),
    };
    let outer = OuterHessian {
        value: Point3::new(0.0, 0.0, 0.0),
        d_x: Point3::new(1.0, 0.0, 0.0),
        d_y: Point3::new(0.0, 1.0, 0.0),
        d_xx: Point3::new(0.0, 0.0, 1.0),
        d_xy: Point3::new(1.0, 0.0, 0.0),
        d_yy: Point3::new(0.0, 1.0, 0.0),
    };
    let j = compose_chain_rule(&outer, &chi);
    assert!(approx_p(j.d_u, Point3::new(a, 0.0, 0.0)));
    assert!(approx_p(j.d_v, Point3::new(0.0, b, 0.0)));
    assert!(approx_p(j.d_uu, Point3::new(0.0, 0.0, a * a)));
    assert!(approx_p(j.d_uv, Point3::new(a * b, 0.0, 0.0)));
    assert!(approx_p(j.d_vv, Point3::new(0.0, b * b, 0.0)));
}

// ---- outer_hessian_at_natural ----

#[test]
fn outer_hessian_at_natural_inside_simplex() {
    // Construct a sector where every CP equals a constant P. Then
    // g△ is the constant function P on the whole domain, and its
    // OuterHessian at any interior point should be P with zero
    // partials.
    let p = Point3::new(1.5, -2.0, 3.25);
    let s = super::super::kp13_bi5_gtri5::Gtri5Sector { cps: [p; 21] };
    let oh = outer_hessian_at_natural(&s, 0.3, 0.4);
    assert!(approx_p(oh.value, p));
    assert!(approx_p(oh.d_x, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(oh.d_y, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(oh.d_xx, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(oh.d_xy, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(oh.d_yy, Point3::new(0.0, 0.0, 0.0)));
}

#[test]
fn outer_hessian_at_natural_outside_simplex() {
    // At (u, v) = (1, 1) (w = −1, outside natural simplex): constant
    // CPs still give constant value, zero partials.
    let p = Point3::new(7.0, -3.0, 1.0);
    let s = super::super::kp13_bi5_gtri5::Gtri5Sector { cps: [p; 21] };
    let oh = outer_hessian_at_natural(&s, 1.0, 1.0);
    assert!(approx_p(oh.value, p));
    assert!(approx_p(oh.d_x, Point3::new(0.0, 0.0, 0.0)));
    assert!(approx_p(oh.d_y, Point3::new(0.0, 0.0, 0.0)));
}

#[test]
fn outer_hessian_linear_function() {
    // Construct CPs corresponding to the linear function
    // f(u, v) = 2u + 3v + 1. Then f_x = 2, f_y = 3, higher partials = 0.
    let mut s = super::super::kp13_bi5_gtri5::Gtri5Sector::zero();
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            s.set_cp(i, j, Point3::new(2.0 * u + 3.0 * v + 1.0, 0.0, 0.0));
        }
    }
    let oh = outer_hessian_at_natural(&s, 0.3, 0.4);
    // value at (0.3, 0.4) = 2*0.3 + 3*0.4 + 1 = 2.8
    assert!((oh.value.x - 2.8).abs() < 1e-10);
    assert!((oh.d_x.x - 2.0).abs() < 1e-10);
    assert!((oh.d_y.x - 3.0).abs() < 1e-10);
    assert!(oh.d_xx.x.abs() < 1e-10);
    assert!(oh.d_xy.x.abs() < 1e-10);
    assert!(oh.d_yy.x.abs() < 1e-10);
}

// ---- Chain rule against a known composite, by FD check ----

/// Bilinear χ at n=6, sector 0, evaluate chain rule at corner
/// `Star` (= (0, 0)). Compare to a direct FD of the composite.
#[test]
fn chain_rule_matches_fd_on_bilinear_chi_quadratic_outer() {
    // Outer: f(x, y) = (x², x·y + y, 3x + 2y + 1) in Point3.
    // f_x  = (2x, y, 3)
    // f_y  = (0, x + 1, 2)
    // f_xx = (2, 0, 0)
    // f_xy = (0, 1, 0)
    // f_yy = (0, 0, 0)
    let eval_f = |x: f64, y: f64| Point3::new(x * x, x * y + y, 3.0 * x + 2.0 * y + 1.0);

    let n = 6;
    let k = 0;
    let chi = BilinearChi::new(n);
    let corner = ChiCorner::Star;
    let (u0, v0) = corner.uv();
    let chi_jet = chi.jet(k, corner);
    let (x0, y0) = chi_jet.pos;
    let outer = OuterHessian {
        value: eval_f(x0, y0),
        d_x: Point3::new(2.0 * x0, y0, 3.0),
        d_y: Point3::new(0.0, x0 + 1.0, 2.0),
        d_xx: Point3::new(2.0, 0.0, 0.0),
        d_xy: Point3::new(0.0, 1.0, 0.0),
        d_yy: Point3::new(0.0, 0.0, 0.0),
    };
    let j = compose_chain_rule(&outer, &chi_jet);

    // FD check: sample the composite at a few perturbations.
    let h = 1e-5;
    let comp = |u: f64, v: f64| {
        let (x, y) = chi.eval(k, u, v);
        eval_f(x, y)
    };
    // ∂u
    let fd_u_x = (comp(u0 + h, v0).x - comp(u0 - h, v0).x) / (2.0 * h);
    assert!(
        (fd_u_x - j.d_u.x).abs() < 1e-5,
        "fd_u.x {} vs analytic {}",
        fd_u_x,
        j.d_u.x
    );
    // ∂v
    let fd_v_y = (comp(u0, v0 + h).y - comp(u0, v0 - h).y) / (2.0 * h);
    assert!((fd_v_y - j.d_v.y).abs() < 1e-5);
    // ∂uu
    let fd_uu_x = (comp(u0 + h, v0).x - 2.0 * comp(u0, v0).x + comp(u0 - h, v0).x) / (h * h);
    assert!(
        (fd_uu_x - j.d_uu.x).abs() < 1e-3,
        "fd_uu.x {} vs analytic {}",
        fd_uu_x,
        j.d_uu.x
    );
}
