//! Tests for [`super`] degree-elevation helpers.
//!
//! The defining property of degree elevation is that the elevated
//! polynomial evaluates to the **same curve** at every `t`. All
//! tests pin this identity at f64 precision, plus the endpoint-
//! preservation and constant-curve corollaries.

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-12;

/// Evaluate a Bezier curve of arbitrary degree (`cps.len() − 1`)
/// at parameter `t` via direct Bernstein-basis sum.
fn eval_bezier_vec(cps: &[Vec3], t: f64) -> Vec3 {
    let n = cps.len() - 1;
    let omt = 1.0 - t;
    let mut out = Vec3::zero();
    // Stable binomial via Pascal's-row recurrence.
    let mut binom = vec![0.0_f64; n + 1];
    binom[0] = 1.0;
    for k in 1..=n {
        binom[k] = binom[k - 1] * ((n - k + 1) as f64) / (k as f64);
    }
    for i in 0..=n {
        let w = binom[i] * t.powi(i as i32) * omt.powi((n - i) as i32);
        out.x += cps[i].x * w;
        out.y += cps[i].y * w;
        out.z += cps[i].z * w;
    }
    out
}

fn eval_bezier_point(cps: &[Point3], t: f64) -> Point3 {
    let n = cps.len() - 1;
    let omt = 1.0 - t;
    let mut out = Point3::new(0.0, 0.0, 0.0);
    let mut binom = vec![0.0_f64; n + 1];
    binom[0] = 1.0;
    for k in 1..=n {
        binom[k] = binom[k - 1] * ((n - k + 1) as f64) / (k as f64);
    }
    for i in 0..=n {
        let w = binom[i] * t.powi(i as i32) * omt.powi((n - i) as i32);
        out.x += cps[i].x * w;
        out.y += cps[i].y * w;
        out.z += cps[i].z * w;
    }
    out
}

fn vec3_max(a: Vec3, b: Vec3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

// ---- elevate_cubic_vec_to_nine ----

#[test]
fn elevate_cubic_vec_preserves_endpoints() {
    let cubic = [
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(4.0, 5.0, 6.0),
        Vec3::new(7.0, 8.0, 9.0),
        Vec3::new(10.0, 11.0, 12.0),
    ];
    let elevated = elevate_cubic_vec_to_nine(cubic);
    assert!(vec3_max(elevated[0], cubic[0]) < EPS, "CP 0 preserved");
    assert!(vec3_max(elevated[9], cubic[3]) < EPS, "CP 9 preserved");
}

#[test]
fn elevate_cubic_vec_preserves_curve_pointwise() {
    // A curvy cubic with all components varying.
    let cubic = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 3.0, -1.0),
        Vec3::new(2.0, -2.0, 4.0),
        Vec3::new(3.0, 1.0, 2.0),
    ];
    let elevated = elevate_cubic_vec_to_nine(cubic);
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let cubic_eval = eval_bezier_vec(&cubic, t);
        let elevated_eval = eval_bezier_vec(&elevated, t);
        assert!(
            vec3_max(cubic_eval, elevated_eval) < EPS,
            "t = {t}: cubic {cubic_eval:?} vs elevated {elevated_eval:?}"
        );
    }
}

#[test]
fn elevate_cubic_vec_constant_curve_stays_constant() {
    let c = Vec3::new(7.0, -3.5, 2.25);
    let elevated = elevate_cubic_vec_to_nine([c, c, c, c]);
    for i in 0..10 {
        assert!(
            vec3_max(elevated[i], c) < EPS,
            "CP {i} drifted from constant {c:?}: {:?}",
            elevated[i]
        );
    }
}

#[test]
fn elevate_cubic_vec_zero_curve_is_zero() {
    let zero = Vec3::zero();
    let elevated = elevate_cubic_vec_to_nine([zero; 4]);
    for i in 0..10 {
        assert!(vec3_max(elevated[i], zero) < EPS);
    }
}

#[test]
fn elevate_cubic_vec_linear_curve_stays_linear() {
    // A linear curve expressed as a degenerate cubic: CPs on a
    // straight line with equal spacing. Degree elevation should
    // keep every CP on that line.
    let a = Vec3::new(0.0, 0.0, 0.0);
    let b = Vec3::new(1.0, 2.0, 3.0);
    let cubic = [
        a,
        Vec3::new(b.x / 3.0, b.y / 3.0, b.z / 3.0),
        Vec3::new(2.0 * b.x / 3.0, 2.0 * b.y / 3.0, 2.0 * b.z / 3.0),
        b,
    ];
    let elevated = elevate_cubic_vec_to_nine(cubic);
    for i in 0..10 {
        let t = (i as f64) / 9.0;
        let expected = Vec3::new(t * b.x, t * b.y, t * b.z);
        assert!(
            vec3_max(elevated[i], expected) < EPS,
            "CP {i} at t={t}: expected {expected:?}, got {:?}",
            elevated[i]
        );
    }
}

// ---- elevate_cubic_point_to_nine ----

#[test]
fn elevate_cubic_point_preserves_curve_pointwise() {
    let cubic = [
        Point3::new(0.0, 1.0, -1.0),
        Point3::new(1.5, 2.5, 3.0),
        Point3::new(2.5, -1.0, 1.0),
        Point3::new(4.0, 2.0, 0.5),
    ];
    let elevated = elevate_cubic_point_to_nine(cubic);
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let c_eval = eval_bezier_point(&cubic, t);
        let e_eval = eval_bezier_point(&elevated, t);
        assert!(
            point_max(c_eval, e_eval) < EPS,
            "t = {t}: cubic {c_eval:?} vs elevated {e_eval:?}"
        );
    }
}

#[test]
fn elevate_cubic_point_preserves_endpoints() {
    let a = Point3::new(1.0, 2.0, 3.0);
    let d = Point3::new(10.0, 11.0, 12.0);
    let cubic = [a, Point3::new(4.0, 5.0, 6.0), Point3::new(7.0, 8.0, 9.0), d];
    let elevated = elevate_cubic_point_to_nine(cubic);
    assert!(point_max(elevated[0], a) < EPS);
    assert!(point_max(elevated[9], d) < EPS);
}

// ---- elevate_point_curve_to_nine ----

#[test]
fn elevate_point_curve_degree_7_matches_pointwise() {
    // A degree-7 Bezier curve (8 CPs). Elevate to degree 9 and
    // verify pointwise identity.
    let cps: Vec<Point3> = (0..8)
        .map(|i| {
            let t = i as f64;
            Point3::new(t, (t * 0.5).sin(), (t * 0.3).cos())
        })
        .collect();
    let elevated = elevate_point_curve_to_nine(&cps).expect("should succeed for 8 CPs");
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let orig_eval = eval_bezier_point(&cps, t);
        let elev_eval = eval_bezier_point(&elevated, t);
        assert!(
            point_max(orig_eval, elev_eval) < EPS,
            "t = {t}: orig {orig_eval:?} vs elevated {elev_eval:?}"
        );
    }
    assert!(point_max(elevated[0], cps[0]) < EPS, "endpoint 0 preserved");
    assert!(
        point_max(elevated[9], cps[7]) < EPS,
        "endpoint last preserved"
    );
}

#[test]
fn elevate_point_curve_already_degree_9_is_passthrough() {
    let cps: Vec<Point3> = (0..10)
        .map(|i| Point3::new(i as f64, 2.0 * i as f64, 0.0))
        .collect();
    let elevated = elevate_point_curve_to_nine(&cps).expect("10 CPs is valid");
    for i in 0..10 {
        assert!(
            point_max(elevated[i], cps[i]) < EPS,
            "CP {i}: expected {:?}, got {:?}",
            cps[i],
            elevated[i]
        );
    }
}

#[test]
fn elevate_point_curve_degree_1_line_evaluates_as_line() {
    // A degree-1 curve (2 CPs) is a straight line. Elevated to
    // degree 9, every CP should lie on that line at parameter j/9.
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(9.0, 18.0, -9.0);
    let elevated = elevate_point_curve_to_nine(&[a, b]).expect("2 CPs valid");
    for i in 0..10 {
        let t = (i as f64) / 9.0;
        let expected = Point3::new(t * b.x, t * b.y, t * b.z);
        assert!(
            point_max(elevated[i], expected) < EPS,
            "CP {i} at t={t}: expected {expected:?}, got {:?}",
            elevated[i]
        );
    }
}

#[test]
fn elevate_point_curve_empty_returns_none() {
    assert!(elevate_point_curve_to_nine(&[]).is_none());
}

#[test]
fn elevate_point_curve_too_many_cps_returns_none() {
    let cps: Vec<Point3> = (0..11).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
    assert!(elevate_point_curve_to_nine(&cps).is_none());
}

#[test]
fn elevate_point_curve_degree_2_quadratic_matches_pointwise() {
    // Quadratic (3 CPs) → degree 9.
    let cps = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 2.0, -1.0),
        Point3::new(2.0, 0.0, 0.5),
    ];
    let elevated = elevate_point_curve_to_nine(&cps).expect("3 CPs valid");
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let orig_eval = eval_bezier_point(&cps, t);
        let elev_eval = eval_bezier_point(&elevated, t);
        assert!(point_max(orig_eval, elev_eval) < EPS);
    }
}
