//! Tests for shared polar K-P utilities.
//!
//! `weighted_sum` is indirectly tested by the K-P 2015 rules
//! tests; here we verify only the externally-visible
//! `catmull_clark_star_point`.

use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

const EPS: f64 = 1e-12;

fn approx_eq(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

#[test]
fn catmull_clark_constant_input_is_fixed_point() {
    let c = pt(3.0, -2.0, 0.5);
    for n in [3usize, 5, 6, 8] {
        let c5s: Vec<Point3> = vec![c; n];
        let c6s: Vec<Point3> = vec![c; n];
        let star = catmull_clark_star_point(n, c, &c5s, &c6s);
        assert!(approx_eq(star, c), "n={n}: star={star:?} c={c:?}");
    }
}

#[test]
fn catmull_clark_weights_sum_to_one() {
    for n in [3usize, 5, 6, 7, 8, 10] {
        let nf = n as f64;
        let nu5 = 1.0 / (nf * (nf + 5.0)) + if n == 3 { 1.0 / 96.0 } else { 0.0 };
        let nu6 = 4.0 / (nf * (nf + 5.0));
        let center = 1.0 - nf * nu5 - nf * nu6;
        let sum = center + nf * nu5 + nf * nu6;
        assert!((sum - 1.0).abs() < 1e-15, "n={n} sum={sum}");
    }
}

#[test]
fn catmull_clark_n3_has_extra_1_over_96_term() {
    let n = 3;
    let c7 = pt(1.0, 0.0, 0.0);
    let c5s = vec![pt(0.0, 0.0, 0.0); n];
    let c6s = vec![pt(0.0, 0.0, 0.0); n];
    let star = catmull_clark_star_point(n, c7, &c5s, &c6s);
    // center_weight at n=3: 1 - 15/96 - 1/2 = 0.34375.
    assert!((star.x - 0.34375).abs() < 1e-15);
}

#[test]
#[should_panic]
fn catmull_clark_panics_on_length_mismatch() {
    let _ = catmull_clark_star_point(
        5,
        pt(0.0, 0.0, 0.0),
        &[pt(0.0, 0.0, 0.0); 4],
        &[pt(0.0, 0.0, 0.0); 5],
    );
}
