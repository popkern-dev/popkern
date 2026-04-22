//! Tests for bi-3 ↔ bi-5 degree raise + bi-5 subdivision helpers.

#![allow(clippy::needless_range_loop)]

use super::super::bernstein_bi5::{bernstein_bi3_all, bernstein_bi5_all};
use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

fn eval_bi3_1d(cps: &[Point3; 4], t: f64) -> Point3 {
    let b = [
        (1.0 - t).powi(3),
        3.0 * t * (1.0 - t).powi(2),
        3.0 * t * t * (1.0 - t),
        t * t * t,
    ];
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        out.x += b[i] * cps[i].x;
        out.y += b[i] * cps[i].y;
        out.z += b[i] * cps[i].z;
    }
    out
}

fn eval_bi5_1d(cps: &[Point3; 6], t: f64) -> Point3 {
    let b = bernstein_bi5_all(t);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        out.x += b[i] * cps[i].x;
        out.y += b[i] * cps[i].y;
        out.z += b[i] * cps[i].z;
    }
    out
}

fn eval_bi3_grid(cps: &[[Point3; 4]; 4], u: f64, v: f64) -> Point3 {
    let bu = bernstein_bi3_all(u);
    let bv = bernstein_bi3_all(v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        for j in 0..4 {
            let w = bu[i] * bv[j];
            out.x += w * cps[i][j].x;
            out.y += w * cps[i][j].y;
            out.z += w * cps[i][j].z;
        }
    }
    out
}

fn eval_bi5_grid(cps: &[[Point3; 6]; 6], u: f64, v: f64) -> Point3 {
    let bu = bernstein_bi5_all(u);
    let bv = bernstein_bi5_all(v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        for j in 0..6 {
            let w = bu[i] * bv[j];
            out.x += w * cps[i][j].x;
            out.y += w * cps[i][j].y;
            out.z += w * cps[i][j].z;
        }
    }
    out
}

// ---- 1-D degree raise ----

#[test]
fn degree_raise_1d_preserves_polynomial_3_to_4() {
    let cps: [Point3; 4] = [
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(4.0, -1.0, 0.5),
        Point3::new(-2.0, 3.0, 1.0),
        Point3::new(0.5, 0.5, -2.0),
    ];
    let raised = degree_raise_1d_3_to_4(&cps);
    for t in [0.0_f64, 0.25, 0.5, 0.75, 1.0] {
        // Evaluate bi-4 polynomial at t.
        let b4 = [
            (1.0 - t).powi(4),
            4.0 * t * (1.0 - t).powi(3),
            6.0 * t * t * (1.0 - t).powi(2),
            4.0 * t * t * t * (1.0 - t),
            t.powi(4),
        ];
        let mut v4 = Point3::new(0.0, 0.0, 0.0);
        for i in 0..5 {
            v4.x += b4[i] * raised[i].x;
            v4.y += b4[i] * raised[i].y;
            v4.z += b4[i] * raised[i].z;
        }
        let v3 = eval_bi3_1d(&cps, t);
        assert!(approx_p(v3, v4, 1e-10), "t={t}: {v3:?} vs {v4:?}");
    }
}

#[test]
fn degree_raise_1d_preserves_polynomial_3_to_5() {
    let cps: [Point3; 4] = [
        Point3::new(0.1, 0.2, 0.3),
        Point3::new(0.4, 0.5, 0.6),
        Point3::new(0.7, 0.8, 0.9),
        Point3::new(1.0, 1.1, 1.2),
    ];
    let raised = degree_raise_1d_3_to_5(&cps);
    for t in [0.0_f64, 0.17, 0.33, 0.5, 0.77, 1.0] {
        let v3 = eval_bi3_1d(&cps, t);
        let v5 = eval_bi5_1d(&raised, t);
        assert!(approx_p(v3, v5, 1e-10), "t={t}: {v3:?} vs {v5:?}");
    }
}

#[test]
fn degree_raise_preserves_endpoints() {
    let cps: [Point3; 4] = [
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(4.0, 5.0, 6.0),
        Point3::new(7.0, 8.0, 9.0),
        Point3::new(10.0, 11.0, 12.0),
    ];
    let raised = degree_raise_1d_3_to_5(&cps);
    assert_eq!(raised[0], cps[0]);
    assert_eq!(raised[5], cps[3]);
}

// ---- 2-D degree raise ----

#[test]
fn bi3_to_bi5_grid_preserves_polynomial() {
    let mut bi3 = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            let t = (i + j) as f64 * 0.1;
            bi3[i][j] = Point3::new(t.sin(), t.cos(), 0.5 - t);
        }
    }
    let bi5 = bi3_to_bi5_grid(&bi3);
    for &(u, v) in &[
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (1.0, 1.0),
        (0.3, 0.7),
        (0.5, 0.5),
    ] {
        let v3 = eval_bi3_grid(&bi3, u, v);
        let v5 = eval_bi5_grid(&bi5, u, v);
        assert!(approx_p(v3, v5, 1e-10), "({u}, {v}): {v3:?} vs {v5:?}");
    }
}

#[test]
fn bi3_to_bi5_grid_preserves_corners() {
    let bi3: [[Point3; 4]; 4] = [
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        [
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
        ],
        [
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(3.0, 2.0, 0.0),
        ],
        [
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
            Point3::new(3.0, 3.0, 0.0),
        ],
    ];
    let bi5 = bi3_to_bi5_grid(&bi3);
    assert_eq!(bi5[0][0], bi3[0][0]);
    assert_eq!(bi5[5][0], bi3[3][0]);
    assert_eq!(bi5[0][5], bi3[0][3]);
    assert_eq!(bi5[5][5], bi3[3][3]);
}

// ---- de Casteljau subdivision ----

#[test]
fn decasteljau_half_reproduces_polynomial() {
    // Evaluate on left and right sub-curves at various t, compare
    // with the original curve.
    let cps: [Point3; 6] = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(4.0, 5.0, 6.0),
        Point3::new(7.0, 8.0, 9.0),
        Point3::new(10.0, 11.0, 12.0),
        Point3::new(13.0, 14.0, 15.0),
    ];
    let (left, right) = decasteljau_bi5_at_half(&cps);
    for t in [0.0_f64, 0.3, 0.7, 1.0] {
        let v_left = eval_bi5_1d(&left, t);
        let v_orig_left = eval_bi5_1d(&cps, t / 2.0);
        assert!(
            approx_p(v_left, v_orig_left, 1e-10),
            "left at t={t}: {v_left:?} vs orig {v_orig_left:?}"
        );

        let v_right = eval_bi5_1d(&right, t);
        let v_orig_right = eval_bi5_1d(&cps, 0.5 + t / 2.0);
        assert!(
            approx_p(v_right, v_orig_right, 1e-10),
            "right at t={t}: {v_right:?} vs orig {v_orig_right:?}"
        );
    }
}

#[test]
fn decasteljau_half_endpoints() {
    let cps: [Point3; 6] = [
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(20.0, 0.0, 0.0),
        Point3::new(30.0, 0.0, 0.0),
        Point3::new(40.0, 0.0, 0.0),
        Point3::new(50.0, 0.0, 0.0),
        Point3::new(60.0, 0.0, 0.0),
    ];
    let (left, right) = decasteljau_bi5_at_half(&cps);
    assert_eq!(left[0], cps[0]);
    assert_eq!(right[5], cps[5]);
    // left[5] and right[0] both equal the curve's value at t=1/2.
    assert!(approx_p(left[5], right[0], 1e-10));
}

// ---- 2-D subdivision ----

#[test]
fn subdivide_bi5_2x2_reproduces_bi5_at_sub_points() {
    let mut cps = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let t = (i * 6 + j) as f64 * 0.17;
            cps[i][j] = Point3::new(t.sin(), t.cos(), t * 0.3);
        }
    }
    let subs = subdivide_bi5_2x2(&cps);

    for &(u, v) in &[(0.0, 0.0), (0.25, 0.75), (0.5, 0.5), (0.75, 0.25)] {
        let orig = eval_bi5_grid(&cps, u, v);
        // Which sub-patch covers (u, v)?
        let iu = if u > 0.5 { 1 } else { 0 };
        let jv = if v > 0.5 { 1 } else { 0 };
        let sub_u = if iu == 0 { u * 2.0 } else { (u - 0.5) * 2.0 };
        let sub_v = if jv == 0 { v * 2.0 } else { (v - 0.5) * 2.0 };
        let sub_val = eval_bi5_grid(&subs[iu][jv], sub_u, sub_v);
        assert!(
            approx_p(orig, sub_val, 1e-9),
            "({u}, {v}) in sub [{iu}][{jv}] at ({sub_u}, {sub_v}): orig {orig:?} vs sub {sub_val:?}"
        );
    }
}

#[test]
fn subdivide_preserves_outer_corners() {
    let cps: [[Point3; 6]; 6] =
        std::array::from_fn(|i| std::array::from_fn(|j| Point3::new(i as f64, j as f64, 0.0)));
    let subs = subdivide_bi5_2x2(&cps);
    // sub[0][0] bottom-left corner = original (0, 0) corner.
    assert_eq!(subs[0][0][0][0], cps[0][0]);
    // sub[1][1] top-right corner = original (1, 1) corner.
    assert_eq!(subs[1][1][5][5], cps[5][5]);
    // sub[1][0] bottom-right = original (1, 0).
    assert_eq!(subs[1][0][5][0], cps[5][0]);
    // sub[0][1] top-left = original (0, 1).
    assert_eq!(subs[0][1][0][5], cps[0][5]);
}

#[test]
fn subdivide_shared_boundaries_agree() {
    // Subdivided sub-patches share edges; the shared CPs should
    // agree between neighboring sub-patches.
    let mut cps = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            cps[i][j] = Point3::new(i as f64 * 0.3, j as f64 * 0.4, (i + j) as f64 * 0.1);
        }
    }
    let subs = subdivide_bi5_2x2(&cps);
    // sub[0][0]'s u=1 edge = sub[1][0]'s u=0 edge.
    for j in 0..6 {
        assert!(
            approx_p(subs[0][0][5][j], subs[1][0][0][j], 1e-10),
            "u-split mismatch at j={j}"
        );
    }
    // sub[0][0]'s v=1 edge = sub[0][1]'s v=0 edge.
    for i in 0..6 {
        assert!(
            approx_p(subs[0][0][i][5], subs[0][1][i][0], 1e-10),
            "v-split mismatch at i={i}"
        );
    }
}
