//! Tests for [`pns_center`]: eq (4) `p_{22}` + eq (5) central point.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

#[test]
fn p22_zero_inputs_give_zero() {
    let z = pt(0.0, 0.0, 0.0);
    assert!(approx_pt(
        p22_symmetric_average(z, z, z, z, z, z, z, z),
        z,
        1e-12
    ));
}

#[test]
fn p22_constant_inputs_preserve_value() {
    // If every input is the same point P, both the horizontal and
    // vertical reductions give P: (2/3)(P+P) − (1/6)(P+P) = P.
    // Averaging gives P.
    let p = pt(3.0, 2.5, 1.5);
    assert!(approx_pt(
        p22_symmetric_average(p, p, p, p, p, p, p, p),
        p,
        1e-12
    ));
}

#[test]
fn p22_degree_reduction_identity_on_linear_row() {
    // If column 2 is a linear sampling (p_{i2} on a line), then
    // the horizontal reduction `(2/3)(p12+p32) − (1/6)(p02+p42)`
    // equals the value at the midpoint (degree-reduction preserves
    // straight lines).
    //
    // Line: p_{i2} = (i, 0, 0) for i = 0..4. Midpoint at i=2 is (2, 0, 0).
    // Check (2/3)((1,0,0)+(3,0,0)) − (1/6)((0,0,0)+(4,0,0))
    //     = (2/3)(4,0,0) − (1/6)(4,0,0) = (8/3 − 4/6, 0, 0) = (8/3 − 2/3, 0, 0)
    //     = (2, 0, 0) ✓
    //
    // Full test with symmetric inputs (row 2 ≡ column 2):
    let p02 = pt(0.0, 0.0, 0.0);
    let p12 = pt(1.0, 0.0, 0.0);
    let p32 = pt(3.0, 0.0, 0.0);
    let p42 = pt(4.0, 0.0, 0.0);
    let p20 = pt(0.0, 0.0, 0.0);
    let p21 = pt(1.0, 0.0, 0.0);
    let p23 = pt(3.0, 0.0, 0.0);
    let p24 = pt(4.0, 0.0, 0.0);
    let got = p22_symmetric_average(p02, p12, p32, p42, p20, p21, p23, p24);
    assert!(approx_pt(got, pt(2.0, 0.0, 0.0), 1e-12));
}

#[test]
fn central_point_zero_inputs() {
    let z = pt(0.0, 0.0, 0.0);
    let d1 = vec![z, z, z, z, z];
    let d24 = vec![z, z, z, z, z];
    for &e in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        assert!(approx_pt(central_point(&d1, &d24, e), z, 1e-12));
    }
}

#[test]
fn central_point_constant_inputs() {
    // All d1 = P, all (d2+d4)/2 = Q → eop = e·P + (1-e)·Q.
    let p = pt(1.0, 0.0, 0.0);
    let q = pt(0.0, 1.0, 0.0);
    let d1 = vec![p; 6];
    let d24 = vec![q; 6];
    for &e in &[0.0, 0.1, 0.5, 0.9, 1.0] {
        let expected = pt(e * p.x + (1.0 - e) * q.x, e * p.y + (1.0 - e) * q.y, 0.0);
        assert!(approx_pt(central_point(&d1, &d24, e), expected, 1e-12));
    }
}

#[test]
fn central_point_e_zero_averages_d24() {
    // e = 0 → output is average of d24 across k.
    let d1 = vec![pt(99.0, 99.0, 99.0); 5];
    let d24 = vec![
        pt(1.0, 0.0, 0.0),
        pt(2.0, 0.0, 0.0),
        pt(3.0, 0.0, 0.0),
        pt(4.0, 0.0, 0.0),
        pt(5.0, 0.0, 0.0),
    ];
    let got = central_point(&d1, &d24, 0.0);
    assert!(approx_pt(got, pt(3.0, 0.0, 0.0), 1e-12));
}

#[test]
fn central_point_e_one_averages_d1() {
    let d1 = vec![
        pt(1.0, 0.0, 0.0),
        pt(2.0, 0.0, 0.0),
        pt(3.0, 0.0, 0.0),
        pt(4.0, 0.0, 0.0),
        pt(5.0, 0.0, 0.0),
    ];
    let d24 = vec![pt(99.0, 99.0, 99.0); 5];
    let got = central_point(&d1, &d24, 1.0);
    assert!(approx_pt(got, pt(3.0, 0.0, 0.0), 1e-12));
}

#[test]
fn e_def_bi4_paper_values() {
    // Paper: ε = 0 for n ∈ {5, 6, 7}, ε = n/100 for n ∈ {8, 9, 10}.
    //
    // e_def(5) = 1.5 / 6.6
    // e_def(6) = 1.5 / 7.6
    // e_def(7) = 1.5 / 8.6
    // e_def(8) = (1.5 + 0.08) / 9.6
    // e_def(9) = (1.5 + 0.09) / 10.6
    // e_def(10) = (1.5 + 0.10) / 11.6
    assert!((e_def_bi4(5) - 1.5 / 6.6).abs() < 1e-12);
    assert!((e_def_bi4(6) - 1.5 / 7.6).abs() < 1e-12);
    assert!((e_def_bi4(7) - 1.5 / 8.6).abs() < 1e-12);
    assert!((e_def_bi4(8) - (1.5 + 0.08) / 9.6).abs() < 1e-12);
    assert!((e_def_bi4(9) - (1.5 + 0.09) / 10.6).abs() < 1e-12);
    assert!((e_def_bi4(10) - (1.5 + 0.10) / 11.6).abs() < 1e-12);
}

#[test]
fn e_def_bi3_n3_n5_is_022() {
    assert_eq!(E_DEF_BI3_N3_N5, 0.22);
}

#[test]
fn e_def_bi3_2x2_paper_values() {
    // e_def(n) = (0.02n + 1.5) / (n + 1.75).
    // n=6 → (0.12 + 1.5) / 7.75 = 1.62 / 7.75 ≈ 0.209032...
    let n6 = e_def_bi3_2x2(6);
    assert!((n6 - 1.62 / 7.75).abs() < 1e-12);
    let n8 = e_def_bi3_2x2(8);
    assert!((n8 - (0.02 * 8.0 + 1.5) / (8.0 + 1.75)).abs() < 1e-12);
}
