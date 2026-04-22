//! Tests for F_1 functional on the bi-5 cap.

#![allow(clippy::needless_range_loop)]

use super::super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, CapPieceId};
use super::*;
use crate::geometry::Point3;

const EPS: f64 = 1e-9;

#[test]
fn f1_of_constant_piece_is_zero() {
    let mut piece = Bi5CapPiece::zero();
    for i in 0..6 {
        for j in 0..6 {
            piece.cps[i][j] = Point3::new(5.0, -2.7, 1.4);
        }
    }
    assert!(f1_piece(&piece) < 1e-12);
}

#[test]
fn f1_of_zero_piece_is_zero() {
    let piece = Bi5CapPiece::zero();
    assert_eq!(f1_piece(&piece), 0.0);
}

#[test]
fn f1_scales_quadratically_with_cps() {
    // F_1 is a quadratic form → F_1(2·P) = 4·F_1(P).
    let mut p1 = Bi5CapPiece::zero();
    for i in 0..6 {
        for j in 0..6 {
            p1.cps[i][j] = Point3::new(i as f64 * 0.1, j as f64 * 0.1, (i + j) as f64 * 0.05);
        }
    }
    let mut p2 = p1;
    for i in 0..6 {
        for j in 0..6 {
            p2.cps[i][j].x *= 2.0;
            p2.cps[i][j].y *= 2.0;
            p2.cps[i][j].z *= 2.0;
        }
    }
    let f1_a = f1_piece(&p1);
    let f1_b = f1_piece(&p2);
    assert!(
        (f1_b - 4.0 * f1_a).abs() < 1e-9 * f1_a.abs().max(1.0),
        "F_1(2P) = {f1_b}, expected 4·F_1(P) = {}",
        4.0 * f1_a
    );
}

#[test]
fn f1_of_linear_piece_matches_analytic() {
    // Linear function f(u, v) = a u + b v (scalar).
    // ∂u f = a, ∂v f = b. (∂u f)² + (∂v f)² = a² + b².
    // F_1 = ∫∫ (a² + b²) du dv = a² + b² over [0,1]² = a² + b².
    //
    // In bi-5 BB form, the CPs for f(u, v) = a u + b v are CP[i][j]
    // = a (i/5) + b (j/5) (linear has CP = function value at nodes).
    let mut piece = Bi5CapPiece::zero();
    let a = 2.0;
    let b = 3.0;
    for i in 0..6 {
        for j in 0..6 {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            piece.cps[i][j] = Point3::new(a * u + b * v, 0.0, 0.0);
        }
    }
    // Only x-axis non-zero, so F_1 = (a² + b²) · 1 for x-axis.
    // Other axes contribute 0. So total = a² + b² = 4 + 9 = 13.
    let f1 = f1_piece(&piece);
    assert!(
        (f1 - 13.0).abs() < 1e-10,
        "linear a={a}, b={b}: F_1 = {f1}, expected 13"
    );
}

#[test]
fn f1_cap_sums_across_pieces_and_sectors() {
    // Construct a cap where only sector 0's O piece has a known
    // F_1, rest zero. Verify f1_cap sums correctly.
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    for i in 0..6 {
        for j in 0..6 {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            o.cps[i][j] = Point3::new(u + v, 0.0, 0.0);
        }
    }
    let expected_o = 2.0; // a=1, b=1 → 1² + 1² = 2.

    let f1 = f1_cap(&cap);
    assert!(
        (f1 - expected_o).abs() < 1e-10,
        "F_1(cap) = {f1}, expected {expected_o}"
    );
}

#[test]
fn f1_difference_is_zero_when_p_equals_a() {
    // F_1(p − a) = 0 when p = a.
    let mut a = Bi5Cap::zero(5);
    // Non-zero a.
    for k in 0..5 {
        let o = a.sector_mut(k).piece_mut(CapPieceId::O);
        for i in 0..6 {
            for j in 0..6 {
                o.cps[i][j] = Point3::new(k as f64 + i as f64, j as f64, 0.0);
            }
        }
    }
    let p = a.clone();
    let f1d = f1_difference(&p, &a);
    assert!(f1d.abs() < 1e-12, "F_1(a - a) = {f1d} ≠ 0");
}

#[test]
fn f1_difference_measures_cp_difference() {
    // If p and a differ only by a linear function f(u, v) = c_u u + c_v v
    // on one piece, F_1(p − a) should equal c_u² + c_v².
    let a = Bi5Cap::zero(4);
    let mut p = Bi5Cap::zero(4);
    let c_u = 0.5;
    let c_v = 0.7;
    let o_p = p.sector_mut(0).piece_mut(CapPieceId::O);
    for i in 0..6 {
        for j in 0..6 {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            o_p.cps[i][j] = Point3::new(c_u * u + c_v * v, 0.0, 0.0);
        }
    }
    let _ = a.clone(); // ensure a stays zero
    let f1d = f1_difference(&p, &a);
    let expected = c_u * c_u + c_v * c_v;
    assert!(
        (f1d - expected).abs() < 1e-10,
        "F_1(p - a) = {f1d}, expected {expected}"
    );
}

#[test]
fn f1_cap_deterministic() {
    let mut cap = Bi5Cap::zero(5);
    for k in 0..5 {
        for id in CapPieceId::all() {
            let p = cap.sector_mut(k).piece_mut(id);
            for i in 0..6 {
                for j in 0..6 {
                    p.cps[i][j] = Point3::new((k + i) as f64 * 0.1, j as f64 * 0.2, 0.0);
                }
            }
        }
    }
    let a = f1_cap(&cap);
    let b = f1_cap(&cap);
    assert_eq!(a, b);
}
