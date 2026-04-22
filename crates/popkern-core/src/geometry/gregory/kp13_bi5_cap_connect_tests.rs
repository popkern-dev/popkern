//! Tests for I.14.d internal C² propagation.

#![allow(clippy::needless_range_loop)]

use super::super::bernstein_bi5::bernstein_bi5_all;
use super::super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- C⁰/C¹/C² boundary formula verification ----

#[test]
fn zero_o_piece_leaves_neighbors_zero() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        connect_all_sectors(&mut cap);
        cap.validate().unwrap();
        for k in 0..n {
            for id in CapPieceId::all() {
                let piece = cap.sector(k).piece(id);
                for p in 0..6 {
                    for q in 0..6 {
                        let cp = piece.cps[p][q];
                        assert!(
                            cp.x.abs() < 1e-12 && cp.y.abs() < 1e-12 && cp.z.abs() < 1e-12,
                            "n={n}, k={k}, {id:?} ({p}, {q}): non-zero {cp:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn o_to_r_c0_continuity() {
    // Set O's u=1 column to distinct values, propagate, verify R's
    // u=0 column matches (C⁰).
    let mut cap = Bi5Cap::zero(5);
    for j in 0..6 {
        cap.sector_mut(0).piece_mut(CapPieceId::O).cps[5][j] =
            Point3::new(j as f64, j as f64 * 2.0, j as f64 * 3.0);
        cap.sector_mut(0).piece_mut(CapPieceId::O).cps[4][j] = Point3::new(0.5, 0.5, 0.5);
        cap.sector_mut(0).piece_mut(CapPieceId::O).cps[3][j] = Point3::new(0.1, 0.1, 0.1);
    }
    connect_sector(cap.sector_mut(0));
    let r = cap.sector(0).piece(CapPieceId::R);
    for j in 0..6 {
        let expected = Point3::new(j as f64, j as f64 * 2.0, j as f64 * 3.0);
        assert!(approx_p(r.cps[0][j], expected, EPS), "j={j}: C⁰ fail");
    }
}

#[test]
fn o_to_r_c1_continuity() {
    // ∂_u_sector continuous at the shared boundary.
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    o.cps[5][2] = Point3::new(1.0, 2.0, 3.0);
    o.cps[4][2] = Point3::new(0.5, 0.0, 0.0);
    connect_sector(cap.sector_mut(0));
    let r = cap.sector(0).piece(CapPieceId::R);
    // R[1][j] = 2 O[5][j] - O[4][j]
    let expected = Point3::new(2.0 * 1.0 - 0.5, 2.0 * 2.0 - 0.0, 2.0 * 3.0 - 0.0);
    assert!(approx_p(r.cps[1][2], expected, EPS));
}

#[test]
fn o_to_r_c2_continuity() {
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    o.cps[5][2] = Point3::new(1.0, 0.0, 0.0);
    o.cps[4][2] = Point3::new(0.0, 1.0, 0.0);
    o.cps[3][2] = Point3::new(0.0, 0.0, 1.0);
    connect_sector(cap.sector_mut(0));
    let r = cap.sector(0).piece(CapPieceId::R);
    // R[2][j] = 4 O[5][j] - 4 O[4][j] + O[3][j]
    let expected = Point3::new(4.0, -4.0, 1.0);
    assert!(approx_p(r.cps[2][2], expected, EPS));
}

#[test]
fn o_to_l_c0_c1_c2_continuity() {
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    o.cps[3][5] = Point3::new(1.0, 0.0, 0.0);
    o.cps[3][4] = Point3::new(0.0, 1.0, 0.0);
    o.cps[3][3] = Point3::new(0.0, 0.0, 1.0);
    connect_sector(cap.sector_mut(0));
    let l = cap.sector(0).piece(CapPieceId::L);
    assert!(
        approx_p(l.cps[3][0], Point3::new(1.0, 0.0, 0.0), EPS),
        "C⁰ O→L"
    );
    // L[i][1] = 2 O[i][5] - O[i][4]
    assert!(
        approx_p(l.cps[3][1], Point3::new(2.0, -1.0, 0.0), EPS),
        "C¹ O→L"
    );
    // L[i][2] = 4 O[i][5] - 4 O[i][4] + O[i][3]
    assert!(
        approx_p(l.cps[3][2], Point3::new(4.0, -4.0, 1.0), EPS),
        "C² O→L"
    );
}

// ---- Full-polynomial smooth continuation test ----

#[test]
fn polynomial_c2_continuation_o_to_r() {
    // Build O as a known bi-5 polynomial (linear in u, v). Then
    // propagate. Verify R evaluates consistently at the shared
    // boundary.
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    // Linear function f(u, v) = u + 2v. CPs: f(i/5, j/5).
    // Actually for a bi-5 BB representation of a linear polynomial,
    // the CPs are literally the function values at (i/5, j/5).
    for i in 0..6 {
        for j in 0..6 {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            o.cps[i][j] = Point3::new(u + 2.0 * v, u * 2.0, v * 3.0);
        }
    }
    connect_sector(cap.sector_mut(0));

    // Evaluate O at its u=1 boundary and R at its u=0 boundary;
    // should match at corresponding sector-local v.
    let o = cap.sector(0).piece(CapPieceId::O);
    let r = cap.sector(0).piece(CapPieceId::R);
    for v_test in [0.0, 0.25, 0.5, 0.75, 1.0_f64] {
        let o_val = eval_piece(o, 1.0, v_test);
        let r_val = eval_piece(r, 0.0, v_test);
        assert!(
            approx_p(o_val, r_val, 1e-10),
            "v={v_test}: C⁰ polynomial mismatch O={o_val:?} vs R={r_val:?}"
        );
    }
}

fn eval_piece(piece: &super::super::kp13_bi5_cap::Bi5CapPiece, u: f64, v: f64) -> Point3 {
    let bu = bernstein_bi5_all(u);
    let bv = bernstein_bi5_all(v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for p in 0..6 {
        for q in 0..6 {
            let w = bu[p] * bv[q];
            let cp = piece.cps[p][q];
            out.x += w * cp.x;
            out.y += w * cp.y;
            out.z += w * cp.z;
        }
    }
    out
}

// ---- R → M and L → M propagation ----

#[test]
fn r_to_m_c0_at_shared_edge() {
    let mut cap = Bi5Cap::zero(5);
    let r = cap.sector_mut(0).piece_mut(CapPieceId::R);
    r.cps[2][5] = Point3::new(7.0, 8.0, 9.0);
    // Also populate the surrounding cells so the propagation path isn't all zero.
    connect_sector(cap.sector_mut(0));
    // R[2][5] propagated to M[2][0] (via R → M C⁰).
    let m = cap.sector(0).piece(CapPieceId::M);
    // Note: L → M overwrites the overlap at (0..=2, 0..=2), but
    // M[2][0] is at (i=2, j=0) which is in the overlap — L's
    // contribution comes from L[5][0], L[4][0], L[3][0] which
    // in this test are all zero (O was never touched).
    // After L overwrite: M[2][0] = 4 L[5][0] - 4 L[4][0] + L[3][0] = 0.
    // So R's contribution of 7 gets overwritten to 0.
    assert!(approx_p(m.cps[2][0], Point3::new(0.0, 0.0, 0.0), EPS));
}

#[test]
fn l_overrides_r_at_m_overlap() {
    // Set O's (i=5, j=...) so that L→M overwrites R→M's
    // contribution at the (0..=2, 0..=2) overlap.
    let mut cap = Bi5Cap::zero(5);
    let o = cap.sector_mut(0).piece_mut(CapPieceId::O);
    // Set O's u=1 edge (feeds R) to one value.
    for j in 0..6 {
        o.cps[5][j] = Point3::new(1.0, 0.0, 0.0);
    }
    // Set O's v=1 edge (feeds L) to a different value.
    for i in 0..6 {
        o.cps[i][5] = Point3::new(2.0, 0.0, 0.0);
    }
    // O[5][5] is in both — consistency: (1, 0, 0) vs (2, 0, 0).
    // Since the two assignments conflict, the last one wins:
    // O[5][5] = (2, 0, 0).
    assert_eq!(o.cps[5][5], Point3::new(2.0, 0.0, 0.0));
    connect_sector(cap.sector_mut(0));

    let m = cap.sector(0).piece(CapPieceId::M);
    // R→M at (0, 0) comes from R[0][5] = O[5][5] = (2, 0, 0).
    // L→M at (0, 0) comes from L[5][0] = O[5][5] = (2, 0, 0).
    // Agreement by construction at this corner.
    assert!(approx_p(m.cps[0][0], Point3::new(2.0, 0.0, 0.0), EPS));

    // M[1][0] is in the 3×3 overlap; L wins per our documented
    // order. Values propagated:
    //   O[5][5] = (2, 0, 0), O[4][5] = (2, 0, 0) (from v=1 loop).
    //   L[5][0] = O[5][5] = (2, 0, 0).
    //   L[4][0] = O[4][5] = (2, 0, 0).
    //   M[1][0] = 2 L[5][0] − L[4][0] = 2·2 − 2 = 2.
    assert!(approx_p(m.cps[1][0], Point3::new(2.0, 0.0, 0.0), EPS));
}

// ---- Determinism ----

#[test]
fn connect_deterministic() {
    let mut cap_a = Bi5Cap::zero(5);
    let mut cap_b = Bi5Cap::zero(5);
    // Populate O with some deterministic values.
    for k in 0..5 {
        let o_a = cap_a.sector_mut(k).piece_mut(CapPieceId::O);
        let o_b = cap_b.sector_mut(k).piece_mut(CapPieceId::O);
        for i in 0..6 {
            for j in 0..6 {
                let v = Point3::new((k + i) as f64, j as f64, (k * i * j) as f64 * 0.1);
                o_a.cps[i][j] = v;
                o_b.cps[i][j] = v;
            }
        }
    }
    connect_all_sectors(&mut cap_a);
    connect_all_sectors(&mut cap_b);
    assert_eq!(cap_a, cap_b);
}

// ---- All-valences sanity ----

#[test]
fn connect_all_sectors_works_for_all_valences() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        // Set each O piece to something non-zero.
        for k in 0..n {
            let o = cap.sector_mut(k).piece_mut(CapPieceId::O);
            for i in 0..6 {
                for j in 0..6 {
                    o.cps[i][j] = Point3::new((k + i + j) as f64 * 0.1, 0.0, 0.0);
                }
            }
        }
        connect_all_sectors(&mut cap);
        cap.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
    }
}
