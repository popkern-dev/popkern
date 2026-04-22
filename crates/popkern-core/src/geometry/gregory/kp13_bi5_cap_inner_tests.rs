//! Tests for the I.14.c `O` piece populator from the guide.

#![allow(clippy::needless_range_loop)]

use super::super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use super::super::kp13_bi5_chi::cc_subdominant_eigenvalue;
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_guide::Bi5Guide;
use super::*;

const EPS: f64 = 1e-9;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- Structural ----

#[test]
fn zero_guide_leaves_o_piece_zero() {
    for n in [5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
        let guide = Bi5Guide::from_gtri5(&g_tri);

        let mut cap = Bi5Cap::zero(n);
        populate_o_piece_from_guide_simplified(&mut cap, &guide);
        for k in 0..n {
            let o = cap.sector(k).piece(CapPieceId::O);
            for p in 0..6 {
                for q in 0..6 {
                    let cp = o.cps[p][q];
                    assert!(
                        cp.x.abs() < 1e-9 && cp.y.abs() < 1e-9 && cp.z.abs() < 1e-9,
                        "n={n}, k={k}, O ({p}, {q}): non-zero {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn populate_does_not_touch_other_pieces() {
    let n = 5;
    let params = Bi5FreeParams::zero(n);
    let g_tri = apply_free_params(&params, Point3::new(1.0, 2.0, 3.0)).unwrap();
    let guide = Bi5Guide::from_gtri5(&g_tri);

    let mut cap = Bi5Cap::zero(n);
    populate_o_piece_from_guide_simplified(&mut cap, &guide);

    // L, R, M pieces should all remain zero.
    for k in 0..n {
        for id in [CapPieceId::L, CapPieceId::R, CapPieceId::M] {
            let piece = cap.sector(k).piece(id);
            for p in 0..6 {
                for q in 0..6 {
                    let cp = piece.cps[p][q];
                    assert!(
                        cp.x.abs() < 1e-9 && cp.y.abs() < 1e-9 && cp.z.abs() < 1e-9,
                        "n={n}, k={k}, {id:?} ({p}, {q}): non-zero {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn populate_sets_only_corner_block_cps() {
    // Of the 36 CPs in each O piece, exactly 24 (the four 6-CP
    // corner blocks) are set; the 12 non-corner-block CPs remain
    // at their pre-call values (zero from Bi5Cap::zero).
    let n = 6;
    let params = Bi5FreeParams::zero(n);
    let mut p = params.clone();
    p.q.alpha = Point3::new(0.5, 0.5, 0.5);
    p.q.delta = Point3::new(0.3, 0.3, 0.3);
    let g_tri = apply_free_params(&p, Point3::new(0.1, 0.2, 0.3)).unwrap();
    let guide = Bi5Guide::from_gtri5(&g_tri);

    let mut cap = Bi5Cap::zero(n);
    populate_o_piece_from_guide_simplified(&mut cap, &guide);

    // Non-corner-block CPs that should stay zero. For a 6x6 bi-5
    // piece with 4 corner blocks at {[0..=2]×[0..=2], [3..=5]×[0..=2],
    // [0..=2]×[3..=5], [3..=5]×[3..=5]}, the corner-block-jet
    // setters touch exactly 6 CPs per block. The 3 untouched CPs
    // per block are at positions:
    //   (0, 0) block: P[2][1], P[1][2], P[2][2]
    //   (5, 0) block: P[3][1], P[4][2], P[3][2]
    //   (0, 5) block: P[1][3], P[2][4], P[2][3]
    //   (5, 5) block: P[3][4], P[4][3], P[3][3]
    // These 12 CPs stay zero.
    let untouched = [
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (4, 2),
        (3, 2),
        (1, 3),
        (2, 4),
        (2, 3),
        (3, 4),
        (4, 3),
        (3, 3),
    ];

    for k in 0..n {
        let o = cap.sector(k).piece(CapPieceId::O);
        for &(p, q) in &untouched {
            let cp = o.cps[p][q];
            assert!(
                cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                "n={n}, k={k}, O ({p}, {q}): should be zero, got {cp:?}"
            );
        }
    }
}

// ---- Correctness: corner value matches guide eval at (u_p/(2λ)) ----

#[test]
fn corner_values_match_guide_at_scaled_position() {
    let n = 5;
    let params = Bi5FreeParams::zero(n);
    let mut p = params.clone();
    p.q.alpha = Point3::new(0.1, 0.2, 0.3);
    p.q.beta = Point3::new(0.4, 0.5, 0.6);
    let star = Point3::new(0.7, 0.8, 0.9);
    let g_tri = apply_free_params(&p, star).unwrap();
    let guide = Bi5Guide::from_gtri5(&g_tri);

    let mut cap = Bi5Cap::zero(n);
    populate_o_piece_from_guide_simplified(&mut cap, &guide);

    let lambda = cc_subdominant_eigenvalue(n);
    let inv_2lambda = 0.5 / lambda;

    for k in 0..n {
        let o = cap.sector(k).piece(CapPieceId::O);
        let g = guide.sector(k);

        // Corner (0, 0) piece-local → (0, 0) guide; value = P[0][0].
        assert!(
            approx_p(o.cps[0][0], g.eval(0.0, 0.0), EPS),
            "k={k}, corner (0, 0): o{:?} vs g{:?}",
            o.cps[0][0],
            g.eval(0.0, 0.0)
        );
        // Corner (1, 0) → (inv_2lambda, 0); value = P[5][0].
        let expected_10 = g.eval(inv_2lambda, 0.0);
        assert!(
            approx_p(o.cps[5][0], expected_10, EPS),
            "k={k}, corner (1, 0): o{:?} vs g{:?}",
            o.cps[5][0],
            expected_10
        );
        // Corner (0, 1) → (0, inv_2lambda); value = P[0][5].
        let expected_01 = g.eval(0.0, inv_2lambda);
        assert!(
            approx_p(o.cps[0][5], expected_01, EPS),
            "k={k}, corner (0, 1): o{:?} vs g{:?}",
            o.cps[0][5],
            expected_01
        );
        // Corner (1, 1) → (inv_2lambda, inv_2lambda); value = P[5][5].
        let expected_11 = g.eval(inv_2lambda, inv_2lambda);
        assert!(
            approx_p(o.cps[5][5], expected_11, EPS),
            "k={k}, corner (1, 1): o{:?} vs g{:?}",
            o.cps[5][5],
            expected_11
        );
    }
}

// ---- Star-point invariant: P[0][0] of every sector's O piece = star ----

#[test]
fn o_piece_star_corner_equals_cc_limit_via_guide() {
    let n = 7;
    let params = Bi5FreeParams::zero(n);
    let star = Point3::new(1.5, -0.5, 0.25);
    let g_tri = apply_free_params(&params, star).unwrap();
    let guide = Bi5Guide::from_gtri5(&g_tri);

    let mut cap = Bi5Cap::zero(n);
    populate_o_piece_from_guide_simplified(&mut cap, &guide);

    // Every sector's O[0][0] equals the guide's value at (0, 0),
    // which is the CC limit point (star) by I.12's degree-raise.
    for k in 0..n {
        let o = cap.sector(k).piece(CapPieceId::O);
        assert!(
            approx_p(o.cps[0][0], star, EPS),
            "k={k}: O[0][0]={:?} vs star={:?}",
            o.cps[0][0],
            star
        );
    }
}

// ---- All valences covered ----

#[test]
fn works_for_all_valences() {
    for n in [3, 4, 5, 6, 7, 8] {
        let params = Bi5FreeParams::zero(n);
        let g_tri = apply_free_params(&params, Point3::new(0.2, 0.3, 0.4)).unwrap();
        let guide = Bi5Guide::from_gtri5(&g_tri);
        let mut cap = Bi5Cap::zero(n);
        populate_o_piece_from_guide_simplified(&mut cap, &guide);
        cap.validate().unwrap_or_else(|e| panic!("n={n}: {e}"));
    }
}

// ---- Determinism ----

#[test]
fn deterministic() {
    let n = 5;
    let params = Bi5FreeParams::zero(n);
    let mut p = params.clone();
    p.q.alpha = Point3::new(0.1, 0.0, 0.0);
    let g_tri = apply_free_params(&p, Point3::new(0.5, 0.5, 0.5)).unwrap();
    let guide = Bi5Guide::from_gtri5(&g_tri);
    let mut cap_a = Bi5Cap::zero(n);
    let mut cap_b = Bi5Cap::zero(n);
    populate_o_piece_from_guide_simplified(&mut cap_a, &guide);
    populate_o_piece_from_guide_simplified(&mut cap_b, &guide);
    assert_eq!(cap_a, cap_b);
}

// ---- n-mismatch guard ----

#[test]
#[should_panic]
fn populate_rejects_n_mismatch() {
    let guide = {
        let params = Bi5FreeParams::zero(5);
        let g_tri = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
        Bi5Guide::from_gtri5(&g_tri)
    };
    let mut cap = Bi5Cap::zero(6);
    populate_o_piece_from_guide_simplified(&mut cap, &guide);
}
