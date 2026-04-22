//! Tests for the bi-5 cap scaffolding (KP13 § 5, I.14.a).

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- CapPieceId ----

#[test]
fn piece_ids_canonical_order() {
    let all = CapPieceId::all();
    assert_eq!(all[0], CapPieceId::O);
    assert_eq!(all[1], CapPieceId::L);
    assert_eq!(all[2], CapPieceId::R);
    assert_eq!(all[3], CapPieceId::M);
}

#[test]
fn piece_id_indices_unique_and_in_range() {
    let indices: Vec<usize> = CapPieceId::all().iter().map(|p| p.index()).collect();
    assert_eq!(indices, vec![0, 1, 2, 3]);
}

// ---- Bi5CapPiece ----

#[test]
fn zero_piece_evaluates_to_zero() {
    let p = Bi5CapPiece::zero();
    for &(u, v) in &[(0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (0.3, 0.7)] {
        let v_out = p.eval(u, v);
        assert!(
            approx_p(v_out, Point3::new(0.0, 0.0, 0.0), 1e-12),
            "({u}, {v}): {v_out:?}"
        );
    }
}

#[test]
fn piece_eval_at_corners_equals_corner_cps() {
    let mut p = Bi5CapPiece::zero();
    p.cps[0][0] = Point3::new(1.0, 2.0, 3.0);
    p.cps[5][0] = Point3::new(4.0, 5.0, 6.0);
    p.cps[0][5] = Point3::new(7.0, 8.0, 9.0);
    p.cps[5][5] = Point3::new(10.0, 11.0, 12.0);
    assert!(approx_p(p.eval(0.0, 0.0), p.cps[0][0], 1e-12));
    assert!(approx_p(p.eval(1.0, 0.0), p.cps[5][0], 1e-12));
    assert!(approx_p(p.eval(0.0, 1.0), p.cps[0][5], 1e-12));
    assert!(approx_p(p.eval(1.0, 1.0), p.cps[5][5], 1e-12));
}

// ---- dispatch_piece ----

#[test]
fn dispatch_piece_corners() {
    // (0, 0) → O with local coords (0, 0).
    let (id, u, v) = dispatch_piece(0.0, 0.0);
    assert_eq!(id, CapPieceId::O);
    assert!(u.abs() < EPS && v.abs() < EPS);
    // (0.25, 0.25) → O with local (0.5, 0.5).
    let (id, u, v) = dispatch_piece(0.25, 0.25);
    assert_eq!(id, CapPieceId::O);
    assert!((u - 0.5).abs() < EPS && (v - 0.5).abs() < EPS);
    // (0.75, 0.25) → R with local (0.5, 0.5).
    let (id, u, v) = dispatch_piece(0.75, 0.25);
    assert_eq!(id, CapPieceId::R);
    assert!((u - 0.5).abs() < EPS && (v - 0.5).abs() < EPS);
    // (0.25, 0.75) → L with local (0.5, 0.5).
    let (id, u, v) = dispatch_piece(0.25, 0.75);
    assert_eq!(id, CapPieceId::L);
    assert!((u - 0.5).abs() < EPS && (v - 0.5).abs() < EPS);
    // (0.75, 0.75) → M with local (0.5, 0.5).
    let (id, u, v) = dispatch_piece(0.75, 0.75);
    assert_eq!(id, CapPieceId::M);
    assert!((u - 0.5).abs() < EPS && (v - 0.5).abs() < EPS);
    // (1.0, 1.0) → M with local (1, 1).
    let (id, u, v) = dispatch_piece(1.0, 1.0);
    assert_eq!(id, CapPieceId::M);
    assert!((u - 1.0).abs() < EPS && (v - 1.0).abs() < EPS);
}

#[test]
fn dispatch_piece_boundary_resolves_to_lower() {
    // On the internal boundary (u = 0.5), resolve to left-side.
    let (id, _, _) = dispatch_piece(0.5, 0.3);
    assert_eq!(id, CapPieceId::O);
    // On (v = 0.5), resolve to lower-side.
    let (id, _, _) = dispatch_piece(0.3, 0.5);
    assert_eq!(id, CapPieceId::O);
    // At (0.5, 0.5) → O.
    let (id, _, _) = dispatch_piece(0.5, 0.5);
    assert_eq!(id, CapPieceId::O);
}

// ---- Bi5CapSector ----

#[test]
fn sector_eval_dispatches_correctly() {
    // Set each piece to a distinct constant and verify eval
    // dispatches to the right piece.
    let mut sector = Bi5CapSector::zero();
    let values = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    for (id, &val) in CapPieceId::all().iter().zip(values.iter()) {
        let piece = sector.piece_mut(*id);
        for p in 0..6 {
            for q in 0..6 {
                piece.cps[p][q] = val;
            }
        }
    }
    assert!(approx_p(sector.eval(0.25, 0.25), values[0], EPS));
    assert!(approx_p(sector.eval(0.25, 0.75), values[1], EPS));
    assert!(approx_p(sector.eval(0.75, 0.25), values[2], EPS));
    assert!(approx_p(sector.eval(0.75, 0.75), values[3], EPS));
}

// ---- Bi5Cap ----

#[test]
fn zero_cap_validates() {
    for n in [3, 4, 5, 6, 7, 10] {
        let cap = Bi5Cap::zero(n);
        cap.validate().unwrap();
        assert_eq!(cap.n, n);
        assert_eq!(cap.sectors.len(), n);
    }
}

#[test]
#[should_panic]
fn zero_cap_below_three_panics() {
    let _ = Bi5Cap::zero(2);
}

#[test]
fn validate_detects_sector_count_mismatch() {
    let mut cap = Bi5Cap::zero(5);
    cap.sectors.pop();
    assert!(cap.validate().is_err());
}

#[test]
fn validate_detects_non_finite_cp() {
    let mut cap = Bi5Cap::zero(5);
    cap.sectors[3].piece_mut(CapPieceId::M).cps[2][4].y = f64::NAN;
    assert!(cap.validate().is_err());
}

#[test]
fn sector_mut_round_trips() {
    let mut cap = Bi5Cap::zero(4);
    let star = Point3::new(42.0, -7.0, 9.0);
    cap.sector_mut(2).piece_mut(CapPieceId::O).cps[0][0] = star;
    assert_eq!(cap.sector(2).piece(CapPieceId::O).cps[0][0], star);
    // Other sectors untouched.
    assert_eq!(
        cap.sector(0).piece(CapPieceId::O).cps[0][0],
        Point3::new(0.0, 0.0, 0.0)
    );
}

#[test]
fn sector_access_wraps_mod_n() {
    let cap = Bi5Cap::zero(5);
    let s0 = cap.sector(0);
    let s5 = cap.sector(5);
    assert_eq!(s0, s5);
    let s_big = cap.sector(5 * 100 + 3);
    let s3 = cap.sector(3);
    assert_eq!(s_big, s3);
}

// ---- Piece counts ----

#[test]
fn constants_match_bi5_tensor_product() {
    assert_eq!(BI5_PIECE_CP_COUNT, 6 * 6);
    assert_eq!(PIECES_PER_SECTOR, 4);
}
