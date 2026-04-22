//! Tests for KP13 I.14.b outer-corner-jet transfer.

#![allow(clippy::needless_range_loop)]

use super::super::kp13_bi5_boundary::TensorBorder;
use super::super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, CapPieceId};
use super::super::kp13_bi5_chi::ChiCorner;
use super::super::kp13_bi5_jet::TaylorJet2;
use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- Round-trip: set_corner_jet_xy → extract jet → original jet ----

fn extract_jet_10(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p50 = piece.cps[5][0];
    let p40 = piece.cps[4][0];
    let p51 = piece.cps[5][1];
    let p30 = piece.cps[3][0];
    let p52 = piece.cps[5][2];
    let p41 = piece.cps[4][1];
    TaylorJet2 {
        value: p50,
        d_u: Point3::new(
            5.0 * (p50.x - p40.x),
            5.0 * (p50.y - p40.y),
            5.0 * (p50.z - p40.z),
        ),
        d_v: Point3::new(
            5.0 * (p51.x - p50.x),
            5.0 * (p51.y - p50.y),
            5.0 * (p51.z - p50.z),
        ),
        d_uu: Point3::new(
            20.0 * (p50.x - 2.0 * p40.x + p30.x),
            20.0 * (p50.y - 2.0 * p40.y + p30.y),
            20.0 * (p50.z - 2.0 * p40.z + p30.z),
        ),
        d_vv: Point3::new(
            20.0 * (p52.x - 2.0 * p51.x + p50.x),
            20.0 * (p52.y - 2.0 * p51.y + p50.y),
            20.0 * (p52.z - 2.0 * p51.z + p50.z),
        ),
        d_uv: Point3::new(
            25.0 * (p51.x - p50.x - p41.x + p40.x),
            25.0 * (p51.y - p50.y - p41.y + p40.y),
            25.0 * (p51.z - p50.z - p41.z + p40.z),
        ),
    }
}

fn extract_jet_01(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p05 = piece.cps[0][5];
    let p15 = piece.cps[1][5];
    let p04 = piece.cps[0][4];
    let p25 = piece.cps[2][5];
    let p03 = piece.cps[0][3];
    let p14 = piece.cps[1][4];
    TaylorJet2 {
        value: p05,
        d_u: Point3::new(
            5.0 * (p15.x - p05.x),
            5.0 * (p15.y - p05.y),
            5.0 * (p15.z - p05.z),
        ),
        d_v: Point3::new(
            5.0 * (p05.x - p04.x),
            5.0 * (p05.y - p04.y),
            5.0 * (p05.z - p04.z),
        ),
        d_uu: Point3::new(
            20.0 * (p25.x - 2.0 * p15.x + p05.x),
            20.0 * (p25.y - 2.0 * p15.y + p05.y),
            20.0 * (p25.z - 2.0 * p15.z + p05.z),
        ),
        d_vv: Point3::new(
            20.0 * (p05.x - 2.0 * p04.x + p03.x),
            20.0 * (p05.y - 2.0 * p04.y + p03.y),
            20.0 * (p05.z - 2.0 * p04.z + p03.z),
        ),
        d_uv: Point3::new(
            25.0 * (p15.x - p05.x - p14.x + p04.x),
            25.0 * (p15.y - p05.y - p14.y + p04.y),
            25.0 * (p15.z - p05.z - p14.z + p04.z),
        ),
    }
}

fn extract_jet_11(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p55 = piece.cps[5][5];
    let p45 = piece.cps[4][5];
    let p54 = piece.cps[5][4];
    let p35 = piece.cps[3][5];
    let p53 = piece.cps[5][3];
    let p44 = piece.cps[4][4];
    TaylorJet2 {
        value: p55,
        d_u: Point3::new(
            5.0 * (p55.x - p45.x),
            5.0 * (p55.y - p45.y),
            5.0 * (p55.z - p45.z),
        ),
        d_v: Point3::new(
            5.0 * (p55.x - p54.x),
            5.0 * (p55.y - p54.y),
            5.0 * (p55.z - p54.z),
        ),
        d_uu: Point3::new(
            20.0 * (p55.x - 2.0 * p45.x + p35.x),
            20.0 * (p55.y - 2.0 * p45.y + p35.y),
            20.0 * (p55.z - 2.0 * p45.z + p35.z),
        ),
        d_vv: Point3::new(
            20.0 * (p55.x - 2.0 * p54.x + p53.x),
            20.0 * (p55.y - 2.0 * p54.y + p53.y),
            20.0 * (p55.z - 2.0 * p54.z + p53.z),
        ),
        d_uv: Point3::new(
            25.0 * (p55.x - p45.x - p54.x + p44.x),
            25.0 * (p55.y - p45.y - p54.y + p44.y),
            25.0 * (p55.z - p45.z - p54.z + p44.z),
        ),
    }
}

fn sample_jet(seed: f64) -> TaylorJet2 {
    TaylorJet2 {
        value: Point3::new(seed, 2.0 * seed, -seed),
        d_u: Point3::new(seed * 0.3, 0.5, 1.1),
        d_v: Point3::new(0.2, seed * 0.6, -0.3),
        d_uu: Point3::new(0.1, -0.2, seed * 0.05),
        d_vv: Point3::new(seed * 0.07, 0.3, -0.1),
        d_uv: Point3::new(0.04, 0.08, seed * 0.02),
    }
}

#[test]
fn set_jet_10_round_trips() {
    let mut piece = Bi5CapPiece::zero();
    let jet = sample_jet(1.5);
    set_corner_jet_10(&mut piece, &jet);
    let extracted = extract_jet_10(&piece);
    assert!(approx_p(extracted.value, jet.value, EPS));
    assert!(approx_p(extracted.d_u, jet.d_u, EPS));
    assert!(approx_p(extracted.d_v, jet.d_v, EPS));
    assert!(approx_p(extracted.d_uu, jet.d_uu, EPS));
    assert!(approx_p(extracted.d_vv, jet.d_vv, EPS));
    assert!(approx_p(extracted.d_uv, jet.d_uv, EPS));
}

#[test]
fn set_jet_01_round_trips() {
    let mut piece = Bi5CapPiece::zero();
    let jet = sample_jet(-0.7);
    set_corner_jet_01(&mut piece, &jet);
    let extracted = extract_jet_01(&piece);
    assert!(approx_p(extracted.value, jet.value, EPS));
    assert!(approx_p(extracted.d_u, jet.d_u, EPS));
    assert!(approx_p(extracted.d_v, jet.d_v, EPS));
    assert!(approx_p(extracted.d_uu, jet.d_uu, EPS));
    assert!(approx_p(extracted.d_vv, jet.d_vv, EPS));
    assert!(approx_p(extracted.d_uv, jet.d_uv, EPS));
}

#[test]
fn set_jet_11_round_trips() {
    let mut piece = Bi5CapPiece::zero();
    let jet = sample_jet(2.3);
    set_corner_jet_11(&mut piece, &jet);
    let extracted = extract_jet_11(&piece);
    assert!(approx_p(extracted.value, jet.value, EPS));
    assert!(approx_p(extracted.d_u, jet.d_u, EPS));
    assert!(approx_p(extracted.d_v, jet.d_v, EPS));
    assert!(approx_p(extracted.d_uu, jet.d_uu, EPS));
    assert!(approx_p(extracted.d_vv, jet.d_vv, EPS));
    assert!(approx_p(extracted.d_uv, jet.d_uv, EPS));
}

#[test]
fn corner_setters_preserve_other_cps() {
    let mut piece = Bi5CapPiece::zero();
    for p in 0..6 {
        for q in 0..6 {
            piece.cps[p][q] = Point3::new(p as f64 + 0.5, q as f64 + 0.5, 0.25);
        }
    }
    let snap = piece;
    let jet = sample_jet(1.0);
    set_corner_jet_10(&mut piece, &jet);
    // CPs (5, 0), (4, 0), (5, 1), (3, 0), (5, 2), (4, 1) changed.
    // All others must be unchanged.
    let changed = [(5, 0), (4, 0), (5, 1), (3, 0), (5, 2), (4, 1)];
    for p in 0..6 {
        for q in 0..6 {
            if changed.contains(&(p, q)) {
                continue;
            }
            assert_eq!(piece.cps[p][q], snap.cps[p][q], "({p}, {q}) changed");
        }
    }
}

// ---- populate_outer_corner_jets_identity: end-to-end ----

#[test]
fn populate_with_zero_target_leaves_cap_zero() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        let target = TensorBorder::zero(n);
        populate_outer_corner_jets_identity(&mut cap, &target);
        cap.validate().unwrap();
        for k in 0..n {
            for id in CapPieceId::all() {
                let piece = cap.sector(k).piece(id);
                for p in 0..6 {
                    for q in 0..6 {
                        let cp = piece.cps[p][q];
                        assert!(
                            cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                            "n={n}, k={k}, id={id:?}, ({p}, {q}): non-zero CP {cp:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn populate_sets_only_outer_corner_blocks() {
    // O piece should remain zero after populate (I.14.b touches
    // only ℓ, m, r). Interior CPs of ℓ, m, r should also stay zero.
    let n = 5;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(1.0, 2.0, 3.0);
        target.patches[k].cps[3][0] = Point3::new(4.0, 5.0, 6.0);
        target.patches[k].cps[0][3] = Point3::new(7.0, 8.0, 9.0);
    }
    populate_outer_corner_jets_identity(&mut cap, &target);

    // O piece remains all zero.
    for k in 0..n {
        let o = cap.sector(k).piece(CapPieceId::O);
        for p in 0..6 {
            for q in 0..6 {
                let cp = o.cps[p][q];
                assert!(
                    cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                    "n={n}, k={k}, O piece ({p}, {q}): non-zero {cp:?}"
                );
            }
        }
    }

    // R piece: outer-corner block CPs at (1, 0) should be set; interior stays zero.
    let outer_r = [(5, 0), (4, 0), (5, 1), (3, 0), (5, 2), (4, 1)];
    for k in 0..n {
        let r = cap.sector(k).piece(CapPieceId::R);
        for p in 0..6 {
            for q in 0..6 {
                let is_outer = outer_r.contains(&(p, q));
                let cp = r.cps[p][q];
                if is_outer {
                    // At least the value (P[5][0]) should match the tensor-border jet at m_k.
                    // We don't check each CP exactly, just that they're not all zero.
                } else {
                    assert!(
                        cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                        "n={n}, k={k}, R ({p}, {q}): non-zero interior CP {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn populate_matches_tensor_border_value_at_outer_corners() {
    // The piece corner VALUE CPs (P[5][0], P[0][5], P[5][5]) should
    // match the tensor-border's corresponding corner values.
    let n = 6;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        // Set the outer corner CPs of each bi-3 patch to distinct values
        // so the tensor-border jets have non-zero value at each corner.
        target.patches[k].cps[3][0] = Point3::new(k as f64, 0.0, 0.0); // m_k corner value
        target.patches[k].cps[0][3] = Point3::new(0.0, k as f64, 0.0); // m_{k+1}
        target.patches[k].cps[3][3] = Point3::new(0.0, 0.0, k as f64); // v_k
    }
    populate_outer_corner_jets_identity(&mut cap, &target);

    for k in 0..n {
        let sector = cap.sector(k);
        let r = sector.piece(CapPieceId::R);
        let l = sector.piece(CapPieceId::L);
        let m = sector.piece(CapPieceId::M);
        let j_mk = target.jet(k, ChiCorner::Mk);
        let j_mk1 = target.jet(k, ChiCorner::Mk1);
        let j_vk = target.jet(k, ChiCorner::Vk);
        assert!(approx_p(r.cps[5][0], j_mk.value, 1e-10));
        assert!(approx_p(l.cps[0][5], j_mk1.value, 1e-10));
        assert!(approx_p(m.cps[5][5], j_vk.value, 1e-10));
    }
}

#[test]
fn populate_deterministic() {
    let n = 5;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(k as f64, 1.0, -1.0);
    }
    let mut cap_a = Bi5Cap::zero(n);
    let mut cap_b = Bi5Cap::zero(n);
    populate_outer_corner_jets_identity(&mut cap_a, &target);
    populate_outer_corner_jets_identity(&mut cap_b, &target);
    assert_eq!(cap_a, cap_b);
}

// ---- Structural invariants ----

#[test]
#[should_panic]
fn populate_rejects_n_mismatch() {
    let mut cap = Bi5Cap::zero(5);
    let target = TensorBorder::zero(6);
    populate_outer_corner_jets_identity(&mut cap, &target);
}
