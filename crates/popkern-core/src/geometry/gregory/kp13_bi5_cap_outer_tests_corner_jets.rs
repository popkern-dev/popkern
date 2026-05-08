use super::tests_support::*;
use super::*;

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
