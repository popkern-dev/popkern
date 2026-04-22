//! Tests for the tri-3 `g̃^△` star-jet extractor.

use super::*;

const EPS: f64 = 1e-10;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

fn approx_jet(a: &TaylorJet2, b: &TaylorJet2, eps: f64) -> bool {
    approx_p(a.value, b.value, eps)
        && approx_p(a.d_u, b.d_u, eps)
        && approx_p(a.d_v, b.d_v, eps)
        && approx_p(a.d_uu, b.d_uu, eps)
        && approx_p(a.d_uv, b.d_uv, eps)
        && approx_p(a.d_vv, b.d_vv, eps)
}

// ---- tri3_star_jet ----

#[test]
fn zero_cps_gives_zero_jet() {
    let s = Gtri3SixCps::zero();
    let jet = tri3_star_jet(&s);
    let zero = Point3::new(0.0, 0.0, 0.0);
    assert_eq!(jet.value, zero);
    assert_eq!(jet.d_u, zero);
    assert_eq!(jet.d_v, zero);
    assert_eq!(jet.d_uu, zero);
    assert_eq!(jet.d_uv, zero);
    assert_eq!(jet.d_vv, zero);
}

#[test]
fn star_value_equals_c_300() {
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(7.0, 8.0, 9.0);
    let jet = tri3_star_jet(&s);
    assert_eq!(jet.value, s.cps[0]);
}

#[test]
fn partial_u_formula_matches() {
    // ∂_u = 3 (c_210 - c_300).
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(1.0, 0.0, 0.0); // c_300
    s.cps[3] = Point3::new(2.0, 1.0, 0.0); // c_210
    let jet = tri3_star_jet(&s);
    let expected = Point3::new(3.0 * (2.0 - 1.0), 3.0 * 1.0, 0.0);
    assert_eq!(jet.d_u, expected);
}

#[test]
fn partial_v_formula_matches() {
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(0.0, 0.0, 0.0);
    s.cps[1] = Point3::new(0.0, 1.0, 0.0); // c_201
    let jet = tri3_star_jet(&s);
    let expected = Point3::new(0.0, 3.0, 0.0);
    assert_eq!(jet.d_v, expected);
}

#[test]
fn second_partials_formulas_match() {
    // ∂²_u = 6 (c_120 - 2 c_210 + c_300)
    // ∂²_v = 6 (c_102 - 2 c_201 + c_300)
    // ∂_u∂_v = 6 (c_111 - c_210 - c_201 + c_300)
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(1.0, 0.0, 0.0); // c_300
    s.cps[3] = Point3::new(0.0, 0.0, 0.0); // c_210
    s.cps[5] = Point3::new(0.0, 0.0, 0.0); // c_120
    let jet = tri3_star_jet(&s);
    // ∂²_u = 6 (0 - 0 + 1) = 6 in x direction.
    assert_eq!(jet.d_uu, Point3::new(6.0, 0.0, 0.0));

    let mut s2 = Gtri3SixCps::zero();
    s2.cps[0] = Point3::new(0.0, 1.0, 0.0); // c_300
    s2.cps[3] = Point3::new(0.0, 0.0, 0.0); // c_210
    s2.cps[1] = Point3::new(0.0, 0.0, 0.0); // c_201
    s2.cps[4] = Point3::new(0.0, 2.0, 0.0); // c_111
    let jet2 = tri3_star_jet(&s2);
    // ∂_u∂_v = 6 (2 - 0 - 0 + 1) = 18 in y direction.
    assert_eq!(jet2.d_uv, Point3::new(0.0, 18.0, 0.0));
}

// ---- rotate_star_jet_to_sector ----

#[test]
fn rotation_by_zero_is_identity() {
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(1.0, 2.0, 3.0);
    s.cps[3] = Point3::new(4.0, 5.0, 6.0);
    let jet = tri3_star_jet(&s);
    let rotated = rotate_star_jet_to_sector(&jet, 0, 5);
    assert!(approx_jet(&jet, &rotated, 1e-14));
}

#[test]
fn rotation_preserves_value() {
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(5.0, 6.0, 7.0);
    s.cps[3] = Point3::new(1.0, 0.0, 0.0);
    let jet = tri3_star_jet(&s);
    for k in 0..7 {
        let rotated = rotate_star_jet_to_sector(&jet, k, 7);
        assert_eq!(rotated.value, jet.value);
    }
}

#[test]
fn rotation_preserves_partial_norms() {
    // |∂_u|² + |∂_v|² is rotation-invariant.
    let mut s = Gtri3SixCps::zero();
    s.cps[3] = Point3::new(1.0, 2.0, 0.0);
    s.cps[1] = Point3::new(3.0, 0.0, 1.0);
    let jet = tri3_star_jet(&s);
    let norm2 = |p: Point3, q: Point3| -> f64 {
        p.x * p.x + p.y * p.y + p.z * p.z + q.x * q.x + q.y * q.y + q.z * q.z
    };
    let original = norm2(jet.d_u, jet.d_v);
    for k in 1..5 {
        let rotated = rotate_star_jet_to_sector(&jet, k, 5);
        let rotated_norm = norm2(rotated.d_u, rotated.d_v);
        assert!(
            (rotated_norm - original).abs() < 1e-10,
            "k={k}: |∂|² not preserved: {original} vs {rotated_norm}"
        );
    }
}

#[test]
fn rotation_by_full_cycle_returns_identity() {
    // n identical rotations of 2π/n each = 2π = identity.
    let mut s = Gtri3SixCps::zero();
    s.cps[0] = Point3::new(1.0, 2.0, 3.0);
    s.cps[3] = Point3::new(0.5, 0.0, 0.0);
    s.cps[1] = Point3::new(0.0, 0.5, 0.0);
    let jet = tri3_star_jet(&s);
    let rotated_n = rotate_star_jet_to_sector(&jet, 5, 5);
    assert!(approx_jet(&jet, &rotated_n, 1e-10));
}

#[test]
fn rotation_by_half_on_n_equal_2_inverts_partials() {
    // At n=2, k=1, θ = π. ∂_u → -∂_u, ∂_v → -∂_v.
    // (cos π = -1, sin π = 0.)
    let mut s = Gtri3SixCps::zero();
    s.cps[3] = Point3::new(1.0, 0.0, 0.0);
    s.cps[1] = Point3::new(0.0, 1.0, 0.0);
    let jet = tri3_star_jet(&s);
    let rotated = rotate_star_jet_to_sector(&jet, 1, 2);
    let expected_du = Point3::new(-jet.d_u.x, -jet.d_u.y, -jet.d_u.z);
    let expected_dv = Point3::new(-jet.d_v.x, -jet.d_v.y, -jet.d_v.z);
    assert!(approx_p(rotated.d_u, expected_du, 1e-10));
    assert!(approx_p(rotated.d_v, expected_dv, 1e-10));
    // ∂²'s: c² = 1, s² = 0, cs = 0 → unchanged.
    assert!(approx_p(rotated.d_uu, jet.d_uu, 1e-10));
    assert!(approx_p(rotated.d_vv, jet.d_vv, 1e-10));
    assert!(approx_p(rotated.d_uv, jet.d_uv, 1e-10));
}

#[test]
fn rotation_second_partials_trace_invariant() {
    // trace(H) = ∂²_u + ∂²_v is rotation-invariant.
    let mut s = Gtri3SixCps::zero();
    s.cps[5] = Point3::new(1.0, 0.0, 0.0);
    s.cps[3] = Point3::new(0.0, 1.0, 0.0);
    s.cps[4] = Point3::new(0.0, 0.0, 1.0);
    let jet = tri3_star_jet(&s);
    let trace_orig = Point3::new(
        jet.d_uu.x + jet.d_vv.x,
        jet.d_uu.y + jet.d_vv.y,
        jet.d_uu.z + jet.d_vv.z,
    );
    for k in 0..5 {
        let rotated = rotate_star_jet_to_sector(&jet, k, 5);
        let trace_rot = Point3::new(
            rotated.d_uu.x + rotated.d_vv.x,
            rotated.d_uu.y + rotated.d_vv.y,
            rotated.d_uu.z + rotated.d_vv.z,
        );
        assert!(
            approx_p(trace_orig, trace_rot, 1e-10),
            "k={k}: trace not invariant"
        );
    }
}

// ---- scale_jet_for_piece_local ----

#[test]
fn scale_jet_identity_at_s_one() {
    let jet = TaylorJet2 {
        value: Point3::new(1.0, 2.0, 3.0),
        d_u: Point3::new(4.0, 5.0, 6.0),
        d_v: Point3::new(7.0, 8.0, 9.0),
        d_uu: Point3::new(10.0, 11.0, 12.0),
        d_uv: Point3::new(13.0, 14.0, 15.0),
        d_vv: Point3::new(16.0, 17.0, 18.0),
    };
    let scaled = scale_jet_for_piece_local(&jet, 1.0);
    assert!(approx_jet(&jet, &scaled, 1e-14));
}

#[test]
fn scale_jet_halves_partials_at_s_half() {
    let jet = TaylorJet2 {
        value: Point3::new(0.0, 0.0, 0.0),
        d_u: Point3::new(1.0, 0.0, 0.0),
        d_v: Point3::new(0.0, 1.0, 0.0),
        d_uu: Point3::new(0.0, 0.0, 1.0),
        d_uv: Point3::new(1.0, 1.0, 0.0),
        d_vv: Point3::new(1.0, 0.0, 1.0),
    };
    let scaled = scale_jet_for_piece_local(&jet, 0.5);
    assert_eq!(scaled.value, jet.value);
    assert_eq!(scaled.d_u.x, 0.5);
    assert_eq!(scaled.d_v.y, 0.5);
    assert_eq!(scaled.d_uu.z, 0.25);
    assert_eq!(scaled.d_uv.x, 0.25);
    assert_eq!(scaled.d_vv.x, 0.25);
}
