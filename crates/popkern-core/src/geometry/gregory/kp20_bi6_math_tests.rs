use super::*;

const EPS: f64 = 1e-12;

#[test]
fn tau_table_matches_kp20_eq23() {
    let expected = [
        (3, 185.0),
        (5, 1033.0),
        (6, 352.0),
        (7, -145.0),
        (8, -469.0),
        (9, -687.0),
        (10, -840.0),
    ];
    for (n, kappa) in expected {
        let params = Kp20TauParams::for_valence(n).unwrap();
        assert_eq!(params.kappa, kappa);
        assert!((params.mu0 - (1.0 + kappa / 105.0)).abs() < EPS);
        assert!((params.mu1 - 1.0 / 3.0).abs() < EPS);
        assert!(params.mu2.is_finite());
        assert!(params.mu3.is_finite());
    }
}

#[test]
fn tau_rejects_unpublished_valences() {
    assert!(Kp20TauParams::for_valence(4).is_err());
    assert!(Kp20TauParams::for_valence(11).is_err());
}

#[test]
fn n3_limit_weight_expands_kp15_correction() {
    let (gamma5, gamma6) = kp20_limit_weights(3).unwrap();
    assert!((gamma5 - (1.0 / 24.0 + 1.0 / 96.0)).abs() < EPS);
    assert!((gamma5 - 5.0 / 96.0).abs() < EPS);
    assert!((gamma6 - 1.0 / 6.0).abs() < EPS);
}

#[test]
fn general_limit_weights_match_kp18_eq7() {
    let (gamma5, gamma6) = kp20_limit_weights(5).unwrap();
    assert!((gamma5 - 1.0 / 50.0).abs() < EPS);
    assert!((gamma6 - 4.0 / 50.0).abs() < EPS);
}

#[test]
fn b_coefficients_encode_n3_special_pinning() {
    let c = (core::f64::consts::TAU / 3.0).cos();
    let b = kp20_b_coeffs_near_o(c);
    assert!((b[1] + 1.0).abs() < EPS);
}

#[test]
fn g2_strip_weights_follow_eq13() {
    let (w1, w2, w3) = kp20_g2_strip_weights(2.0);
    assert!((w1 + 0.02).abs() < EPS);
    assert!((w2 - 0.02).abs() < EPS);
    assert!((w3 - 2.0 / 30.0).abs() < EPS);
}

#[test]
fn hermite_to_bi6_reproduces_doc_coefficients() {
    let mut h = [[0.0_f64; 4]; 4];
    for (p, row) in h.iter_mut().enumerate() {
        for (q, hpq) in row.iter_mut().enumerate() {
            *hpq = (10 * p + q) as f64;
        }
    }
    let f = hermite_to_bi6_corner_block(h);
    let c21 =
        h[0][0] + h[1][0] / 3.0 + h[2][0] / 30.0 + h[0][1] / 6.0 + h[1][1] / 18.0 + h[2][1] / 180.0;
    let c12 =
        h[0][0] + h[1][0] / 6.0 + h[0][1] / 3.0 + h[0][2] / 30.0 + h[1][1] / 18.0 + h[1][2] / 180.0;
    assert!((f[2][1] - c21).abs() < EPS);
    assert!((f[1][2] - c12).abs() < EPS);
}
