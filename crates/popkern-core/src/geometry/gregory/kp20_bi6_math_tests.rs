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
fn characteristic_scalars_match_kp18_pr08_formulas() {
    let n4 = Kp20CharacteristicScalars::for_valence(4).unwrap();
    assert!((n4.c - 0.0).abs() < EPS);
    assert!((n4.lambda - 0.5).abs() < EPS);
    assert!((n4.a - 0.25).abs() < EPS);
    assert!((n4.b - 0.125).abs() < EPS);

    let n5 = Kp20CharacteristicScalars::for_valence(5).unwrap();
    let c = (core::f64::consts::TAU / 5.0).cos();
    let lambda = (c + 5.0 + ((c + 1.0) * (c + 9.0)).sqrt()) / 16.0;
    assert!((n5.lambda - lambda).abs() < EPS);
    assert_eq!(kp20_cc_subdominant_eigenvalue(5).unwrap(), n5.lambda);
    assert!(Kp20CharacteristicScalars::for_valence(2).is_err());
}

#[test]
fn stam_characteristic_eigenvector_uses_kp18_symmetry() {
    let scalars = Kp20CharacteristicScalars::for_valence(5).unwrap();
    let e = scalars.eigenvector;
    assert_eq!(e[0], e[4]);
    assert_eq!(e[1], e[3]);
    assert_eq!(e[2], e[5]);

    let expected = kp20_stam_characteristic_eigenvector(scalars.lambda, scalars.a, scalars.b);
    for i in 0..6 {
        assert!((e[i] - expected[i]).abs() < EPS);
    }
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
fn refined_s0_free_counts_match_kp20_dimension_claim() {
    for n in [3, 5, 6, 7, 8, 9, 10, 12] {
        assert_eq!(kp20_s0_free_count(n).unwrap(), 7 * n + 10);
        assert_eq!(kp20_refined_s0_free_count(n).unwrap(), 16 * n + 10);
        assert_eq!(kp20_refinement_added_free_count(n).unwrap(), 9 * n);
    }
}

#[test]
fn refined_s0_free_counts_reject_invalid_valence() {
    assert!(kp20_s0_free_count(2).is_err());
    assert!(kp20_refined_s0_free_count(2).is_err());
    assert!(kp20_refinement_added_free_count(2).is_err());
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
fn kp18_c2_join_rule_matches_second_difference_identity() {
    let a_d_minus_2 = 2.5;
    let a_d = -1.25;
    let b_0 = a_d;
    let b_2 = 4.75;
    let a_d_minus_1 = kp18_c2_join_left(a_d_minus_2, a_d, b_2);
    let b_1 = kp18_c2_join_right(a_d_minus_2, b_0, b_2);

    assert!((b_1 - (2.0 * a_d - a_d_minus_1)).abs() < EPS);
    let left_second = a_d - 2.0 * a_d_minus_1 + a_d_minus_2;
    let right_second = b_2 - 2.0 * b_1 + b_0;
    assert!((left_second - right_second).abs() < EPS);
}

#[test]
fn hermite_to_bi5_reproduces_kp18_corner_jet_coefficients() {
    let mut h = [[0.0_f64; 3]; 3];
    for (p, row) in h.iter_mut().enumerate() {
        for (q, hpq) in row.iter_mut().enumerate() {
            *hpq = (10 * p + q) as f64;
        }
    }
    let f = hermite_to_bi5_corner_block(h);
    let c21 = h[0][0]
        + 2.0 * h[1][0] / 5.0
        + h[2][0] / 20.0
        + h[0][1] / 5.0
        + 2.0 * h[1][1] / 25.0
        + h[2][1] / 100.0;
    let c12 = h[0][0]
        + h[1][0] / 5.0
        + 2.0 * h[0][1] / 5.0
        + h[0][2] / 20.0
        + 2.0 * h[1][1] / 25.0
        + h[1][2] / 100.0;
    assert!((f[2][1] - c21).abs() < EPS);
    assert!((f[1][2] - c12).abs() < EPS);
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
