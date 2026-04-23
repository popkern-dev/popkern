//! Tests for [`pns_sigma`] — paper-default σ and construction
//! invariants.

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-12;

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

#[test]
fn d1_d2_paper_formulas() {
    // d_1(n) = (0.24n − 0.75) / (n − 3.16)
    for &n in &[5_usize, 6, 7, 8, 9, 10] {
        let expected = (0.24 * (n as f64) - 0.75) / ((n as f64) - 3.16);
        assert!(approx(d1_def(n), expected, EPS), "d_1({n})");
    }
    // d_2(n) = (0.44n − 0.75) / ((n − 4) + 2.02√2)
    for &n in &[5_usize, 6, 7, 8, 9, 10] {
        let expected =
            (0.44 * (n as f64) - 0.75) / ((n as f64 - 4.0) + 2.02 * std::f64::consts::SQRT_2);
        assert!(approx(d2_def(n), expected, EPS), "d_2({n})");
    }
}

#[test]
fn sigma_00_is_origin() {
    for &n in &[5_usize, 6, 7, 8, 10] {
        let jet = sigma_center_jet_default(n);
        assert!(approx(jet.s00[0], 0.0, EPS));
        assert!(approx(jet.s00[1], 0.0, EPS));
    }
}

#[test]
fn sigma_10_x_is_4_times_d1() {
    // ∂_u σ(0, 0) = 4·σ_{10} = 4·(d_1, 0).
    for &n in &[5_usize, 6, 7, 8, 10] {
        let jet = sigma_center_jet_default(n);
        let d1 = d1_def(n);
        assert!(
            approx(jet.ds_u[0], 4.0 * d1, EPS),
            "∂_u σ_x at n={n}: got {}, expected {}",
            jet.ds_u[0],
            4.0 * d1
        );
        assert!(approx(jet.ds_u[1], 0.0, EPS));
    }
}

#[test]
fn sigma_01_is_d1_rotated_by_alpha() {
    // σ_{01} = (d_1 cos α, d_1 sin α), and ∂_v σ(0, 0) = 4·σ_{01}.
    for &n in &[5_usize, 6, 7, 8, 10] {
        let jet = sigma_center_jet_default(n);
        let d1 = d1_def(n);
        let alpha = 2.0 * std::f64::consts::PI / (n as f64);
        let expected_x = 4.0 * d1 * alpha.cos();
        let expected_y = 4.0 * d1 * alpha.sin();
        assert!(approx(jet.ds_v[0], expected_x, 1e-12), "∂_v σ_x at n={n}");
        assert!(approx(jet.ds_v[1], expected_y, 1e-12), "∂_v σ_y at n={n}");
    }
}

#[test]
fn sigma_jet_has_rotational_symmetry_in_u_v() {
    // Under flip u ↔ v, the σ jet should transform according to
    // reflection across the bisector (angle α/2). Equivalently:
    //   |ds_u| = |ds_v|  (same magnitude, different direction)
    //   |ds_uu| = |ds_vv|
    //   angle(ds_v) − angle(ds_u) = α  (ds_v = ds_u rotated by α)
    for &n in &[5_usize, 6, 7, 8, 10] {
        let jet = sigma_center_jet_default(n);
        let mag_u = (jet.ds_u[0].powi(2) + jet.ds_u[1].powi(2)).sqrt();
        let mag_v = (jet.ds_v[0].powi(2) + jet.ds_v[1].powi(2)).sqrt();
        assert!(
            approx(mag_u, mag_v, 1e-12),
            "magnitudes at n={n}: |ds_u|={mag_u}, |ds_v|={mag_v}"
        );
        let mag_uu = (jet.ds_uu[0].powi(2) + jet.ds_uu[1].powi(2)).sqrt();
        let mag_vv = (jet.ds_vv[0].powi(2) + jet.ds_vv[1].powi(2)).sqrt();
        assert!(
            approx(mag_uu, mag_vv, 1e-12),
            "2nd-order magnitudes at n={n}"
        );
    }
}

#[test]
fn sigma_uv_lies_on_bisector_for_symmetric_construction() {
    // σ_{11} is on the bisector by construction. The mixed
    // derivative ∂²_{uv} σ(0, 0) = 16·(σ_{11} − σ_{10} − σ_{01}).
    // σ_{10} + σ_{01} points "symmetric pair" which by the
    // bisector-symmetry reflects into itself — its direction IS
    // the bisector. And σ_{11} is on the bisector. So the
    // difference is also on the bisector (angle α/2).
    for &n in &[5_usize, 6, 7, 8, 10] {
        let jet = sigma_center_jet_default(n);
        let alpha = 2.0 * std::f64::consts::PI / (n as f64);
        let bisector_angle = alpha * 0.5;
        let mixed_angle = jet.ds_uv[1].atan2(jet.ds_uv[0]);
        // Should match bisector_angle up to sign (depending on which
        // sign d_1, d_2 pan out to).
        let diff_mod_pi = ((mixed_angle - bisector_angle).abs())
            .min(((mixed_angle - bisector_angle).abs() - std::f64::consts::PI).abs());
        assert!(
            diff_mod_pi < 1e-10,
            "ds_uv should lie on bisector at n={n}: mixed_angle={mixed_angle}, bisector={bisector_angle}"
        );
    }
}

#[test]
fn explicit_n6_numerics() {
    // Sanity-check numeric output at n=6 against explicit hand calc:
    //   α = 2π/6 = π/3, c = cos(π/3) = 1/2, α/2 = π/6, cos(α/2) = √3/2
    //   d_1 = (0.24·6 − 0.75)/(6 − 3.16) = (1.44 − 0.75)/2.84 = 0.69/2.84 ≈ 0.24296
    //   d_2 = (0.44·6 − 0.75)/((6 − 4) + 2.02·√2) = (2.64 − 0.75)/(2 + 2.857)
    //       ≈ 1.89 / 4.857 ≈ 0.3892
    let jet = sigma_center_jet_default(6);
    let d1 = 0.69 / 2.84;
    let d2 = 1.89 / (2.0 + 2.02 * std::f64::consts::SQRT_2);
    assert!(approx(d1_def(6), d1, 1e-8));
    assert!(approx(d2_def(6), d2, 1e-8));

    // σ_{10} = (d_1, 0), so ∂_u σ(0, 0) = 4·(d_1, 0)
    assert!(approx(jet.ds_u[0], 4.0 * d1, 1e-10));
    assert!(approx(jet.ds_u[1], 0.0, 1e-10));
}
