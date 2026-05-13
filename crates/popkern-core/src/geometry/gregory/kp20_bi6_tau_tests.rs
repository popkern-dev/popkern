use super::*;
use crate::geometry::gregory::{Kp20AffineGuideMap, kp20_sample_corner_block, tri6_flat_index};

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn assert_close(a: Point3, b: Point3) {
    assert!(
        a.approx_eq(&b, 1e-10),
        "expected {a:?} to match {b:?}, distance {}",
        a.distance_to(&b)
    );
}

fn affine_guide() -> Kp20Tri6GuideSector {
    let mut guide = Kp20Tri6GuideSector::zero();
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let idx = tri6_flat_index(6, i, j);
            let u = i as f64 / 6.0;
            let v = j as f64 / 6.0;
            guide.cps[idx] = p(2.0 * u + 3.0 * v + 1.0, -u + 0.5 * v, 4.0);
        }
    }
    guide
}

#[test]
fn tau_map_matches_kp20_eq22_sector_coefficients() {
    let tau = Kp20TauMap::for_sector(5, 0).unwrap();
    let params = Kp20TauParams::for_valence(5).unwrap();
    let alpha = core::f64::consts::TAU / 5.0;
    let next = [alpha.cos(), alpha.sin()];

    assert_eq!(tau.cps[flat_index(3, 0, 0)], [0.0, 0.0]);
    assert_eq!(tau.cps[flat_index(3, 0, 1)], [1.0 / 3.0, 0.0]);
    assert_eq!(tau.cps[flat_index(3, 0, 2)], [2.0 / 3.0, 0.0]);
    assert_eq!(tau.cps[flat_index(3, 0, 3)], [params.mu0, 0.0]);
    assert!((tau.cps[flat_index(3, 3, 0)][0] - params.mu0 * next[0]).abs() < 1e-12);
    assert!((tau.cps[flat_index(3, 3, 0)][1] - params.mu0 * next[1]).abs() < 1e-12);
    assert!((tau.cps[flat_index(3, 1, 1)][0] - params.mu1 * (1.0 + next[0])).abs() < 1e-12);
    assert!((tau.cps[flat_index(3, 1, 1)][1] - params.mu1 * next[1]).abs() < 1e-12);
}

#[test]
fn affine_tau_jet_has_constant_first_partials_and_zero_higher_partials() {
    let tau = Kp20TauMap::affine([0.1, 0.2], [0.5, -0.25], [0.125, 0.75]);

    let jet = tau.jet3(0.3, 0.4);

    assert!((jet.d[0][0][0] - (0.1 + 0.5 * 0.3 + 0.125 * 0.4)).abs() < 1e-12);
    assert!((jet.d[0][0][1] - (0.2 - 0.25 * 0.3 + 0.75 * 0.4)).abs() < 1e-12);
    assert!((jet.d[1][0][0] - 0.5).abs() < 1e-12);
    assert!((jet.d[1][0][1] + 0.25).abs() < 1e-12);
    assert!((jet.d[0][1][0] - 0.125).abs() < 1e-12);
    assert!((jet.d[0][1][1] - 0.75).abs() < 1e-12);
    assert!(jet.d[2][0][0].abs() < 1e-12);
    assert!(jet.d[2][0][1].abs() < 1e-12);
    assert!(jet.d[1][1][0].abs() < 1e-12);
    assert!(jet.d[1][1][1].abs() < 1e-12);
    assert!(jet.d[0][2][0].abs() < 1e-12);
    assert!(jet.d[0][2][1].abs() < 1e-12);
    assert!(jet.d[3][0][0].abs() < 1e-12);
    assert!(jet.d[3][0][1].abs() < 1e-12);
}

#[test]
fn cubic_tau_composition_reduces_to_affine_sampler_for_affine_tau() {
    let guide = affine_guide();
    let tau = Kp20TauMap::affine([0.1, 0.2], [0.5, 0.0], [0.0, 0.25]);
    let affine = Kp20AffineGuideMap::new([0.1, 0.2], [0.5, 0.0], [0.0, 0.25]);

    let via_tau = kp20_sample_corner_block_with_tau(&guide, &tau, 0.0, 0.0);
    let via_affine = kp20_sample_corner_block(&guide, affine);

    for i in 0..4 {
        for j in 0..4 {
            assert_close(via_tau[i][j], via_affine[i][j]);
        }
    }
}

#[test]
fn tau_affine_corner_sampling_precomposes_piece_map() {
    let guide = affine_guide();
    let tau_origin = [0.1, -0.2];
    let tau_du = [0.5, 0.125];
    let tau_dv = [-0.25, 0.75];
    let tau = Kp20TauMap::affine(tau_origin, tau_du, tau_dv);
    let map = Kp20AffineGuideMap::new([0.2, 0.4], [0.5, -0.25], [0.125, 0.25]);

    let composed = Kp20AffineGuideMap::new(
        [
            tau_origin[0] + tau_du[0] * map.origin[0] + tau_dv[0] * map.origin[1],
            tau_origin[1] + tau_du[1] * map.origin[0] + tau_dv[1] * map.origin[1],
        ],
        [
            tau_du[0] * map.du[0] + tau_dv[0] * map.du[1],
            tau_du[1] * map.du[0] + tau_dv[1] * map.du[1],
        ],
        [
            tau_du[0] * map.dv[0] + tau_dv[0] * map.dv[1],
            tau_du[1] * map.dv[0] + tau_dv[1] * map.dv[1],
        ],
    );

    let via_tau = kp20_sample_corner_block_with_tau_affine(&guide, &tau, map);
    let via_affine = kp20_sample_corner_block(&guide, composed);

    for i in 0..4 {
        for j in 0..4 {
            assert_close(via_tau[i][j], via_affine[i][j]);
        }
    }
}

#[test]
fn affine_guide_composed_with_paper_tau_uses_tau_derivatives() {
    let guide = affine_guide();
    let tau = Kp20TauMap::for_sector(5, 1).unwrap();
    let tau_jet = tau.jet3(0.0, 0.0);

    let composed = kp20_compose_guide_tau_jet(&guide, &tau, 0.0, 0.0);

    assert_close(
        composed.h[0][0],
        p(
            2.0 * tau_jet.d[0][0][0] + 3.0 * tau_jet.d[0][0][1] + 1.0,
            -tau_jet.d[0][0][0] + 0.5 * tau_jet.d[0][0][1],
            4.0,
        ),
    );
    assert_close(
        composed.h[1][0],
        p(
            2.0 * tau_jet.d[1][0][0] + 3.0 * tau_jet.d[1][0][1],
            -tau_jet.d[1][0][0] + 0.5 * tau_jet.d[1][0][1],
            0.0,
        ),
    );
    assert_close(
        composed.h[2][1],
        p(
            2.0 * tau_jet.d[2][1][0] + 3.0 * tau_jet.d[2][1][1],
            -tau_jet.d[2][1][0] + 0.5 * tau_jet.d[2][1][1],
            0.0,
        ),
    );
}
