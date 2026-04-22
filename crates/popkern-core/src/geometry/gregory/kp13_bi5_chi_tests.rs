//! Tests for the bilinear χ characteristic-map approximation.

use super::*;
use std::f64::consts::PI;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

#[test]
#[should_panic]
fn new_panics_below_three() {
    let _ = BilinearChi::new(2);
}

// ---- Corner values ----

#[test]
fn star_corner_maps_to_origin_every_sector() {
    for n in [3, 5, 6, 7, 8] {
        let chi = BilinearChi::new(n);
        for k in 0..n {
            let jet = chi.jet(k, ChiCorner::Star);
            assert!(
                approx_eq(jet.pos.0, 0.0) && approx_eq(jet.pos.1, 0.0),
                "n={n} k={k}: star != origin"
            );
        }
    }
}

#[test]
fn mk_corner_equals_sector_corner_at_k() {
    for n in [3, 5, 6, 7, 8] {
        let chi = BilinearChi::new(n);
        for k in 0..n {
            let (expected_x, expected_y) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let jet = chi.jet(k, ChiCorner::Mk);
            assert!(
                approx_eq(jet.pos.0, expected_x) && approx_eq(jet.pos.1, expected_y),
                "n={n} k={k}: Mk corner"
            );
        }
    }
}

#[test]
fn mk1_corner_equals_sector_corner_at_k_plus_1() {
    for n in [3, 5, 6, 7, 8] {
        let chi = BilinearChi::new(n);
        for k in 0..n {
            let (expected_x, expected_y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
            let jet = chi.jet(k, ChiCorner::Mk1);
            assert!(
                approx_eq(jet.pos.0, expected_x) && approx_eq(jet.pos.1, expected_y),
                "n={n} k={k}: Mk+1 corner"
            );
        }
    }
}

#[test]
fn vk_corner_equals_mk_plus_mk1() {
    let n = 6;
    let chi = BilinearChi::new(n);
    for k in 0..n {
        let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
        let (mk1x, mk1y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
        let jet = chi.jet(k, ChiCorner::Vk);
        assert!(approx_eq(jet.pos.0, mkx + mk1x));
        assert!(approx_eq(jet.pos.1, mky + mk1y));
    }
}

// ---- Derivatives ----

#[test]
fn partials_are_corner_constants() {
    // Bilinear: ∂_u χ = m_k (constant in u, v); ∂_v χ = m_{k+1}.
    for n in [3, 5, 6, 7] {
        let chi = BilinearChi::new(n);
        for k in 0..n {
            let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let (mk1x, mk1y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
            for corner in [
                ChiCorner::Star,
                ChiCorner::Mk,
                ChiCorner::Mk1,
                ChiCorner::Vk,
            ] {
                let jet = chi.jet(k, corner);
                assert!(
                    approx_eq(jet.d_u.0, mkx) && approx_eq(jet.d_u.1, mky),
                    "n={n} k={k} corner={corner:?}: d_u != m_k"
                );
                assert!(
                    approx_eq(jet.d_v.0, mk1x) && approx_eq(jet.d_v.1, mk1y),
                    "n={n} k={k} corner={corner:?}: d_v != m_{{k+1}}"
                );
            }
        }
    }
}

#[test]
fn second_partials_are_zero_bilinear() {
    let chi = BilinearChi::new(7);
    for k in 0..7 {
        for corner in [
            ChiCorner::Star,
            ChiCorner::Mk,
            ChiCorner::Mk1,
            ChiCorner::Vk,
        ] {
            let jet = chi.jet(k, corner);
            assert!(approx_eq(jet.d_uu.0, 0.0) && approx_eq(jet.d_uu.1, 0.0));
            assert!(approx_eq(jet.d_vv.0, 0.0) && approx_eq(jet.d_vv.1, 0.0));
            assert!(approx_eq(jet.d_uv.0, 0.0) && approx_eq(jet.d_uv.1, 0.0));
        }
    }
}

#[test]
fn bilinear_agrees_with_eval() {
    // Sanity: jet.pos at a corner should equal eval at that corner.
    let chi = BilinearChi::new(8);
    for k in 0..8 {
        for corner in [
            ChiCorner::Star,
            ChiCorner::Mk,
            ChiCorner::Mk1,
            ChiCorner::Vk,
        ] {
            let (u, v) = corner.uv();
            let direct = chi.eval(k, u, v);
            let jet = chi.jet(k, corner);
            assert!(approx_eq(direct.0, jet.pos.0));
            assert!(approx_eq(direct.1, jet.pos.1));
        }
    }
}

#[test]
fn n6_regular_chi_is_identity_on_edge_length() {
    // For n = 6, sector_corner gives the regular hexagon vertices,
    // each at distance 1 from origin. Bilinear χ^k at (1, 0) must
    // land on m_k with unit distance.
    let chi = BilinearChi::new(6);
    for k in 0..6 {
        let jet = chi.jet(k, ChiCorner::Mk);
        let dist = (jet.pos.0 * jet.pos.0 + jet.pos.1 * jet.pos.1).sqrt();
        assert!(approx_eq(dist, 1.0), "n=6 k={k}: m_k not at unit radius");
    }
}

#[test]
fn angle_between_adjacent_mk_matches_valence() {
    // Angle between consecutive m_k vertices should be 2π/n.
    let n = 7;
    let chi = BilinearChi::new(n);
    let j0 = chi.jet(0, ChiCorner::Mk).pos;
    let j1 = chi.jet(1, ChiCorner::Mk).pos;
    let dot = j0.0 * j1.0 + j0.1 * j1.1;
    let expected_cos = (2.0 * PI / n as f64).cos();
    assert!(
        (dot - expected_cos).abs() < EPS,
        "angle between m_0 and m_1 doesn't match n={n}"
    );
}
