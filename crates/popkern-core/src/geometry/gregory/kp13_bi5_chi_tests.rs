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

// ---- CatmullClarkChi: subdominant eigenvalue ----

#[test]
fn lambda_n4_equals_one_half() {
    // Regular quad CC at n=4: λ = 1/2.
    let l = cc_subdominant_eigenvalue(4);
    assert!((l - 0.5).abs() < EPS);
}

#[test]
fn lambda_monotonic_in_valence() {
    // λ increases monotonically toward 1 as n grows past regular.
    let l5 = cc_subdominant_eigenvalue(5);
    let l6 = cc_subdominant_eigenvalue(6);
    let l7 = cc_subdominant_eigenvalue(7);
    let l8 = cc_subdominant_eigenvalue(8);
    let l12 = cc_subdominant_eigenvalue(12);
    assert!(l5 < l6);
    assert!(l6 < l7);
    assert!(l7 < l8);
    assert!(l8 < l12);
    assert!(l12 < 1.0);
}

#[test]
fn lambda_below_3_panics() {
    for bad in [0usize, 1, 2] {
        let r = std::panic::catch_unwind(|| cc_subdominant_eigenvalue(bad));
        assert!(r.is_err(), "expected panic for n={bad}");
    }
}

// ---- CatmullClarkChi: corner values ----

#[test]
fn cc_corner_positions_match_bilinear() {
    // By construction, CC corner positions equal bilinear at all 4.
    for n in [3, 5, 6, 7, 8] {
        let bil = BilinearChi::new(n);
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            for corner in [
                ChiCorner::Star,
                ChiCorner::Mk,
                ChiCorner::Mk1,
                ChiCorner::Vk,
            ] {
                let b = bil.jet(k, corner).pos;
                let c = cc.jet(k, corner).pos;
                assert!(
                    approx_eq(b.0, c.0) && approx_eq(b.1, c.1),
                    "n={n} k={k} corner={corner:?}: CC vs bilinear position mismatch"
                );
            }
        }
    }
}

// ---- CatmullClarkChi: partials at the star ----

#[test]
fn cc_star_first_partials_scale_by_lambda() {
    for n in [3, 5, 6, 7, 8, 12] {
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let (mk1x, mk1y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
            let j = cc.jet(k, ChiCorner::Star);
            assert!(approx_eq(j.d_u.0, cc.lambda * mkx));
            assert!(approx_eq(j.d_u.1, cc.lambda * mky));
            assert!(approx_eq(j.d_v.0, cc.lambda * mk1x));
            assert!(approx_eq(j.d_v.1, cc.lambda * mk1y));
        }
    }
}

#[test]
fn cc_star_second_partials_scale_by_1_minus_lambda() {
    let n = 7;
    let cc = CatmullClarkChi::new(n);
    let s1 = 1.0 - cc.lambda;
    for k in 0..n {
        let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
        let (mk1x, mk1y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
        let j = cc.jet(k, ChiCorner::Star);
        assert!(approx_eq(j.d_uu.0, 4.0 * s1 * mkx));
        assert!(approx_eq(j.d_uu.1, 4.0 * s1 * mky));
        assert!(approx_eq(j.d_vv.0, 4.0 * s1 * mk1x));
        assert!(approx_eq(j.d_vv.1, 4.0 * s1 * mk1y));
        assert!(approx_eq(j.d_uv.0, 3.0 * s1 * (mkx + mk1x)));
        assert!(approx_eq(j.d_uv.1, 3.0 * s1 * (mky + mk1y)));
    }
}

#[test]
fn cc_non_star_partials_match_bilinear() {
    // At (1, 0), (0, 1), (1, 1), CC's partials equal bilinear's by
    // construction — only the star corner picks up CC behavior.
    let n = 6;
    let bil = BilinearChi::new(n);
    let cc = CatmullClarkChi::new(n);
    for k in 0..n {
        for corner in [ChiCorner::Mk, ChiCorner::Mk1, ChiCorner::Vk] {
            let jb = bil.jet(k, corner);
            let jc = cc.jet(k, corner);
            assert!(approx_eq(jb.d_u.0, jc.d_u.0));
            assert!(approx_eq(jb.d_u.1, jc.d_u.1));
            assert!(approx_eq(jb.d_v.0, jc.d_v.0));
            assert!(approx_eq(jb.d_v.1, jc.d_v.1));
            // Bilinear has zero 2nd partials; CC inherits that.
            assert!(approx_eq(jc.d_uu.0, 0.0));
            assert!(approx_eq(jc.d_vv.0, 0.0));
            assert!(approx_eq(jc.d_uv.0, 0.0));
        }
    }
}

// ---- CatmullClarkChi: eval at corners ----

#[test]
fn cc_eval_at_corners_matches_jet_pos() {
    let n = 7;
    let cc = CatmullClarkChi::new(n);
    for k in 0..n {
        for corner in [
            ChiCorner::Star,
            ChiCorner::Mk,
            ChiCorner::Mk1,
            ChiCorner::Vk,
        ] {
            let (u, v) = corner.uv();
            let direct = cc.eval(k, u, v);
            let jet_pos = cc.jet(k, corner).pos;
            assert!(
                (direct.0 - jet_pos.0).abs() < 1e-12,
                "n={n} k={k} corner={corner:?}: eval vs jet.pos x"
            );
            assert!((direct.1 - jet_pos.1).abs() < 1e-12);
        }
    }
}

#[test]
fn cc_eval_interior_differs_from_bilinear_when_lambda_lt_one() {
    // At n > 4, λ < 1, so the CC-modified CP P_{1,0} changes the
    // bi-cubic. Evaluating at u=0.1 (well near the star) should give
    // a different position from bilinear.
    let n = 6;
    let bil = BilinearChi::new(n);
    let cc = CatmullClarkChi::new(n);
    let diff = {
        let (bx, by) = bil.eval(0, 0.1, 0.0);
        let (cx, cy) = cc.eval(0, 0.1, 0.0);
        (cx - bx).hypot(cy - by)
    };
    assert!(
        diff > 1e-4,
        "CC eval should diverge from bilinear near star; diff = {diff}"
    );
}

// ---- CatmullClarkChi: FD partial check ----

#[test]
fn cc_star_partials_match_fd_of_eval() {
    // FD check: ∂u χ(0, 0) computed via central difference should
    // match the analytic jet value (up to FD error).
    let n = 7;
    let cc = CatmullClarkChi::new(n);
    let h = 1e-5;
    for k in 0..n {
        let jet = cc.jet(k, ChiCorner::Star);
        let (fwd_x, fwd_y) = cc.eval(k, h, 0.0);
        let (bwd_x, bwd_y) = cc.eval(k, -h, 0.0);
        let fd_u_x = (fwd_x - bwd_x) / (2.0 * h);
        let fd_u_y = (fwd_y - bwd_y) / (2.0 * h);
        assert!(
            (fd_u_x - jet.d_u.0).abs() < 1e-6,
            "k={k} ∂u_x fd={fd_u_x} vs analytic={}",
            jet.d_u.0
        );
        assert!((fd_u_y - jet.d_u.1).abs() < 1e-6);
    }
}

// ---- Smoothness: CC reduces to bilinear when λ → 1 ----

#[test]
fn cc_reduces_to_bilinear_in_lambda_one_limit() {
    // Manually construct a CC with λ = 1 — should give bilinear
    // positions AND partials at every corner (2nd partials → 0).
    let cc_limit = CatmullClarkChi { n: 6, lambda: 1.0 };
    let bil = BilinearChi::new(6);
    for k in 0..6 {
        for corner in [
            ChiCorner::Star,
            ChiCorner::Mk,
            ChiCorner::Mk1,
            ChiCorner::Vk,
        ] {
            let jc = cc_limit.jet(k, corner);
            let jb = bil.jet(k, corner);
            assert!(approx_eq(jc.d_u.0, jb.d_u.0));
            assert!(approx_eq(jc.d_v.0, jb.d_v.0));
            assert!(approx_eq(jc.d_uu.0, 0.0));
            assert!(approx_eq(jc.d_vv.0, 0.0));
            assert!(approx_eq(jc.d_uv.0, 0.0));
        }
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
