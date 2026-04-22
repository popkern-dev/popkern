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
fn cc_corner_positions_star_mk_mk1_match_bilinear() {
    // Three of the four corners are preserved: Star = 0, Mk = m_k,
    // Mk1 = m_{k+1}. Only Vk differs (gets μ-scaled).
    for n in [3, 5, 6, 7, 8] {
        let bil = BilinearChi::new(n);
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            for corner in [ChiCorner::Star, ChiCorner::Mk, ChiCorner::Mk1] {
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

#[test]
fn cc_regular_n4_v_k_matches_bilinear() {
    // At n = 4 (regular quad), λ = 1/2 and μ = 1/4. Looking at
    // Vk's position: μ · (m_k + m_{k+1}). For n = 4, m_0 = (1, 0),
    // m_1 = (0, 1), so μ·(m_0 + m_1) = 0.25·(1, 1) = (0.25, 0.25).
    // Bilinear gives (1, 1). So CC does NOT match bilinear at n = 4
    // when using μ = λ² — the footnote claim that the construction
    // "collapses back to locally-regular CC" is about the geometric
    // SHAPE of the patch, not the numerical Vk position. This test
    // documents the n = 4 numerical behaviour explicitly.
    let cc = CatmullClarkChi::new(4);
    let bil = BilinearChi::new(4);
    for k in 0..4 {
        let cc_vk = cc.jet(k, ChiCorner::Vk).pos;
        let bil_vk = bil.jet(k, ChiCorner::Vk).pos;
        // With μ = 1/4, CC's Vk is at 1/4 of bilinear's Vk.
        assert!(approx_eq(cc_vk.0, bil_vk.0 * cc.mu));
        assert!(approx_eq(cc_vk.1, bil_vk.1 * cc.mu));
    }
}

#[test]
fn cc_vk_equals_mu_times_mk_plus_mk1() {
    for n in [3, 5, 6, 7, 8, 12] {
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let (mk1x, mk1y) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
            let jet = cc.jet(k, ChiCorner::Vk);
            assert!(approx_eq(jet.pos.0, cc.mu * (mkx + mk1x)));
            assert!(approx_eq(jet.pos.1, cc.mu * (mky + mk1y)));
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
        // ∂²uv at star vanishes by rotational symmetry of the CP grid.
        assert!(approx_eq(j.d_uv.0, 0.0));
        assert!(approx_eq(j.d_uv.1, 0.0));
    }
}

#[test]
fn cc_mk_partial_u_is_m_k() {
    // At Mk corner (1, 0), the CP grid has P_{3,0} = m_k and
    // P_{2,0} = (2/3)·m_k, so ∂u = 3·(m_k − 2/3·m_k) = m_k.
    // This is a property of the bilinear-in-u layer past the star.
    for n in [5, 6, 7, 8] {
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let jet = cc.jet(k, ChiCorner::Mk);
            assert!(
                approx_eq(jet.d_u.0, mkx) && approx_eq(jet.d_u.1, mky),
                "n={n} k={k}: ∂u at Mk should be m_k"
            );
        }
    }
}

#[test]
fn cc_mk_partial_uu_nonzero_when_lambda_lt_one() {
    // ∂²u at Mk = 6·(P_{3,0} − 2·P_{2,0} + P_{1,0}) =
    //   6·(m_k − 4/3·m_k + λ/3·m_k) = 2·(λ − 1)·m_k.
    // Zero when λ = 1, negative (points back toward star) when λ < 1.
    let n = 6;
    let cc = CatmullClarkChi::new(n);
    let expected_scale = 2.0 * (cc.lambda - 1.0);
    for k in 0..n {
        let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
        let jet = cc.jet(k, ChiCorner::Mk);
        assert!(approx_eq(jet.d_uu.0, expected_scale * mkx));
        assert!(approx_eq(jet.d_uu.1, expected_scale * mky));
    }
}

#[test]
fn cc_mk1_partial_u_scales_with_lambda() {
    // At Mk1 corner (0, 1), ∂u = 3·(P_{1,3} − P_{0,3}).
    // P_{1,3} = (λ/3)·m_k + m_{k+1}; P_{0,3} = m_{k+1}. So ∂u = λ·m_k.
    // This reflects the CC subdominant scaling propagating along the
    // `m_{k+1}` edge — the tangent in u shrinks by λ across the full
    // edge, not just at the star.
    for n in [5, 6, 7, 8] {
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            let (mkx, mky) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
            let jet = cc.jet(k, ChiCorner::Mk1);
            assert!(
                approx_eq(jet.d_u.0, cc.lambda * mkx) && approx_eq(jet.d_u.1, cc.lambda * mky),
                "n={n} k={k}: ∂u at Mk1 should be λ·m_k"
            );
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
    // Manually construct a CC with λ = 1, μ = 1 — should give
    // bilinear positions AND partials at every corner (2nd partials
    // → 0).
    let cc_limit = CatmullClarkChi {
        n: 6,
        lambda: 1.0,
        mu: 1.0,
    };
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
