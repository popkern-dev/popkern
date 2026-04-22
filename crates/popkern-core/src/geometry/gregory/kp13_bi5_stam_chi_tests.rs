//! Tests for the Stam-1998 CC characteristic-map evaluator.

use super::super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner};
use super::*;

const EPS: f64 = 1e-10;

fn approx_pair(a: (f64, f64), b: (f64, f64), eps: f64) -> bool {
    (a.0 - b.0).abs() < eps && (a.1 - b.1).abs() < eps
}

// ---- Construction ----

#[test]
fn new_shares_lambda_with_catmull_clark() {
    for n in [3, 4, 5, 6, 8] {
        let stam = StamChi::new(n);
        let cc = CatmullClarkChi::new(n);
        assert!(
            (stam.lambda - cc.lambda).abs() < 1e-15,
            "n={n}: stam.lambda {} vs cc.lambda {}",
            stam.lambda,
            cc.lambda
        );
        assert_eq!(stam.n, n);
    }
}

// ---- Corner values (should match CC exactly — corners are in far region) ----

#[test]
fn eval_at_corners_matches_cc() {
    for n in [3, 4, 5, 6, 7, 8] {
        let stam = StamChi::new(n);
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            for &(u, v) in &[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
                let s = stam.eval(k, u, v);
                let c = cc.eval(k, u, v);
                assert!(
                    approx_pair(s, c, EPS),
                    "n={n}, k={k}, (u,v)=({u},{v}): Stam {s:?} vs CC {c:?}"
                );
            }
        }
    }
}

#[test]
fn jet_at_corners_matches_cc() {
    for n in [3, 4, 5, 6, 7, 8] {
        let stam = StamChi::new(n);
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            for corner in [
                ChiCorner::Star,
                ChiCorner::Mk,
                ChiCorner::Mk1,
                ChiCorner::Vk,
            ] {
                let js = stam.jet(k, corner);
                let jc = cc.jet(k, corner);
                assert!(approx_pair(js.pos, jc.pos, EPS));
                assert!(approx_pair(js.d_u, jc.d_u, EPS));
                assert!(approx_pair(js.d_v, jc.d_v, EPS));
                assert!(approx_pair(js.d_uu, jc.d_uu, EPS));
                assert!(approx_pair(js.d_vv, jc.d_vv, EPS));
                assert!(approx_pair(js.d_uv, jc.d_uv, EPS));
            }
        }
    }
}

// ---- Far-region agreement with CC ----

#[test]
fn far_region_eval_matches_cc() {
    // Any (u, v) with max(u, v) > STAR_BOUNDARY should bypass the
    // recursion entirely — StamChi delegates directly.
    let far_points = [
        (0.51, 0.0),
        (0.0, 0.51),
        (0.6, 0.2),
        (0.2, 0.6),
        (0.7, 0.7),
        (0.9, 0.9),
        (0.51, 0.51),
    ];
    for n in [3, 4, 5, 6] {
        let stam = StamChi::new(n);
        let cc = CatmullClarkChi::new(n);
        for k in 0..n {
            for &(u, v) in &far_points {
                let s = stam.eval(k, u, v);
                let c = cc.eval(k, u, v);
                assert!(
                    approx_pair(s, c, EPS),
                    "n={n}, k={k}, (u,v)=({u},{v}): Stam {s:?} vs CC {c:?}"
                );
            }
        }
    }
}

// ---- Self-similarity law ----

#[test]
fn self_similarity_law_holds_exactly() {
    // χ(u, v) = λ · χ(2u, 2v) for max(u, v) ≤ 1/2.
    // This is the Stam 1998 self-similarity identity.
    for n in [3, 5, 6, 8] {
        let stam = StamChi::new(n);
        let near_star = [(0.1, 0.1), (0.25, 0.25), (0.4, 0.3), (0.1, 0.4)];
        for k in 0..n {
            for &(u, v) in &near_star {
                let small = stam.eval(k, u, v);
                let big = stam.eval(k, 2.0 * u, 2.0 * v);
                let scaled = (stam.lambda * big.0, stam.lambda * big.1);
                assert!(
                    approx_pair(small, scaled, EPS),
                    "n={n}, k={k}, (u,v)=({u},{v}): χ(u,v) = {small:?}, λ·χ(2u,2v) = {scaled:?}"
                );
            }
        }
    }
}

#[test]
fn self_similarity_ratio_is_lambda() {
    // At (u, v) = (ε, ε), χ ≈ λ^depth · χ(ε·2^depth, ε·2^depth).
    // So ‖χ(u, v)‖ / ‖χ(2u, 2v)‖ ≈ λ exactly.
    for n in [3, 5, 6, 7] {
        let stam = StamChi::new(n);
        let k = 0;
        for &u in &[0.02_f64, 0.05, 0.1, 0.25, 0.3] {
            let small = stam.eval(k, u, u);
            let big = stam.eval(k, 2.0 * u, 2.0 * u);
            let n_small = (small.0.powi(2) + small.1.powi(2)).sqrt();
            let n_big = (big.0.powi(2) + big.1.powi(2)).sqrt();
            if n_big > 1e-12 {
                let ratio = n_small / n_big;
                assert!(
                    (ratio - stam.lambda).abs() < 1e-9,
                    "n={n}, u={u}: ratio = {ratio}, λ = {}",
                    stam.lambda
                );
            }
        }
    }
}

// ---- Near-star limit ----

#[test]
fn eval_at_exact_origin_is_zero() {
    for n in [3, 4, 5, 6, 7, 8] {
        let stam = StamChi::new(n);
        for k in 0..n {
            let v = stam.eval(k, 0.0, 0.0);
            // After MAX_DEPTH recursions at (0, 0), result is
            // λ^MAX_DEPTH · (far-region eval at (0, 0)), which is
            // below machine epsilon.
            assert!(
                v.0.abs() < 1e-12 && v.1.abs() < 1e-12,
                "n={n}, k={k}: {v:?}"
            );
        }
    }
}

#[test]
fn eval_converges_to_zero_near_origin() {
    // For a sequence (u_k, u_k) → (0, 0), ‖χ‖ → 0. (Not necessarily
    // monotone: the far-region bi-cubic is not monotone in the
    // radial direction, so intermediate values can fluctuate by
    // O(λ^k) between adjacent scales. Only the limit matters.)
    let stam = StamChi::new(5);
    let k = 0;
    for scale in [1e-2_f64, 1e-4, 1e-8, 1e-16] {
        let (x, y) = stam.eval(k, scale, scale);
        let norm = (x * x + y * y).sqrt();
        // Upper bound: ‖far-region χ‖ ≤ ~2 (sum of two unit-circle
        // points), scaled by λ^depth where depth = ceil(log₂(0.5/scale)).
        let depth = StamChi::depth_for(scale, scale);
        let expected_order = stam.lambda.powi(depth as i32) * 2.0;
        assert!(
            norm <= expected_order + 1e-12,
            "scale {scale}: norm {norm} > λ^{depth} · 2 = {expected_order}"
        );
    }
    // At (0, 0) exactly, result is < 1e-12 (by eval_at_exact_origin_is_zero).
}

// ---- depth_for ----

#[test]
fn depth_zero_for_far_region() {
    assert_eq!(StamChi::depth_for(0.6, 0.6), 0);
    assert_eq!(StamChi::depth_for(0.51, 0.1), 0);
    assert_eq!(StamChi::depth_for(1.0, 1.0), 0);
}

#[test]
fn depth_grows_as_uv_shrinks() {
    // max(u, v) = 1/2: depth 0 (exactly on boundary, not in near-star).
    // max(u, v) = 1/4: depth 1.
    // max(u, v) = 1/8: depth 2.
    // max(u, v) = 1/16: depth 3.
    assert_eq!(StamChi::depth_for(0.5, 0.5), 0);
    assert_eq!(StamChi::depth_for(0.25, 0.25), 1);
    assert_eq!(StamChi::depth_for(0.125, 0.125), 2);
    assert_eq!(StamChi::depth_for(0.0625, 0.0625), 3);
}

#[test]
fn depth_capped_at_max() {
    assert_eq!(StamChi::depth_for(0.0, 0.0), MAX_DEPTH);
    // Extremely small but non-zero:
    assert_eq!(StamChi::depth_for(1e-300, 1e-300), MAX_DEPTH);
}

// ---- Determinism ----

#[test]
fn eval_is_deterministic() {
    let stam = StamChi::new(5);
    let a = stam.eval(2, 0.17, 0.31);
    let b = stam.eval(2, 0.17, 0.31);
    assert_eq!(a, b);
}

// ---- n=4 regular: self-similarity combined with CC bi-cubic ----

#[test]
fn n4_near_star_self_consistent() {
    // At n=4, λ=1/2 exactly (standard bi-cubic B-spline subdivision).
    let stam = StamChi::new(4);
    assert!((stam.lambda - 0.5).abs() < 1e-15);
    // χ(u, u) in near-star should satisfy
    //   χ(u, u) = (1/2) · χ(2u, 2u) = (1/2)^k · χ(2^k u, 2^k u)
    // for any k such that 2^k u is still in near-star.
    for &u in &[0.05_f64, 0.1, 0.15, 0.2, 0.24] {
        let small = stam.eval(0, u, u);
        let big = stam.eval(0, 2.0 * u, 2.0 * u);
        let scaled = (0.5 * big.0, 0.5 * big.1);
        assert!(
            approx_pair(small, scaled, EPS),
            "u={u}: small={small:?}, scaled={scaled:?}"
        );
    }
}

// ---- Structural sanity: finite, bounded ----

#[test]
fn eval_always_finite() {
    // Across a sweep of (u, v), values should be finite. Covers
    // near-star recursion, boundary crossings, and far region.
    for n in [3, 4, 5, 6, 7, 8] {
        let stam = StamChi::new(n);
        for k in 0..n {
            for i in 0..=20 {
                for j in 0..=20 {
                    let u = i as f64 / 20.0;
                    let v = j as f64 / 20.0;
                    let (x, y) = stam.eval(k, u, v);
                    assert!(
                        x.is_finite() && y.is_finite(),
                        "n={n}, k={k}, (u,v)=({u},{v}): ({x}, {y}) not finite"
                    );
                    // Rough bound: |χ| should be on the order of
                    // the unit disk radius (m_k at radius 1). We
                    // allow 5x headroom for the bi-cubic overshoot.
                    assert!(
                        x.abs() < 5.0 && y.abs() < 5.0,
                        "n={n}, k={k}, (u,v)=({u},{v}): ({x}, {y}) too large"
                    );
                }
            }
        }
    }
}
