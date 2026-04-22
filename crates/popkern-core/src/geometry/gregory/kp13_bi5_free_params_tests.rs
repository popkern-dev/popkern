//! Tests for Bi5FreeParams + apply_free_params.

use super::super::kp13_bi5_circulant::circulant_residual;
use super::super::kp13_bi5_connect::{ConnectWeights, c1_join, paper_flat};
use super::super::kp13_bi5_gtri5::{IDX_005, IDX_050, IDX_500};
use super::super::kp13_bi5_star_quadratic::{NEAR_STAR_FLAT, StarQuadratic};
use super::*;

const EPS: f64 = 1e-9;

fn approx_p(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- Bi5FreeParams ----

#[test]
fn zero_has_correct_shape() {
    for n in [3, 5, 6, 7, 8, 10] {
        let p = Bi5FreeParams::zero(n);
        assert_eq!(p.valence(), n);
        p.validate().unwrap();
    }
}

#[test]
#[should_panic]
fn zero_below_three_panics() {
    let _ = Bi5FreeParams::zero(2);
}

#[test]
fn validate_detects_length_mismatch() {
    let mut p = Bi5FreeParams::zero(5);
    p.g_30.pop();
    let err = p.validate().unwrap_err();
    assert!(err.contains("g_30"));
}

// ---- apply_free_params: zero params ----

#[test]
fn apply_zero_params_with_zero_star_gives_valid_gtri5() {
    for n in [5, 6, 7, 8, 10] {
        let params = Bi5FreeParams::zero(n);
        let star = Point3::new(0.0, 0.0, 0.0);
        let g = apply_free_params(&params, star).unwrap();
        g.validate().unwrap();
    }
}

#[test]
fn apply_zero_params_rejects_n3_n4() {
    for n in [3, 4] {
        let params = Bi5FreeParams::zero(n);
        let star = Point3::new(0.0, 0.0, 0.0);
        let err = apply_free_params(&params, star).unwrap_err();
        assert!(matches!(err, ApplyError::Circulant(_)));
    }
}

#[test]
fn apply_pins_star_to_cc_limit() {
    // Any params + any cc_limit should produce a g△ whose flat 0 CP
    // is the limit.
    let params = Bi5FreeParams {
        q: StarQuadratic {
            alpha: Point3::new(0.3, 0.0, 0.0),
            beta: Point3::new(0.0, 0.2, 0.0),
            gamma: Point3::new(-0.1, 0.0, 0.0),
            delta: Point3::new(0.5, 0.0, 0.0),
            epsilon: Point3::new(0.0, -0.3, 0.0),
            zeta: Point3::new(999.0, 999.0, 999.0), // overridden
        },
        ..Bi5FreeParams::zero(6)
    };
    let cc = Point3::new(1.5, -0.5, 0.25);
    let g = apply_free_params(&params, cc).unwrap();
    for k in 0..6 {
        assert!(approx_p(g.sector(k).cps[IDX_500], cc));
    }
}

// ---- apply_free_params: Q propagation ----

#[test]
fn apply_propagates_q_to_near_star_across_sectors() {
    let q = StarQuadratic {
        alpha: Point3::new(0.3, 0.0, 0.0),
        beta: Point3::new(0.0, 0.2, 0.0),
        gamma: Point3::new(0.1, 0.0, 0.0),
        delta: Point3::new(0.5, 0.5, 0.0),
        epsilon: Point3::new(-0.1, 0.0, 0.3),
        zeta: Point3::new(0.0, 0.0, 0.0), // will be overridden
    };
    let n = 7;
    let params = Bi5FreeParams {
        q,
        ..Bi5FreeParams::zero(n)
    };
    let cc = Point3::new(1.0, 0.0, 0.5);
    let g = apply_free_params(&params, cc).unwrap();
    // The star-quadratic identity says all sectors' 6 near-star CPs
    // are consistent with a single global Q. We check that sector 0's
    // near-star recovers Q with ζ = cc_limit (not the user's zeta).
    let near0 = [
        g.sector(0).cps[NEAR_STAR_FLAT[0]],
        g.sector(0).cps[NEAR_STAR_FLAT[1]],
        g.sector(0).cps[NEAR_STAR_FLAT[2]],
        g.sector(0).cps[NEAR_STAR_FLAT[3]],
        g.sector(0).cps[NEAR_STAR_FLAT[4]],
        g.sector(0).cps[NEAR_STAR_FLAT[5]],
    ];
    let q_rec = StarQuadratic::from_sector0_near_star(&near0, n);
    assert!(approx_p(q_rec.zeta, cc));
    assert!(approx_p(q_rec.alpha, q.alpha));
    assert!(approx_p(q_rec.beta, q.beta));
    assert!(approx_p(q_rec.gamma, q.gamma));
    assert!(approx_p(q_rec.delta, q.delta));
    assert!(approx_p(q_rec.epsilon, q.epsilon));
}

// ---- apply_free_params: mid-edge and shared-edge invariants ----

#[test]
fn apply_sets_mid_edge_invariant() {
    let n = 6;
    let mk: Vec<Point3> = (0..n)
        .map(|k| {
            let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            Point3::new(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let params = Bi5FreeParams {
        m_k: mk.clone(),
        ..Bi5FreeParams::zero(n)
    };
    let g = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    for (k, &m_k) in mk.iter().enumerate() {
        // Sector k's IDX_050 holds m_k.
        assert!(approx_p(g.sector(k).cps[IDX_050], m_k));
        // Sector (k-1)'s IDX_005 also holds m_k.
        let prev = (k + n - 1) % n;
        assert!(approx_p(g.sector(prev).cps[IDX_005], m_k));
    }
}

#[test]
fn apply_shared_edge_g30_equals_next_sector_g03() {
    // g^i_{30} (flat 15) of sector i should match g^{i+1}_{03} (flat
    // 3) of sector i+1 — the sync_shared_edge_all_depths call inside
    // complete_derivative_layers_with_higher_j enforces this.
    let n = 7;
    let g_30: Vec<Point3> = (0..n)
        .map(|k| Point3::new(k as f64 * 0.1, 0.0, 0.0))
        .collect();
    let params = Bi5FreeParams {
        g_30: g_30.clone(),
        ..Bi5FreeParams::zero(n)
    };
    let g = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    for i in 0..n {
        let next = (i + 1) % n;
        let idx_03 = paper_flat(0, 3);
        let idx_30 = paper_flat(3, 0);
        assert!(approx_p(
            g.sector(i).cps[idx_03],
            g.sector(next).cps[idx_30]
        ));
    }
}

// ---- apply_free_params: circulant + higher-j invariants ----

#[test]
fn apply_satisfies_circulant_three_term_recurrence() {
    // Every valence supported by the circulant (n ≥ 5, pseudo-inverse
    // at n = 6).
    let cc = Point3::new(0.0, 0.0, 0.0);
    for n in [5, 6, 7, 8, 10] {
        let params = {
            let mut p = Bi5FreeParams::zero(n);
            p.q.alpha = Point3::new(0.2, 0.0, 0.0);
            p.q.delta = Point3::new(0.3, 0.0, 0.0);
            for k in 0..n {
                let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                p.g_30[k] = Point3::new(0.1 * a.cos(), 0.0, 0.0);
                p.m_k[k] = Point3::new(a.cos(), a.sin(), 0.0);
            }
            p
        };
        let g = apply_free_params(&params, cc).unwrap();
        for i in 0..n {
            let r = circulant_residual(&g, i);
            // For n != 6, residual should be zero. For n == 6, non-
            // null-mode components should be zero; the null-mode
            // component may be non-zero but here the input RHS is
            // constructed from symmetric rotational data which has
            // no (−1)^k component, so residual should be zero.
            assert!(
                r.x.abs() < 1e-7 && r.y.abs() < 1e-7 && r.z.abs() < 1e-7,
                "n={n} sector {i}: residual {r:?}"
            );
        }
    }
}

#[test]
fn apply_satisfies_higher_j_c1_joins() {
    // After apply, eq (3) at j ∈ {0..=4} should hold across every
    // shared edge.
    let n = 7;
    let params = {
        let mut p = Bi5FreeParams::zero(n);
        p.q.alpha = Point3::new(0.3, 0.0, 0.0);
        p.q.beta = Point3::new(-0.2, 0.0, 0.0);
        for k in 0..n {
            let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            p.g_30[k] = Point3::new(0.1 * a.cos(), 0.0, 0.0);
            p.g_40[k] = Point3::new(0.2 * a.sin(), 0.0, 0.0);
            p.m_k[k] = Point3::new(a.cos(), a.sin(), 0.0);
        }
        p
    };
    let g = apply_free_params(&params, Point3::new(0.0, 0.0, 0.0)).unwrap();
    let w = ConnectWeights::for_valence(n);
    for i in 0..n {
        let next = (i + 1) % n;
        for j in 0..=4 {
            let predicted = c1_join(&w, g.sector(i), j);
            let actual = g.sector(next).cp(j, 1);
            let dx = (predicted.x - actual.x).abs();
            let dy = (predicted.y - actual.y).abs();
            let dz = (predicted.z - actual.z).abs();
            assert!(
                dx < 1e-8 && dy < 1e-8 && dz < 1e-8,
                "C¹ fails at i={i} → {next}, j={j}: predicted {predicted:?} actual {actual:?}"
            );
        }
    }
}

// ---- End-to-end: full free-param set ----

#[test]
fn end_to_end_apply_produces_valid_gtri5_all_valences() {
    // Feed a realistic set of free params and check validate() passes.
    let cc = Point3::new(0.0, 0.0, 1.0);
    for n in [5, 6, 7, 8, 10] {
        let mut params = Bi5FreeParams::zero(n);
        params.q = StarQuadratic {
            alpha: Point3::new(0.3, 0.1, 0.0),
            beta: Point3::new(-0.2, 0.0, 0.0),
            gamma: Point3::new(0.15, 0.0, -0.1),
            delta: Point3::new(0.5, 0.3, 0.0),
            epsilon: Point3::new(-0.1, 0.3, 0.15),
            zeta: Point3::new(999.0, 999.0, 999.0), // overridden
        };
        for k in 0..n {
            let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            params.m_k[k] = Point3::new(a.cos(), a.sin(), 0.0);
            params.g_30[k] = Point3::new(0.2 * a.cos(), 0.0, 0.0);
            params.g_40[k] = Point3::new(0.3 * a.sin(), 0.0, 0.0);
            params.g_122[k] = Point3::new(0.15, 0.1, 0.05);
            params.g_023[k] = Point3::new(0.0, 0.2, 0.0);
            params.g_032[k] = Point3::new(0.1, 0.0, 0.1);
        }
        let g = apply_free_params(&params, cc).unwrap();
        g.validate().unwrap();
        assert_eq!(g.n, n);
    }
}
