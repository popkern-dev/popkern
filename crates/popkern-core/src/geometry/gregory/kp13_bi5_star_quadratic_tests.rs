//! Tests for the KP09a §2.3(a) star-point quadratic identity.
//!
//! # Property tests
//!
//! 1. **Polar form at P, P** equals `Q(P)`.
//! 2. **Identity round-trip:** `from_sector0_near_star(Q.apply_to_sector0(...)) == Q`.
//! 3. **Constant quadratic** — when `Q ≡ c`, every near-star coef
//!    across every sector equals `c`.
//! 4. **Linear Q + C¹/C² join** — a linear `Q` applied to all sectors
//!    should produce near-star tri-5 coefs that satisfy KP09a eq (3)
//!    across adjacent sectors.
//! 5. **Degree-raise 6 × 6 invertibility.**
//! 6. **Cross-sector consistency** — the shared star `IDX_500 = 0`
//!    lies in the near-star region for every sector, and all sectors
//!    must agree on it.

use super::super::kp13_bi5_connect::{ConnectWeights, c1_join, c2_join};
use super::super::kp13_bi5_gtri5::{Gtri5, IDX_500};
use super::*;

const EPS: f64 = 1e-10;

fn approx_point(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

fn approx_points(a: &[Point3], b: &[Point3]) -> bool {
    a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| approx_point(*x, *y))
}

// ---- polar form properties ----

#[test]
fn polar_diagonal_equals_q() {
    let q = StarQuadratic {
        alpha: Point3::new(1.0, 0.0, 0.5),
        beta: Point3::new(0.0, 1.0, 0.0),
        gamma: Point3::new(0.5, 0.0, 1.0),
        delta: Point3::new(1.0, 0.0, 0.0),
        epsilon: Point3::new(0.0, 1.0, 0.0),
        zeta: Point3::new(0.0, 0.0, 1.0),
    };
    for &(x, y) in &[(0.1, 0.2), (0.5, -0.3), (1.0, 0.0), (0.0, 1.0)] {
        let b = q.polar((x, y), (x, y));
        let direct = q.evaluate(x, y);
        assert!(approx_point(b, direct), "B(P, P) != Q(P) at ({x}, {y})");
    }
}

#[test]
fn polar_is_symmetric() {
    let q = StarQuadratic {
        alpha: Point3::new(0.5, 0.0, 0.0),
        beta: Point3::new(1.0, 0.0, 0.0),
        gamma: Point3::new(-0.5, 0.0, 0.0),
        delta: Point3::new(0.25, 0.0, 0.0),
        epsilon: Point3::new(-0.25, 0.0, 0.0),
        zeta: Point3::new(1.0, 0.0, 0.0),
    };
    let p1 = (0.3, 0.4);
    let p2 = (-0.2, 0.7);
    assert!(approx_point(q.polar(p1, p2), q.polar(p2, p1)));
}

// ---- degree-raise + its inverse ----

#[test]
fn degree_raise_roundtrip() {
    let q = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 0.3, 0.0),
        Point3::new(-0.2, 0.0, 0.8),
        Point3::new(0.0, 0.6, -0.4),
        Point3::new(0.9, -0.1, 0.2),
        Point3::new(0.1, 0.4, 1.1),
    ];
    let b = tri2_to_tri5_near_star(&q);
    let q_recovered = tri5_near_star_to_tri2(&b);
    assert!(approx_points(&q, &q_recovered));
}

#[test]
fn degree_raise_constant_preserves_constant() {
    let c = Point3::new(2.0, -3.0, 5.0);
    let q = [c; 6];
    let b = tri2_to_tri5_near_star(&q);
    for (slot, bb) in b.iter().enumerate() {
        assert!(
            approx_point(*bb, c),
            "slot {slot}: expected constant {c:?}, got {bb:?}"
        );
    }
}

#[test]
fn degree_raise_q_002_pure() {
    // Only q_{0,0,2} = 1, rest zero.
    // Expected b: b_{0,0,5}=1, b_{0,1,4}=6/10=0.6, b_{0,2,3}=3/10=0.3,
    // b_{1,0,4}=0.6, b_{1,1,3}=0.3, b_{2,0,3}=0.3.
    let one = Point3::new(1.0, 0.0, 0.0);
    let zero = Point3::new(0.0, 0.0, 0.0);
    let q = [one, zero, zero, zero, zero, zero];
    let b = tri2_to_tri5_near_star(&q);
    assert!(approx_point(b[0], one));
    assert!(approx_point(b[1], Point3::new(0.6, 0.0, 0.0)));
    assert!(approx_point(b[2], Point3::new(0.3, 0.0, 0.0)));
    assert!(approx_point(b[3], Point3::new(0.6, 0.0, 0.0)));
    assert!(approx_point(b[4], Point3::new(0.3, 0.0, 0.0)));
    assert!(approx_point(b[5], Point3::new(0.3, 0.0, 0.0)));
}

// ---- StarQuadratic round-trip ----

#[test]
fn quadratic_round_trip_via_sector0() {
    let q = StarQuadratic {
        alpha: Point3::new(0.7, 0.0, 0.0),
        beta: Point3::new(-0.4, 0.0, 0.0),
        gamma: Point3::new(0.3, 0.0, 0.0),
        delta: Point3::new(0.9, 0.0, 0.0),
        epsilon: Point3::new(-0.2, 0.0, 0.0),
        zeta: Point3::new(1.1, 0.0, 0.0),
    };
    for n in [3, 4, 5, 6, 7, 8, 10, 12] {
        // Sector 0 corners are (1, 0) and (cos 2π/n, sin 2π/n).
        let v_a = (1.0_f64, 0.0);
        let v_b = (
            (2.0 * std::f64::consts::PI / n as f64).cos(),
            (2.0 * std::f64::consts::PI / n as f64).sin(),
        );
        let near = q.tri5_near_star_on_wedge(v_a, v_b);
        let q_rec = StarQuadratic::from_sector0_near_star(&near, n);
        assert!(approx_point(q.alpha, q_rec.alpha), "n = {n}: α mismatch");
        assert!(approx_point(q.beta, q_rec.beta), "n = {n}: β mismatch");
        assert!(approx_point(q.gamma, q_rec.gamma), "n = {n}: γ mismatch");
        assert!(approx_point(q.delta, q_rec.delta), "n = {n}: δ mismatch");
        assert!(
            approx_point(q.epsilon, q_rec.epsilon),
            "n = {n}: ε mismatch"
        );
        assert!(approx_point(q.zeta, q_rec.zeta), "n = {n}: ζ mismatch");
    }
}

// ---- apply_to_sector / apply_to_all_sectors ----

#[test]
fn constant_quadratic_fills_near_star_with_c() {
    let c = Point3::new(1.5, -2.5, 3.5);
    let q = StarQuadratic::constant(c);
    for n in [3, 4, 5, 6, 7, 10] {
        let mut g = Gtri5::new(n);
        q.apply_to_all_sectors(&mut g);
        for k in 0..n {
            for &flat in &NEAR_STAR_FLAT {
                assert!(
                    approx_point(g.sector(k).cps[flat], c),
                    "n = {n}, sector {k}, flat {flat}: expected {c:?}, got {:?}",
                    g.sector(k).cps[flat]
                );
            }
        }
    }
}

#[test]
fn all_sectors_share_star_point_cp() {
    // The star point (flat IDX_500 = 0) is the first near-star slot
    // and equals Q(o) = ζ, independent of sector.
    let q = StarQuadratic {
        alpha: Point3::new(1.0, 0.0, 0.0),
        beta: Point3::new(0.0, 1.0, 0.0),
        gamma: Point3::new(0.0, 0.0, 1.0),
        delta: Point3::new(0.5, 0.5, 0.5),
        epsilon: Point3::new(-0.5, 0.5, 0.5),
        zeta: Point3::new(0.1, 0.2, 0.3),
    };
    let n = 7;
    let mut g = Gtri5::new(n);
    q.apply_to_all_sectors(&mut g);
    let expected = q.zeta;
    for k in 0..n {
        assert!(
            approx_point(g.sector(k).cps[IDX_500], expected),
            "sector {k}: star != ζ"
        );
    }
}

// ---- C¹/C² compatibility (the purpose of the identity) ----

#[test]
fn quadratic_near_star_satisfies_c1_join() {
    // If a single global quadratic populates all sectors' near-star
    // tri-5 coefs, then the C¹ join equation (KP09a eq 3) should hold
    // exactly across adjacent-sector boundaries. Only `j ∈ {0, 1}` sits
    // entirely in the near-star block of sector i (needs (0, j), (0, j+1),
    // (1, j)) and entirely in near-star of sector i+1 (needs (j, 1)).
    let q = StarQuadratic {
        alpha: Point3::new(0.3, 0.0, 0.0),
        beta: Point3::new(-0.2, 0.0, 0.0),
        gamma: Point3::new(0.1, 0.0, 0.0),
        delta: Point3::new(0.5, 0.0, 0.0),
        epsilon: Point3::new(-0.3, 0.0, 0.0),
        zeta: Point3::new(0.7, 0.0, 0.0),
    };
    for n in [5, 6, 7, 8] {
        let mut g = Gtri5::new(n);
        q.apply_to_all_sectors(&mut g);
        let w = ConnectWeights::for_valence(n);
        for k in 0..n {
            let kp1 = (k + 1) % n;
            let s_k = g.sector(k);
            let s_kp1 = g.sector(kp1);
            for j in 0..=1 {
                let predicted = c1_join(&w, s_k, j);
                let actual = s_kp1.cp(j, 1);
                assert!(
                    (predicted.x - actual.x).abs() < 1e-9,
                    "n = {n}, k = {k}, j = {j}: predicted {predicted:?} != actual {actual:?}"
                );
            }
        }
    }
}

#[test]
fn quadratic_near_star_satisfies_c2_join_at_j0() {
    // C² needs p^i_{0,j}, p^i_{0,j+1}, p^i_{0,j+2}, p^i_{1,j},
    // p^i_{1,j+1}, p^i_{2,j} — all in near-star only when j = 0.
    let q = StarQuadratic {
        alpha: Point3::new(0.25, 0.0, 0.0),
        beta: Point3::new(0.4, 0.0, 0.0),
        gamma: Point3::new(-0.3, 0.0, 0.0),
        delta: Point3::new(0.1, 0.0, 0.0),
        epsilon: Point3::new(-0.6, 0.0, 0.0),
        zeta: Point3::new(1.2, 0.0, 0.0),
    };
    for n in [5, 6, 7, 8] {
        let mut g = Gtri5::new(n);
        q.apply_to_all_sectors(&mut g);
        let w = ConnectWeights::for_valence(n);
        for k in 0..n {
            let kp1 = (k + 1) % n;
            let predicted = c2_join(&w, g.sector(k), 0);
            let actual = g.sector(kp1).cp(0, 2);
            assert!(
                (predicted.x - actual.x).abs() < 1e-9,
                "n = {n}, k = {k}: predicted {predicted:?} != actual {actual:?}"
            );
        }
    }
}

// ---- sanity ----

#[test]
fn near_star_flat_indices_no_duplicates() {
    let mut xs = NEAR_STAR_FLAT;
    xs.sort();
    for w in xs.windows(2) {
        assert_ne!(w[0], w[1]);
    }
}

#[test]
fn near_star_flat_includes_star_point_only_once() {
    let count = NEAR_STAR_FLAT.iter().filter(|&&f| f == IDX_500).count();
    assert_eq!(count, 1);
}
