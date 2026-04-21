//! Tests for K-P 2015 bi-6 circulant system solvers.
//!
//! Strategy: build deterministic per-sector input data, solve
//! the circulant system, then independently apply `rule_s2_p21`
//! / `rule_s2_p12` / `rule_s2_p22` at every sector pair and
//! verify the solver's output *satisfies* the rule assignment
//! to machine precision. This is the strongest correctness
//! check short of cross-referencing a reference implementation.

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::kp15_bi6_mu::{mu_ddot, mu_dot, mu_tdot};
use crate::geometry::gregory::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::gregory::kp15_bi6_rules::{rule_s2_p12, rule_s2_p21, rule_s2_p22};

const EPS: f64 = 1e-9;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_eq(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

/// Deterministic non-trivial per-sector inputs. Varies with k to
/// avoid the constant-input null-space edge cases.
fn make_sectors(n: usize) -> Vec<CirculantSectorRhs> {
    let c_star = pt(0.0, 0.0, 1.0); // shared star point
    (0..n)
        .map(|k| {
            let theta = (k as f64) * 2.0 * std::f64::consts::PI / (n as f64);
            let dx = theta.cos();
            let dy = theta.sin();
            let r = 1.0 + 0.1 * (k as f64); // Radial scale differs per sector.
            CirculantSectorRhs {
                p_00: c_star,
                p_10: pt(dx * 0.1 * r, dy * 0.1 * r, 0.9),
                p_20: pt(dx * 0.2 * r, dy * 0.2 * r, 0.7),
                p_40: pt(dx * 0.8 * r, dy * 0.8 * r, 0.3),
                p_50: pt(dx * 0.9 * r, dy * 0.9 * r, 0.2),
                p_60: pt(dx * 1.0 * r, dy * 1.0 * r, 0.1),
                p_01: pt(-dy * 0.1 * r, dx * 0.1 * r, 0.9),
                p_11: pt(dx * 0.1 * r - dy * 0.1, dy * 0.1 * r + dx * 0.1, 0.75),
                p_02: pt(-dy * 0.2 * r, dx * 0.2 * r, 0.7),
                p_42: pt(-dy * 0.5, dx * 0.5, 0.05),
                p_52: pt(-dy * 0.6, dx * 0.6, 0.03),
                p_62: pt(-dy * 0.7, dx * 0.7, 0.01),
            }
        })
        .collect()
}

// ---- 2n×2n system (S2) eq (7) ----

/// Cross-sector round-trip verification at a single k.
///
/// Per the paper's u-v swap convention [KP15 §3.2 note after eq
/// (7)], in the local (p̀, ṕ) frame the previous sector's BB-CPs
/// map as:
///   p̀_{21}  ↔  global q^{k-1}_{12}
///   p̀_{12}  ↔  global q^{k-1}_{21}
/// (diagonal / row-0 indices are unswapped).
///
/// So the round-trip check feeds `rule_s2_p21` / `rule_s2_p12`
/// with p_prev_21 := q^{k-1}_{12} (= res.p_12[km1]) and
/// p_prev_12 := q^{k-1}_{21} (= res.p_21[km1]).
fn assert_eq7_round_trip(
    n: usize,
    md: [f64; 5],
    mdd: [f64; 8],
    sectors: &[CirculantSectorRhs],
    res: &S2Eq7Result,
) {
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let s = &sectors[km1];
        let e21 = rule_s2_p21(
            md,
            s.p_00,
            s.p_10,
            s.p_20,
            s.p_40,
            s.p_50,
            s.p_60,
            res.p_12[km1],
        );
        assert!(
            approx_eq(e21, res.p_21[k]),
            "p_21 round-trip fail at n={n} k={k}: expected {:?}, got {:?}",
            e21,
            res.p_21[k]
        );
        let e12 = rule_s2_p12(
            mdd,
            s.p_00,
            s.p_01,
            s.p_10,
            s.p_11,
            s.p_20,
            res.p_12[km1], // p_prev_21 (local) = q^{k-1}_{12}
            s.p_40,
            s.p_50,
            s.p_60,
            res.p_21[km1], // p_prev_12 (local) = q^{k-1}_{21}
        );
        assert!(
            approx_eq(e12, res.p_12[k]),
            "p_12 round-trip fail at n={n} k={k}: expected {:?}, got {:?}",
            e12,
            res.p_12[k]
        );
    }
}

#[test]
fn s2_eq7_round_trip_at_n5() {
    let params = Bi6Reparams::for_valence(5);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(5);
    let res = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).expect("n=5 solvable");
    assert_eq!(res.p_21.len(), 5);
    assert_eq!(res.p_12.len(), 5);
    assert_eq7_round_trip(5, md, mdd, &sectors, &res);
}

#[test]
fn s2_eq7_round_trip_at_n7() {
    let params = Bi6Reparams::for_valence(7);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(7);
    let res = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).expect("n=7 solvable");
    assert_eq7_round_trip(7, md, mdd, &sectors, &res);
}

#[test]
fn s2_eq7_round_trip_at_n8() {
    let params = Bi6Reparams::for_valence(8);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(8);
    let res = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).expect("n=8 solvable");
    assert_eq7_round_trip(8, md, mdd, &sectors, &res);
}

#[test]
fn s2_eq7_requires_pinned_p0_21_at_n6() {
    let params = Bi6Reparams::for_valence(6);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(6);
    // Without pinning: should error.
    let err = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).unwrap_err();
    assert!(err.contains("pinned"));
}

#[test]
fn s2_eq7_at_n6_honours_pinned_p0_21() {
    let params = Bi6Reparams::for_valence(6);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(6);
    let pinned = pt(0.42, -0.13, 0.05);
    let res =
        solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, Some(pinned)).expect("n=6 with pin");
    // p⁰_{21} exactly equals the pinned value.
    assert!(approx_eq(res.p_21[0], pinned));
    // Rule (S2) consistency on every k ∈ 1..n (skip k=0 since
    // that's where we pinned, replacing the eq). Uses the
    // paper's u-v swap convention.
    for k in 1..6 {
        let km1 = k - 1;
        let s = &sectors[km1];
        let e21 = rule_s2_p21(
            md,
            s.p_00,
            s.p_10,
            s.p_20,
            s.p_40,
            s.p_50,
            s.p_60,
            res.p_12[km1],
        );
        assert!(approx_eq(e21, res.p_21[k]), "n=6 p_21 ring fails at k={k}");
    }
}

#[test]
fn s2_eq7_requires_pin_at_n3() {
    let params = Bi6Reparams::for_valence(3);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let sectors = make_sectors(3);
    let err = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).unwrap_err();
    assert!(err.contains("pinned"));
    // With pin, solvable.
    let pinned = pt(0.0, 0.0, 0.0);
    let _res =
        solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, Some(pinned)).expect("n=3 pinned");
}

// ---- n×n system (S2) eq (8) ----

#[test]
fn s2_eq8_round_trip_at_n5() {
    let params = Bi6Reparams::for_valence(5);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let mt = mu_tdot(&params);
    let sectors = make_sectors(5);
    let res21 = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).unwrap();
    let p0_22_pinned = pt(0.0, 0.0, 0.8);
    let p22 = solve_bi6_circulant_s2_eq8(&params, mt, &sectors, &res21.p_12, p0_22_pinned)
        .expect("n=5 p22 solvable");

    assert!(approx_eq(p22[0], p0_22_pinned));

    // Round-trip: rule_s2_p22 with p̀ inputs and p̀_{22} should
    // produce ṕ_{22} = p22[k]. Per the swap convention, the
    // rule's `p_prev_21` parameter gets the global (1, 2) CP:
    // res21.p_12[km1].
    for k in 1..5 {
        let km1 = k - 1;
        let s_prev = &sectors[km1];
        let s_cur = &sectors[k];
        let diff_42 = s_prev.p_42 - s_cur.p_42;
        let diff_52 = s_prev.p_52 - s_cur.p_52;
        let diff_62 = s_prev.p_62 - s_cur.p_62;
        let expected = rule_s2_p22(
            mt,
            s_prev.p_00,
            s_prev.p_01,
            s_prev.p_10,
            s_prev.p_11,
            s_prev.p_20,
            res21.p_12[km1], // p_prev_21 (local) = q^{k-1}_{12}
            s_prev.p_40,
            s_prev.p_50,
            s_prev.p_60,
            p22[km1],
            Point3::new(diff_42.x, diff_42.y, diff_42.z),
            Point3::new(diff_52.x, diff_52.y, diff_52.z),
            Point3::new(diff_62.x, diff_62.y, diff_62.z),
        );
        assert!(
            approx_eq(expected, p22[k]),
            "p22 round-trip fail at k={k}: expected {:?}, got {:?}",
            expected,
            p22[k]
        );
    }
}

#[test]
fn s2_eq8_rejected_at_n3() {
    let params = Bi6Reparams::for_valence(3);
    let mt = mu_tdot(&params);
    let sectors = make_sectors(3);
    let p_12 = vec![pt(0.0, 0.0, 0.0); 3];
    let err =
        solve_bi6_circulant_s2_eq8(&params, mt, &sectors, &p_12, pt(0.0, 0.0, 0.0)).unwrap_err();
    assert!(err.contains("n = 3"));
}

#[test]
fn s2_eq8_p0_22_pin_sets_first_entry() {
    // Regardless of other inputs, solved p22[0] equals the pin.
    let params = Bi6Reparams::for_valence(7);
    let md = mu_dot(&params);
    let mdd = mu_ddot(&params);
    let mt = mu_tdot(&params);
    let sectors = make_sectors(7);
    let res21 = solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).unwrap();
    for &pin in &[pt(0.0, 0.0, 0.0), pt(1.0, 2.0, 3.0), pt(-0.5, 0.7, 0.1)] {
        let p22 = solve_bi6_circulant_s2_eq8(&params, mt, &sectors, &res21.p_12, pin).unwrap();
        assert!(approx_eq(p22[0], pin));
    }
}

// ---- Cross-valence smoke ----

#[test]
fn cross_valence_smoke() {
    for n in [5usize, 7, 8, 9, 10, 12] {
        let params = Bi6Reparams::for_valence(n);
        let md = mu_dot(&params);
        let mdd = mu_ddot(&params);
        let mt = mu_tdot(&params);
        let sectors = make_sectors(n);
        let res21 =
            solve_bi6_circulant_s2_eq7(&params, md, mdd, &sectors, None).expect("eq7 solvable");
        assert_eq!(res21.p_21.len(), n);
        for cp in &res21.p_21 {
            assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
        }
        let p22 = solve_bi6_circulant_s2_eq8(&params, mt, &sectors, &res21.p_12, pt(0.0, 0.0, 0.0))
            .expect("eq8 solvable");
        assert_eq!(p22.len(), n);
    }
}
