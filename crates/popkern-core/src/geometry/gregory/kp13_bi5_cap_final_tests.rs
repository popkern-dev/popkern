//! Tests for I.14.e final cap `p` construction.

#![allow(clippy::needless_range_loop)]

use super::super::kp13_bi5_cap::{Bi5Cap, CapPieceId};
use super::*;

const EPS: f64 = 1e-9;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- Free-param struct ----

#[test]
fn zero_params_has_correct_shape() {
    for n in [3, 4, 5, 6, 7, 8] {
        let p = Bi5CapFreeParams::zero(n);
        assert_eq!(p.n, n);
        assert_eq!(p.g3.len(), n + 6);
        assert_eq!(p.p_o_i2.len(), n);
        assert_eq!(p.num_scalar_dofs_per_axis(), 5 * n + 6);
    }
}

#[test]
#[should_panic]
fn zero_params_below_three_panics() {
    let _ = Bi5CapFreeParams::zero(2);
}

// ---- Eq (18) strip transfer ----

#[test]
fn eq18_copies_correct_cps_from_a_to_p() {
    let n = 5;
    let mut a = Bi5Cap::zero(n);
    // Populate a's r-piece and l-piece with distinct patterns.
    for k in 0..n {
        let r = a.sector_mut(k).piece_mut(CapPieceId::R);
        for i in 0..6 {
            for j in 0..6 {
                r.cps[i][j] = Point3::new((k * 100 + i * 10 + j) as f64, 0.0, 0.0);
            }
        }
        let l = a.sector_mut(k).piece_mut(CapPieceId::L);
        for i in 0..6 {
            for j in 0..6 {
                l.cps[i][j] = Point3::new(0.0, (k * 100 + i * 10 + j) as f64, 0.0);
            }
        }
    }
    let mut p = Bi5Cap::zero(n);
    apply_eq18_strips(&mut p, &a);

    // Check specific CPs transferred.
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let p_r_prev = p.sector(k_prev).piece(CapPieceId::R);
        let a_r_prev = a.sector(k_prev).piece(CapPieceId::R);
        for i in 3..=5 {
            for j in 0..=2 {
                assert_eq!(p_r_prev.cps[i][j], a_r_prev.cps[i][j]);
            }
        }
        let p_l = p.sector(k).piece(CapPieceId::L);
        let a_l = a.sector(k).piece(CapPieceId::L);
        for i in 3..=5 {
            for j in 0..=2 {
                assert_eq!(p_l.cps[i][j], a_l.cps[i][j]);
            }
        }
    }

    // CPs outside the (i ∈ [3, 5], j ∈ [0, 2]) strip should still be zero.
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let p_r_prev = p.sector(k_prev).piece(CapPieceId::R);
        for i in 0..3 {
            for j in 0..6 {
                assert_eq!(p_r_prev.cps[i][j], Point3::new(0.0, 0.0, 0.0));
            }
        }
        for i in 3..=5 {
            for j in 3..6 {
                assert_eq!(p_r_prev.cps[i][j], Point3::new(0.0, 0.0, 0.0));
            }
        }
    }
}

// ---- build_p_baseline ----

#[test]
fn build_p_baseline_starts_from_a_with_eq18_override() {
    let n = 4;
    let mut a = Bi5Cap::zero(n);
    for k in 0..n {
        for id in CapPieceId::all() {
            let piece = a.sector_mut(k).piece_mut(id);
            for i in 0..6 {
                for j in 0..6 {
                    piece.cps[i][j] = Point3::new(k as f64, id as usize as f64, (i * 7 + j) as f64);
                }
            }
        }
    }
    let p = build_p_baseline(&a);
    p.validate().unwrap();
    // Eq (18) transfers specific CPs that match a's. Since all of
    // a's CPs are set, p's corresponding CPs should also match a's.
    // Most of p's CPs just retain a's values (from the clone).
    for k in 0..n {
        for id in [CapPieceId::O, CapPieceId::M] {
            // O and M are unchanged by eq (18).
            assert_eq!(
                p.sector(k).piece(id).cps,
                a.sector(k).piece(id).cps,
                "k={k}, {id:?}: baseline should copy from a"
            );
        }
    }
}

// ---- apply_free_params_p ----

#[test]
fn zero_free_params_equals_baseline() {
    let n = 5;
    let mut a = Bi5Cap::zero(n);
    for k in 0..n {
        let o = a.sector_mut(k).piece_mut(CapPieceId::O);
        for i in 0..6 {
            for j in 0..6 {
                o.cps[i][j] = Point3::new((k + i + j) as f64 * 0.1, 0.0, 0.0);
            }
        }
    }
    let params = Bi5CapFreeParams::zero(n);
    let p = apply_free_params_p(&params, &a);
    let baseline = build_p_baseline(&a);

    // With zero free params, apply_free_params_p STILL writes zeros
    // at the star-region 6 CPs of sector 0 and p^{o,k}_{3,0} +
    // p^{o,k}_{i, 2}. So if baseline had non-zero values there,
    // apply zeros them out.
    //
    // For this test, we check that NON-overridden CPs match
    // baseline.
    for k in 0..n {
        let p_o = p.sector(k).piece(CapPieceId::O);
        let b_o = baseline.sector(k).piece(CapPieceId::O);
        // Check the untouched interior (not in any override position).
        for i in 0..6 {
            for j in 3..6 {
                // Skip positions that overlap with overrides.
                if (i == 3 || i == 2 || i == 4 || i == 5) && j == 2 {
                    continue;
                }
                assert_eq!(
                    p_o.cps[i][j], b_o.cps[i][j],
                    "k={k}, O ({i}, {j}): p differs from baseline"
                );
            }
        }
    }
}

#[test]
fn free_params_write_to_designated_positions() {
    let n = 5;
    let a = Bi5Cap::zero(n);
    let mut params = Bi5CapFreeParams::zero(n);
    // Set G_3[0] to a distinctive value.
    params.g3[0] = Point3::new(42.0, 0.0, 0.0);
    // Set G_3[6 + 0] (sector 0's mid-edge CP).
    params.g3[6] = Point3::new(0.0, 77.0, 0.0);
    // Set p_o_i2[1][2] (sector 1's p^{o,1}_{4, 2} analog).
    params.p_o_i2[1][2] = Point3::new(0.0, 0.0, 99.0);

    let p = apply_free_params_p(&params, &a);
    // Sector 0's O piece cps[0][0] should have G_3[0]'s value.
    assert_eq!(
        p.sector(0).piece(CapPieceId::O).cps[0][0],
        Point3::new(42.0, 0.0, 0.0)
    );
    // Sector 0's O piece cps[3][0] should have G_3[6+0]'s value.
    assert_eq!(
        p.sector(0).piece(CapPieceId::O).cps[3][0],
        Point3::new(0.0, 77.0, 0.0)
    );
    // Sector 1's O piece cps[4][2] should have p_o_i2[1][2]'s value
    // (idx=2 means i = 2+2 = 4).
    assert_eq!(
        p.sector(1).piece(CapPieceId::O).cps[4][2],
        Point3::new(0.0, 0.0, 99.0)
    );
}

// ---- F_1 minimisation ----

#[test]
fn solve_f1_min_from_zero_a_gives_zero_p() {
    // Trivial: a = 0 → residual min is p = 0.
    for n in [3, 4, 5, 6, 7, 8] {
        let a = Bi5Cap::zero(n);
        let p = solve_f1_min(&a).unwrap();
        for k in 0..n {
            for id in CapPieceId::all() {
                let piece = p.sector(k).piece(id);
                for i in 0..6 {
                    for j in 0..6 {
                        let cp = piece.cps[i][j];
                        assert!(
                            cp.x.abs() < 1e-8 && cp.y.abs() < 1e-8 && cp.z.abs() < 1e-8,
                            "n={n}, k={k}, {id:?} ({i}, {j}): non-zero {cp:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn solve_f1_min_deterministic() {
    let n = 5;
    let mut a = Bi5Cap::zero(n);
    for k in 0..n {
        let o = a.sector_mut(k).piece_mut(CapPieceId::O);
        for i in 0..6 {
            for j in 0..6 {
                o.cps[i][j] = Point3::new((k + i * j) as f64 * 0.1, 0.0, 0.0);
            }
        }
    }
    let p_a = solve_f1_min(&a).unwrap();
    let p_b = solve_f1_min(&a).unwrap();
    assert_eq!(p_a, p_b);
}

#[test]
fn solve_f1_min_preserves_eq18_strips() {
    // After solve, eq (18) strip CPs of p must still equal a's.
    let n = 5;
    let mut a = Bi5Cap::zero(n);
    for k in 0..n {
        let r = a.sector_mut(k).piece_mut(CapPieceId::R);
        for i in 3..=5 {
            for j in 0..=2 {
                r.cps[i][j] = Point3::new((k + 1) as f64, i as f64, j as f64);
            }
        }
        let l = a.sector_mut(k).piece_mut(CapPieceId::L);
        for i in 3..=5 {
            for j in 0..=2 {
                l.cps[i][j] = Point3::new(-(k as f64 + 1.0), -(i as f64), -(j as f64));
            }
        }
    }

    let p = solve_f1_min(&a).unwrap();
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let p_r = p.sector(k_prev).piece(CapPieceId::R);
        let a_r = a.sector(k_prev).piece(CapPieceId::R);
        for i in 3..=5 {
            for j in 0..=2 {
                assert!(
                    approx_p(p_r.cps[i][j], a_r.cps[i][j], 1e-8),
                    "k={k_prev}, R ({i}, {j}) eq (18): p={:?} vs a={:?}",
                    p_r.cps[i][j],
                    a_r.cps[i][j]
                );
            }
        }
        let p_l = p.sector(k).piece(CapPieceId::L);
        let a_l = a.sector(k).piece(CapPieceId::L);
        for i in 3..=5 {
            for j in 0..=2 {
                assert!(
                    approx_p(p_l.cps[i][j], a_l.cps[i][j], 1e-8),
                    "k={k}, L ({i}, {j}) eq (18): p={:?} vs a={:?}",
                    p_l.cps[i][j],
                    a_l.cps[i][j]
                );
            }
        }
    }
}

#[test]
fn solve_f1_min_works_for_all_valences() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut a = Bi5Cap::zero(n);
        for k in 0..n {
            let o = a.sector_mut(k).piece_mut(CapPieceId::O);
            o.cps[1][1] = Point3::new(k as f64 * 0.5, 0.1, 0.2);
        }
        let p = solve_f1_min(&a).unwrap_or_else(|e| panic!("n={n}: {e}"));
        p.validate().unwrap();
    }
}

#[test]
fn solve_f1_min_reduces_f1_residual_vs_zero_params() {
    // The solver minimises over the free-param space; since
    // zero free params are in that space, F_1(solved − a) must
    // be ≤ F_1(apply(zero, a) − a).
    use super::super::kp13_bi5_cap_f1::f1_difference;

    let n = 5;
    let mut a = Bi5Cap::zero(n);
    for k in 0..n {
        let o = a.sector_mut(k).piece_mut(CapPieceId::O);
        for i in 0..6 {
            for j in 0..6 {
                let u = i as f64 / 5.0;
                let v = j as f64 / 5.0;
                o.cps[i][j] = Point3::new(u + 0.3 * v, 0.5 * v - 0.2 * u, u * v * 0.4);
            }
        }
    }
    let zero_params = Bi5CapFreeParams::zero(n);
    let p_zero = apply_free_params_p(&zero_params, &a);
    let f1_zero = f1_difference(&p_zero, &a);
    let p_solved = solve_f1_min(&a).unwrap();
    let f1_solved = f1_difference(&p_solved, &a);
    assert!(
        f1_solved <= f1_zero + 1e-9,
        "F_1(solved − a) = {f1_solved} should be ≤ zero-params baseline {f1_zero}"
    );
}

#[test]
fn f1_cap_finite_on_all_valences() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        for k in 0..n {
            for id in CapPieceId::all() {
                let p = cap.sector_mut(k).piece_mut(id);
                for i in 0..6 {
                    for j in 0..6 {
                        p.cps[i][j] = Point3::new((k + i) as f64 * 0.1, j as f64 * 0.2, 0.0);
                    }
                }
            }
        }
        let f = cap_f1(&cap);
        assert!(f.is_finite(), "n={n}: F_1 = {f} not finite");
    }
}
