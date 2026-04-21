//! Tests for the KP09a §2.3(b) circulant solver.
//!
//! # Strategy
//!
//! 1. **RHS constant check**: if every input CP equals the same
//!    point `P`, `R^i = (2c + 1)·P` (sum of paper coefficients on a
//!    constant).
//! 2. **Circulant identity**: after solving, the three-term
//!    recurrence `g^{i-1}_{21} + 4c·g^i_{21} + g^{i+1}_{21} = 2·R^i`
//!    holds for every sector.
//! 3. **Constant invariance**: constant RHS → constant g^i_{21} = R^i.
//! 4. **End-to-end**: quadratic Q populates near-star, then set free
//!    `g^i_{30}`, then run `complete_second_derivative_layer` — the
//!    circulant identity holds and `g^i_{12}` satisfies its derived
//!    formula.
//! 5. **Excluded valences**: error on `n ∈ {3, 4, 6}`.

use super::super::kp13_bi5_connect::paper_flat;
use super::super::kp13_bi5_gtri5::Gtri5;
use super::super::kp13_bi5_star_quadratic::StarQuadratic;
use super::*;

const EPS: f64 = 1e-9;

fn approx_point(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- RHS ----

#[test]
fn rhs_constant_matches_expected_sum() {
    // All input CPs = P → R^i = (2c + 1)·P.
    for n in [5, 7, 8, 10] {
        let p = Point3::new(1.0, 0.5, -0.25);
        let mut g = Gtri5::new(n);
        for k in 0..n {
            g.sector_mut(k).cps = [p; 21];
        }
        let c = (2.0 * std::f64::consts::PI / n as f64).cos();
        let expected_scalar = 2.0 * c + 1.0;
        let expected = Point3::new(
            expected_scalar * p.x,
            expected_scalar * p.y,
            expected_scalar * p.z,
        );
        for i in 0..n {
            let r = compute_r_i(g.sector(i), n);
            assert!(
                approx_point(r, expected),
                "n = {n}, sector {i}: R = {r:?}, expected {expected:?}"
            );
        }
    }
}

// ---- Excluded valences ----

#[test]
fn solve_rejects_n3() {
    let mut g = Gtri5::new(3);
    assert_eq!(
        solve_g21_circulant(&mut g),
        Err(CirculantError::UnsupportedValence(3))
    );
}

#[test]
fn solve_rejects_n4() {
    let mut g = Gtri5::new(4);
    assert_eq!(
        solve_g21_circulant(&mut g),
        Err(CirculantError::UnsupportedValence(4))
    );
}

#[test]
fn solve_rejects_n6() {
    let mut g = Gtri5::new(6);
    assert_eq!(
        solve_g21_circulant(&mut g),
        Err(CirculantError::UnsupportedValence(6))
    );
}

// ---- Circulant identity holds after solve ----

/// Verify `g^{i-1}_{21} + 4c·g^i_{21} + g^{i+1}_{21} = 2·R^i` for
/// every sector.
fn assert_circulant_identity(gtri5: &Gtri5) {
    let n = gtri5.n;
    let c = (2.0 * std::f64::consts::PI / n as f64).cos();
    let idx_21 = paper_flat(2, 1);
    for i in 0..n {
        let prev = (i + n - 1) % n;
        let next = (i + 1) % n;
        let gm = gtri5.sector(prev).cps[idx_21];
        let gi = gtri5.sector(i).cps[idx_21];
        let gp = gtri5.sector(next).cps[idx_21];
        let lhs = Point3::new(
            gm.x + 4.0 * c * gi.x + gp.x,
            gm.y + 4.0 * c * gi.y + gp.y,
            gm.z + 4.0 * c * gi.z + gp.z,
        );
        let r = compute_r_i(gtri5.sector(i), n);
        let rhs = Point3::new(2.0 * r.x, 2.0 * r.y, 2.0 * r.z);
        let dx = (lhs.x - rhs.x).abs();
        let dy = (lhs.y - rhs.y).abs();
        let dz = (lhs.z - rhs.z).abs();
        assert!(
            dx < 1e-8 && dy < 1e-8 && dz < 1e-8,
            "sector {i}: circulant identity fails, lhs = {lhs:?}, rhs = {rhs:?}"
        );
    }
}

#[test]
fn circulant_identity_after_solve_constant() {
    let p = Point3::new(1.0, 2.0, 3.0);
    for n in [5, 7, 8, 9, 10] {
        let mut g = Gtri5::new(n);
        for k in 0..n {
            g.sector_mut(k).cps = [p; 21];
        }
        solve_g21_circulant(&mut g).unwrap();
        assert_circulant_identity(&g);
    }
}

#[test]
fn circulant_identity_after_solve_from_quadratic_plus_g30() {
    let q = StarQuadratic {
        alpha: Point3::new(0.3, 0.0, 0.0),
        beta: Point3::new(-0.2, 0.0, 0.0),
        gamma: Point3::new(0.1, 0.0, 0.0),
        delta: Point3::new(0.5, 0.0, 0.0),
        epsilon: Point3::new(-0.3, 0.0, 0.0),
        zeta: Point3::new(0.7, 0.0, 0.0),
    };
    for n in [5, 7, 8, 9, 10] {
        let mut g = Gtri5::new(n);
        q.apply_to_all_sectors(&mut g);
        // Set each sector's g^i_{30} to some sector-dependent value.
        let idx_30 = paper_flat(3, 0);
        for k in 0..n {
            let angle = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            g.sector_mut(k).cps[idx_30] = Point3::new(angle.cos() * 0.5, 0.0, angle.sin() * 0.5);
        }
        // Sync edges, then solve.
        sync_edge_layer_03_and_30(&mut g);
        solve_g21_circulant(&mut g).unwrap();
        assert_circulant_identity(&g);
    }
}

// ---- Constant invariance ----

#[test]
fn solve_constant_gives_constant_g21() {
    // Constant R^i = R across sectors → circulant has constant
    // solution `G` with `(2 + 4c)·G = 2·R`, so `G = R / (1 + 2c)`.
    // With all CPs = P: R = (1 + 2c)·P (checked in
    // `rhs_constant_matches_expected_sum`), so `G = P`.
    let p = Point3::new(1.0, -0.5, 0.25);
    for n in [5, 7, 8] {
        let mut g = Gtri5::new(n);
        for k in 0..n {
            g.sector_mut(k).cps = [p; 21];
        }
        solve_g21_circulant(&mut g).unwrap();
        let idx_21 = paper_flat(2, 1);
        for i in 0..n {
            assert!(
                approx_point(g.sector(i).cps[idx_21], p),
                "n = {n}, sector {i}: g_21 = {:?}, expected {p:?}",
                g.sector(i).cps[idx_21]
            );
        }
    }
}

// ---- edge sync ----

#[test]
fn sync_edge_03_30_makes_shared_consistent() {
    let n = 7;
    let mut g = Gtri5::new(n);
    let idx_30 = paper_flat(3, 0);
    let idx_03 = paper_flat(0, 3);
    // Populate each sector's (3, 0) with a unique value.
    for k in 0..n {
        g.sector_mut(k).cps[idx_30] = Point3::new(k as f64, 0.0, 0.0);
    }
    sync_edge_layer_03_and_30(&mut g);
    // After sync, sector i's (0, 3) should equal sector (i+1)'s (3, 0).
    for i in 0..n {
        let next = (i + 1) % n;
        assert!(
            approx_point(g.sector(i).cps[idx_03], g.sector(next).cps[idx_30]),
            "sector {i} (0,3) != sector {next} (3,0)"
        );
    }
}

#[test]
fn sync_shared_edge_all_depths_covers_j3_and_j4() {
    let n = 6;
    let mut g = Gtri5::new(n);
    let idx_30 = paper_flat(3, 0);
    let idx_40 = paper_flat(4, 0);
    let idx_03 = paper_flat(0, 3);
    let idx_04 = paper_flat(0, 4);
    // Unique values per sector on (3, 0) and (4, 0).
    for k in 0..n {
        g.sector_mut(k).cps[idx_30] = Point3::new(10.0 + k as f64, 0.0, 0.0);
        g.sector_mut(k).cps[idx_40] = Point3::new(20.0 + k as f64, 0.0, 0.0);
    }
    sync_shared_edge_all_depths(&mut g);
    for i in 0..n {
        let next = (i + 1) % n;
        assert!(
            approx_point(g.sector(i).cps[idx_03], g.sector(next).cps[idx_30]),
            "j=3 mismatch at sector {i}"
        );
        assert!(
            approx_point(g.sector(i).cps[idx_04], g.sector(next).cps[idx_40]),
            "j=4 mismatch at sector {i}"
        );
    }
}

// ---- g^i_{12} back-fill ----

#[test]
fn fill_g12_matches_eq3_formula() {
    // After solving the circulant, g^i_{12} should satisfy
    // g^i_{12} = ν₀·g^i_{02} + ν₁·g^i_{03} − g^{i+1}_{21}.
    let q = StarQuadratic {
        alpha: Point3::new(0.2, 0.0, 0.0),
        beta: Point3::new(0.0, 0.0, 0.0),
        gamma: Point3::new(-0.3, 0.0, 0.0),
        delta: Point3::new(0.7, 0.0, 0.0),
        epsilon: Point3::new(0.4, 0.0, 0.0),
        zeta: Point3::new(1.1, 0.0, 0.0),
    };
    let n = 7;
    let mut g = Gtri5::new(n);
    q.apply_to_all_sectors(&mut g);
    let idx_30 = paper_flat(3, 0);
    for k in 0..n {
        g.sector_mut(k).cps[idx_30] = Point3::new(0.1 * k as f64, 0.0, 0.0);
    }
    complete_second_derivative_layer(&mut g).unwrap();

    let c = (2.0 * std::f64::consts::PI / n as f64).cos();
    let c_bar = 1.0 - c;
    let nu0 = 2.0 * c_bar;
    let nu1 = 2.0 * c;
    let idx_02 = paper_flat(0, 2);
    let idx_03 = paper_flat(0, 3);
    let idx_12 = paper_flat(1, 2);
    let idx_21 = paper_flat(2, 1);

    for i in 0..n {
        let s = g.sector(i);
        let g_02 = s.cps[idx_02];
        let g_03 = s.cps[idx_03];
        let g_next_21 = g.sector((i + 1) % n).cps[idx_21];
        let expected = Point3::new(
            nu0 * g_02.x + nu1 * g_03.x - g_next_21.x,
            nu0 * g_02.y + nu1 * g_03.y - g_next_21.y,
            nu0 * g_02.z + nu1 * g_03.z - g_next_21.z,
        );
        assert!(
            approx_point(s.cps[idx_12], expected),
            "sector {i}: g_12 = {:?}, expected {expected:?}",
            s.cps[idx_12]
        );
    }
}
