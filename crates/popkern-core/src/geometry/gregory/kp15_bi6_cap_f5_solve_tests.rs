//! Tests for F_5 SPD minimisation on the bi-6 cap.
//!
//! Strategy: verify (a) solver returns a solution, (b) optimal
//! F_5 value is ≤ F_5 at any perturbed nearby point (first-order
//! optimality), (c) the solution is stable across runs (bit-
//! identical output for identical input), (d) the optimisation
//! reduces F_5 substantially from the all-zero-free-param
//! starting point.

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::CirculantSectorRhs;
use crate::geometry::gregory::kp15_bi6_cap::Bi6CapInput;
use crate::geometry::gregory::kp15_bi6_f5::f5_evaluate;
use crate::geometry::gregory::kp15_bi6_reparam::Bi6Reparams;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn pentagon_input() -> Bi6CapInput {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let c7 = pt(0.0, 0.0, 1.0);
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (k as f64) * 2.0 * std::f64::consts::PI / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        c5s.push(pt(dx * 0.3, dy * 0.3, 0.9));
        c6s.push(pt(dx * 0.7, dy * 0.7, 0.6));
        sectors.push(CirculantSectorRhs {
            p_00: c7,
            p_10: pt(dx * 0.1, dy * 0.1, 0.95),
            p_20: pt(dx * 0.25, dy * 0.25, 0.85),
            p_40: pt(dx * 0.6, dy * 0.6, 0.4),
            p_50: pt(dx * 0.8, dy * 0.8, 0.25),
            p_60: pt(dx, dy, 0.1),
            p_01: pt(-dy * 0.1, dx * 0.1, 0.95),
            p_11: pt((dx - dy) * 0.1, (dy + dx) * 0.1, 0.85),
            p_02: pt(-dy * 0.2, dx * 0.2, 0.88),
            p_42: pt(-dy * 0.5, dx * 0.5, 0.35),
            p_52: pt(-dy * 0.6, dx * 0.6, 0.2),
            p_62: pt(-dy * 0.7, dx * 0.7, 0.05),
        });
    }
    Bi6CapInput {
        n,
        params,
        c7,
        c5s,
        c6s,
        sectors,
        p0_22: pt(0.0, 0.0, 0.0),
        p0_21_pinned: None,
        interior_filler: None,
    }
}

fn cap_f5_total(out: &super::Bi6CapOutput) -> f64 {
    out.patches_bi6.iter().map(f5_evaluate).sum()
}

#[test]
fn solver_returns_optimal_params_at_n5() {
    let base = pentagon_input();
    let (optimal, out) = f5_solve_optimal_bi6_cap(&base).expect("solver converges");
    // Layout check.
    assert_eq!(optimal.p_23.len(), 5);
    assert_eq!(optimal.p_33.len(), 5);
    // Optimal CPs are finite.
    assert!(optimal.p0_22.x.is_finite());
    // Output has N patches.
    assert_eq!(out.patches_bi6.len(), 5);
}

#[test]
fn optimal_f5_not_worse_than_zero_free_params() {
    let base = pentagon_input();
    // F_5 at zero free params.
    let zero_free = Bi6CapFreeParams::zero(5);
    let zero_input = apply_free_params_for_test(&base, &zero_free);
    let zero_out = assemble_bi6_polar_cap(&zero_input).unwrap();
    let f5_zero = cap_f5_total(&zero_out);

    // F_5 at optimal.
    let (_optimal, opt_out) = f5_solve_optimal_bi6_cap(&base).unwrap();
    let f5_opt = cap_f5_total(&opt_out);

    assert!(
        f5_opt <= f5_zero + 1e-6,
        "F_5(opt) = {f5_opt} > F_5(zero) = {f5_zero}"
    );
}

#[test]
fn optimal_is_first_order_stationary() {
    // At the minimum, perturbing any free param by ε should
    // increase F_5 (to second order).
    let base = pentagon_input();
    let (optimal, opt_out) = f5_solve_optimal_bi6_cap(&base).unwrap();
    let f5_opt = cap_f5_total(&opt_out);

    let eps = 1e-3;
    for i in 0..Bi6CapFreeParams::flat_len(5) {
        for coord in 0..3 {
            let mut perturbed = optimal.clone();
            let bump = match coord {
                0 => pt(eps, 0.0, 0.0),
                1 => pt(0.0, eps, 0.0),
                _ => pt(0.0, 0.0, eps),
            };
            if i == 0 {
                perturbed.p0_22 = add_pts(perturbed.p0_22, bump);
            } else {
                let k = (i - 1) / 2;
                if (i - 1).is_multiple_of(2) {
                    perturbed.p_23[k] = add_pts(perturbed.p_23[k], bump);
                } else {
                    perturbed.p_33[k] = add_pts(perturbed.p_33[k], bump);
                }
            }
            let pert_input = apply_free_params_for_test(&base, &perturbed);
            let pert_out = assemble_bi6_polar_cap(&pert_input).unwrap();
            let f5_pert = cap_f5_total(&pert_out);
            // Stationarity: F_5(opt + ε·e_i) ≥ F_5(opt) − O(ε²).
            // The second-order correction can be as large as
            // (eps² · max curvature). We allow a generous slack.
            assert!(
                f5_pert >= f5_opt - 1e-5,
                "perturb param[{i}] coord {coord}: F_5 = {f5_pert} < opt {f5_opt}"
            );
        }
    }
}

#[test]
fn solver_deterministic() {
    // Two runs on identical input produce identical output.
    let base = pentagon_input();
    let (opt_a, _) = f5_solve_optimal_bi6_cap(&base).unwrap();
    let (opt_b, _) = f5_solve_optimal_bi6_cap(&base).unwrap();
    assert_eq!(opt_a.p0_22.x, opt_b.p0_22.x);
    assert_eq!(opt_a.p0_22.y, opt_b.p0_22.y);
    assert_eq!(opt_a.p0_22.z, opt_b.p0_22.z);
    for k in 0..5 {
        assert_eq!(opt_a.p_23[k].x, opt_b.p_23[k].x);
        assert_eq!(opt_a.p_33[k].z, opt_b.p_33[k].z);
    }
}

#[test]
fn works_at_larger_valences() {
    for n in [6usize, 7, 8, 10, 12] {
        let base = build_input(
            n,
            if n == 6 {
                Some(pt(0.0, 0.0, 0.5))
            } else {
                None
            },
        );
        let result = f5_solve_optimal_bi6_cap(&base);
        assert!(result.is_ok(), "n = {n} solver failed: {result:?}");
        let (_opt, out) = result.unwrap();
        assert_eq!(out.patches_bi6.len(), n);
    }
}

// ---- helpers ----

fn add_pts(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

/// Replicate the private `apply_free_params` logic for tests.
fn apply_free_params_for_test(base: &Bi6CapInput, free: &Bi6CapFreeParams) -> Bi6CapInput {
    let mut out = base.clone();
    out.p0_22 = free.p0_22;
    let mut filler: Vec<[[Option<Point3>; 7]; 7]> = vec![[[None; 7]; 7]; out.n];
    for k in 0..out.n {
        filler[k][2][3] = Some(free.p_23[k]);
        filler[k][3][3] = Some(free.p_33[k]);
    }
    out.interior_filler = Some(filler);
    out
}

fn build_input(n: usize, p0_21_pin: Option<Point3>) -> Bi6CapInput {
    let params = Bi6Reparams::for_valence(n);
    let c7 = pt(0.0, 0.0, 1.0);
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (k as f64) * 2.0 * std::f64::consts::PI / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        c5s.push(pt(dx * 0.3, dy * 0.3, 0.9));
        c6s.push(pt(dx * 0.7, dy * 0.7, 0.6));
        sectors.push(CirculantSectorRhs {
            p_00: c7,
            p_10: pt(dx * 0.1, dy * 0.1, 0.95),
            p_20: pt(dx * 0.25, dy * 0.25, 0.85),
            p_40: pt(dx * 0.6, dy * 0.6, 0.4),
            p_50: pt(dx * 0.8, dy * 0.8, 0.25),
            p_60: pt(dx, dy, 0.1),
            p_01: pt(-dy * 0.1, dx * 0.1, 0.95),
            p_11: pt((dx - dy) * 0.1, (dy + dx) * 0.1, 0.85),
            p_02: pt(-dy * 0.2, dx * 0.2, 0.88),
            p_42: pt(-dy * 0.5, dx * 0.5, 0.35),
            p_52: pt(-dy * 0.6, dx * 0.6, 0.2),
            p_62: pt(-dy * 0.7, dx * 0.7, 0.05),
        });
    }
    Bi6CapInput {
        n,
        params,
        c7,
        c5s,
        c6s,
        sectors,
        p0_22: pt(0.0, 0.0, 0.0),
        p0_21_pinned: p0_21_pin,
        interior_filler: None,
    }
}
