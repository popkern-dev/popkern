//! MMA against Svanberg 1987 §6, test problem 1: the five-element
//! cantilever of least weight C₁ Σ x_j under a tip-deflection limit
//! 61/x₁³ + 37/x₂³ + 19/x₃³ + 7/x₄³ + 1/x₅³ ≤ C₂, C₁ = 0.0624, C₂ = 1,
//! from x = 5 (feasible, weight 1.560). The optimum is
//! x = (6.016, 5.309, 4.494, 3.502, 2.153), weight 1.340 (eq 22).

use super::*;

const C1: f64 = 0.0624;
const C2: f64 = 1.0;
const A: [f64; 5] = [61.0, 37.0, 19.0, 7.0, 1.0];
const OPTIMUM: [f64; 5] = [6.016, 5.309, 4.494, 3.502, 2.153];

fn weight(x: &[f64]) -> f64 {
    C1 * x.iter().sum::<f64>()
}

/// f₁ = Σ a/x³ − C₂ and its gradient.
fn deflection(x: &[f64]) -> (f64, Vec<f64>) {
    let f = A.iter().zip(x).map(|(a, x)| a / (x * x * x)).sum::<f64>() - C2;
    let g = A
        .iter()
        .zip(x)
        .map(|(a, x)| -3.0 * a / (x * x * x * x))
        .collect();
    (f, g)
}

/// Eq 25.
fn infeasibility(x: &[f64]) -> f64 {
    (deflection(x).0 / C2).max(0.0)
}

/// Iterates with eq 23's asymptotes L = t x, U = x / t and eq 24's move
/// limits; `t = 0` is the 'traditional' method (L = 0, U → ∞).
fn run(t: f64, updates: usize) -> Vec<(f64, f64)> {
    let mut x = vec![5.0; 5];
    let mut out = Vec::new();
    for _ in 0..updates {
        let lower: Vec<f64> = x.iter().map(|x| t * x).collect();
        let upper: Vec<f64> = x
            .iter()
            .map(|x| if t > 0.0 { x / t } else { 1e12 * x })
            .collect();
        let alpha: Vec<f64> = x
            .iter()
            .zip(&lower)
            .map(|(x, l)| (0.5 * x).max(1.01 * l))
            .collect();
        let beta: Vec<f64> = x
            .iter()
            .zip(&upper)
            .map(|(x, u)| (2.0 * x).min(0.99 * u))
            .collect();
        let (f1, df1) = deflection(&x);
        let df0 = vec![C1; 5];
        let (next, _) = Subproblem {
            x: &x,
            lower: &lower,
            upper: &upper,
            alpha: &alpha,
            beta: &beta,
            df0: &df0,
            f1,
            df1: &df1,
            exact: None,
        }
        .solve();
        x = next;
        out.push((weight(&x), infeasibility(&x)));
    }
    out
}

/// Table I, t = 1/4: the iterates' weight and infeasibility to the
/// table's digits, the optimum reached at the third.
#[test]
fn mma_reproduces_svanbergs_table_one() {
    let got = run(0.25, 3);
    let table = [(1.309, 0.10), (1.335, 0.01), (1.340, 0.0005)];
    for (k, ((w, i), (tw, ti))) in got.iter().zip(table).enumerate() {
        assert!((w - tw).abs() < 6e-4, "iteration {}: {w} vs {tw}", k + 1);
        assert!((i - ti).abs() < 0.006, "iteration {}: {i} vs {ti}", k + 1);
    }
    // Every t of the table (1/16 to 3/4) reaches the optimum.
    for t in [1.0 / 16.0, 0.125, 0.25, 1.0 / 3.0, 0.5, 2.0 / 3.0, 0.75] {
        let &(w, i) = run(t, 20).last().unwrap_or(&(0.0, 1.0));
        assert!((w / 1.340 - 1.0).abs() < 1e-3 && i < 1e-3, "t {t}: {w} {i}");
    }
}

/// The negative control: the 'traditional' method (t = 0) swings
/// between two infeasible points without converging (Table I: weight
/// 1.25–1.26, infeasibility 0.42–0.44 from the second iteration on).
#[test]
fn the_traditional_method_oscillates() {
    let got = run(0.0, 13);
    for &(w, i) in &got[1..] {
        assert!(i > 0.35 && w < 1.3, "{w} {i}");
    }
}

/// The moving asymptotes of eqs 11–13 on the box [1, 10] reach the
/// optimum point from x = 5.
#[test]
fn moving_asymptotes_reach_the_optimum() {
    let mut mma = Mma::new(vec![1.0; 5], vec![10.0; 5], 1.0);
    let mut x = vec![5.0; 5];
    for _ in 0..40 {
        let (f1, df1) = deflection(&x);
        x = mma.update(&x, &[C1; 5], (f1, &df1), None);
    }
    for (x, o) in x.iter().zip(OPTIMUM) {
        assert!((x - o).abs() < 2e-3, "{x} vs {o}");
    }
    assert!(infeasibility(&x) < 1e-6);
}

/// When no point of the move box meets the constraint, the update takes
/// the box's point closest to meeting it (§5 with an unbounded cost):
/// every height at the top of its box, whatever the weight.
#[test]
fn an_unreachable_constraint_gives_the_closest_point() {
    // Deflection at most 0.01 of C₂: far out of reach of a 0.1 move.
    let mut mma = Mma::new(vec![1.0; 5], vec![10.0; 5], 0.1);
    let x = vec![5.0; 5];
    let (f1, df1) = deflection(&x);
    let tight = f1 + C2 - 0.01;
    let next = mma.update(&x, &[C1; 5], (tight, &df1), None);
    assert!(next.iter().all(|v| v.is_finite()), "{next:?}");
    for (n, x) in next.iter().zip(&x) {
        assert!((n - (x + 0.1 * 9.0)).abs() < 1e-9, "{n}");
    }
    let exact = |z: &[f64]| deflection(z).0 + C2 - 0.01;
    let next = mma.update(&x, &[C1; 5], (tight, &df1), Some(&exact));
    assert!(next.iter().all(|v| v.is_finite()), "{next:?}");
}
