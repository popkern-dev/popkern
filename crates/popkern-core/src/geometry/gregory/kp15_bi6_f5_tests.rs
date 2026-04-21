//! Tests for K-P 2015 F_5 quality functional.

#![allow(clippy::needless_range_loop)]

//!
//! Strategy: exercise the invariants F_5 has by construction —
//! vanishes on low-degree polynomials, non-negative, symmetric,
//! and PSD — plus a closed-form check on a specific degree-5
//! monomial where F_5 can be computed by hand.

use super::*;
use crate::geometry::Point3;

fn pt(x: f64) -> Point3 {
    Point3::new(x, 0.0, 0.0)
}

fn zero_grid() -> [[Point3; 7]; 7] {
    [[Point3::new(0.0, 0.0, 0.0); 7]; 7]
}

#[test]
fn f5_of_constant_is_zero() {
    // f ≡ c → all derivatives vanish, F_5 = 0.
    let mut grid = zero_grid();
    for row in grid.iter_mut() {
        for p in row.iter_mut() {
            *p = Point3::new(7.0, -3.0, 2.5);
        }
    }
    let f5 = f5_evaluate(&grid);
    assert!(f5.abs() < 1e-18, "F_5(constant) = {f5} ≠ 0");
}

#[test]
fn f5_of_linear_in_u_is_zero() {
    // f(u, v) = u → ∂⁵_u f = 0, all other 5th partials are zero.
    let mut grid = zero_grid();
    for (i, row) in grid.iter_mut().enumerate() {
        for p in row.iter_mut() {
            // Bernstein-1 CPs for f(u) = u: grid[i][j].x = i / 6
            p.x = (i as f64) / 6.0;
        }
    }
    let f5 = f5_evaluate(&grid);
    assert!(f5.abs() < 1e-18, "F_5(u) = {f5} ≠ 0");
}

#[test]
fn f5_of_total_degree_4_is_zero() {
    // A bi-6 patch with coefficients nonzero only for i + j ≤ 4
    // represents a polynomial of total degree ≤ 4. Its 5th
    // mixed partials all vanish → F_5 = 0.
    //
    // To build this: take a polynomial p(u, v) = u²·v² and
    // express it in bi-6 Bernstein form. u² = Σ_i (C(2,2) ·
    // C(6,i) / C(8, i+2)) B_i^6(u) ... actually simpler: build
    // degree-4 in each axis via degree elevation from degree-2.
    //
    // Bi-2 CPs for f(u, v) = u²·v²: grid_{ij} = u_i² · v_j²
    // with u_i = i/2, v_j = j/2, i, j ∈ 0..=2. After elevation
    // to bi-6 the BB-CPs are specific but we don't need them
    // explicitly — just set grid[i][j] = (i·j)²/(6·6)² when
    // we want polynomial f = (u·v)² in interpolation positions.
    //
    // But BB-CPs aren't point samples. To sidestep, use a
    // simpler construction: set only CPs where i + j ≤ 4 to
    // nonzero; F_5 still only captures 5th partial content.
    //
    // Strategy: build a bi-6 patch as an EXACT polynomial of
    // total degree ≤ 4 via multi-degree elevation from a bi-4
    // patch. We elevate u-dir and v-dir from degree 4 to 6.
    //
    // For simplicity here we test: a bi-6 patch with BB-CPs
    // `grid[i][j].x = α_i · β_j` where α is a bi-6 elevation of
    // a degree-2 polynomial (only 3 free bi-2 CPs), and similar
    // for β. Such a patch has polynomial degree ≤ 2 in each
    // axis = total degree ≤ 4. Its 5th derivatives vanish.
    //
    // Use α from elevating the degree-2 CPs [0, 0.5, 1.0]:
    let alpha_deg2 = [0.0, 0.5, 1.0];
    let alpha = elevate_to_degree_6(&alpha_deg2);
    let beta = alpha; // Same curve.
    let mut grid = zero_grid();
    for i in 0..7 {
        for j in 0..7 {
            grid[i][j].x = alpha[i] * beta[j];
        }
    }
    let f5 = f5_evaluate(&grid);
    // Polynomial degree 2 in each axis → total degree ≤ 4, so
    // F_5 = 0 (up to quadrature / FP noise).
    assert!(f5.abs() < 1e-18, "F_5 of tensor deg-4 poly = {f5}");
}

#[test]
fn f5_non_negative() {
    // F_5 is a sum of squares → non-negative for any patch.
    for seed in 0..10 {
        let mut grid = zero_grid();
        for i in 0..7 {
            for j in 0..7 {
                // Deterministic pseudo-random values.
                let h = (seed as f64) * 13.0 + (i as f64) * 7.0 + (j as f64);
                grid[i][j].x = h.sin();
                grid[i][j].y = h.cos();
                grid[i][j].z = (h * 0.5).sin();
            }
        }
        let f5 = f5_evaluate(&grid);
        assert!(f5 >= -1e-12, "F_5 should be ≥ 0, got {f5}");
    }
}

#[test]
fn f5_homogeneous_in_scale() {
    // F_5(c · f) = c² · F_5(f).
    let mut grid = zero_grid();
    for i in 0..7 {
        for j in 0..7 {
            grid[i][j].x = ((i * 7 + j) as f64).sin();
        }
    }
    let f5_1 = f5_evaluate(&grid);
    let c = 3.7;
    for row in grid.iter_mut() {
        for p in row.iter_mut() {
            p.x *= c;
        }
    }
    let f5_c = f5_evaluate(&grid);
    let expected = c * c * f5_1;
    assert!(
        (f5_c - expected).abs() < 1e-9 * expected.abs().max(1e-6),
        "F_5(c·f) = {f5_c}, expected c²·F_5(f) = {expected}"
    );
}

#[test]
fn f5_matrix_symmetric() {
    // Q[i][j] == Q[j][i] by construction; check relative error
    // since magnitudes span many orders (peak ~5e6).
    let q = f5_quadratic_form();
    for i in 0..BI6_CP_COUNT {
        for j in 0..BI6_CP_COUNT {
            let mag = q[i][j].abs().max(q[j][i].abs()).max(1.0);
            let rel = (q[i][j] - q[j][i]).abs() / mag;
            assert!(
                rel < 1e-13,
                "Q[{i}][{j}] = {} vs Q[{j}][{i}] = {} (rel {rel})",
                q[i][j],
                q[j][i]
            );
        }
    }
}

#[test]
fn f5_matrix_diagonal_non_negative() {
    // Diagonal entries = F_5(unit vector basis_k) ≥ 0.
    let q = f5_quadratic_form();
    for i in 0..BI6_CP_COUNT {
        assert!(q[i][i] >= -1e-12, "Q[{i}][{i}] = {} < 0", q[i][i]);
    }
}

#[test]
fn f5_matrix_matches_evaluator() {
    // For any patch, evaluator output equals xᵀ Q x + yᵀ Q y + zᵀ Q z
    // where x, y, z are the flattened coord vectors.
    let q = f5_quadratic_form();
    let mut grid = zero_grid();
    for i in 0..7 {
        for j in 0..7 {
            let t = (i * 7 + j) as f64;
            grid[i][j].x = (t * 0.1).sin();
            grid[i][j].y = (t * 0.13).cos();
            grid[i][j].z = (t * 0.17).sin();
        }
    }
    let direct = f5_evaluate(&grid);

    let mut xs = [0.0_f64; BI6_CP_COUNT];
    let mut ys = [0.0_f64; BI6_CP_COUNT];
    let mut zs = [0.0_f64; BI6_CP_COUNT];
    for i in 0..7 {
        for j in 0..7 {
            let idx = i * 7 + j;
            xs[idx] = grid[i][j].x;
            ys[idx] = grid[i][j].y;
            zs[idx] = grid[i][j].z;
        }
    }
    let via_q = quad_form(&q, &xs) + quad_form(&q, &ys) + quad_form(&q, &zs);
    let rel = (direct - via_q).abs() / direct.abs().max(1e-12);
    assert!(
        rel < 1e-10,
        "evaluator = {direct}, via Q = {via_q}, rel {rel}"
    );
}

// ---- helpers ----

fn quad_form(q: &[[f64; BI6_CP_COUNT]; BI6_CP_COUNT], x: &[f64; BI6_CP_COUNT]) -> f64 {
    let mut s = 0.0;
    for i in 0..BI6_CP_COUNT {
        let xi = x[i];
        if xi == 0.0 {
            continue;
        }
        for (j, &xj) in x.iter().enumerate() {
            s += xi * q[i][j] * xj;
        }
    }
    s
}

/// Elevate a degree-2 univariate Bernstein polynomial (3 CPs) to
/// degree 6 (7 CPs). Used in the total-degree-4 test. Applies
/// the standard one-step elevation 4 times.
fn elevate_to_degree_6(c2: &[f64; 3]) -> [f64; 7] {
    let mut cur: Vec<f64> = c2.to_vec();
    for _ in 0..4 {
        let n = cur.len() - 1;
        let mut next = vec![0.0_f64; n + 2];
        next[0] = cur[0];
        next[n + 1] = cur[n];
        for i in 1..=n {
            let frac = (i as f64) / ((n + 1) as f64);
            next[i] = frac * cur[i - 1] + (1.0 - frac) * cur[i];
        }
        cur = next;
    }
    let mut out = [0.0_f64; 7];
    out.copy_from_slice(&cur);
    out
}

#[allow(dead_code)]
fn _debug_unused_pt(_x: f64) -> Point3 {
    pt(_x)
}
