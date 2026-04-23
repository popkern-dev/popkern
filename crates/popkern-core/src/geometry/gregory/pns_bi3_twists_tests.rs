//! Tests for [`pns_bi3_twists`] — the eq (7) pair-sum circulant
//! twist solver for PnS bi-3 caps at n ∈ {3, 5}.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

#[test]
fn zero_inputs_give_zero_twists() {
    let zero = pt(0.0, 0.0, 0.0);
    for &n in &[3_usize, 5] {
        let p01 = vec![zero; n];
        let p21 = vec![zero; n];
        let p31 = vec![zero; n];
        let result = solve_bi3_twists(zero, &p01, &p21, &p31, n);
        assert_eq!(result.len(), n);
        for &r in &result {
            assert!(approx_pt(r, zero, 1e-12), "n={n}");
        }
    }
}

/// The solved twists must satisfy the defining eq (7) exactly.
/// This is the round-trip check: substitute the outputs back into
/// the equation and confirm residuals are zero.
#[test]
fn solution_satisfies_eq7_identity() {
    // Random-ish inputs to exercise the solver.
    for &n in &[3_usize, 5] {
        let p00 = pt(0.1, 0.2, 0.3);
        let p01: Vec<Point3> = (0..n)
            .map(|k| pt(0.5 * k as f64, 0.3 + 0.1 * k as f64, 0.7 - 0.2 * k as f64))
            .collect();
        let p21: Vec<Point3> = (0..n)
            .map(|k| {
                pt(
                    1.1 + 0.2 * k as f64,
                    2.3 + 0.3 * k as f64,
                    0.5 + 0.1 * k as f64,
                )
            })
            .collect();
        let p31: Vec<Point3> = (0..n)
            .map(|k| {
                pt(
                    0.3 + 0.05 * k as f64,
                    0.4 - 0.1 * k as f64,
                    0.6 + 0.05 * k as f64,
                )
            })
            .collect();
        let p11 = solve_bi3_twists(p00, &p01, &p21, &p31, n);

        let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
        for k in 0..n {
            let kp1 = (k + 1) % n;
            // Eq (7):
            //   (3-c)·[3c·(p11[k] + p11[kp1]) − (3-2c)·(p01[k]+p01[kp1] − 2(1-c)·p00)]
            //   = −2c²·[c·(p31[k]+p31[kp1]) − 3·(p21[k]+p21[kp1])]
            let lhs_term1 = 3.0 * c;
            let pa = pt(
                p11[k].x + p11[kp1].x,
                p11[k].y + p11[kp1].y,
                p11[k].z + p11[kp1].z,
            );
            let two_1mc = 2.0 * (1.0 - c);
            let a_k = pt(
                p01[k].x + p01[kp1].x - two_1mc * p00.x,
                p01[k].y + p01[kp1].y - two_1mc * p00.y,
                p01[k].z + p01[kp1].z - two_1mc * p00.z,
            );
            let three_minus_2c = 3.0 - 2.0 * c;
            let lhs = pt(
                (3.0 - c) * (lhs_term1 * pa.x - three_minus_2c * a_k.x),
                (3.0 - c) * (lhs_term1 * pa.y - three_minus_2c * a_k.y),
                (3.0 - c) * (lhs_term1 * pa.z - three_minus_2c * a_k.z),
            );
            let b_k = pt(
                c * (p31[k].x + p31[kp1].x) - 3.0 * (p21[k].x + p21[kp1].x),
                c * (p31[k].y + p31[kp1].y) - 3.0 * (p21[k].y + p21[kp1].y),
                c * (p31[k].z + p31[kp1].z) - 3.0 * (p21[k].z + p21[kp1].z),
            );
            let rhs = pt(
                -2.0 * c * c * b_k.x,
                -2.0 * c * c * b_k.y,
                -2.0 * c * c * b_k.z,
            );
            let residual = pt(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
            let mag = (residual.x.powi(2) + residual.y.powi(2) + residual.z.powi(2)).sqrt();
            assert!(
                mag < 1e-10,
                "Eq (7) residual at n={n}, k={k}: {mag}, residual={residual:?}"
            );
        }
    }
}

#[test]
fn n3_constant_inputs_give_closed_form() {
    // If all per-sector inputs are the same constant Q (and p00
    // is also a specific value), the system reduces to a single
    // equation: x_k + x_{k+1} = rhs (same for all k), solved by
    // x_k = rhs / 2 (constant across k).
    let p00 = pt(1.0, 0.0, 0.0);
    let q01 = pt(2.0, 1.0, 0.5);
    let q21 = pt(3.0, 2.0, 0.5);
    let q31 = pt(4.0, 1.0, 0.0);
    let p01 = vec![q01; 3];
    let p21 = vec![q21; 3];
    let p31 = vec![q31; 3];
    let twists = solve_bi3_twists(p00, &p01, &p21, &p31, 3);
    // All three twists should match.
    for k in 1..3 {
        assert!(approx_pt(twists[0], twists[k], 1e-10), "k={k}");
    }
}
