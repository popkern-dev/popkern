//! PK-572: the banded fit/interpolation solves must be BIT-IDENTICAL
//! to the dense construction they replaced — every fitted SSI record
//! downstream (seam anchoring, pcurve registration, validation
//! verdicts) is pinned to those bits. Each test rebuilds the dense
//! system from the returned curve's own knots, solves it with the
//! dense reference solver, and compares control points bitwise.

use super::numerics::{basis_function, solve_linear_system};
use super::{fit_curve, interpolate_curve};
use crate::geometry::Point3;

/// Deterministic LCG in [0, 1).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A noisy helix with chord-length parameters normalized to [0, 1].
fn noisy_helix(m: usize, seed: u64) -> (Vec<Point3>, Vec<f64>) {
    let mut rng = Lcg(seed);
    let points: Vec<Point3> = (0..m)
        .map(|i| {
            let t = i as f64 / (m - 1) as f64 * 9.0;
            Point3::new(
                t.cos() + 1e-3 * rng.next(),
                t.sin() + 1e-3 * rng.next(),
                0.1 * t + 1e-3 * rng.next(),
            )
        })
        .collect();
    let mut params = vec![0.0];
    let mut total = 0.0;
    for w in points.windows(2) {
        total += w[0].distance_to(&w[1]);
        params.push(total);
    }
    for p in &mut params {
        *p /= total;
    }
    (points, params)
}

fn bits(p: &Point3) -> [u64; 3] {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

#[test]
fn banded_interpolation_is_bit_identical_to_dense() {
    for (m, degree, seed) in [(4, 3, 1), (9, 3, 2), (64, 3, 3), (65, 5, 4), (200, 2, 5)] {
        let (points, params) = noisy_helix(m, seed);
        let curve = interpolate_curve(&points, &params, degree).expect("interpolate");
        let knots = &curve.knots;
        let a: Vec<Vec<f64>> = params
            .iter()
            .map(|&t| {
                (0..m)
                    .map(|j| basis_function(j, degree, t, knots))
                    .collect()
            })
            .collect();
        let coord = |f: fn(&Point3) -> f64| points.iter().map(f).collect::<Vec<f64>>();
        let cx = solve_linear_system(&a, &coord(|p| p.x)).expect("dense x");
        let cy = solve_linear_system(&a, &coord(|p| p.y)).expect("dense y");
        let cz = solve_linear_system(&a, &coord(|p| p.z)).expect("dense z");
        for (i, cp) in curve.control_points.iter().enumerate() {
            let dense = Point3::new(cx[i], cy[i], cz[i]);
            assert_eq!(bits(cp), bits(&dense), "m={m} degree={degree} cp {i}");
        }
    }
}

#[test]
fn banded_least_squares_fit_is_bit_identical_to_dense() {
    for (m, degree, n, seed) in [
        (10, 3, 6, 11),
        (100, 3, 20, 12),
        (500, 3, 100, 13),
        (500, 3, 250, 14),
        (500, 3, 500, 17),
        (300, 2, 150, 15),
        (120, 5, 60, 16),
    ] {
        let (points, params) = noisy_helix(m, seed);
        let curve = fit_curve(&points, &params, degree, n)
            .unwrap_or_else(|e| panic!("fit m={m} n={n} degree={degree}: {e}"));
        let knots = &curve.knots;
        let cp = &curve.control_points;
        let (rows, cols) = (m - 2, n - 2);
        // Dense normal equations exactly as the pre-banding code
        // formed them: full basis rows, row-outer accumulation.
        let basis: Vec<Vec<f64>> = (0..rows)
            .map(|i| {
                (0..n)
                    .map(|j| basis_function(j, degree, params[i + 1], knots))
                    .collect()
            })
            .collect();
        let rhs = |f: fn(&Point3) -> f64| -> Vec<f64> {
            (0..rows)
                .map(|i| {
                    f(&points[i + 1]) - basis[i][0] * f(&cp[0]) - basis[i][n - 1] * f(&cp[n - 1])
                })
                .collect()
        };
        let (rx, ry, rz) = (rhs(|p| p.x), rhs(|p| p.y), rhs(|p| p.z));
        let mut ntn = vec![vec![0.0; cols]; cols];
        let (mut nx, mut ny, mut nz) = (vec![0.0; cols], vec![0.0; cols], vec![0.0; cols]);
        for k in 0..rows {
            for i in 0..cols {
                let vi = basis[k][i + 1];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..cols {
                    let vj = basis[k][j + 1];
                    if vj != 0.0 {
                        ntn[i][j] += vi * vj;
                    }
                }
                nx[i] += vi * rx[k];
                ny[i] += vi * ry[k];
                nz[i] += vi * rz[k];
            }
        }
        let sx = solve_linear_system(&ntn, &nx).expect("dense x");
        let sy = solve_linear_system(&ntn, &ny).expect("dense y");
        let sz = solve_linear_system(&ntn, &nz).expect("dense z");
        for j in 0..cols {
            let dense = Point3::new(sx[j], sy[j], sz[j]);
            assert_eq!(bits(&cp[j + 1]), bits(&dense), "m={m} n={n} cp {}", j + 1);
        }
    }
}

/// The OOM class itself: a trace-length fit near full interpolation
/// must run in bounded memory. Dense storage of this system is
/// 20k² × 8 B = 3.2 GB per copy; banded it is a few MB.
#[test]
fn long_trace_fit_near_interpolation_is_memory_bounded() {
    let m = 20_000;
    let (points, params) = noisy_helix(m, 21);
    let curve = fit_curve(&points, &params, 3, m).expect("long fit");
    let worst = points
        .iter()
        .zip(&params)
        .map(|(p, &t)| {
            use crate::geometry::curve::CurveEval as _;
            p.distance_to(&curve.point_at(t))
        })
        .fold(0.0_f64, f64::max);
    assert!(
        worst < 1e-9,
        "full-count fit must interpolate, worst {worst:e}"
    );
}
