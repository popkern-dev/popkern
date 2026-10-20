//! Exact NURBS curve derivatives against closed forms.

use crate::geometry::{NurbsCurve, Point3};

/// The rational quadratic quarter circle of radius `r` from (r, 0) to
/// (0, r). Its parameter is not arc length, so the test checks the
/// geometric invariants: curvature 1/r, and the curvature VECTOR kN
/// turning with the tangent, d(kN)/ds = -k^2 T.
fn quarter(r: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
        ],
        weights: vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0],
    }
}

#[test]
fn a_polynomial_cubics_derivatives_are_its_power_basis_coefficients() {
    // A Bezier cubic is the power-basis polynomial P(t) = a + b t + c t^2
    // + d t^3 with the control points of its Bernstein form.
    let (a, b, c, d) = (
        [1.0, 2.0, 0.5],
        [3.0, -1.0, 2.0],
        [-2.0, 4.0, 1.0],
        [0.5, -3.0, 2.5],
    );
    let at = |t: f64, k: usize| -> [f64; 3] {
        let mut out = [0.0; 3];
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = match k {
                0 => a[i] + b[i] * t + c[i] * t * t + d[i] * t * t * t,
                1 => b[i] + 2.0 * c[i] * t + 3.0 * d[i] * t * t,
                2 => 2.0 * c[i] + 6.0 * d[i] * t,
                3 => 6.0 * d[i],
                _ => 0.0,
            };
        }
        out
    };
    // Bernstein control points of the power-basis cubic.
    let p0 = a;
    let p1: Vec<f64> = (0..3).map(|i| a[i] + b[i] / 3.0).collect();
    let p2: Vec<f64> = (0..3)
        .map(|i| a[i] + 2.0 * b[i] / 3.0 + c[i] / 3.0)
        .collect();
    let p3: Vec<f64> = (0..3).map(|i| a[i] + b[i] + c[i] + d[i]).collect();
    let curve = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(p0[0], p0[1], p0[2]),
            Point3::new(p1[0], p1[1], p1[2]),
            Point3::new(p2[0], p2[1], p2[2]),
            Point3::new(p3[0], p3[1], p3[2]),
        ],
        weights: vec![1.0; 4],
    };
    for t in [0.0, 0.3, 0.71, 1.0] {
        let ders = curve.derivatives(t, 4);
        for (k, &got) in ders.iter().enumerate() {
            let want = at(t, k);
            let err =
                ((got.x - want[0]).powi(2) + (got.y - want[1]).powi(2) + (got.z - want[2]).powi(2))
                    .sqrt();
            assert!(err < 1e-12, "t {t} order {k}: {got:?} vs {want:?}");
        }
    }
}

#[test]
fn a_rational_quarter_circles_curvature_vector_turns_with_its_tangent() {
    let r = 2.5;
    let curve = quarter(r);
    for t in [0.1, 0.4, 0.5, 0.9] {
        let d = curve.derivatives(t, 3);
        let (x1, x2, x3) = (d[1], d[2], d[3]);
        let q = x1.dot(&x1);
        let kappa = x1.cross(&x2).length() / q.powf(1.5);
        assert!((kappa - 1.0 / r).abs() < 1e-12, "t {t}: kappa {kappa}");
        // The curvature vector's arc-length derivative: -k^2 T.
        let (h, hp) = (x1.dot(&x2), x2.dot(&x2) + x1.dot(&x3));
        let dk = (x3 * (q * q) - x2 * (3.0 * q * h) - x1 * (q * hp) + x1 * (4.0 * h * h))
            * (1.0 / (q * q * q * q.sqrt()));
        let want = x1 * (-(1.0 / (r * r)) / q.sqrt());
        assert!(
            (dk - want).length() < 1e-10,
            "t {t}: dk/ds {dk:?} vs {want:?}"
        );
    }
}
