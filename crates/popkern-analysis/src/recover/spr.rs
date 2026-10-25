//! Superconvergent patch recovery: a weighted linear least-squares fit
//! of the six stress components over a node's 2×2×2 cell patch.

use crate::linalg::Cholesky;

use super::tensor::Voigt;

/// Smallest pivot of the scaled normal matrix, relative to its
/// diagonal, for a fit to count as well spread (coplanar samples, as in
/// a wall one cell thick, fall below it).
pub const SPREAD_RATIO: f64 = 1.0e-4;

/// One sample: fragment centroid (physical), weight (fragment volume)
/// and stress there.
#[derive(Debug, Clone, Copy)]
pub struct Sample {
    /// Position.
    pub x: [f64; 3],
    /// Weight.
    pub w: f64,
    /// Stress.
    pub s: Voigt,
}

/// A linear stress field σ(x) = a + B·(x − origin)/h per component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearFit {
    origin: [f64; 3],
    h: f64,
    /// [component][1, x, y, z] coefficients.
    c: [[f64; 4]; 6],
}

impl LinearFit {
    /// The fitted stress at `x`.
    pub fn at(&self, x: [f64; 3]) -> Voigt {
        let d = [0, 1, 2].map(|a| (x[a] - self.origin[a]) / self.h);
        std::array::from_fn(|k| {
            self.c[k][0] + self.c[k][1] * d[0] + self.c[k][2] * d[1] + self.c[k][3] * d[2]
        })
    }
}

/// The weighted linear fit of `samples` about `origin` (coordinates
/// scaled by `h`), or `None` when the samples are too few or too close
/// to coplanar.
pub fn fit(samples: &[Sample], origin: [f64; 3], h: f64) -> Option<LinearFit> {
    if samples.len() < 4 {
        return None;
    }
    let mut a = [0.0; 16];
    let mut rhs = [[0.0; 4]; 6];
    for s in samples {
        let p = [
            1.0,
            (s.x[0] - origin[0]) / h,
            (s.x[1] - origin[1]) / h,
            (s.x[2] - origin[2]) / h,
        ];
        for r in 0..4 {
            for c in 0..4 {
                a[4 * r + c] += s.w * p[r] * p[c];
            }
            for (k, row) in rhs.iter_mut().enumerate() {
                row[r] += s.w * p[r] * s.s[k];
            }
        }
    }
    let chol = Cholesky::factor(&a, 4, SPREAD_RATIO).ok()?;
    let mut c = [[0.0; 4]; 6];
    for (k, row) in rhs.iter().enumerate() {
        let x = chol.solve(row).ok()?;
        c[k].copy_from_slice(&x);
    }
    Some(LinearFit { origin, h, c })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_linear_field_is_fitted_exactly_and_coplanar_samples_refused() {
        let field = |x: [f64; 3]| -> Voigt {
            std::array::from_fn(|k| k as f64 + 0.5 * x[0] - 0.25 * x[1] + 2.0 * x[2])
        };
        let mut samples = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let x = [
                        0.3 + 0.5 * i as f64,
                        -0.2 + 0.5 * j as f64,
                        0.1 + 0.5 * k as f64,
                    ];
                    samples.push(Sample {
                        x,
                        w: 1.0 + i as f64,
                        s: field(x),
                    });
                }
            }
        }
        let f = fit(&samples, [0.5, 0.0, 0.3], 0.5).unwrap_or_else(|| panic!("fit"));
        let at = [0.9, 0.4, -0.1];
        for (a, b) in f.at(at).iter().zip(field(at)) {
            assert!((a - b).abs() < 1e-12);
        }
        // All samples in one plane: refused.
        let flat: Vec<Sample> = samples.iter().filter(|s| s.x[2] < 0.2).copied().collect();
        assert!(fit(&flat, [0.5, 0.0, 0.3], 0.5).is_none());
    }
}
