//! Weighted normal equations for the ring-wide cross-field fit
//! (`kp20_paper_border_ring_fit`): sparse residual rows over the
//! scalar unknowns `3·index + axis`, accumulated into `AᵀA` / `Aᵀb`,
//! solved by Cholesky with a scale-free ridge.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use popkern_core::geometry::Vec3;

use super::ring_fit::Affine;

/// Relative Tikhonov ridge on the normal equations (scale-free).
const RIDGE: f64 = 1.0e-12;

/// Sparse residual row over the scalar unknowns (`3·index + axis`).
struct Residual {
    entries: Vec<(usize, f64)>,
    rhs: f64,
    weight: f64,
}

pub(super) struct System {
    unknowns: usize,
    ata: Vec<Vec<f64>>,
    atb: Vec<f64>,
}

impl System {
    pub(super) fn new(vec_unknowns: usize) -> Self {
        let n = 3 * vec_unknowns;
        Self {
            unknowns: n,
            ata: vec![vec![0.0; n]; n],
            atb: vec![0.0; n],
        }
    }
    fn add(&mut self, r: &Residual) {
        let w2 = r.weight * r.weight;
        for &(i, a) in &r.entries {
            for &(j, b) in &r.entries {
                self.ata[i][j] += w2 * a * b;
            }
            self.atb[i] += w2 * a * r.rhs;
        }
    }
    /// Scalar-component rows `value[axis] = target[axis]` (three rows).
    pub(super) fn add_vector(&mut self, value: &Affine, target: Vec3, weight: f64) {
        for axis in 0..3 {
            let entries: Vec<(usize, f64)> = value
                .terms
                .iter()
                .map(|&(i, w)| (3 * i + axis, w))
                .collect();
            let rhs = component(target, axis) - component(value.constant, axis);
            self.add(&Residual {
                entries,
                rhs,
                weight,
            });
        }
    }
    /// One row `n · value = 0`.
    pub(super) fn add_normal(&mut self, value: &Affine, n: Vec3, weight: f64) {
        let mut entries = Vec::with_capacity(3 * value.terms.len());
        for &(i, w) in &value.terms {
            for axis in 0..3 {
                entries.push((3 * i + axis, w * component(n, axis)));
            }
        }
        let rhs = -n.dot(&value.constant);
        self.add(&Residual {
            entries,
            rhs,
            weight,
        });
    }
    /// Cholesky solve of the ridge-regularised normal equations.
    pub(super) fn solve(mut self) -> Option<Vec<Vec3>> {
        let n = self.unknowns;
        let trace: f64 = (0..n).map(|i| self.ata[i][i]).sum();
        let ridge = RIDGE * trace / n.max(1) as f64;
        for i in 0..n {
            self.ata[i][i] += ridge;
        }
        // In-place Cholesky (lower).
        for j in 0..n {
            let mut d = self.ata[j][j];
            for k in 0..j {
                d -= self.ata[j][k] * self.ata[j][k];
            }
            if d <= 0.0 || !d.is_finite() {
                return None;
            }
            let d = d.sqrt();
            self.ata[j][j] = d;
            for i in (j + 1)..n {
                let mut s = self.ata[i][j];
                for k in 0..j {
                    s -= self.ata[i][k] * self.ata[j][k];
                }
                self.ata[i][j] = s / d;
            }
        }
        // Forward then back substitution.
        let mut y = vec![0.0; n];
        for i in 0..n {
            let mut s = self.atb[i];
            for (k, yk) in y.iter().enumerate().take(i) {
                s -= self.ata[i][k] * yk;
            }
            y[i] = s / self.ata[i][i];
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut s = y[i];
            for (k, xk) in x.iter().enumerate().skip(i + 1) {
                s -= self.ata[k][i] * xk;
            }
            x[i] = s / self.ata[i][i];
        }
        if x.iter().any(|v| !v.is_finite()) {
            return None;
        }
        Some(
            (0..n / 3)
                .map(|i| Vec3::new(x[3 * i], x[3 * i + 1], x[3 * i + 2]))
                .collect(),
        )
    }
}

pub(super) fn component(v: Vec3, axis: usize) -> f64 {
    match axis {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}
