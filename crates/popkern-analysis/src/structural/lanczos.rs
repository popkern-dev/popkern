//! Shift-invert Lanczos for K x = λ M x, grown a few steps at a time
//! until the wanted modes converge (map Stage 2, item 4).
//!
//! The operator is A = (K − σM)⁻¹M, self-adjoint in the M inner product;
//! each step is one multigrid PCG solve, and the basis is kept whole for
//! full reorthogonalisation (classical Gram–Schmidt, twice). The
//! residual is Lanczos's own, ‖(A − θ)x‖_M / θ = β|y_last|/θ: the
//! Euclidean ‖Kx − λMx‖ is dominated by stiff rows (Nitsche, extended
//! dofs) far below the mode's accuracy.
//!
//! One start vector spans one direction per distinct eigenvalue, so a
//! run alone misses the partners of repeated modes (a free part's six
//! rigid modes, the pairs of symmetric parts). A run may be deflated
//! against locked modes: its start and every step are kept
//! M-orthogonal to them, so it finds what earlier runs missed.

use crate::error::AnalysisError;
use crate::linalg::dense_eigen::symmetric_eigen_dense;
use crate::linalg::sparse::Csr;
use crate::solve::Hierarchy;
use crate::solve::pcg::pcg;

use super::internal;

/// Ritz pairs ((λ, residual), x).
pub(super) type Ritz = (Vec<(f64, f64)>, Vec<Vec<f64>>);

/// A Lanczos run in progress.
pub(super) struct Lanczos<'a> {
    m: &'a Csr,
    hierarchy: &'a Hierarchy,
    sigma: f64,
    max_iters: usize,
    /// The M-orthonormal basis and M times each basis vector.
    basis: Vec<Vec<f64>>,
    m_basis: Vec<Vec<f64>>,
    alpha: Vec<f64>,
    beta: Vec<f64>,
    /// The β after the last step (the residual's scale).
    last_beta: f64,
    /// The next basis vector, once its β is known.
    next: Option<(Vec<f64>, Vec<f64>)>,
    /// An invariant subspace was found: no further steps.
    exhausted: bool,
    /// Locked modes x and M x, kept out of the basis.
    locked: &'a [(Vec<f64>, Vec<f64>)],
    /// PCG iterations spent.
    pub iterations: usize,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

impl<'a> Lanczos<'a> {
    /// A run about `sigma` (the hierarchy is K − σM) from the fixed
    /// start `seed` (every component), deflated against `locked` (each
    /// x M-normalised, with M x).
    pub(super) fn new(
        m: &'a Csr,
        hierarchy: &'a Hierarchy,
        (sigma, max_iters): (f64, usize),
        seed: usize,
        locked: &'a [(Vec<f64>, Vec<f64>)],
    ) -> Result<Self, AnalysisError> {
        let n = m.rows();
        let mut q: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.37 * ((i * 7919 + seed * 104_729) % 101) as f64 / 101.0)
            .collect();
        for _ in 0..2 {
            deflate(&mut q, locked);
        }
        let mq = mul(m, &q)?;
        let norm = dot(&q, &mq).sqrt();
        let start = (
            q.iter().map(|v| v / norm).collect(),
            mq.iter().map(|v| v / norm).collect(),
        );
        Ok(Self {
            m,
            hierarchy,
            sigma,
            max_iters,
            basis: Vec::new(),
            m_basis: Vec::new(),
            alpha: Vec::new(),
            beta: Vec::new(),
            last_beta: 0.0,
            next: Some(start),
            exhausted: norm.is_nan() || norm <= 0.0,
            locked,
            iterations: 0,
        })
    }

    /// Steps taken.
    pub(super) fn steps(&self) -> usize {
        self.alpha.len()
    }

    /// Whether no further step can be taken.
    pub(super) fn exhausted(&self) -> bool {
        self.exhausted
    }

    /// Steps until `steps` are taken (or the space is exhausted).
    pub(super) fn grow(&mut self, steps: usize) -> Result<(), AnalysisError> {
        while self.alpha.len() < steps && !self.exhausted {
            let Some((q, mq)) = self.next.take() else {
                self.exhausted = true;
                break;
            };
            self.basis.push(q);
            self.m_basis.push(mq);
            let j = self.basis.len() - 1;
            let solution = pcg(self.hierarchy, &self.m_basis[j], 1e-12, self.max_iters)
                .map_err(|e| internal("pcg", e))?;
            if !solution.converged {
                return Err(AnalysisError::SolverNotConverged {
                    history: solution.history,
                });
            }
            self.iterations += solution.iterations;
            let mut w = solution.x;
            let a = dot(&self.m_basis[j], &w);
            self.alpha.push(a);
            for _ in 0..2 {
                deflate(&mut w, self.locked);
                let c: Vec<f64> = self.m_basis.iter().map(|mv| dot(mv, &w)).collect();
                for (v, cv) in self.basis.iter().zip(&c) {
                    w.iter_mut().zip(v).for_each(|(x, y)| *x -= cv * y);
                }
            }
            let mw = mul(self.m, &w)?;
            let b = dot(&w, &mw).max(0.0).sqrt();
            self.last_beta = b;
            if b <= 1e-14 * a.abs().max(f64::MIN_POSITIVE) || self.basis.len() >= self.m.rows() {
                self.exhausted = true;
                break;
            }
            self.beta.push(b);
            self.next = Some((
                w.iter().map(|x| x / b).collect(),
                mw.iter().map(|x| x / b).collect(),
            ));
        }
        Ok(())
    }

    /// The `count` lowest ((λ, residual), x) of the steps so far, x
    /// M-normalised (empty unless `vectors`); λ = σ + 1/θ for the
    /// largest Ritz values θ of A.
    pub(super) fn ritz(&self, count: usize, vectors: bool) -> Ritz {
        let size = self.alpha.len();
        let mut t = vec![0.0; size * size];
        for i in 0..size {
            t[i * size + i] = self.alpha[i];
            if i + 1 < size
                && let Some(&b) = self.beta.get(i)
            {
                t[i * size + i + 1] = b;
                t[(i + 1) * size + i] = b;
            }
        }
        let (thetas, ys) = symmetric_eigen_dense(t, size);
        let n = self.m.rows();
        let (mut values, mut xs) = (Vec::new(), Vec::new());
        for col in (0..size).rev().take(count) {
            let theta = thetas[col];
            if theta <= 0.0 {
                break;
            }
            if vectors {
                let mut x = vec![0.0; n];
                for (r, v) in self.basis.iter().enumerate().take(size) {
                    let y = ys[r * size + col];
                    x.iter_mut().zip(v).for_each(|(a, b)| *a += y * b);
                }
                xs.push(x);
            }
            let residual = self.last_beta * ys[(size - 1) * size + col].abs() / theta;
            values.push((self.sigma + 1.0 / theta, residual));
        }
        (values, xs)
    }
}

/// Removes the M-components of `v` along each locked mode.
fn deflate(v: &mut [f64], locked: &[(Vec<f64>, Vec<f64>)]) {
    for (x, mx) in locked {
        let c = dot(mx, v);
        v.iter_mut().zip(x).for_each(|(a, b)| *a -= c * b);
    }
}

fn mul(m: &Csr, v: &[f64]) -> Result<Vec<f64>, AnalysisError> {
    let mut out = vec![0.0; m.rows()];
    m.mul_vec(v, &mut out).map_err(|e| internal("modal", e))?;
    Ok(out)
}
