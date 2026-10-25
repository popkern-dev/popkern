//! Conjugate gradients preconditioned by one multigrid V-cycle.

use crate::linalg::sparse::Csr;

use super::{Hierarchy, MultigridError};

/// The outcome of a PCG solve.
#[derive(Debug, Clone, PartialEq)]
pub struct PcgResult {
    /// The solution, or the best iterate when not converged (untrusted).
    pub x: Vec<f64>,
    /// ‖r‖/‖b‖ after each iteration (the first entry is the initial
    /// residual).
    pub history: Vec<f64>,
    /// Iterations taken.
    pub iterations: usize,
    /// True when the relative residual reached the tolerance.
    pub converged: bool,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solves A x = b from x = 0 with `hierarchy` (whose finest operator is
/// `A`) as preconditioner, to relative residual `rel_tol` within
/// `max_iters` iterations.
pub fn pcg(
    hierarchy: &Hierarchy,
    b: &[f64],
    rel_tol: f64,
    max_iters: usize,
) -> Result<PcgResult, MultigridError> {
    let a: &Csr = &hierarchy.levels[0].a;
    let n = b.len();
    let b_norm = dot(b, b).sqrt();
    let mut x = vec![0.0; n];
    if b_norm == 0.0 {
        return Ok(PcgResult {
            x,
            history: vec![0.0],
            iterations: 0,
            converged: true,
        });
    }
    let mut r = b.to_vec();
    let mut z = hierarchy.v_cycle(0, &r)?;
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut history = vec![1.0];
    let mut best = (1.0, x.clone());
    let mut q = vec![0.0; n];
    for iteration in 1..=max_iters {
        a.mul_vec(&p, &mut q)?;
        let alpha = rz / dot(&p, &q);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        }
        let rel = dot(&r, &r).sqrt() / b_norm;
        history.push(rel);
        if rel < best.0 {
            best = (rel, x.clone());
        }
        if rel <= rel_tol {
            return Ok(PcgResult {
                x,
                history,
                iterations: iteration,
                converged: true,
            });
        }
        z = hierarchy.v_cycle(0, &r)?;
        let rz_next = dot(&r, &z);
        let beta = rz_next / rz;
        rz = rz_next;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    Ok(PcgResult {
        x: best.1,
        iterations: max_iters,
        history,
        converged: false,
    })
}
