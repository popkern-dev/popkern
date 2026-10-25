//! The lowest modes, repeated ones included: Lanczos runs deflated
//! against the modes already locked, restarted from new fixed starts
//! until a run finds nothing below the current highest wanted mode
//! (PH-093).
//!
//! A run grows until its lowest Ritz values are converged: all of the
//! wanted ones while fewer than `count` are locked, otherwise those
//! below the cut (the count-th locked λ), with the first one above it
//! bounded clear of the cut by its residual (`lower_bound`; PH-096: a
//! verifying run need not converge it fully). Ritz values approach from
//! above, so that bounds every missing mode of the run's space. A free part's
//! rigid modes come in locked (`modal::rigid_modes`): they are six-fold
//! repeated per solid and known exactly.

use crate::error::AnalysisError;
use crate::linalg::sparse::Csr;
use crate::solve::Hierarchy;

use super::lanczos::Lanczos;
use super::modal::MODE_TOLERANCE;

/// A locked mode: λ, residual, x (M-normalised) and M x.
pub(super) type Locked = (f64, f64, Vec<f64>, Vec<f64>);

/// Most restarts after the first run; a spectrum still gaining modes
/// past them declines.
pub const MAX_RESTARTS: usize = 12;

/// What the runs spent.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Spent {
    /// Lanczos runs.
    pub runs: usize,
    /// Lanczos steps over every run.
    pub steps: usize,
    /// PCG iterations over every run.
    pub iterations: usize,
}

/// Steps before a verifying run's first check, and at least between
/// checks.
const VERIFY_FIRST: usize = 4;
/// See [`VERIFY_FIRST`].
const VERIFY_EVERY: usize = 2;

/// A mode counts as new below the cut by more than this, relative to
/// the cut's magnitude (a partner of the cut's own repeated mode is not
/// new).
const NEW_BELOW: f64 = 1e-9;

/// The `count` lowest modes about `sigma`, starting from `locked`; at
/// most `max_steps` steps per run, and `restarts` restarts (0: the
/// first run's modes, unchecked, as before restarts existed).
pub(super) fn lowest_modes(
    m: &Csr,
    hierarchy: &Hierarchy,
    (sigma, max_iters): (f64, usize),
    (count, max_steps, restarts): (usize, usize, usize),
    mut locked: Vec<Locked>,
) -> Result<(Vec<Locked>, Spent), AnalysisError> {
    let n = m.rows();
    let count = count.min(n);
    let mut spent = Spent::default();
    sort(&mut locked);
    for seed in 0..=restarts {
        let cut = (locked.len() >= count).then(|| locked[count - 1].0);
        let pairs: Vec<(Vec<f64>, Vec<f64>)> = locked
            .iter()
            .map(|(_, _, x, mx)| (x.clone(), mx.clone()))
            .collect();
        let mut run = Lanczos::new(m, hierarchy, (sigma, max_iters), seed, &pairs)?;
        let want = match cut {
            Some(_) => count,
            None => count - locked.len(),
        };
        // A verifying run checks early and often (each check is a small
        // tridiagonal eigenproblem); a finding run first after want + 10.
        let (first, every) = match cut {
            Some(_) => (VERIFY_FIRST, VERIFY_EVERY),
            None => (want + 10, 10),
        };
        let mut target = first.min(max_steps);
        let found = loop {
            run.grow(target)?;
            let (values, _) = run.ritz(want, false);
            let (needed, settled) = match cut {
                Some(c) => {
                    // The values below the cut converged, and the first
                    // above it bounded clear of the cut.
                    let below = values.iter().filter(|v| v.0 < c).count();
                    let clear = values.get(below).is_none_or(|&(lambda, r)| {
                        r <= MODE_TOLERANCE || lower_bound(lambda, r, sigma) >= c
                    });
                    let settled = values[..below].iter().all(|v| v.1 <= MODE_TOLERANCE) && clear;
                    ((below + 1).min(values.len()), settled)
                }
                None => {
                    let needed = want.min(n - locked.len());
                    let settled = values.len() >= needed
                        && values[..needed].iter().all(|v| v.1 <= MODE_TOLERANCE);
                    (needed, settled)
                }
            };
            if settled || run.exhausted() || run.steps() >= max_steps {
                if !settled && !run.exhausted() {
                    return Err(AnalysisError::SolverNotConverged {
                        history: values.iter().map(|v| v.1).collect(),
                    });
                }
                let (values, vectors) = run.ritz(needed, true);
                break values.into_iter().zip(vectors).collect::<Vec<_>>();
            }
            target = (run.steps() + (run.steps() / 4).max(every)).min(max_steps);
        };
        spent.runs += 1;
        spent.steps += run.steps();
        spent.iterations += run.iterations;
        let mut added = 0;
        for ((lambda, residual), x) in found {
            let new = match cut {
                Some(c) => lambda < c - NEW_BELOW * c.abs().max(f64::MIN_POSITIVE),
                None => true,
            };
            if new && residual <= MODE_TOLERANCE {
                let mut mx = vec![0.0; n];
                m.mul_vec(&x, &mut mx)
                    .map_err(|e| super::internal("modal", e))?;
                locked.push((lambda, residual, x, mx));
                added += 1;
            }
        }
        sort(&mut locked);
        locked.truncate(count);
        if added == 0 || restarts == 0 {
            return Ok((locked, spent));
        }
    }
    Err(AnalysisError::SolverNotConverged {
        history: locked.iter().map(|l| l.1).collect(),
    })
}

/// A lower bound on the eigenvalue a Ritz value λ (residual `r`, about
/// `sigma`) approximates: θ = 1/(λ − σ) has an eigenvalue of A within
/// rθ, so, taking it for the largest as convergence does, the pencil's
/// lowest eigenvalue in the run's space is at least σ + 1/(θ(1 + r)).
fn lower_bound(lambda: f64, r: f64, sigma: f64) -> f64 {
    let theta = 1.0 / (lambda - sigma);
    sigma + 1.0 / (theta * (1.0 + r))
}

/// Ascending λ, ties in their order.
fn sort(locked: &mut [Locked]) {
    locked.sort_by(|a, b| a.0.total_cmp(&b.0));
}
