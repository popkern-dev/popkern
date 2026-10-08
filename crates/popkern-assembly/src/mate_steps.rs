//! The solve steps of one cluster (split from `mate_solve`): the
//! Gauss–Newton iteration and the diagnostics at its solution.

use crate::error::AssemblyError;
use crate::mate_linalg::{Svd, rank_of};
use crate::mate_solve::{MATE_TOLERANCE, Problem, RANK_REL, max_abs};

/// Gauss–Newton with the truncated-SVD step and a backtracking line
/// search, from the problem's base placements: the increment and the
/// iterations taken.
pub(crate) fn gauss_newton(problem: &Problem<'_>) -> Result<(Vec<f64>, usize), AssemblyError> {
    let n = 6 * problem.free.len();
    let mut x = vec![0.0; n];
    let mut r = problem.residual(&x)?;
    let mut iterations = 0;
    while iterations < 80 && max_abs(&r) > MATE_TOLERANCE * 1e-2 && n > 0 {
        iterations += 1;
        let jacobian = problem.jacobian(&x)?;
        let svd = Svd::new(&jacobian, n);
        let minus_r: Vec<f64> = r.iter().map(|v| -v).collect();
        let step = svd.solve(&minus_r, RANK_REL);
        let norm = |v: &[f64]| v.iter().map(|e| e * e).sum::<f64>();
        let mut alpha = 1.0;
        let mut accepted = None;
        while alpha > 1e-8 {
            let trial: Vec<f64> = x.iter().zip(&step).map(|(a, s)| a + alpha * s).collect();
            let tr = problem.residual(&trial)?;
            if norm(&tr) < norm(&r) {
                accepted = Some((trial, tr));
                break;
            }
            alpha *= 0.5;
        }
        let Some((trial, tr)) = accepted else { break };
        x = trial;
        r = tr;
    }
    Ok((x, iterations))
}

/// One cluster's diagnostics at its solution (indices local to the
/// problem's mates and free roots).
pub(crate) struct Diagnosis {
    pub(crate) residual: f64,
    pub(crate) free_dof: Vec<usize>,
    pub(crate) redundant: Vec<usize>,
    pub(crate) conflicting: Vec<usize>,
}

pub(crate) fn diagnose(problem: &Problem<'_>, x: &[f64]) -> Result<Diagnosis, AssemblyError> {
    let n = 6 * problem.free.len();
    let r = problem.residual(x)?;
    let rows = problem.rows(x)?;
    let jacobian = if n > 0 {
        problem.jacobian(x)?
    } else {
        Vec::new()
    };
    let mut free_dof = Vec::new();
    if n > 0 {
        let svd = Svd::new(&jacobian, n);
        let null = svd.null_space(RANK_REL);
        for k in 0..problem.free.len() {
            let block: Vec<Vec<f64>> = null.iter().map(|v| v[6 * k..6 * k + 6].to_vec()).collect();
            free_dof.push(rank_of(&block, 6, 1e-6));
        }
    }
    let mut redundant = Vec::new();
    let mut conflicting = Vec::new();
    let mut so_far: Vec<Vec<f64>> = Vec::new();
    for (index, range) in rows.iter().enumerate() {
        let mate_rows: Vec<Vec<f64>> = jacobian
            .get(range.clone())
            .map(<[Vec<f64>]>::to_vec)
            .unwrap_or_default();
        let own = rank_of(&mate_rows, n, RANK_REL);
        let before = rank_of(&so_far, n, RANK_REL);
        so_far.extend(mate_rows);
        let after = rank_of(&so_far, n, RANK_REL);
        if after - before < own || own == 0 {
            redundant.push(index);
        }
        if r[range.clone()].iter().any(|v| v.abs() > MATE_TOLERANCE) {
            conflicting.push(index);
        }
    }
    Ok(Diagnosis {
        residual: max_abs(&r),
        free_dof,
        redundant,
        conflicting,
    })
}
