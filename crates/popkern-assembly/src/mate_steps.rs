//! The solve steps of one cluster (split from `mate_solve`): the
//! Gauss–Newton iteration and the diagnostics at its solution.

use crate::error::AssemblyError;
use crate::mate_linalg::{Svd, rank_of};
use crate::mate_solve::{MATE_TOLERANCE, Problem, RANK_REL};
use crate::mates::Primitive;

/// Gauss–Newton with the truncated-SVD step and a backtracking line
/// search, from the problem's base placements, retried by
/// Levenberg–Marquardt when it stalls short of the tolerance: the
/// increment and the iterations taken.
pub(crate) fn gauss_newton(problem: &Problem<'_>) -> Result<(Vec<f64>, usize), AssemblyError> {
    let n = 6 * problem.free.len();
    let mut x = vec![0.0; n];
    let mut r = problem.residual(&x)?;
    let r0 = r.clone();
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
    if max_abs(&r) > MATE_TOLERANCE && n > 0 {
        let (lx, lr, li) = levenberg_marquardt(problem, vec![0.0; n], r0)?;
        iterations += li;
        if max_abs(&lr) < max_abs(&r) {
            x = lx;
            r = lr;
        }
    }
    // Both full-coordinate solves missed: a revolute chain's joint-space
    // solve (PK-763), verified on the cluster's own residual.
    if max_abs(&r) > MATE_TOLERANCE
        && n > 0
        && let Some(rx) = crate::mate_reduced::rescue(
            problem,
            max_abs(&problem.residual(&vec![0.0; n])?) > crate::mate_reduced::BRANCH_JUMP,
        )?
    {
        let rr = problem.residual(&rx)?;
        if max_abs(&rr) < max_abs(&r) {
            x = rx;
        }
    }
    Ok((x, iterations))
}

/// Levenberg–Marquardt from the start, after Gauss–Newton stalled at a
/// stationary point of the residual that is not a solution (the
/// minimum-norm step stops descending). Damped steps take a different
/// path from the same start: on the arm's target-posing skeleton
/// (six-axis-arm §11 E4) single-jump convergence went 85 → 90 of 100
/// reachable targets; continuing from the stall point instead
/// recovered one. The caller keeps whichever end has the smaller
/// residual. Gives up when no damping finds descent.
fn levenberg_marquardt(
    problem: &Problem<'_>,
    mut x: Vec<f64>,
    mut r: Vec<f64>,
) -> Result<(Vec<f64>, Vec<f64>, usize), AssemblyError> {
    let n = x.len();
    let norm = |v: &[f64]| v.iter().map(|e| e * e).sum::<f64>();
    let mut lambda = None;
    let mut iterations = 0;
    while iterations < 400 && max_abs(&r) > MATE_TOLERANCE * 1e-2 {
        iterations += 1;
        let jacobian = problem.jacobian(&x)?;
        let svd = Svd::new(&jacobian, n);
        let top = svd.max_sigma();
        if top == 0.0 {
            break;
        }
        let mut damping = lambda.unwrap_or(1e-3 * top);
        let minus_r: Vec<f64> = r.iter().map(|v| -v).collect();
        let mut accepted = false;
        while damping < 1e8 * top {
            let step = svd.damped_solve(&minus_r, damping);
            let trial: Vec<f64> = x.iter().zip(&step).map(|(a, s)| a + s).collect();
            let tr = problem.residual(&trial)?;
            if norm(&tr) < norm(&r) {
                x = trial;
                r = tr;
                accepted = true;
                break;
            }
            damping *= 4.0;
        }
        if !accepted {
            break;
        }
        lambda = Some(damping / 3.0);
    }
    Ok((x, r, iterations))
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

/// The largest magnitude among residual rows.
pub(crate) fn max_abs(r: &[f64]) -> f64 {
    r.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

/// The largest coordinate magnitude of a primitive's point: the
/// assembly's size for scaling linear residuals.
pub(crate) fn coords_max(p: &Primitive) -> f64 {
    let o = match *p {
        Primitive::Plane(o, _) | Primitive::Axis(o, _) | Primitive::Point(o) => o,
    };
    o.x.abs().max(o.y.abs()).max(o.z.abs())
}
