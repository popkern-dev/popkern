//! The block-diagonal system's multigrid solve.

use crate::error::AnalysisError;
use crate::linalg::sparse::Csr;
use crate::options::AnalysisOptions;
use crate::report::SolveReport;
use crate::solve::pcg::pcg;
use crate::solve::{Hierarchy, HierarchyOptions, MultigridError};

use super::internal;

/// Multigrid PCG on the block-diagonal system of `blocks`: the
/// hierarchy, its report fields, and a decline when the solve does not
/// converge.
pub(super) fn solve_system(
    blocks: Vec<Csr>,
    positions: Vec<[usize; 3]>,
    labels: Vec<u8>,
    b: &[f64],
    options: &AnalysisOptions,
    report: &mut SolveReport,
) -> Result<crate::solve::pcg::PcgResult, AnalysisError> {
    solve_system_reusing(blocks, (positions, labels), b, options, report, None).map(|s| s.0)
}

/// [`solve_system`] reusing `previous`'s prolongations where they hold
/// (`Hierarchy::rebuild`), returning the hierarchy for the next solve.
pub(super) fn solve_system_reusing(
    blocks: Vec<Csr>,
    (positions, labels): (Vec<[usize; 3]>, Vec<u8>),
    b: &[f64],
    options: &AnalysisOptions,
    report: &mut SolveReport,
    previous: Option<&Hierarchy>,
) -> Result<(crate::solve::pcg::PcgResult, Hierarchy), AnalysisError> {
    report.dof = b.len();
    let a = Csr::block_diagonal(blocks).map_err(|e| internal("assemble", e))?;
    let hierarchy = Hierarchy::rebuild(
        a,
        positions,
        labels,
        &HierarchyOptions {
            coarsest_dofs: options.coarsest_dofs,
            split_components: true,
        },
        previous,
    )
    .map_err(|e| match e {
        MultigridError::CoarsestTooLarge { dofs, cap } => {
            AnalysisError::CoarseLevelTooLarge { dofs, cap }
        }
        MultigridError::Coarsest(c) => report.singular(c.to_string()),
        other => internal("multigrid", other),
    })?;
    report.levels = hierarchy.levels.len();
    report.level_dofs = hierarchy.levels.iter().map(|l| l.a.rows()).collect();
    let solution = pcg(&hierarchy, b, options.cg_rel_tol, options.cg_max_iters)
        .map_err(|e| internal("pcg", e))?;
    report.cg_iterations = solution.iterations;
    report.residual_history = solution.history.clone();
    if !solution.converged {
        return Err(AnalysisError::SolverNotConverged {
            history: solution.history,
        });
    }
    Ok((solution, hierarchy))
}
