//! The unscaled element blocks the sensitivities read, built once per
//! design: u_eᵀ k_e u_e with k_e the entry's material block plus its
//! Nitsche support terms (both scale with the density).

use crate::error::AnalysisError;

use super::super::assemble::{FaceMechanics, SolidInput, cell_operators, unscaled_terms};
use super::super::solve::Prepared;
use super::super::stability::entry_betas;
use super::super::{ElasticProblem, internal};

/// A 24 × 24 element block.
pub(super) type Block = [[f64; 24]; 24];

/// Each entry's unscaled block for the sensitivities, built once: full
/// cells without a support share one; cut cells, fragments and
/// supported entries keep their own.
pub(super) struct Blocks {
    pub(super) full: Box<Block>,
    pub(super) own: Vec<Vec<Option<Box<Block>>>>,
}

/// The unscaled blocks of `prepared` (see [`Blocks`]).
pub(super) fn sensitivity_blocks(
    prepared: &Prepared,
    problem: &ElasticProblem,
) -> Result<Blocks, AnalysisError> {
    let grid = &prepared.grid;
    let operators = cell_operators(&prepared.vt, grid.cell_size(), &prepared.elastic)
        .map_err(|e| internal("generative", e))?;
    let full = operators
        .of(
            &crate::moments::CellMoments {
                cell: 0,
                kind: crate::moments::CellKind::Full,
                m: crate::element::operators::full_cell_moments(grid.cell_size()),
                fragment: 0,
            },
            false,
        )
        .map_err(|e| internal("generative", e))?;
    let full = Box::new(full.get().0.k);
    let mut own = Vec::with_capacity(prepared.solids.len());
    for (s, (_, clipped, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let input = SolidInput {
            grid,
            clipped,
            moments,
            aggregation,
            faces: &problem.faces,
            prescribed: None,
            thermal: None,
            cores: prepared.cores[s].as_deref(),
            scales: None,
        };
        let betas = entry_betas(&operators, &prepared.bt, &input, &prepared.elastic)
            .map_err(|e| internal("generative", e))?;
        let mut supported = vec![false; moments.cells.len()];
        for (index, piece) in clipped.pieces.iter().enumerate() {
            let held = matches!(
                problem.faces.get(piece.face.0 as usize),
                Some(FaceMechanics::Support(_))
            );
            if held
                && let Some(entry) = moments.entry_of_piece(index, piece.cell)
                && let Some(slot) = supported.get_mut(entry)
            {
                *slot = true;
            }
        }
        let block = |entry: usize| -> Result<Option<Box<Block>>, AnalysisError> {
            let cell = &moments.cells[entry];
            let plain = cell.kind == crate::moments::CellKind::Full
                && !moments.is_fragment(entry)
                && !supported[entry];
            if plain {
                return Ok(None);
            }
            unscaled_terms(
                &operators,
                &prepared.bt,
                &input,
                &prepared.elastic,
                &betas,
                entry,
            )
            .map(|(k, _, _)| Some(Box::new(k)))
            .map_err(|e| internal("generative", e))
        };
        let mut blocks = Vec::with_capacity(moments.cells.len());
        for b in crate::parallel::map_range(0..moments.cells.len(), 256, block) {
            blocks.push(b?);
        }
        own.push(blocks);
    }
    Ok(Blocks { full, own })
}
