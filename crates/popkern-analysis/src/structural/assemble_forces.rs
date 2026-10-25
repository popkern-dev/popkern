//! Face forces of a solved system: the reaction on every support and
//! the applied load per face.

use crate::element::surface_moments;
use crate::element::vector::VectorTables;
use crate::element::vector_boundary::BoundaryTables;
use crate::linalg::sparse::SparseError;

use super::{
    Elastic, FaceMechanics, SolidInput, cell_operators, cell_positions, core_of, enhanced_flux,
    piece_support, piece_traction,
};

/// Force each face exerts on its supports (N, by `FaceId`, summed into
/// `reaction`) and the applied load per face (into `applied`), for
/// nodal displacements `u` (3 per active node).
pub fn face_forces(
    vt: &VectorTables,
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    elastic: &Elastic,
    u: &[f64],
    reaction: &mut [[f64; 3]],
    applied: &mut [[f64; 3]],
) -> Result<(), SparseError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let e = &input.aggregation.extension;
    let pieces = &input.clipped.pieces;
    let operators =
        cell_operators(vt, h, elastic).map_err(|_| SparseError::Shape("cell operator"))?;
    let betas = super::super::stability::entry_betas(&operators, bt, input, elastic)
        .map_err(|_| SparseError::Shape("cell operator"))?;
    for (entry, cell) in input.moments.cells.iter().enumerate() {
        let elastic = &Elastic {
            beta: betas.get(entry).copied().unwrap_or(elastic.beta),
            ..*elastic
        };
        // Only cells with surface pieces carry face forces.
        let start = pieces.partition_point(|p| p.cell < cell.cell);
        let end = pieces.partition_point(|p| p.cell <= cell.cell);
        if start == end {
            continue;
        }
        let condition_of = |index: usize| {
            let face = pieces[index].face.0 as usize;
            let mine = input.moments.entry_of_piece(index, cell.cell) == Some(entry);
            let condition = input
                .faces
                .get(face)
                .copied()
                .unwrap_or(FaceMechanics::Free);
            mine.then_some((face, condition))
        };
        // The cell's field and flux are needed only for support pieces
        // (the reaction); loads need their pieces' areas alone.
        let supported = (start..end).any(|i| {
            matches!(
                condition_of(i),
                Some((_, FaceMechanics::Support(_) | FaceMechanics::Prescribed))
            )
        });
        let lo = grid.cell_min(grid.cell_ijk(cell.cell).unwrap_or_default());
        let operator = if supported {
            let nodes = cell_positions(grid, input.moments, e, entry)?;
            let operator = operators
                .of_core(
                    cell,
                    input.moments.is_fragment(entry),
                    core_of(input, entry),
                )
                .map_err(|_| SparseError::Shape("cell operator"))?;
            let ue: [f64; 24] = std::array::from_fn(|r| u[3 * nodes[r / 3] + r % 3]);
            Some((operator, ue))
        } else {
            None
        };
        let flux = operator
            .as_ref()
            .map(|(operator, ue)| (enhanced_flux(&operator.get().0, cell, elastic), *ue));
        let thermal = operator.as_ref().and_then(|(operator, _)| {
            let op = &operator.get().0;
            super::super::thermoelastic::cell_terms(input, entry, op, elastic)
                .map(|t| (t, elastic.enhanced_flux && op.recovery.is_some()))
        });
        for (index, piece) in pieces[start..end].iter().enumerate() {
            let Some((face, condition)) = condition_of(start + index) else {
                continue;
            };
            match condition {
                FaceMechanics::Support(_) | FaceMechanics::Prescribed => {
                    let Some((enhanced, ue)) = flux else {
                        continue;
                    };
                    let Some(support) = piece_support(condition, input, piece) else {
                        continue;
                    };
                    let s = surface_moments(std::iter::once(piece), lo, h);
                    let (k, f) = bt.support(
                        support,
                        &s,
                        piece.normal,
                        lo,
                        h,
                        elastic.lambda,
                        elastic.mu,
                        elastic.beta,
                        enhanced,
                    );
                    let fth = thermal.as_ref().map_or([0.0; 24], |(t, enhanced)| {
                        super::super::thermoelastic::piece_flux(
                            bt, input, entry, t, elastic, support, piece, *enhanced,
                        )
                    });
                    // Σ_a (K u − f)_a per component: force on the support.
                    for r in 0..24 {
                        let v: f64 = (0..24).map(|c| k[r][c] * ue[c]).sum::<f64>() - f[r] - fth[r];
                        if let Some(slot) = reaction.get_mut(face) {
                            slot[r % 3] += v;
                        }
                    }
                }
                FaceMechanics::Traction(_) | FaceMechanics::Pressure(_) => {
                    let t = match condition {
                        FaceMechanics::Traction(t) => t,
                        FaceMechanics::Pressure(p) => piece.normal.map(|n| -p * n),
                        _ => [0.0; 3],
                    };
                    if let Some(slot) = applied.get_mut(face) {
                        for i in 0..3 {
                            slot[i] += t[i] * piece.area;
                        }
                    }
                }
                FaceMechanics::PrescribedTraction => {
                    // The cap's traction is an applied load (its weak
                    // spring carries next to nothing).
                    if let (Some((t, _, _)), Some(slot)) =
                        (piece_traction(input, elastic, piece), applied.get_mut(face))
                    {
                        for i in 0..3 {
                            slot[i] += t[i] * piece.area;
                        }
                    }
                }
                FaceMechanics::Free => {}
            }
        }
    }
    Ok(())
}
