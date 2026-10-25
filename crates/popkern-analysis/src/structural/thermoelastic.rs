//! Thermoelastic loads (map Stage 2, item 9): a temperature rise θ over
//! a reference expands the material by α θ in every direction. A staggered
//! solve: the temperature comes first (a thermal solve's field, or
//! uniform), then the elastic solve carries ε_th as loads on every cell
//! (`element::thermal_strain`) and as the known part of the Nitsche flux
//! on supports (Zander 2012 eqn 6b), and recovered stresses subtract it.

use crate::element::surface_moments;
use crate::element::thermal_strain::{CellThermal, cell_thermal_parts};
use crate::element::vector::{NODAL_DOFS, VectorCell, fragment_mean};
use crate::element::vector_boundary::{BoundaryTables, PieceSupport};
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::moments::SolidMoments;

use super::assemble::{Elastic, FaceMechanics, SolidInput, piece_support};

/// A thermal load: expansion α (1/K) of the temperature over a
/// reference, °C. `temperature` gives °C at a point (`None` outside).
pub struct ThermalLoad<'a> {
    /// Linear expansion coefficient, 1/K.
    pub alpha: f64,
    /// Stress-free temperature, °C.
    pub reference_c: f64,
    /// The temperature field.
    pub temperature: &'a dyn Fn([f64; 3]) -> Option<f64>,
}

/// Per entry of one solid: the temperature rise over the reference at
/// its 8 nodes (`CELL_CORNERS` order), K, with the expansion α, 1/K.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalStrains {
    /// Linear expansion coefficient, 1/K.
    pub alpha: f64,
    /// Nodal rises per entry.
    pub theta: Vec<[f64; 8]>,
}

impl ThermalStrains {
    /// The rises of `moments`' entries from `temperature` (°C at a
    /// point, `None` outside the part) over `reference_c`. A node is
    /// read a hair toward its entry's material centroid, so a split
    /// cell's fragments read their own side; a node no reading reaches
    /// takes the entry's centroid value.
    pub fn sample(
        grid: &Grid,
        moments: &SolidMoments,
        alpha: f64,
        reference_c: f64,
        temperature: &dyn Fn([f64; 3]) -> Option<f64>,
    ) -> Result<ThermalStrains, AnalysisError> {
        let h = grid.cell_size();
        let mut theta = Vec::with_capacity(moments.cells.len());
        for c in &moments.cells {
            let ijk = grid
                .cell_ijk(c.cell)
                .ok_or_else(|| AnalysisError::Internal {
                    stage: "thermal strain",
                    detail: "cell outside grid".to_string(),
                })?;
            let lo = grid.cell_min(ijk);
            let local = crate::recover::cell::centroid(&c.m);
            let centre = [0, 1, 2].map(|a| lo[a] + 0.5 * h * (local[a] + 1.0));
            let at_centre = temperature(centre);
            let mut t = [0.0; 8];
            for (slot, corner) in t.iter_mut().zip(crate::grid::CELL_CORNERS) {
                let x = [0, 1, 2].map(|a| lo[a] + h * corner[a] as f64);
                let p = [0, 1, 2].map(|a| x[a] + 1.0e-6 * (centre[a] - x[a]));
                let value = temperature(p).or(at_centre).unwrap_or(reference_c);
                *slot = value - reference_c;
            }
            theta.push(t);
        }
        Ok(ThermalStrains { alpha, theta })
    }
}

/// Entry `entry`'s cell thermal terms, when the input carries a strain.
pub(super) fn cell_terms(
    input: &SolidInput<'_>,
    entry: usize,
    op: &VectorCell,
    elastic: &Elastic,
) -> Option<(CellThermal, [f64; 8], f64)> {
    let strains = input.thermal?;
    let theta = *strains.theta.get(entry)?;
    let cell = input.moments.cells.get(entry)?;
    let h = input.grid.cell_size();
    // Walls and an infill core expand alike (one polymer), each with its
    // own stiffness.
    let parts = super::assemble::core_of(input, entry).parts(
        &cell.m,
        (elastic.lambda, elastic.mu),
        elastic.core,
    );
    let parts: Vec<(&[f64; 27], f64, f64)> = parts.iter().map(|(m, l, u)| (m, *l, *u)).collect();
    Some((
        cell_thermal_parts(op, &cell.m, &parts, h, strains.alpha, &theta),
        theta,
        strains.alpha,
    ))
}

/// The known thermal flux of one support piece of entry `entry` (moved
/// to the right-hand side), given the entry's cell terms.
#[allow(clippy::too_many_arguments)]
pub(super) fn piece_flux(
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    entry: usize,
    terms: &(CellThermal, [f64; 8], f64),
    elastic: &Elastic,
    support: PieceSupport,
    piece: &crate::clip::FacetPiece,
    enhanced: bool,
) -> [f64; NODAL_DOFS] {
    let (cell_terms, theta, alpha) = terms;
    let Some(cell) = input.moments.cells.get(entry) else {
        return [0.0; NODAL_DOFS];
    };
    let grid = input.grid;
    let h = grid.cell_size();
    let lo = grid.cell_min(grid.cell_ijk(cell.cell).unwrap_or_default());
    let s = surface_moments(std::iter::once(piece), lo, h);
    let mean = fragment_mean(&cell.m);
    bt.thermal_flux(
        support,
        &s,
        piece.normal,
        h,
        (elastic.lambda, elastic.mu),
        (*alpha, theta),
        enhanced.then_some((mean, &cell_terms.alpha0)),
    )
}

/// The total thermal load of entry `entry`: the cell load plus the
/// known flux of its support pieces.
pub(super) fn entry_load(
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    entry: usize,
    op: &VectorCell,
    elastic: &Elastic,
) -> Option<[f64; NODAL_DOFS]> {
    let terms = cell_terms(input, entry, op, elastic)?;
    let mut load = terms.0.load;
    let cell = input.moments.cells.get(entry)?;
    let pieces = &input.clipped.pieces;
    let start = pieces.partition_point(|p| p.cell < cell.cell);
    let end = pieces.partition_point(|p| p.cell <= cell.cell);
    let enhanced = elastic.enhanced_flux && op.recovery.is_some();
    for (index, piece) in pieces[start..end].iter().enumerate() {
        if input.moments.entry_of_piece(start + index, cell.cell) != Some(entry) {
            continue;
        }
        let condition = input
            .faces
            .get(piece.face.0 as usize)
            .copied()
            .unwrap_or(FaceMechanics::Free);
        let Some(support) = piece_support(condition, input, piece) else {
            continue;
        };
        let f = piece_flux(bt, input, entry, &terms, elastic, support, piece, enhanced);
        for (x, y) in load.iter_mut().zip(f) {
            *x += y;
        }
    }
    Some(load)
}

#[cfg(test)]
#[path = "thermoelastic_tests.rs"]
mod tests;
