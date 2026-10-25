//! Nitsche's β on long aggregates.
//!
//! β = factor·k/h stabilises Nitsche when each aggregate is a few cells
//! across: a small cell's field is its root's trilinear extended a cell
//! or two. Material thinner than the lattice spacing along its whole
//! length (a fin, PH-067) has no large cell near, and its cells extend a
//! root up to tens of cells away. There the root's trilinear's gradient
//! on the fixed faces is no longer bounded by its energy in the
//! aggregate's volume, and the system is indefinite.
//!
//! On a root cell and the members whose nodes are all extended from it
//! alone (`Aggregation::root_weights`: the field there is the root's
//! trilinear), β is raised to
//! 4·k·λ, λ the largest ratio ∫_Γ|∇v|²dA / ∫_Ω|∇v|²dV over that group
//! (Γ its fixed pieces): since |∇v·n| ≤ |∇v|, that bounds the Nitsche
//! consistency term by half the energy plus 2kλ∫_Γv², so the group is
//! coercive. Nitsche is consistent for any β; linear fields stay exact.

use crate::aggregate::cluster::{clusters, gather};
use crate::element::{ScalarTables, surface_moments};
use crate::linalg::dense_eigen::max_generalized_eigenvalue_dense;

use super::assemble::{MmCondition, Physics, SolidInput, piece_condition};

/// Power iterations for a cluster's largest ratio.
const ITERATIONS: usize = 300;

/// Nitsche's β for each fragment entry of the solid: `physics.beta`,
/// raised on each cluster of extended cells to 4·k·λ (above).
pub(super) fn entry_betas(t: &ScalarTables, input: &SolidInput<'_>, physics: &Physics) -> Vec<f64> {
    let grid = input.grid;
    let h = grid.cell_size();
    let moments = input.moments;
    let pieces = &input.clipped.pieces;
    let mut betas = vec![physics.beta; moments.cells.len()];
    let flat = |m: [[f64; 8]; 8]| m.iter().flatten().copied().collect::<Vec<f64>>();
    for cluster in clusters(input.aggregation, grid, moments) {
        let mut volume = gather(&cluster, 1, &mut |entry| {
            moments
                .cells
                .get(entry)
                .map(|c| flat(t.conduction(&c.m, 1.0, h)))
        });
        let surface = gather(&cluster, 1, &mut |entry| {
            let cell = moments.cells.get(entry)?;
            let lo = grid.cell_min(grid.cell_ijk(cell.cell).unwrap_or_default());
            let start = pieces.partition_point(|p| p.cell < cell.cell);
            let end = pieces.partition_point(|p| p.cell <= cell.cell);
            let mut ks = [[0.0; 8]; 8];
            let mut fixed = false;
            for (index, piece) in pieces[start..end].iter().enumerate() {
                if moments.entry_of_piece(start + index, cell.cell) != Some(entry) {
                    continue;
                }
                let condition = input
                    .conditions
                    .get(piece.face.0 as usize)
                    .copied()
                    .unwrap_or(MmCondition::Insulated);
                if !matches!(
                    piece_condition(condition, input, piece),
                    MmCondition::Fixed { .. }
                ) {
                    continue;
                }
                fixed = true;
                let s = surface_moments(std::iter::once(piece), lo, h);
                let k = t.conduction(&s.area, 1.0, h);
                for a in 0..8 {
                    for b in 0..8 {
                        ks[a][b] += k[a][b];
                    }
                }
            }
            fixed.then(|| flat(ks))
        });
        if surface.iter().all(|&x| x == 0.0) {
            continue;
        }
        // Constants have no gradient: lift them off the null space so the
        // volume form factors without changing the largest ratio.
        let n = cluster.dofs.len();
        let lift = (0..n).map(|i| volume[i * n + i]).sum::<f64>() / n as f64;
        volume.iter_mut().for_each(|x| *x += lift);
        let Some(lambda) = max_generalized_eigenvalue_dense(&surface, &volume, n, ITERATIONS)
        else {
            continue;
        };
        let beta = physics.beta.max(4.0 * physics.k * lambda);
        for (entry, _) in &cluster.cells {
            if let Some(slot) = betas.get_mut(*entry) {
                *slot = beta;
            }
        }
    }
    betas
}
