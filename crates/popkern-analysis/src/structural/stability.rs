//! Nitsche's β on long aggregates, structural (PH-068; the thermal
//! counterpart is `thermal::stability`).
//!
//! β = γ_D·(λ + 2μ)/h stabilises the supports when each aggregate is a
//! few cells across. Material thinner than the lattice spacing along
//! its whole length extends roots up to tens of cells away, and its
//! traction σ(v)n on the supports is no longer bounded by its strain
//! energy: the system goes indefinite.
//!
//! On each cluster of extended cells (`aggregate::cluster`), β is raised
//! to 4λ, λ the largest ratio of ∫_Γ|P σ(v) n|² (the flux each support
//! imposes, with the cells' enhanced modes) to vᵀKv (the condensed cell
//! stiffness) over the cluster. Then 2∫(Pσn)·(Pv) ≤ ½vᵀKv + 2λ∫|Pv|²,
//! so the cluster is coercive; Nitsche stays consistent for any β.

use crate::aggregate::cluster::{clusters, gather};
use crate::element::operators::CellOperators;
use crate::element::surface_moments;
use crate::element::vector_boundary::BoundaryTables;
use crate::linalg::dense_eigen::max_generalized_eigenvalue_dense;

use super::assemble::{
    AssemblyError, Elastic, FaceMechanics, SolidInput, core_of, enhanced_flux, piece_support,
};

/// Power iterations for a cluster's largest ratio.
const ITERATIONS: usize = 300;

/// Nitsche's β for each fragment entry of the solid: `elastic.beta`,
/// raised on each cluster of extended cells to 4λ (above).
pub(super) fn entry_betas(
    operators: &CellOperators<'_>,
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    elastic: &Elastic,
) -> Result<Vec<f64>, AssemblyError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let moments = input.moments;
    let pieces = &input.clipped.pieces;
    let mut betas = vec![elastic.beta; moments.cells.len()];
    let flat = |m: &[[f64; 24]; 24]| m.iter().flatten().copied().collect::<Vec<f64>>();
    for cluster in clusters(input.aggregation, grid, moments) {
        let mut failed = None;
        let mut volume = gather(&cluster, 3, &mut |entry| {
            let cell = moments.cells.get(entry)?;
            match operators.of_core(cell, moments.is_fragment(entry), core_of(input, entry)) {
                Ok(op) => Some(flat(&op.get().0.k)),
                Err(e) => {
                    failed = Some(AssemblyError::Enhanced {
                        cell: cell.cell,
                        detail: e.to_string(),
                    });
                    None
                }
            }
        });
        if let Some(e) = failed {
            return Err(e);
        }
        let surface = gather(&cluster, 3, &mut |entry| {
            let cell = moments.cells.get(entry)?;
            let operator = operators
                .of_core(cell, moments.is_fragment(entry), core_of(input, entry))
                .ok()?;
            let op = &operator.get().0;
            let enhanced = enhanced_flux(op, cell, elastic);
            let lo = grid.cell_min(grid.cell_ijk(cell.cell).unwrap_or_default());
            let start = pieces.partition_point(|p| p.cell < cell.cell);
            let end = pieces.partition_point(|p| p.cell <= cell.cell);
            let mut ks = [[0.0; 24]; 24];
            let mut supported = false;
            for (index, piece) in pieces[start..end].iter().enumerate() {
                if moments.entry_of_piece(start + index, cell.cell) != Some(entry) {
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
                supported = true;
                let s = surface_moments(std::iter::once(piece), lo, h);
                let k = bt.flux_gram(
                    support,
                    &s,
                    piece.normal,
                    h,
                    elastic.lambda,
                    elastic.mu,
                    enhanced,
                );
                for r in 0..24 {
                    for c in 0..24 {
                        ks[r][c] += k[r][c];
                    }
                }
            }
            supported.then(|| flat(&ks))
        });
        if surface.iter().all(|&x| x == 0.0) {
            continue;
        }
        // Rigid motions have no strain: lift them off the null space so
        // the stiffness factors without changing the largest ratio.
        let n = 3 * cluster.dofs.len();
        let positions = cluster.positions(grid, input.aggregation);
        let count = positions.len().max(1) as f64;
        let centre: [f64; 3] =
            std::array::from_fn(|a| positions.iter().map(|p| p[a]).sum::<f64>() / count);
        let lift = (0..n).map(|i| volume[i * n + i]).sum::<f64>() / n as f64;
        for mode in 0..6 {
            let v: Vec<f64> = positions
                .iter()
                .flat_map(|p| {
                    let x = [0, 1, 2].map(|a| p[a] - centre[a]);
                    if mode < 3 {
                        std::array::from_fn(|i| if i == mode { 1.0 } else { 0.0 })
                    } else {
                        // e_k × x.
                        let k = mode - 3;
                        let e: [f64; 3] = std::array::from_fn(|i| if i == k { 1.0 } else { 0.0 });
                        [
                            e[1] * x[2] - e[2] * x[1],
                            e[2] * x[0] - e[0] * x[2],
                            e[0] * x[1] - e[1] * x[0],
                        ]
                    }
                })
                .collect();
            let norm2: f64 = v.iter().map(|x| x * x).sum();
            if norm2 == 0.0 {
                continue;
            }
            for i in 0..n {
                for j in 0..n {
                    volume[i * n + j] += lift * v[i] * v[j] / norm2;
                }
            }
        }
        let Some(lambda) = max_generalized_eigenvalue_dense(&surface, &volume, n, ITERATIONS)
        else {
            continue;
        };
        let beta = elastic.beta.max(4.0 * lambda);
        for (entry, _) in &cluster.cells {
            if let Some(slot) = betas.get_mut(*entry) {
                *slot = beta;
            }
        }
    }
    Ok(betas)
}
