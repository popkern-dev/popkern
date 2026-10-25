//! The true error of σ_h against an exact field (verification only).

use crate::error::AnalysisError;

use super::cell::{StressPoly, stress_at};
use super::tensor::{self, Voigt};
use super::{ExactStress, RecoveryInput};

/// Relative true error of σ_h against `exact`, by the midpoint rule on
/// the 8³ sub-cell lattice points inside the solid.
pub(super) fn true_error(
    input: &RecoveryInput<'_>,
    polys: &[StressPoly],
    exact: ExactStress<'_>,
) -> Result<Option<f64>, AnalysisError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let (e_mod, nu) = input.material;
    let n = 8usize;
    let w = (h / n as f64) * (h / n as f64) * (h / n as f64);
    let (mut err, mut norm) = (0.0, 0.0);
    for (k, c) in input.moments.cells.iter().enumerate() {
        let ijk = grid
            .cell_ijk(c.cell)
            .ok_or_else(|| super::internal("cell outside grid".into()))?;
        let lo = grid.cell_min(ijk);
        for kk in 0..n {
            for jj in 0..n {
                for ii in 0..n {
                    let xi = [ii, jj, kk].map(|q| -1.0 + (2.0 * q as f64 + 1.0) / n as f64);
                    let x = [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xi[a] + 1.0));
                    if !crate::moments::lattice::contains(grid, input.clipped, ijk, x) {
                        continue;
                    }
                    let ex = exact(x);
                    let sh = stress_at(&polys[k], xi);
                    let d: Voigt = std::array::from_fn(|q| ex[q] - sh[q]);
                    err += w * tensor::energy_density(&d, e_mod, nu);
                    norm += w * tensor::energy_density(&ex, e_mod, nu);
                }
            }
        }
    }
    Ok((norm > 0.0).then(|| (err.max(0.0) / norm).sqrt()))
}
