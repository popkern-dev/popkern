//! A thermal strain ε_th = α θ I in the vector cells (map Stage 2, item
//! 9; Zander 2012).
//!
//! The cell energy ½∫(ε − ε_th):D:(ε − ε_th) gives loads
//! f = ∫ Bᵀ D ε_th on the nodal and the enhanced dofs alike; D ε_th =
//! (3λ + 2μ) α θ I for an isotropic material. Condensing the enhanced
//! dofs as the stiffness does (α_enh = R u + K_αα⁻¹ f_α, R = −K_αα⁻¹K_αu)
//! leaves the nodal load f_u + Rᵀ f_α and an affine part α₀ = K_αα⁻¹ f_α
//! of the enhanced dofs, which the stress and the Nitsche flux carry.
//! θ is trilinear in the cell from its 8 nodal values, so every integral
//! is exact in the fragment's 27 moments.

use super::poly::{Poly27, integrate, mul, shape, shape_derivative};
use super::vector::{ENHANCED_DOFS, NODAL_DOFS, VectorCell, fragment_mean, xi};

/// θ(ξ) = Σ_b N_b(ξ) θ_b.
pub fn theta_poly(theta: &[f64; 8]) -> Poly27 {
    let mut p = [0.0; 27];
    for (b, t) in theta.iter().enumerate() {
        for (x, y) in p.iter_mut().zip(shape(b)) {
            *x += t * y;
        }
    }
    p
}

/// One cell's thermal load.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellThermal {
    /// The condensed nodal load, N.
    pub load: [f64; NODAL_DOFS],
    /// The enhanced dofs' affine part α₀ (zero without enhanced modes).
    pub alpha0: [f64; ENHANCED_DOFS],
}

/// The thermal load of expansion `alpha` (1/K) and nodal temperature
/// rises `theta` (K) on a cell with moments `m` (physical volume), size
/// `h`, Lamé (λ, μ) and condensed operator `op`.
pub fn cell_thermal(
    op: &VectorCell,
    m: &[f64; 27],
    h: f64,
    (lambda, mu): (f64, f64),
    alpha: f64,
    theta: &[f64; 8],
) -> CellThermal {
    cell_thermal_parts(op, m, &[(m, lambda, mu)], h, alpha, theta)
}

/// [`cell_thermal`] for a cell of several materials sharing `alpha`
/// (walls and an infill core of the same polymer): `parts` gives each
/// material's moments (summing to `m`) and Lamé parameters, as
/// `VectorTables::cell_parts` takes them.
pub fn cell_thermal_parts(
    op: &VectorCell,
    m: &[f64; 27],
    parts: &[(&[f64; 27], f64, f64)],
    h: f64,
    alpha: f64,
    theta: &[f64; 8],
) -> CellThermal {
    let g = 2.0 / h;
    let tp = theta_poly(theta);
    let mean = fragment_mean(m);
    let mut load = [0.0; NODAL_DOFS];
    let mut f_alpha = [0.0; ENHANCED_DOFS];
    for &(pm, lambda, mu) in parts {
        let c = (3.0 * lambda + 2.0 * mu) * alpha;
        for a in 0..8 {
            for i in 0..3 {
                load[3 * a + i] += c * g * integrate(&mul(&shape_derivative(a, i), &tp), pm);
            }
        }
        // Bubble m's gradient, physical: g · (−2)(ξ_m − ξ̄_m), along m only.
        for bubble in 0..3 {
            let moment = integrate(&mul(&xi(bubble), &tp), pm) - mean[bubble] * integrate(&tp, pm);
            f_alpha[3 * bubble + bubble] += c * g * -2.0 * moment;
        }
    }
    let mut alpha0 = [0.0; ENHANCED_DOFS];
    if let (Some(recovery), Some(kaa_inv)) = (&op.recovery, &op.kaa_inv) {
        for (r, row) in kaa_inv.iter().enumerate() {
            alpha0[r] = row.iter().zip(&f_alpha).map(|(k, f)| k * f).sum();
        }
        for (col, slot) in load.iter_mut().enumerate() {
            *slot += (0..ENHANCED_DOFS)
                .map(|r| recovery[r][col] * f_alpha[r])
                .sum::<f64>();
        }
    }
    CellThermal { load, alpha0 }
}

#[cfg(test)]
#[path = "thermal_strain_tests.rs"]
mod tests;
