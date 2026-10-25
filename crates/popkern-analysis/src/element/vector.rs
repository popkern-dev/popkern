//! Vector Q1 cells with incompatible (Wilson) modes, condensed per cell
//! (map §1.2, notes 03).
//!
//! The cell's functions are the 8 trilinear N_a and 3 bubbles
//! 1 − ξ_m², each in 3 components: 24 nodal and 9 enhanced dofs. All
//! integrands have exponents ≤ 2 per axis, so every cell matrix is an
//! exact combination of the fragment's 27 moments — full and cut cells
//! share one path. The bubble gradients are taken relative to their
//! fragment mean (Simo–Rifai eq. 62), which is zero on full cells; the
//! 9 enhanced dofs are then condensed out.

use crate::linalg::{Cholesky, CholeskyError};
use crate::moments::moment_index;

use super::poly::{Poly27, integrate, mul, shape_derivative};

/// Nodal dofs per cell (8 nodes × 3 components), ordered 3a + i.
pub const NODAL_DOFS: usize = 24;
/// Enhanced dofs per cell (3 bubbles × 3 components), ordered 3m + i.
pub const ENHANCED_DOFS: usize = 9;

/// A 24×24 cell matrix.
pub type Matrix24 = [[f64; NODAL_DOFS]; NODAL_DOFS];

/// Lamé parameters from Young's modulus and Poisson's ratio.
pub fn lame(youngs: f64, poisson: f64) -> (f64, f64) {
    let lambda = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let mu = youngs / (2.0 * (1.0 + poisson));
    (lambda, mu)
}

/// Polynomial tables: ∫ ∂_k f ∂_l g for the 8 nodal functions, and the
/// pieces needed for the bubbles.
#[derive(Debug, Clone)]
pub struct VectorTables {
    /// ∂_k N_a ∂_l N_b, indexed [a][b][k][l].
    nn: Vec<[[Poly27; 3]; 3]>,
    /// ∂_k N_a · ξ_m, indexed [a][k][m].
    n_xi: Vec<[[Poly27; 3]; 3]>,
    /// ∂_k N_a, indexed [a][k].
    dn: [[Poly27; 3]; 8],
}

impl Default for VectorTables {
    fn default() -> Self {
        Self::new()
    }
}

pub(crate) fn xi(m: usize) -> Poly27 {
    let mut p = [0.0; 27];
    let mut e = [0; 3];
    e[m] = 1;
    p[moment_index(e[0], e[1], e[2])] = 1.0;
    p
}

/// Volume fraction below which a cut cell skips the enhanced modes
/// (plain Q1 with the fragment's exact moments; map §1.2).
pub const ETA_EAS: f64 = 0.05;

/// How a cell is enhanced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Enhancement {
    /// Plain Q1.
    None,
    /// Wilson modes with gradients relative to their fragment mean
    /// (Simo–Rifai eq. 62): the production element.
    FragmentMean,
    /// Wilson modes without the mean shift. Fails the cut-cell patch
    /// test; exists only to show what the correction prevents.
    Raw,
}

/// What to keep from a cell computation.
#[derive(Debug, Clone, PartialEq)]
pub struct VectorCell {
    /// The condensed 24×24 stiffness.
    pub k: Matrix24,
    /// α = `recovery` · u: the enhanced dofs of nodal displacements u.
    pub recovery: Option<[[f64; NODAL_DOFS]; ENHANCED_DOFS]>,
    /// K_αα⁻¹, for loads on the enhanced dofs (a thermal strain).
    pub kaa_inv: Option<[[f64; ENHANCED_DOFS]; ENHANCED_DOFS]>,
}

/// The fragment mean ξ̄ (local coordinates) the enhanced modes are
/// centred on: m₁₀₀/m₀₀₀ etc., or 0 for an empty fragment.
pub fn fragment_mean(m: &[f64; 27]) -> [f64; 3] {
    std::array::from_fn(|q| {
        let mut e = [0; 3];
        e[q] = 1;
        if m[0] > 0.0 {
            m[moment_index(e[0], e[1], e[2])] / m[0]
        } else {
            0.0
        }
    })
}

impl VectorTables {
    /// Builds the tables.
    pub fn new() -> VectorTables {
        let dn: [[Poly27; 3]; 8] =
            std::array::from_fn(|a| std::array::from_fn(|k| shape_derivative(a, k)));
        let mut nn = Vec::with_capacity(64);
        let mut n_xi = Vec::with_capacity(8);
        for a in 0..8 {
            n_xi.push(std::array::from_fn(|k| {
                std::array::from_fn(|m| mul(&dn[a][k], &xi(m)))
            }));
            for b in 0..8 {
                nn.push(std::array::from_fn(|k| {
                    std::array::from_fn(|l| mul(&dn[a][k], &dn[b][l]))
                }));
            }
        }
        VectorTables { nn, n_xi, dn }
    }

    /// The production cell operator: fragment-mean EAS, or plain Q1
    /// when the volume fraction is below [`ETA_EAS`]. Returns the cell
    /// and whether EAS was used.
    pub fn cell_operator(
        &self,
        m: &[f64; 27],
        h: f64,
        lambda: f64,
        mu: f64,
    ) -> Result<(VectorCell, bool), CholeskyError> {
        let use_eas = m[0] / (h * h * h) >= ETA_EAS;
        let enhancement = if use_eas {
            Enhancement::FragmentMean
        } else {
            Enhancement::None
        };
        Ok((self.cell(m, h, lambda, mu, enhancement)?, use_eas))
    }

    /// [`Self::cell_operator`] for an entry of a cell. A fragment of a
    /// split cell (`fragment` true) has lattice moments, which need not
    /// be consistent enough for the enhanced modes: when its enhanced
    /// block is not positive definite it falls back to plain Q1.
    pub fn entry_operator(
        &self,
        m: &[f64; 27],
        h: f64,
        lambda: f64,
        mu: f64,
        fragment: bool,
    ) -> Result<(VectorCell, bool), CholeskyError> {
        match self.cell_operator(m, h, lambda, mu) {
            Err(_) if fragment => Ok((self.cell(m, h, lambda, mu, Enhancement::None)?, false)),
            other => other,
        }
    }

    /// The cell stiffness for fragment moments `m` on a cell of size
    /// `h` with Lamé parameters (λ, μ) and the given enhancement; the 9
    /// bubble dofs are condensed out.
    pub fn cell(
        &self,
        m: &[f64; 27],
        h: f64,
        lambda: f64,
        mu: f64,
        enhancement: Enhancement,
    ) -> Result<VectorCell, CholeskyError> {
        self.cell_parts(m, &[(m, lambda, mu)], h, enhancement)
    }

    /// [`Self::cell`] for a cell of several materials: `parts` gives each
    /// material's moments (summing to `m`, the fragment's) and Lamé
    /// parameters. The blocks are summed before condensation, and the
    /// bubbles are centred on the whole fragment's mean.
    pub fn cell_parts(
        &self,
        m: &[f64; 27],
        parts: &[(&[f64; 27], f64, f64)],
        h: f64,
        enhancement: Enhancement,
    ) -> Result<VectorCell, CholeskyError> {
        let enhanced = enhancement != Enhancement::None;
        let mean = if enhancement == Enhancement::FragmentMean {
            fragment_mean(m)
        } else {
            [0.0; 3]
        };
        let mut kuu = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        let mut kua = [[0.0; ENHANCED_DOFS]; NODAL_DOFS];
        let mut kaa = [0.0; ENHANCED_DOFS * ENHANCED_DOFS];
        for &(pm, lambda, mu) in parts {
            let (uu, ua, aa) = self.blocks(pm, h, (lambda, mu), enhanced.then_some(mean));
            for r in 0..NODAL_DOFS {
                for c in 0..NODAL_DOFS {
                    kuu[r][c] += uu[r][c];
                }
                for c in 0..ENHANCED_DOFS {
                    kua[r][c] += ua[r][c];
                }
            }
            for (x, y) in kaa.iter_mut().zip(aa) {
                *x += y;
            }
        }
        if !enhanced {
            return Ok(VectorCell {
                k: kuu,
                recovery: None,
                kaa_inv: None,
            });
        }
        let chol = Cholesky::factor(&kaa, ENHANCED_DOFS, 0.0)?;
        // recovery = −Kαα⁻¹ Kαu; condensed K = Kuu + Kuα · recovery.
        let mut recovery = [[0.0; NODAL_DOFS]; ENHANCED_DOFS];
        for col in 0..NODAL_DOFS {
            let rhs: Vec<f64> = (0..ENHANCED_DOFS).map(|r| kua[col][r]).collect();
            let x = chol.solve(&rhs)?;
            for r in 0..ENHANCED_DOFS {
                recovery[r][col] = -x[r];
            }
        }
        let mut kaa_inv = [[0.0; ENHANCED_DOFS]; ENHANCED_DOFS];
        for col in 0..ENHANCED_DOFS {
            let mut e = vec![0.0; ENHANCED_DOFS];
            e[col] = 1.0;
            let x = chol.solve(&e)?;
            for r in 0..ENHANCED_DOFS {
                kaa_inv[r][col] = x[r];
            }
        }
        let mut k = kuu;
        for (row, krow) in k.iter_mut().enumerate() {
            for (col, v) in krow.iter_mut().enumerate() {
                *v += (0..ENHANCED_DOFS)
                    .map(|r| kua[row][r] * recovery[r][col])
                    .sum::<f64>();
            }
        }
        Ok(VectorCell {
            k,
            recovery: Some(recovery),
            kaa_inv: Some(kaa_inv),
        })
    }
    /// The uncondensed blocks (K_uu, K_uα, K_αα) of one material with
    /// moments `m` and Lamé (λ, μ); the bubbles centred on
    /// `enhanced_mean`, or none.
    #[allow(clippy::type_complexity)]
    fn blocks(
        &self,
        m: &[f64; 27],
        h: f64,
        (lambda, mu): (f64, f64),
        enhanced_mean: Option<[f64; 3]>,
    ) -> (
        Matrix24,
        [[f64; ENHANCED_DOFS]; NODAL_DOFS],
        [f64; ENHANCED_DOFS * ENHANCED_DOFS],
    ) {
        let scale = 4.0 / (h * h);
        // I[a][b][k][l] = ∫ ∂_k N_a ∂_l N_b.
        let int_nn =
            |a: usize, b: usize, k: usize, l: usize| integrate(&self.nn[a * 8 + b][k][l], m);
        let pair = |dfk: &dyn Fn(usize, usize) -> f64| -> [[f64; 3]; 3] {
            std::array::from_fn(|i| std::array::from_fn(|j| dfk(i, j)))
        };
        let mut kuu = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        for a in 0..8 {
            for b in 0..8 {
                let g = pair(&|k, l| int_nn(a, b, k, l));
                let trace = g[0][0] + g[1][1] + g[2][2];
                for i in 0..3 {
                    for j in 0..3 {
                        let diag = if i == j { mu * trace } else { 0.0 };
                        kuu[3 * a + i][3 * b + j] =
                            scale * (lambda * g[i][j] + mu * g[j][i] + diag);
                    }
                }
            }
        }
        let Some(mean) = enhanced_mean else {
            return (
                kuu,
                [[0.0; ENHANCED_DOFS]; NODAL_DOFS],
                [0.0; ENHANCED_DOFS * ENHANCED_DOFS],
            );
        };
        // Bubble m: ∂_l φ_m = −2 (ξ_m − ξ̄_m) δ_lm, ξ̄_m the fragment mean.
        let volume = m[0];
        // ∫ ∂_k N_a ∂_l φ_m = −2 δ_lm (∫ ∂_k N_a ξ_m − ξ̄_m ∫ ∂_k N_a).
        let int_n_phi = |a: usize, k: usize, mm: usize| {
            -2.0 * (integrate(&self.n_xi[a][k][mm], m) - mean[mm] * integrate(&self.dn[a][k], m))
        };
        // ∫ ∂φ_m ∂φ_n (only l = m, l' = n): 4 ∫ (ξ_m − ξ̄_m)(ξ_n − ξ̄_n).
        let second = |p: usize, q: usize| {
            let mut e = [0; 3];
            e[p] += 1;
            e[q] += 1;
            m[moment_index(e[0], e[1], e[2])]
        };
        let first = |p: usize| mean[p] * volume;
        let int_phi_phi = |p: usize, q: usize| {
            4.0 * (second(p, q) - mean[q] * first(p) - mean[p] * first(q)
                + mean[p] * mean[q] * volume)
        };
        let mut kua = [[0.0; ENHANCED_DOFS]; NODAL_DOFS];
        for a in 0..8 {
            for bubble in 0..3 {
                // gradient of the bubble lives only in direction `bubble`
                let g: [[f64; 3]; 3] = std::array::from_fn(|k| {
                    std::array::from_fn(|l| {
                        if l == bubble {
                            int_n_phi(a, k, bubble)
                        } else {
                            0.0
                        }
                    })
                });
                let trace = g[0][0] + g[1][1] + g[2][2];
                for i in 0..3 {
                    for j in 0..3 {
                        let diag = if i == j { mu * trace } else { 0.0 };
                        kua[3 * a + i][3 * bubble + j] =
                            scale * (lambda * g[i][j] + mu * g[j][i] + diag);
                    }
                }
            }
        }
        let mut kaa = [0.0; ENHANCED_DOFS * ENHANCED_DOFS];
        for p in 0..3 {
            for q in 0..3 {
                let g: [[f64; 3]; 3] = std::array::from_fn(|k| {
                    std::array::from_fn(|l| {
                        if k == p && l == q {
                            int_phi_phi(p, q)
                        } else {
                            0.0
                        }
                    })
                });
                let trace = g[0][0] + g[1][1] + g[2][2];
                for i in 0..3 {
                    for j in 0..3 {
                        let diag = if i == j { mu * trace } else { 0.0 };
                        kaa[(3 * p + i) * ENHANCED_DOFS + 3 * q + j] =
                            scale * (lambda * g[i][j] + mu * g[j][i] + diag);
                    }
                }
            }
        }
        let mut kaa = [0.0; ENHANCED_DOFS * ENHANCED_DOFS];
        for p in 0..3 {
            for q in 0..3 {
                let g: [[f64; 3]; 3] = std::array::from_fn(|k| {
                    std::array::from_fn(|l| {
                        if k == p && l == q {
                            int_phi_phi(p, q)
                        } else {
                            0.0
                        }
                    })
                });
                let trace = g[0][0] + g[1][1] + g[2][2];
                for i in 0..3 {
                    for j in 0..3 {
                        let diag = if i == j { mu * trace } else { 0.0 };
                        kaa[(3 * p + i) * ENHANCED_DOFS + 3 * q + j] =
                            scale * (lambda * g[i][j] + mu * g[j][i] + diag);
                    }
                }
            }
        }
        (kuu, kua, kaa)
    }
}
