//! Boundary terms of vector cells on one clipped facet piece: Nitsche
//! supports (fixed or roller) and traction loads.
//!
//! The piece's normal n is constant, so every term is an integral over
//! the piece of N_a, ∂_k N_b or their products, computed from the
//! piece's surface moments. Nitsche is consistent: a displacement field
//! the element can represent is reproduced exactly (Thread 2).

use super::poly::{Poly27, integrate, mul, shape, shape_derivative};
use super::scalar::SurfaceMoments;
use super::vector::{ENHANCED_DOFS, Matrix24, NODAL_DOFS, xi};

/// A cell's condensed enhanced modes, for the boundary flux: α =
/// `recovery` · u, and bubble m has ∂_m φ_m = −2 (ξ_m − ξ̄_m) with ξ̄ the
/// fragment `mean`. With it, the Nitsche flux is the enhanced stress
/// D(B u + G̃ α), the stress whose energy the cell's stiffness holds;
/// the plain Q1 flux carries the locked volumetric part that the EAS
/// energy no longer sees, which near incompressibility outgrows β
/// (PH-034).
#[derive(Debug, Clone, Copy)]
pub struct EnhancedFlux<'a> {
    /// α = recovery · u.
    pub recovery: &'a [[f64; NODAL_DOFS]; ENHANCED_DOFS],
    /// Fragment mean ξ̄.
    pub mean: [f64; 3],
}

/// Traction σ n of the unit displacement gradient E_{im} (σ = λ tr I +
/// μ (E + Eᵀ)), component `q`.
fn unit_traction(i: usize, m: usize, q: usize, n: [f64; 3], lambda: f64, mu: f64) -> f64 {
    let d = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
    lambda * d(i, m) * n[q] + mu * (d(q, i) * n[m] + d(q, m) * n[i])
}

/// Polynomials for the boundary terms.
#[derive(Debug, Clone)]
pub struct BoundaryTables {
    /// N_a.
    value: [Poly27; 8],
    /// ∂_k N_a (local).
    derivative: [[Poly27; 3]; 8],
    /// N_a ∂_k N_b (local), indexed [a][b][k].
    value_derivative: Vec<[Poly27; 3]>,
    /// N_a N_b.
    value_value: Vec<Poly27>,
    /// N_a ξ_m, indexed [a][m].
    value_xi: [[Poly27; 3]; 8],
    /// ξ_m.
    xi: [Poly27; 3],
}

impl Default for BoundaryTables {
    fn default() -> Self {
        Self::new()
    }
}

/// A support on a piece.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PieceSupport {
    /// u = u0 + G·x on the piece (all components).
    Fixed {
        /// Displacement at the origin.
        u0: [f64; 3],
        /// Displacement gradient G[i][j] = ∂u_i/∂x_j.
        gradient: [[f64; 3]; 3],
    },
    /// u·n = 0 on the piece; tangential motion free. On a curved face
    /// each piece uses its own normal (radial on a cylinder).
    Roller,
    /// u·d = 0 for a fixed unit direction d; motion normal to d free.
    Directional {
        /// The constrained unit direction.
        direction: [f64; 3],
    },
}

impl BoundaryTables {
    /// Builds the tables.
    pub fn new() -> BoundaryTables {
        let value: [Poly27; 8] = std::array::from_fn(shape);
        let derivative: [[Poly27; 3]; 8] =
            std::array::from_fn(|a| std::array::from_fn(|k| shape_derivative(a, k)));
        let mut value_derivative = Vec::with_capacity(64);
        let mut value_value = Vec::with_capacity(64);
        for a in 0..8 {
            for b in 0..8 {
                value_derivative.push(std::array::from_fn(|k| mul(&value[a], &derivative[b][k])));
                value_value.push(mul(&value[a], &value[b]));
            }
        }
        let xi: [Poly27; 3] = std::array::from_fn(xi);
        let value_xi = std::array::from_fn(|a| std::array::from_fn(|m| mul(&value[a], &xi[m])));
        BoundaryTables {
            value,
            derivative,
            value_derivative,
            value_value,
            value_xi,
            xi,
        }
    }

    /// Traction load ∫ N_a t dA for a constant traction `t`.
    pub fn traction(&self, s: &SurfaceMoments, t: [f64; 3]) -> [f64; NODAL_DOFS] {
        let mut f = [0.0; NODAL_DOFS];
        for a in 0..8 {
            let w = integrate(&self.value[a], &s.area);
            for i in 0..3 {
                f[3 * a + i] = t[i] * w;
            }
        }
        f
    }

    /// Nitsche matrix and load of `support` on a piece with unit normal
    /// `n` and surface moments `s`, on a cell of size `h` with minimum
    /// corner `lo`, Lamé (λ, μ) and β; the flux includes the cell's
    /// enhanced modes when `enhanced` is given.
    #[allow(clippy::too_many_arguments)]
    pub fn support(
        &self,
        support: PieceSupport,
        s: &SurfaceMoments,
        n: [f64; 3],
        lo: [f64; 3],
        h: f64,
        lambda: f64,
        mu: f64,
        beta: f64,
        enhanced: Option<EnhancedFlux<'_>>,
    ) -> (Matrix24, [f64; NODAL_DOFS]) {
        let g = 2.0 / h;
        // c[a][b][k] = ∫ N_a ∂_k N_b dA (physical derivative).
        let c = |a: usize, b: usize, k: usize| {
            g * integrate(&self.value_derivative[a * 8 + b][k], &s.area)
        };
        let mass = |a: usize, b: usize| integrate(&self.value_value[a * 8 + b], &s.area);
        let dn = |a: usize, b: usize| (0..3).map(|k| n[k] * c(a, b, k)).sum::<f64>();
        let mut t = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        let mut k = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        // A roller constrains u·d only: d = n, or the given direction.
        let constrained = match support {
            PieceSupport::Fixed { .. } => None,
            PieceSupport::Roller => Some(n),
            PieceSupport::Directional { direction } => Some(direction),
        };
        let dot_dn = constrained.map_or(0.0, |d| (0..3).map(|k| d[k] * n[k]).sum::<f64>());
        for a in 0..8 {
            for b in 0..8 {
                let m = mass(a, b);
                let dnb = dn(a, b);
                for i in 0..3 {
                    for j in 0..3 {
                        // t: ∫ (σ(N_b e_j) n) · (test N_a e_i), or for a
                        // roller ∫ (d·σ(N_b e_j) n)(d · test N_a e_i).
                        t[3 * a + i][3 * b + j] = match constrained {
                            Some(d) => {
                                let ddb = (0..3).map(|k| d[k] * c(a, b, k)).sum::<f64>();
                                d[i] * (lambda * dot_dn * c(a, b, j)
                                    + mu * (d[j] * dnb + n[j] * ddb))
                            }
                            None => {
                                lambda * n[i] * c(a, b, j)
                                    + mu * (if i == j { dnb } else { 0.0 } + n[j] * c(a, b, i))
                            }
                        };
                        let pen = match constrained {
                            Some(d) => d[i] * d[j],
                            None if i == j => 1.0,
                            None => 0.0,
                        };
                        k[3 * a + i][3 * b + j] = beta * pen * m;
                    }
                }
            }
        }
        if let Some(e) = enhanced {
            // P[a][m] = ∫ N_a ∂_m φ_m dA (physical).
            let p: [[f64; 3]; 8] = std::array::from_fn(|a| {
                std::array::from_fn(|m| {
                    -2.0 * g
                        * (integrate(&self.value_xi[a][m], &s.area)
                            - e.mean[m] * integrate(&self.value[a], &s.area))
                })
            });
            // tr[i][r]: traction component i (its d-d part for a roller)
            // of bubble r = 3m + i' at unit amplitude, per ∫N_a ∂φ_m.
            let tr: [[f64; ENHANCED_DOFS]; 3] = std::array::from_fn(|i| {
                std::array::from_fn(|r| {
                    let (m, ip) = (r / 3, r % 3);
                    match constrained {
                        Some(d) => {
                            d[i] * (0..3)
                                .map(|q| d[q] * unit_traction(ip, m, q, n, lambda, mu))
                                .sum::<f64>()
                        }
                        None => unit_traction(ip, m, i, n, lambda, mu),
                    }
                })
            });
            for (a, pa) in p.iter().enumerate() {
                for (i, tri) in tr.iter().enumerate() {
                    for (col, slot) in t[3 * a + i].iter_mut().enumerate() {
                        // Traction of α(u = e_col), tested with N_a e_i.
                        *slot += (0..ENHANCED_DOFS)
                            .map(|r| e.recovery[r][col] * pa[r / 3] * tri[r])
                            .sum::<f64>();
                    }
                }
            }
        }
        for r in 0..NODAL_DOFS {
            for col in 0..NODAL_DOFS {
                k[r][col] -= t[r][col] + t[col][r];
            }
        }
        let mut f = [0.0; NODAL_DOFS];
        if let PieceSupport::Fixed { u0, gradient } = support {
            let centre = lo.map(|x| x + 0.5 * h);
            // g_k as a polynomial in local coordinates.
            let field: [Poly27; 3] = std::array::from_fn(|i| {
                let mut p = [0.0; 27];
                p[0] = u0[i] + (0..3).map(|j| gradient[i][j] * centre[j]).sum::<f64>();
                p[1] = 0.5 * h * gradient[i][0];
                p[3] = 0.5 * h * gradient[i][1];
                p[9] = 0.5 * h * gradient[i][2];
                p
            });
            for a in 0..8 {
                // d[l][q] = ∫ ∂_l N_a g_q dA (physical derivative).
                let d: [[f64; 3]; 3] = std::array::from_fn(|l| {
                    std::array::from_fn(|q| {
                        g * integrate(&mul(&self.derivative[a][l], &field[q]), &s.area)
                    })
                });
                for i in 0..3 {
                    // σ(N_a e_i) n · g
                    let flux = lambda * (0..3).map(|q| n[q] * d[i][q]).sum::<f64>()
                        + mu * ((0..3).map(|m| n[m] * d[m][i]).sum::<f64>()
                            + n[i] * (0..3).map(|q| d[q][q]).sum::<f64>());
                    let value = integrate(&mul(&self.value[a], &field[i]), &s.area);
                    f[3 * a + i] = -flux + beta * value;
                }
            }
            if let Some(e) = enhanced {
                // W[m][q] = ∫ ∂_m φ_m g_q dA; the test function's
                // enhanced flux against g.
                let w: [[f64; 3]; 3] = std::array::from_fn(|m| {
                    std::array::from_fn(|q| {
                        -2.0 * g
                            * (integrate(&mul(&self.xi[m], &field[q]), &s.area)
                                - e.mean[m] * integrate(&field[q], &s.area))
                    })
                });
                // c[r]: bubble r's traction against g.
                let c: [f64; ENHANCED_DOFS] = std::array::from_fn(|r| {
                    let (m, ip) = (r / 3, r % 3);
                    (0..3)
                        .map(|q| unit_traction(ip, m, q, n, lambda, mu) * w[m][q])
                        .sum()
                });
                for (col, slot) in f.iter_mut().enumerate() {
                    *slot -= (0..ENHANCED_DOFS)
                        .map(|r| e.recovery[r][col] * c[r])
                        .sum::<f64>();
                }
            }
        }
        (k, f)
    }

    /// ∫ (P σ(u) n) · (P σ(w) n) dA over a piece for unit nodal dofs u,
    /// w (the flux `support` imposes, P its constrained part: all of it
    /// for a fixed support, d·(σ n) along d for a roller), with the
    /// cell's enhanced modes when given: the form the Nitsche
    /// consistency term is bounded by (PH-068).
    #[allow(clippy::too_many_arguments)]
    pub fn flux_gram(
        &self,
        support: PieceSupport,
        s: &SurfaceMoments,
        n: [f64; 3],
        h: f64,
        lambda: f64,
        mu: f64,
        enhanced: Option<EnhancedFlux<'_>>,
    ) -> Matrix24 {
        let g = 2.0 / h;
        // Traction component q of each unit dof, a polynomial.
        let traction = |col: usize, q: usize| -> Poly27 {
            let (b, j) = (col / 3, col % 3);
            let mut p = [0.0; 27];
            for m in 0..3 {
                let c = g * unit_traction(j, m, q, n, lambda, mu);
                for (x, y) in p.iter_mut().zip(&self.derivative[b][m]) {
                    *x += c * y;
                }
            }
            if let Some(e) = enhanced {
                for r in 0..ENHANCED_DOFS {
                    let (m, ip) = (r / 3, r % 3);
                    let c = -2.0 * g * e.recovery[r][col] * unit_traction(ip, m, q, n, lambda, mu);
                    for (x, y) in p.iter_mut().zip(&self.xi[m]) {
                        *x += c * y;
                    }
                    p[0] -= c * e.mean[m];
                }
            }
            p
        };
        let parts: Vec<[Poly27; NODAL_DOFS]> = match support {
            PieceSupport::Fixed { .. } => (0..3)
                .map(|q| std::array::from_fn(|col| traction(col, q)))
                .collect(),
            PieceSupport::Roller | PieceSupport::Directional { .. } => {
                let d = match support {
                    PieceSupport::Directional { direction } => direction,
                    _ => n,
                };
                vec![std::array::from_fn(|col| {
                    let mut p = [0.0; 27];
                    for (q, dq) in d.iter().enumerate() {
                        for (x, y) in p.iter_mut().zip(traction(col, q)) {
                            *x += dq * y;
                        }
                    }
                    p
                })]
            }
        };
        let mut t = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        for part in &parts {
            for r in 0..NODAL_DOFS {
                for c in r..NODAL_DOFS {
                    let v = integrate(&mul(&part[r], &part[c]), &s.area);
                    t[r][c] += v;
                    if c != r {
                        t[c][r] += v;
                    }
                }
            }
        }
        t
    }

    /// ∫ (P N_a e_i) · (σ_k n) dA over a piece: the Nitsche flux of the
    /// known stress σ_k = D ε(α₀) − (3λ + 2μ) α θ I of a thermal strain
    /// (`element::thermal_strain`), moved to the right-hand side (the
    /// consistency term −∫ v·σ n with σ = σ(u) + σ_k). P is the support's
    /// constrained part; `theta` the nodal temperature rises, `alpha` the
    /// expansion; `enhanced` the fragment mean and α₀.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_flux(
        &self,
        support: PieceSupport,
        s: &SurfaceMoments,
        n: [f64; 3],
        h: f64,
        (lambda, mu): (f64, f64),
        (alpha, theta): (f64, &[f64; 8]),
        enhanced: Option<([f64; 3], &[f64; ENHANCED_DOFS])>,
    ) -> [f64; NODAL_DOFS] {
        let g = 2.0 / h;
        let c = (3.0 * lambda + 2.0 * mu) * alpha;
        let mut f = [0.0; NODAL_DOFS];
        for a in 0..8 {
            // ∫ N_a θ dA.
            let n_theta: f64 = (0..8)
                .map(|b| theta[b] * integrate(&self.value_value[a * 8 + b], &s.area))
                .sum();
            let mut t = n.map(|nq| -c * nq * n_theta);
            if let Some((mean, alpha0)) = enhanced {
                // ∫ N_a ∂_m φ_m dA (physical).
                let p: [f64; 3] = std::array::from_fn(|m| {
                    -2.0 * g
                        * (integrate(&self.value_xi[a][m], &s.area)
                            - mean[m] * integrate(&self.value[a], &s.area))
                });
                for (r, a0) in alpha0.iter().enumerate() {
                    let (m, ip) = (r / 3, r % 3);
                    for (q, slot) in t.iter_mut().enumerate() {
                        *slot += a0 * unit_traction(ip, m, q, n, lambda, mu) * p[m];
                    }
                }
            }
            let constrained = match support {
                PieceSupport::Fixed { .. } => None,
                PieceSupport::Roller => Some(n),
                PieceSupport::Directional { direction } => Some(direction),
            };
            for i in 0..3 {
                f[3 * a + i] = match constrained {
                    None => t[i],
                    Some(d) => d[i] * (0..3).map(|q| d[q] * t[q]).sum::<f64>(),
                };
            }
        }
        f
    }

    /// A pure penalty spring of `stiffness` (N/mm³) pulling the piece
    /// with surface moments `s` toward u = u0 + G·x, with no flux terms:
    /// a weak anchor (a force-driven submodel cap) rather than a
    /// constraint.
    pub fn spring(
        &self,
        s: &SurfaceMoments,
        stiffness: f64,
        target: ([f64; 3], [[f64; 3]; 3]),
        lo: [f64; 3],
        h: f64,
    ) -> (Matrix24, [f64; NODAL_DOFS]) {
        let (u0, gradient) = target;
        let mut k = [[0.0; NODAL_DOFS]; NODAL_DOFS];
        for a in 0..8 {
            for b in 0..8 {
                let m = stiffness * integrate(&self.value_value[a * 8 + b], &s.area);
                for i in 0..3 {
                    k[3 * a + i][3 * b + i] = m;
                }
            }
        }
        let centre = lo.map(|x| x + 0.5 * h);
        let field: [Poly27; 3] = std::array::from_fn(|i| {
            let mut p = [0.0; 27];
            p[0] = u0[i] + (0..3).map(|j| gradient[i][j] * centre[j]).sum::<f64>();
            p[1] = 0.5 * h * gradient[i][0];
            p[3] = 0.5 * h * gradient[i][1];
            p[9] = 0.5 * h * gradient[i][2];
            p
        });
        let mut f = [0.0; NODAL_DOFS];
        for a in 0..8 {
            for i in 0..3 {
                f[3 * a + i] = stiffness * integrate(&mul(&self.value[a], &field[i]), &s.area);
            }
        }
        (k, f)
    }
}

#[cfg(test)]
#[path = "vector_boundary_tests.rs"]
mod tests;
