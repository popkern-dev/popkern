//! Scalar Q1 cell operators (conduction) from fragment moments.
//!
//! On a cell of size h the trilinear shape functions N_a are
//! polynomials of degree ≤ 1 per local axis, so ∇N_a·∇N_b and N_a N_b
//! have degree ≤ 2 per axis and their integrals over any fragment are
//! exact linear combinations of its 27 moments. The coefficient tables
//! are built once and passed explicitly.

use crate::clip::FacetPiece;
use crate::moments::triangle::monomial_integrals;

use super::poly::{Poly27, integrate, mul, shape, shape_derivative};

/// An 8×8 cell matrix in `CELL_CORNERS` order.
pub type Matrix8 = [[f64; 8]; 8];

/// Polynomial tables for the scalar Q1 cell.
#[derive(Debug, Clone)]
pub struct ScalarTables {
    /// Σ_d ∂_d N_a ∂_d N_b in local derivatives.
    grad_grad: [[Poly27; 8]; 8],
    /// N_a N_b.
    value_value: [[Poly27; 8]; 8],
    /// N_a.
    value: [Poly27; 8],
    /// ∂_d N_a in local derivatives.
    derivative: [[Poly27; 3]; 8],
}

impl Default for ScalarTables {
    fn default() -> Self {
        Self::new()
    }
}

impl ScalarTables {
    /// Builds the tables.
    pub fn new() -> ScalarTables {
        let value: [Poly27; 8] = std::array::from_fn(shape);
        let d: [[Poly27; 3]; 8] =
            std::array::from_fn(|a| std::array::from_fn(|k| shape_derivative(a, k)));
        let grad_grad = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let mut s = [0.0; 27];
                for (da, db) in d[a].iter().zip(&d[b]) {
                    for (x, y) in s.iter_mut().zip(mul(da, db)) {
                        *x += y;
                    }
                }
                s
            })
        });
        let value_value =
            std::array::from_fn(|a| std::array::from_fn(|b| mul(&value[a], &value[b])));
        ScalarTables {
            grad_grad,
            value_value,
            value,
            derivative: d,
        }
    }

    /// Conduction matrix ∫ k ∇N_a·∇N_b dV for fragment moments `m`
    /// (physical volume units) on a cell of size `h`.
    pub fn conduction(&self, m: &[f64; 27], k: f64, h: f64) -> Matrix8 {
        let scale = 4.0 * k / (h * h);
        std::array::from_fn(|a| {
            std::array::from_fn(|b| scale * integrate(&self.grad_grad[a][b], m))
        })
    }

    /// Mass (capacity) matrix ∫ N_a N_b dV.
    pub fn mass(&self, m: &[f64; 27]) -> Matrix8 {
        std::array::from_fn(|a| std::array::from_fn(|b| integrate(&self.value_value[a][b], m)))
    }

    /// Load ∫ N_a dV (multiply by a uniform source density).
    pub fn volume_load(&self, m: &[f64; 27]) -> [f64; 8] {
        std::array::from_fn(|a| integrate(&self.value[a], m))
    }

    /// Facet mass ∫_Γ N_a N_b dA from surface moments (Robin terms).
    pub fn facet_mass(&self, s: &SurfaceMoments) -> Matrix8 {
        std::array::from_fn(|a| {
            std::array::from_fn(|b| integrate(&self.value_value[a][b], &s.area))
        })
    }

    /// Facet load ∫_Γ N_a dA (multiply by a uniform flux density).
    pub fn facet_load(&self, s: &SurfaceMoments) -> [f64; 8] {
        std::array::from_fn(|a| integrate(&self.value[a], &s.area))
    }

    /// ∫_Γ N_a (g · n) dA for a constant vector g: the Neumann load of
    /// a field whose flux density is g.
    pub fn normal_load(&self, s: &SurfaceMoments, g: [f64; 3]) -> [f64; 8] {
        std::array::from_fn(|a| {
            (0..3)
                .map(|d| g[d] * integrate(&self.value[a], &s.normal[d]))
                .sum()
        })
    }
}

impl ScalarTables {
    /// Nitsche matrix for a weakly fixed value on facets with surface
    /// moments `s`: −∫k(∇N_a·n)N_b − ∫k N_a(∇N_b·n) + β∫N_a N_b, on a
    /// cell of size `h`.
    pub fn nitsche(&self, s: &SurfaceMoments, k: f64, beta: f64, h: f64) -> Matrix8 {
        // c[a][b] = ∫ (∇N_a·n) N_b dA, physical gradient = (2/h)∂_ξ.
        let c: Matrix8 = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                (0..3)
                    .map(|d| integrate(&mul(&self.derivative[a][d], &self.value[b]), &s.normal[d]))
                    .sum::<f64>()
                    * 2.0
                    / h
            })
        });
        let mass = self.facet_mass(s);
        std::array::from_fn(|a| {
            std::array::from_fn(|b| -k * (c[a][b] + c[b][a]) + beta * mass[a][b])
        })
    }

    /// Nitsche load for fixed data `g` (a polynomial in local
    /// coordinates): −∫k(∇N_a·n)g + β∫N_a g.
    pub fn nitsche_load(
        &self,
        s: &SurfaceMoments,
        g: &Poly27,
        k: f64,
        beta: f64,
        h: f64,
    ) -> [f64; 8] {
        std::array::from_fn(|a| {
            let flux: f64 = (0..3)
                .map(|d| integrate(&mul(&self.derivative[a][d], g), &s.normal[d]))
                .sum::<f64>()
                * 2.0
                / h;
            -k * flux + beta * integrate(&mul(&self.value[a], g), &s.area)
        })
    }
}

/// Surface moments of the facet pieces in one cell: ∫ ξ^a η^b ζ^c dA
/// and the same weighted by each normal component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceMoments {
    /// ∫ ξ^a η^b ζ^c dA.
    pub area: [f64; 27],
    /// ∫ ξ^a η^b ζ^c n_d dA for d = x, y, z.
    pub normal: [[f64; 27]; 3],
}

/// Surface moments of `pieces` (all in the cell with minimum corner
/// `lo` and size `h`).
pub fn surface_moments<'a>(
    pieces: impl IntoIterator<Item = &'a FacetPiece>,
    lo: [f64; 3],
    h: f64,
) -> SurfaceMoments {
    let mut s = SurfaceMoments {
        area: [0.0; 27],
        normal: [[0.0; 27]; 3],
    };
    let local = |p: [f64; 3]| [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0);
    for piece in pieces {
        let v = piece.polygon.vertices();
        for i in 1..v.len().saturating_sub(1) {
            let tri = [v[0], v[i], v[i + 1]];
            let u = [0, 1, 2].map(|k| tri[1][k] - tri[0][k]);
            let w = [0, 1, 2].map(|k| tri[2][k] - tri[0][k]);
            let x = [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ];
            let area = 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            if area == 0.0 {
                continue;
            }
            let m = monomial_integrals(tri.map(local), area);
            for idx in 0..27 {
                let v = m[idx % 3][(idx / 3) % 3][idx / 9];
                s.area[idx] += v;
                for (normal, n) in s.normal.iter_mut().zip(piece.normal) {
                    normal[idx] += n * v;
                }
            }
        }
    }
    s
}
