//! The fairness energy as a dense Gram matrix per sector, over the
//! sector's coefficients of all three coordinates (local index
//! `coord·m² + i·m + j`).
//!
//! Two integrands, both evaluated by Gauss quadrature on every span
//! pair of `K`, with a metric G = RᵀR and M = R⁻ᵀ (so G⁻¹ = MᵀM):
//!
//! - **Parametric thin plate** per coordinate f:
//!   `tr(G⁻¹HG⁻¹H) + ε tr(G⁻¹∇fᵀ∇f) = ‖M H Mᵀ‖_F² + ε ‖M ∇fᵀ‖²`, H the
//!   parametric Hessian. Coordinate-separable. With the characteristic
//!   metric it is the regular thin plate at n = 4, but it is not
//!   geometric: it penalises the bilinear chart's own distortion, and
//!   over a spherical n-gon the interior sags with n (measured 0.6 % at
//!   n = 4 to 11 % of the radius at n = 8).
//! - **Bending** with a frozen reference normal N and metric G:
//!   `tr(G⁻¹BG⁻¹B)` with `B_ab = N·∂_a∂_b x`, i.e. κ₁² + κ₂² of the
//!   surface when the reference is the surface itself. The spherical
//!   cap is a critical point of this energy under sphere G1 data, so
//!   the fit reproduces it. It sees no tangential motion, so a small
//!   parametric term regularises the parameterisation.

use super::space::{Basis1, SectorSpace, gauss_legendre};
use crate::geometry::Vec3;

/// Per-sector Gram blocks, `(3m²) × (3m²)` row-major.
pub(super) struct EnergyGram {
    /// One dense block per sector.
    pub blocks: Vec<Vec<f64>>,
    /// Block side (`3m²`).
    pub size: usize,
}

/// What the integrand needs at a quadrature point.
pub(super) enum Integrand<'a> {
    /// Parametric thin plate with metric `metric(s, u, v)` and membrane
    /// weight `eps`.
    Parametric {
        metric: &'a dyn Fn(usize, f64, f64) -> [f64; 3],
        eps: f64,
    },
    /// Bending against `reference(s, u, v) -> (normal, metric)`, plus
    /// `regular` times the parametric thin plate and `eps` membrane.
    Bending {
        reference: &'a dyn Fn(usize, f64, f64) -> Option<(Vec3, [f64; 3])>,
        regular: f64,
        eps: f64,
    },
}

/// Assemble the Gram blocks for `space`.
pub(super) fn energy_gram(space: &SectorSpace, integrand: &Integrand<'_>) -> EnergyGram {
    let m = space.m;
    let per = m * m;
    let size = 3 * per;
    let nodes = gauss_legendre(space.degree + 2);
    let spans = space.spans();
    let mut blocks = Vec::with_capacity(space.n);
    for sector in 0..space.n {
        let mut block = vec![0.0; size * size];
        for &(ua, ub) in &spans {
            for &(va, vb) in &spans {
                for &(xu, wu) in &nodes {
                    let u = ua + (ub - ua) * xu;
                    let bu = space.basis(u);
                    for &(xv, wv) in &nodes {
                        let v = va + (vb - va) * xv;
                        let bv = space.basis(v);
                        let area = wu * wv * (ub - ua) * (vb - va);
                        match integrand {
                            Integrand::Parametric { metric, eps } => {
                                let g = metric(sector, u, v);
                                if let Some((mm, det)) = inverse_root(g) {
                                    let w = area * det.sqrt();
                                    parametric(&mut block, m, &bu, &bv, &mm, w, *eps);
                                }
                            }
                            Integrand::Bending {
                                reference,
                                regular,
                                eps,
                            } => {
                                let Some((normal, g)) = reference(sector, u, v) else {
                                    continue;
                                };
                                if let Some((mm, det)) = inverse_root(g) {
                                    let w = area * det.sqrt();
                                    bending(&mut block, m, &bu, &bv, &mm, normal, w);
                                    parametric(
                                        &mut block,
                                        m,
                                        &bu,
                                        &bv,
                                        &mm,
                                        w * regular,
                                        *eps / regular.max(1e-300),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        blocks.push(block);
    }
    EnergyGram { blocks, size }
}

/// `M = R⁻ᵀ` (lower) and `det G` for `G = [[a, b], [b, c]]`; `None`
/// when G is not positive definite.
fn inverse_root([a, b, c]: [f64; 3]) -> Option<([[f64; 2]; 2], f64)> {
    let det = a * c - b * b;
    if !(a > 0.0 && det > 0.0 && det.is_finite()) {
        return None;
    }
    let r11 = a.sqrt();
    let r12 = b / r11;
    let r22 = (c - r12 * r12).sqrt();
    Some(([[1.0 / r11, 0.0], [-r12 / (r11 * r22), 1.0 / r22]], det))
}

/// `‖M H Mᵀ‖_F²` as three squared functionals of the Hessian entries.
fn hessian_rows(mm: &[[f64; 2]; 2], h: [[f64; 2]; 2]) -> [f64; 3] {
    let k = |p: usize, q: usize| {
        let mut acc = 0.0;
        for c in 0..2 {
            for d in 0..2 {
                acc += mm[p][c] * h[c][d] * mm[q][d];
            }
        }
        acc
    };
    [k(0, 0), std::f64::consts::SQRT_2 * k(0, 1), k(1, 1)]
}

/// Basis-product supports: local index and `(value, ∂u, ∂v, ∂uu, ∂uv, ∂vv)`.
fn supports(m: usize, bu: &Basis1, bv: &Basis1) -> Vec<(usize, [f64; 6])> {
    let mut out = Vec::with_capacity(bu.funcs.len() * bv.funcs.len());
    for (ri, fu) in bu.funcs.iter().enumerate() {
        for (rj, fv) in bv.funcs.iter().enumerate() {
            out.push((
                (bu.first + ri) * m + bv.first + rj,
                [
                    fu[0] * fv[0],
                    fu[1] * fv[0],
                    fu[0] * fv[1],
                    fu[2] * fv[0],
                    fu[1] * fv[1],
                    fu[0] * fv[2],
                ],
            ));
        }
    }
    out
}

fn parametric(
    block: &mut [f64],
    m: usize,
    bu: &Basis1,
    bv: &Basis1,
    mm: &[[f64; 2]; 2],
    w: f64,
    eps: f64,
) {
    let per = m * m;
    let size = 3 * per;
    let se = eps.sqrt();
    let rows: Vec<(usize, [f64; 5])> = supports(m, bu, bv)
        .into_iter()
        .map(|(l, d)| {
            let h = hessian_rows(mm, [[d[3], d[4]], [d[4], d[5]]]);
            let g0 = mm[0][0] * d[1] + mm[0][1] * d[2];
            let g1 = mm[1][0] * d[1] + mm[1][1] * d[2];
            (l, [h[0], h[1], h[2], se * g0, se * g1])
        })
        .collect();
    for c in 0..3 {
        let base = c * per;
        for &(i, fi) in &rows {
            for &(j, fj) in &rows {
                let dot: f64 = fi.iter().zip(fj.iter()).map(|(a, b)| a * b).sum();
                block[(base + i) * size + base + j] += w * dot;
            }
        }
    }
}

fn bending(
    block: &mut [f64],
    m: usize,
    bu: &Basis1,
    bv: &Basis1,
    mm: &[[f64; 2]; 2],
    normal: Vec3,
    w: f64,
) {
    let per = m * m;
    let size = 3 * per;
    let nc = [normal.x, normal.y, normal.z];
    // Rows over (coord, local): B_ab = Σ_c N_c ∂_a∂_b x_c.
    let mut rows: Vec<(usize, [f64; 3])> = Vec::new();
    for (l, d) in supports(m, bu, bv) {
        let h = hessian_rows(mm, [[d[3], d[4]], [d[4], d[5]]]);
        for (c, ncc) in nc.iter().enumerate() {
            rows.push((c * per + l, [h[0] * ncc, h[1] * ncc, h[2] * ncc]));
        }
    }
    for &(i, fi) in &rows {
        for &(j, fj) in &rows {
            let dot = fi[0] * fj[0] + fi[1] * fj[1] + fi[2] * fj[2];
            block[i * size + j] += w * dot;
        }
    }
}
