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
//! - **Mean-curvature variation** `∫ |∇_g H|² √det G` with `H = N·Δ_g x`,
//!   N, g and the Laplace–Beltrami coefficients frozen at a reference
//!   surface (Schneider & Kobbelt's ΔH = 0 fairness, linearised). Zero on
//!   every constant-mean-curvature surface, so a sphere is an exact
//!   minimiser; ∇H is taken by central differences of the linear
//!   functional H. It sees no tangential motion, so a parametric thin
//!   plate regularises the parameterisation. Measured on spherical
//!   n-gons (theta 0.6), interior sag 1.8e-4 (n = 4) to 2.0e-3 (n = 6)
//!   against 5e-3..3e-2 for a frozen-normal bending energy
//!   (∫ tr(G⁻¹BG⁻¹B), B = N·∂∂x), which is not stationary at a sphere.
//!   It is sixth order in x, so it needs the curvature rows at the rim:
//!   under G1 rows alone it folded (cube corner 0.05/0.2/0.3, pyramid4
//!   h=3), where the fourth-order bending energy serves.
//! - **Bending** (fourth order), the energy of a G1-only fit.

use super::corner_cells::CornerOwner;
use super::energy_mcv::{Reference, combine, mean_curvature_gradient, outer};
use super::space::{SectorSpace, gauss_legendre};
use crate::geometry::Vec3;

/// Per-sector Gram blocks, `(3m²) × (3m²)` row-major.
pub(super) struct EnergyGram {
    /// One dense block per sector.
    pub blocks: Vec<Vec<f64>>,
    /// Block side (`3m²`).
    pub size: usize,
    /// Linear term per sector: the energy is `xᵀ B x − 2 bᵀ x`.
    pub linear: Vec<Vec<f64>>,
}

/// Tangential anchoring: `weight · ∫ Σ_i (e_i·(x − x_ref))²`, with
/// `(x_ref, [e_1, e_2])` a reference point and orthonormal tangent frame.
/// It holds the parameterisation without resisting normal (shape)
/// motion: the fairness energies act on normal components only, and
/// the G1 rows are homogeneous in the cross derivative, so without it
/// a fit can collapse its parameterisation into a fold (measured), and
/// a parametric regulariser strong enough to prevent that also flattens
/// the shape (measured sag).
pub(super) struct Anchor<'a> {
    pub reference: &'a (dyn Fn(usize, f64, f64) -> Option<(Vec3, [Vec3; 2])> + Sync),
    pub weight: f64,
}

/// What the integrand needs at a quadrature point.
pub(super) enum Integrand<'a> {
    /// Parametric thin plate with metric `metric(s, u, v)` and membrane
    /// weight `eps`.
    Parametric {
        metric: &'a (dyn Fn(usize, f64, f64) -> [f64; 3] + Sync),
        eps: f64,
    },
    /// Bending `∫ tr(G⁻¹BG⁻¹B)`, `B_ab = N·∂_a∂_b x`, against a frozen
    /// reference `(normal, metric)`, plus `regular` times the parametric
    /// thin plate and `eps` membrane. Fourth order: the energy for a fit
    /// constrained at the rim to G1 only.
    Bending {
        reference: &'a (dyn Fn(usize, f64, f64) -> Option<(Vec3, [f64; 3])> + Sync),
        regular: f64,
        eps: f64,
    },
    /// Mean-curvature variation `∫ |∇_g H|²` against a frozen reference
    /// (see [`Reference`]), plus `regular` times the parametric thin
    /// plate and `eps` membrane.
    MeanCurvatureVariation {
        reference: &'a (dyn Fn(usize, f64, f64) -> Option<Reference> + Sync),
        regular: f64,
        eps: f64,
    },
}

/// One sector's energy, regulariser and anchor blocks and anchor
/// linear term.
type SectorBlocks = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

/// Assemble the Gram blocks for `space`.
pub(super) fn energy_gram(
    space: &SectorSpace,
    integrand: &Integrand<'_>,
    anchor: Option<&Anchor<'_>>,
) -> EnergyGram {
    let per = space.per();
    let size = 3 * per;
    let nodes = gauss_legendre(space.degree + 2);
    let spans = space.spans();
    // The mean-curvature-variation part and its parametric regulariser
    // are assembled apart and balanced by trace before weighting: the
    // third-order MCV entries are orders larger than the thin plate's,
    // and under one shared normalisation the regulariser vanished and
    // the tangential directions went singular (measured: the cube corner
    // 0.05/0.2/0.3 declined `Singular`, pyramid4 h=3 folded).
    // Sectors are independent quadratures: they assemble in parallel
    // with an ordered collect, so the blocks are bit-identical to the
    // serial loop (PK-598: 0.09-0.11 s per bending refresh at n = 8).
    let assemble = |sector: usize| {
        let mut block = vec![0.0; size * size];
        let mut regular_block = vec![0.0; size * size];
        let mut anchor_block = vec![0.0; size * size];
        let mut anchor_linear = vec![0.0; size];
        for &(ua, ub) in &spans {
            for &(va, vb) in &spans {
                for &(xu, wu) in &nodes {
                    let u = ua + (ub - ua) * xu;
                    let bu = space.basis(u);
                    for &(xv, wv) in &nodes {
                        let v = va + (vb - va) * xv;
                        let bv = space.basis(v);
                        let area = wu * wv * (ub - ua) * (vb - va);
                        let sup = space.local_supports(&bu, &bv, CornerOwner::Interior);
                        if let Some(a) = anchor
                            && let Some((p, frame)) = (a.reference)(sector, u, v)
                        {
                            anchor_rows(
                                &mut anchor_block,
                                &mut anchor_linear,
                                per,
                                &sup,
                                p,
                                &frame,
                                area,
                            );
                        }
                        match integrand {
                            Integrand::Parametric { metric, eps } => {
                                let g = metric(sector, u, v);
                                if let Some((mm, det)) = inverse_root(g) {
                                    let w = area * det.sqrt();
                                    parametric(&mut block, per, &sup, &mm, w, *eps);
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
                                    bending(&mut block, per, &sup, &mm, normal, w);
                                    parametric(&mut regular_block, per, &sup, &mm, w, *eps);
                                    let _ = regular;
                                }
                            }
                            Integrand::MeanCurvatureVariation {
                                reference,
                                regular,
                                eps,
                            } => {
                                let Some(r) = reference(sector, u, v) else {
                                    continue;
                                };
                                let Some((mm, det)) = inverse_root(r.metric) else {
                                    continue;
                                };
                                let w = area * det.sqrt();
                                let h = 1e-3 * (ub - ua).min(vb - va);
                                if let Some([hu, hv]) =
                                    mean_curvature_gradient(space, sector, u, v, h, *reference)
                                {
                                    for mp in &mm {
                                        let row = combine(&[(mp[0], &hu), (mp[1], &hv)]);
                                        outer(&mut block, size, &row, w);
                                    }
                                }
                                let _ = regular;
                                parametric(&mut regular_block, per, &sup, &mm, w, *eps);
                            }
                        }
                    }
                }
            }
        }
        (block, regular_block, anchor_block, anchor_linear)
    };
    #[cfg(feature = "parallel")]
    let per_sector: Vec<SectorBlocks> = {
        use rayon::prelude::*;
        (0..space.n).into_par_iter().map(assemble).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let per_sector: Vec<SectorBlocks> = (0..space.n).map(assemble).collect();
    let mut blocks = Vec::with_capacity(space.n);
    let mut regular_blocks = Vec::with_capacity(space.n);
    let mut anchor_blocks = Vec::with_capacity(space.n);
    let mut linear = Vec::with_capacity(space.n);
    for (block, regular_block, anchor_block, anchor_linear) in per_sector {
        blocks.push(block);
        regular_blocks.push(regular_block);
        anchor_blocks.push(anchor_block);
        linear.push(anchor_linear);
    }
    if let Integrand::MeanCurvatureVariation { regular, .. } | Integrand::Bending { regular, .. } =
        integrand
    {
        let trace = |bs: &[Vec<f64>]| -> f64 {
            bs.iter()
                .map(|b| (0..size).map(|i| b[i * size + i]).sum::<f64>())
                .sum()
        };
        let (t_mcv, t_reg) = (trace(&blocks), trace(&regular_blocks));
        let scale = if t_mcv > 0.0 { t_reg / t_mcv } else { 0.0 };
        let t_anchor = trace(&anchor_blocks);
        let anchor_scale = match anchor {
            Some(a) if t_anchor > 0.0 => a.weight * t_reg / t_anchor,
            _ => 0.0,
        };
        for ((b, r), a) in blocks.iter_mut().zip(&regular_blocks).zip(&anchor_blocks) {
            for ((x, y), z) in b.iter_mut().zip(r).zip(a) {
                *x = scale * *x + regular * y + anchor_scale * z;
            }
        }
        for l in &mut linear {
            for x in l.iter_mut() {
                *x *= anchor_scale;
            }
        }
    } else {
        for l in &mut linear {
            l.fill(0.0);
        }
    }
    EnergyGram {
        blocks,
        size,
        linear,
    }
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

fn parametric(
    block: &mut [f64],
    per: usize,
    sup: &[(usize, [f64; 6])],
    mm: &[[f64; 2]; 2],
    w: f64,
    eps: f64,
) {
    let size = 3 * per;
    let se = eps.sqrt();
    let rows: Vec<(usize, [f64; 5])> = sup
        .iter()
        .map(|&(l, d)| {
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

/// Anchor rows at one point: for each tangent `e`, the row
/// `Σ_c e_c N_l x_{c,l}` with target `e·p`, accumulated into the block
/// and the linear term (`E = xᵀBx − 2bᵀx`).
#[allow(clippy::too_many_arguments)]
fn anchor_rows(
    block: &mut [f64],
    linear: &mut [f64],
    per: usize,
    sup: &[(usize, [f64; 6])],
    p: Vec3,
    frame: &[Vec3; 2],
    w: f64,
) {
    let size = 3 * per;
    for e in frame {
        let ec = [e.x, e.y, e.z];
        let target = e.dot(&p);
        let mut row: Vec<(usize, f64)> = Vec::new();
        for &(l, d) in sup {
            for (c, ecc) in ec.iter().enumerate() {
                row.push((c * per + l, ecc * d[0]));
            }
        }
        for &(i, a) in &row {
            linear[i] += w * target * a;
            for &(j, b) in &row {
                block[i * size + j] += w * a * b;
            }
        }
    }
}

/// `‖M B Mᵀ‖_F²` with `B_ab = N·∂_a∂_b x`: rows over `(coord, local)`.
fn bending(
    block: &mut [f64],
    per: usize,
    sup: &[(usize, [f64; 6])],
    mm: &[[f64; 2]; 2],
    normal: Vec3,
    w: f64,
) {
    let size = 3 * per;
    let nc = [normal.x, normal.y, normal.z];
    let mut rows: Vec<(usize, [f64; 3])> = Vec::new();
    for &(l, d) in sup {
        let h = hessian_rows(mm, [[d[3], d[4]], [d[4], d[5]]]);
        for (c, ncc) in nc.iter().enumerate() {
            rows.push((c * per + l, [h[0] * ncc, h[1] * ncc, h[2] * ncc]));
        }
    }
    for &(i, fi) in &rows {
        for &(j, fj) in &rows {
            block[i * size + j] += w * (fi[0] * fj[0] + fi[1] * fj[1] + fi[2] * fj[2]);
        }
    }
}
