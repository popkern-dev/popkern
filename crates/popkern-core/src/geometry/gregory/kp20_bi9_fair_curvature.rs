//! The curvature-variation fair interior of the KP20 bi-9 cap (PK-601).
//!
//! `X = T₉ + D` minimises, under the cap's homogeneous relation rows,
//!
//! `CV_{T₉}(X) + μ ∫ |P_T(X − T₉)|² √g + ν TP3(X − S₉)`
//!
//! - `CV_{T₉}(X) = Σ ∫ g^{kl} ⟨∂_k S, ∂_l S⟩ √g` is the variation of the
//!   shape operator `S^i_j = g^{ik} II_kj`, `II_ij = N·(X_ij − Γ^l_ij X_l)`,
//!   with `N`, `g`, `Γ` frozen at `T₉` (Moreton–Séquin's curvature
//!   variation, linearised). On a sphere `S = I/R`, so the energy
//!   vanishes; unlike the mean-curvature variation it also sees the
//!   principal curvatures splitting at constant `H` (PK-601 measured the
//!   mean-curvature variation tripling max κ on a strongly curved n = 3
//!   border: 36 → 75; curvature variation 31).
//! - The tangential anchor holds the parameterisation, which `CV` cannot
//!   see (PK-600's single-freeze attempt folded without it).
//! - `ν TP3(X − S₉)`, PK-600's energy, keeps the exact composition where
//!   curvature variation alone overshoots the rim transition.
//!
//! Each Gauss point adds columns to a factor `Q` of each Gram (`g⁻¹ =
//! LLᵀ` for the variation, the tangent frame for the anchor), and each
//! Gram's upper triangle is one BLAS-3 product `QQᵀ`. `∂_k S` is a
//! central difference of the frozen functional (`h = 1e-4` on
//! polynomial data: O(h²) relative), as PK-574's `energy_mcv` does.

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::kp20_bi9_seams::{Net, b9};
use crate::geometry::Vec3;

/// Unknowns of one piece: 100 control points × 3 coordinates, ordered
/// `3·(10 i + j) + coordinate`.
pub(super) const PIECE_DIM: usize = 300;

/// Central-difference step for `∂_k S` in the unit parameter square.
const DIFFERENCE_STEP: f64 = 1e-4;

/// Weight of the tangential anchor against the variation, each
/// normalised by its mean diagonal (PK-601 sweep with the mean-curvature
/// variation: 1e-1, 1e-2 and 1 moved max κ·r by < 0.07 on every
/// fixture; 1e-2 is the middle).
const ANCHOR_WEIGHT: f64 = 1e-2;

/// Weight of PK-600's third-order pull towards the guide sampling `S₉`
/// (PK-601 sweep: 1e-3 and 1e-2 left the pentagon flagship at 1.96 /
/// 1.95 against its 1.90 pin, 1 lifted the cube to 1.98; 1e-1 is the one
/// setting that met every fairness row).
const GUIDE_WEIGHT: f64 = 1e-1;

/// Gauss–Legendre points per direction for the frozen Grams (PK-601:
/// 6 and 8 agree on every fixture's max κ·r within 0.03; the Gram's
/// BLAS-3 product dominates the step and scales with the point count).
const GAUSS_POINTS: usize = 6;

/// Why the curvature-variation fair interior declined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Kp20CurvatureVariationDecline {
    /// The relation space could not be built.
    RelationSpace,
    /// A Gauss point of the reference cap has a degenerate frame.
    DegenerateFrame,
    /// The coupled solve's pattern or factor could not be built.
    Structure,
    /// The solve produced non-finite values.
    NonFinite,
    /// The refined residual stayed above its gate (relative value).
    Residual(f64),
}

/// The frozen fields at one parameter point.
struct Frame {
    normal: Vec3,
    /// `g⁻¹` as `[g^uu, g^uv, g^vv]`.
    inverse: [f64; 3],
    /// `Γ^·_uu`, `Γ^·_uv`, `Γ^·_vv`, each `[Γ^u, Γ^v]`.
    christoffel: [[f64; 2]; 3],
    sqrt_det: f64,
    su: Vec3,
    bu: [[f64; 3]; 10],
    bv: [[f64; 3]; 10],
}

/// Degree-9 basis values and first two derivatives at one abscissa.
type Basis = [[f64; 3]; 10];

fn basis(t: f64) -> Basis {
    std::array::from_fn(|i| b9(i, t))
}

/// The basis at every Gauss abscissa and its `±h` offsets, shared by
/// every piece: `[t, t + h, t − h]` per abscissa.
fn basis_table(gauss: &[(f64, f64)]) -> Vec<[Basis; 3]> {
    gauss
        .iter()
        .map(|&(t, _)| {
            [
                basis(t),
                basis(t + DIFFERENCE_STEP),
                basis(t - DIFFERENCE_STEP),
            ]
        })
        .collect()
}

fn frame(net: &Net, bu: &Basis, bv: &Basis) -> Option<Frame> {
    let (bu, bv) = (*bu, *bv);
    let partial = |p: usize, q: usize| {
        let mut acc = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            for j in 0..10 {
                let x = net[i][j];
                acc = acc + Vec3::new(x.x, x.y, x.z) * (bu[i][p] * bv[j][q]);
            }
        }
        acc
    };
    let (su, sv) = (partial(1, 0), partial(0, 1));
    let (suu, suv, svv) = (partial(2, 0), partial(1, 1), partial(0, 2));
    let normal = su.cross(&sv).normalized()?;
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f * f;
    if det.is_nan() || det <= 1e-24 * e * g {
        return None;
    }
    let inverse = [g / det, -f / det, e / det];
    // Γ^k_ij = g^{kl} <x_ij, x_l>.
    let gamma = |xij: Vec3| {
        let (p, q) = (xij.dot(&su), xij.dot(&sv));
        [
            inverse[0] * p + inverse[1] * q,
            inverse[1] * p + inverse[2] * q,
        ]
    };
    Some(Frame {
        normal,
        inverse,
        christoffel: [gamma(suu), gamma(suv), gamma(svv)],
        sqrt_det: det.sqrt(),
        su,
        bu,
        bv,
    })
}

/// The shape operator's four components `(S^u_u, S^u_v, S^v_u, S^v_v)`
/// at `(u, v)` as functionals over the piece unknowns, the fields frozen
/// from `net` at the same point.
fn shape_functionals(net: &Net, bu: &Basis, bv: &Basis) -> Option<[[f64; PIECE_DIM]; 4]> {
    let fr = frame(net, bu, bv)?;
    let (inv, [guu, guv, gvv]) = (fr.inverse, fr.christoffel);
    let n = [fr.normal.x, fr.normal.y, fr.normal.z];
    let mut out = [[0.0; PIECE_DIM]; 4];
    for i in 0..10 {
        for j in 0..10 {
            let (bi, bj) = (fr.bu[i], fr.bv[j]);
            let (du, dv) = (bi[1] * bj[0], bi[0] * bj[1]);
            let iiuu = bi[2] * bj[0] - guu[0] * du - guu[1] * dv;
            let iiuv = bi[1] * bj[1] - guv[0] * du - guv[1] * dv;
            let iivv = bi[0] * bj[2] - gvv[0] * du - gvv[1] * dv;
            let s = [
                inv[0] * iiuu + inv[1] * iiuv,
                inv[0] * iiuv + inv[1] * iivv,
                inv[1] * iiuu + inv[2] * iiuv,
                inv[1] * iiuv + inv[2] * iivv,
            ];
            let a = 3 * (10 * i + j);
            for (component, value) in s.into_iter().enumerate() {
                for c in 0..3 {
                    out[component][a + c] = n[c] * value;
                }
            }
        }
    }
    Some(out)
}

/// The piece's curvature-variation and tangential-anchor Grams
/// (`PIECE_DIM²`, upper triangle filled), frozen at `net`; `None` when a
/// Gauss point's frame is degenerate.
pub(super) fn curvature_piece_grams(
    net: &Net,
    gauss: &[(f64, f64)],
    table: &[[Basis; 3]],
) -> Option<(Mat<f64>, Mat<f64>)> {
    let points = gauss.len() * gauss.len();
    let mut q_variation = Mat::<f64>::zeros(PIECE_DIM, 8 * points);
    let mut q_anchor = Mat::<f64>::zeros(PIECE_DIM, 2 * points);
    let h = DIFFERENCE_STEP;
    let mut point = 0;
    for (iu, &(_, wu)) in gauss.iter().enumerate() {
        for (iv, &(_, wv)) in gauss.iter().enumerate() {
            let ([bu, bu_plus, bu_minus], [bv, bv_plus, bv_minus]) = (&table[iu], &table[iv]);
            let fr = frame(net, bu, bv)?;
            let (up, um) = (
                shape_functionals(net, bu_plus, bv)?,
                shape_functionals(net, bu_minus, bv)?,
            );
            let (vp, vm) = (
                shape_functionals(net, bu, bv_plus)?,
                shape_functionals(net, bu, bv_minus)?,
            );
            let root_w = (wu * wv * fr.sqrt_det).sqrt();
            // g⁻¹ = LLᵀ (2×2 Cholesky; g⁻¹ is SPD on a regular frame).
            let l11 = fr.inverse[0].sqrt();
            let l21 = fr.inverse[1] / l11;
            let l22 = (fr.inverse[2] - l21 * l21).max(0.0).sqrt();
            for component in 0..4 {
                let col = 8 * point + 2 * component;
                for a in 0..PIECE_DIM {
                    let ru = (up[component][a] - um[component][a]) / (2.0 * h);
                    let rv = (vp[component][a] - vm[component][a]) / (2.0 * h);
                    q_variation.write(a, col, root_w * (l11 * ru + l21 * rv));
                    q_variation.write(a, col + 1, root_w * l22 * rv);
                }
            }
            // P_T = t₁t₁ᵀ + t₂t₂ᵀ with an orthonormal tangent frame.
            let t1 = fr.su.normalized()?;
            let t2 = fr.normal.cross(&t1);
            for i in 0..10 {
                for j in 0..10 {
                    let basis = root_w * fr.bu[i][0] * fr.bv[j][0];
                    let a = 3 * (10 * i + j);
                    for (c, (x1, x2)) in [(t1.x, t2.x), (t1.y, t2.y), (t1.z, t2.z)]
                        .into_iter()
                        .enumerate()
                    {
                        q_anchor.write(a + c, 2 * point, basis * x1);
                        q_anchor.write(a + c, 2 * point + 1, basis * x2);
                    }
                }
            }
            point += 1;
        }
    }
    // The upper triangle of Q Qᵀ (PK-676: a fixed summation order).
    let gram = |q: &Mat<f64>| {
        let g = popkern_math::linalg::syrk_upper(&crate::geometry::dense::dense(q));
        debug_assert_eq!(g.nrows(), PIECE_DIM);
        crate::geometry::dense::mat(&g)
    };
    Some((gram(&q_variation), gram(&q_anchor)))
}

/// Entry `(a, b)` of a symmetric matrix stored as its upper triangle.
fn upper(m: &Mat<f64>, a: usize, b: usize) -> f64 {
    m.read(a.min(b), a.max(b))
}

/// Relax `nets` (which satisfy every relation) to the curvature-variation
/// fair interior, the reference frozen at the incoming `nets` (one step:
/// PK-601 measured re-frozen steps moving max κ by < 1 %). Returns the
/// refined residual of the coupled solve.
pub(super) fn fair_interior_curvature(
    n: usize,
    nets: &mut [[Net; 4]],
    s9: &[[Net; 4]],
    b: [[f64; 3]; 2],
    tau: &super::kp20_bi6_cap::Kp20MacroPatch,
) -> Result<f64, Kp20CurvatureVariationDecline> {
    use Kp20CurvatureVariationDecline as D;
    let space = super::kp20_bi9_relations::relation_space(n, b).ok_or(D::RelationSpace)?;
    let gauss = crate::geometry::variational_cap::space::gauss_legendre(GAUSS_POINTS);
    let table = basis_table(&gauss);
    // Each piece's frozen Grams read only its own net: under the
    // `parallel` feature they are formed concurrently, collected in
    // piece order (identical results; PK-600: 0.18 s serial at n = 8).
    let pieces: Vec<&Net> = nets.iter().flatten().collect();
    let grams = |net: &&Net| curvature_piece_grams(net, &gauss, &table);
    #[cfg(feature = "parallel")]
    let frozen: Option<Vec<_>> = {
        use rayon::prelude::*;
        pieces.par_iter().map(grams).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let frozen: Option<Vec<_>> = pieces.iter().map(grams).collect();
    let frozen = frozen.ok_or(D::DegenerateFrame)?;
    let guide: Vec<Vec<Vec<f64>>> = tau
        .pieces
        .iter()
        .map(super::kp20_bi9_fair::metric_third_order_gram)
        .collect();
    let trace = |m: &Mat<f64>| (0..PIECE_DIM).map(|a| m.read(a, a)).sum::<f64>();
    let dims = frozen.len() as f64 * PIECE_DIM as f64;
    let variation_scale = dims / frozen.iter().map(|g| trace(&g.0)).sum::<f64>().max(1e-300);
    let anchor_scale =
        ANCHOR_WEIGHT * dims / frozen.iter().map(|g| trace(&g.1)).sum::<f64>().max(1e-300);
    let guide_trace: f64 = guide
        .iter()
        .map(|g| (0..100).map(|a| g[a][a]).sum::<f64>())
        .sum();
    let guide_scale = GUIDE_WEIGHT * 400.0 / guide_trace.max(1e-300);
    let mut grams = Vec::with_capacity(frozen.len());
    let mut rhs = vec![0.0f64; 3 * space.count];
    for s in 0..n {
        for k in 0..4 {
            let (variation, anchor) = &frozen[4 * s + k];
            let flat = |net: &Net| -> Vec<f64> {
                (0..PIECE_DIM)
                    .map(|a| {
                        let p = net[a / 30][(a / 3) % 10];
                        [p.x, p.y, p.z][a % 3]
                    })
                    .collect()
            };
            let (x, target) = (flat(&nets[s][k]), flat(&s9[s][k]));
            let guide_entry = |a: usize, c: usize| {
                if a % 3 == c % 3 {
                    guide[k][a / 3][c / 3] * guide_scale
                } else {
                    0.0
                }
            };
            let combined = Mat::<f64>::from_fn(PIECE_DIM, PIECE_DIM, |a, c| {
                if c < a {
                    return 0.0;
                }
                variation.read(a, c) * variation_scale
                    + anchor.read(a, c) * anchor_scale
                    + guide_entry(a, c)
            });
            for a in 0..PIECE_DIM {
                let unknown = space.index[s][k][a / 30][(a / 3) % 10];
                if unknown == usize::MAX {
                    continue;
                }
                // At D = 0 the anchor's gradient vanishes; the variation
                // pulls towards a constant shape operator, the guide term
                // towards S₉.
                let mut r = 0.0;
                for c in 0..PIECE_DIM {
                    r -= upper(variation, a, c) * variation_scale * x[c]
                        + guide_entry(a, c) * (x[c] - target[c]);
                }
                rhs[3 * unknown + a % 3] = r;
            }
            grams.push(combined);
        }
    }
    let (displacement, residual) = super::kp20_bi9_kkt::solve_coupled(n, b, &space, &grams, &rhs)
        .map_err(|e| match e {
        super::kp20_bi9_kkt::KktDecline::Structure => D::Structure,
        super::kp20_bi9_kkt::KktDecline::NonFinite => D::NonFinite,
        super::kp20_bi9_kkt::KktDecline::Residual(r) => D::Residual(r),
    })?;
    for s in 0..n {
        for k in 0..4 {
            for i in 0..10 {
                for j in 0..10 {
                    let x = space.index[s][k][i][j];
                    if x == usize::MAX {
                        continue;
                    }
                    let p = nets[s][k][i][j];
                    nets[s][k][i][j] = crate::geometry::Point3::new(
                        p.x + displacement[3 * x],
                        p.y + displacement[3 * x + 1],
                        p.z + displacement[3 * x + 2],
                    );
                }
            }
        }
    }
    Ok(residual)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    /// A bi-9 net whose control points sit on a sphere, optionally with a
    /// bump on the interior points.
    fn sphere_net(bump: f64) -> Net {
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let (u, v) = (i as f64 / 9.0 - 0.5, j as f64 / 9.0 - 0.5);
                let (x, y) = (0.4 * u, 0.4 * v);
                let z = (1.0 - x * x - y * y).sqrt();
                let lift = if (3..=6).contains(&i) && (3..=6).contains(&j) {
                    bump
                } else {
                    0.0
                };
                Point3::new(x, y, z + lift)
            })
        })
    }

    /// The frozen curvature-variation energy at its own reference.
    fn energy(net: &Net) -> f64 {
        let gauss = crate::geometry::variational_cap::space::gauss_legendre(GAUSS_POINTS);
        let (g, _) = curvature_piece_grams(net, &gauss, &basis_table(&gauss)).expect("regular net");
        let x: Vec<f64> = (0..PIECE_DIM)
            .map(|a| {
                let p = net[a / 30][(a / 3) % 10];
                [p.x, p.y, p.z][a % 3]
            })
            .collect();
        (0..PIECE_DIM)
            .map(|a| {
                (0..PIECE_DIM)
                    .map(|b| x[a] * upper(&g, a, b) * x[b])
                    .sum::<f64>()
            })
            .sum()
    }

    /// Near zero on a near-sphere, large once the interior is bumped.
    #[test]
    fn curvature_variation_is_small_on_a_sphere_and_sees_a_bump() {
        let sphere = energy(&sphere_net(0.0));
        let bumped = energy(&sphere_net(0.01));
        assert!(
            sphere < 1e-2 * bumped,
            "sphere {sphere:e} bumped {bumped:e}"
        );
    }

    /// The upper-triangular product equals the full `QQᵀ` there, with a
    /// non-negative diagonal.
    #[test]
    fn grams_are_positive_on_the_diagonal() {
        let gauss = crate::geometry::variational_cap::space::gauss_legendre(GAUSS_POINTS);
        let (g, t) = curvature_piece_grams(&sphere_net(0.003), &gauss, &basis_table(&gauss))
            .expect("regular");
        for a in 0..PIECE_DIM {
            assert!(g.read(a, a) >= 0.0 && t.read(a, a) >= 0.0);
        }
    }
}
