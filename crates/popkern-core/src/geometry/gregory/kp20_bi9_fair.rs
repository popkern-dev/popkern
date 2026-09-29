//! Fair interior for the KP20 bi-9 cap (PK-600): the deviation `D` from
//! the resampled cap that minimises the thin-plate energy of the cap's
//! departure from the accurate guide sampling `S₉`, under every
//! relation the cap keeps.
//!
//! The guide misfits the border along the whole rim (PK-599's M2), and
//! placement swaps the border rows in over three coefficient rows, so
//! the cap relaxes from the border to the guide within a few rows. Here
//! the relaxation is chosen, not inherited: `X = T₉ + D` minimises
//! `Σ_pieces TP(X − S₉)` with the rim bands and the o-blocks held and
//! internal C², sector-seam C⁰ and eq (8)/(9) exact. `T₉` satisfies every
//! relation, so the rows are homogeneous in `D`.

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::kp20_bi9_seams::{Net, b9};

fn binom(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |c, i| c * (n - i) as f64 / (i + 1) as f64)
}

/// Third-order (linearised mean-curvature-variation analogue) Gram of a
/// bi-9 piece: `∫∫ Σ_{i+j=3} C(3,i) |∂ᵢᵤ∂ʲᵥ D|² s⁻⁶ √det g`, with `s² =
/// tr g / 2` the isotropic scale of τ's frozen metric. Three held rim
/// rows clamp value, slope and second derivative; a sixth-order problem
/// takes all three, where bending forces a boundary layer.
fn metric_third_order_gram(tau: &super::kp20_bi6_cap::Kp20Bi6Patch) -> Vec<Vec<f64>> {
    let gl = crate::geometry::variational_cap::space::gauss_legendre(12);
    let bern = |m: usize, k: isize, t: f64| {
        if k < 0 || k as usize > m {
            0.0
        } else {
            binom(m, k as usize) * t.powi(k as i32) * (1.0 - t).powi((m - k as usize) as i32)
        }
    };
    // Degree-9 basis derivatives of order 0..=3.
    let d9 = |i: usize, p: usize, t: f64| -> f64 {
        let k = i as isize;
        match p {
            0 => bern(9, k, t),
            1 => 9.0 * (bern(8, k - 1, t) - bern(8, k, t)),
            2 => 72.0 * (bern(7, k - 2, t) - 2.0 * bern(7, k - 1, t) + bern(7, k, t)),
            _ => {
                504.0
                    * (bern(6, k - 3, t) - 3.0 * bern(6, k - 2, t) + 3.0 * bern(6, k - 1, t)
                        - bern(6, k, t))
            }
        }
    };
    let d6 = |i: usize, t: f64| 6.0 * (bern(5, i as isize - 1, t) - bern(5, i as isize, t));
    let mut g = vec![vec![0.0; 100]; 100];
    for &(u, wu) in &gl {
        for &(v, wv) in &gl {
            let (mut tu, mut tv) = ([0.0f64; 2], [0.0f64; 2]);
            for i in 0..7 {
                for j in 0..7 {
                    let p = tau.cps[i][j];
                    let (a, b) = (
                        d6(i, u) * bern(6, j as isize, v),
                        bern(6, i as isize, u) * d6(j, v),
                    );
                    tu[0] += a * p.x;
                    tu[1] += a * p.y;
                    tv[0] += b * p.x;
                    tv[1] += b * p.y;
                }
            }
            let (e, f, gg) = (
                tu[0] * tu[0] + tu[1] * tu[1],
                tu[0] * tv[0] + tu[1] * tv[1],
                tv[0] * tv[0] + tv[1] * tv[1],
            );
            let det = e * gg - f * f;
            if det <= 1e-30 {
                continue;
            }
            let s2 = 0.5 * (e + gg);
            let w = wu * wv * det.sqrt() / (s2 * s2 * s2);
            let terms: Vec<[f64; 4]> = (0..100)
                .map(|a| {
                    let (i, j) = (a / 10, a % 10);
                    [
                        d9(i, 3, u) * d9(j, 0, v),
                        d9(i, 2, u) * d9(j, 1, v),
                        d9(i, 1, u) * d9(j, 2, v),
                        d9(i, 0, u) * d9(j, 3, v),
                    ]
                })
                .collect();
            let c = [1.0, 3.0, 3.0, 1.0];
            for a in 0..100 {
                for b in a..100 {
                    let t: f64 = (0..4).map(|k| c[k] * terms[a][k] * terms[b][k]).sum();
                    g[a][b] += w * t;
                    if b != a {
                        g[b][a] += w * t;
                    }
                }
            }
        }
    }
    g
}

/// Relaxes `nets` (which satisfy every relation) to the fair interior:
/// per coordinate, the sparse saddle-point system
/// `[[G + ρI, Cᵀ], [C, −εI]] [D; λ] = [−G B; 0]` with `G` the
/// block-diagonal third-order Gram, `C` the unit-norm relation rows and
/// `B = nets − s9`. The `−εI` block absorbs redundant rows; the caller
/// re-measures every relation against its gates. Factorised once for
/// the three coordinates, strictly single-threaded.
pub(super) fn fair_interior(
    n: usize,
    nets: &mut [[Net; 4]],
    s9: &[[Net; 4]],
    b: [[f64; 3]; 2],
    tau: &super::kp20_bi6_cap::Kp20MacroPatch,
) -> Option<()> {
    use faer::dyn_stack::{GlobalPodBuffer, PodStack};
    use faer::sparse::SparseColMat;
    use faer::sparse::linalg::lu::{NumericLu, factorize_symbolic_lu};
    use faer::{Conj, Parallelism};
    let space = super::kp20_bi9_relations::relation_space(n, b)?;
    let (count, m) = (space.count, space.rows.len());
    let grams: Vec<Vec<Vec<f64>>> = tau.pieces.iter().map(metric_third_order_gram).collect();
    let trace: f64 = grams
        .iter()
        .map(|g| (0..100).map(|a| g[a][a]).sum::<f64>())
        .sum::<f64>()
        / 400.0;
    // The Gram is normalised by its mean diagonal so the relation block
    // (unit-norm rows) and the −εI regularisation share one O(1) scale.
    let scale = 1.0 / trace.max(1e-300);
    let (ridge, eps) = (1e-12, 1e-12);
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    let mut rhs = Mat::<f64>::zeros(count + m, 3);
    for s in 0..n {
        for k in 0..4 {
            let g = &grams[k];
            for a in 0..100 {
                let xa = space.index[s][k][a / 10][a % 10];
                if xa == usize::MAX {
                    continue;
                }
                let mut gb = [0.0f64; 3];
                for bi in 0..100 {
                    let gv = g[a][bi] * scale;
                    if gv == 0.0 {
                        continue;
                    }
                    let xb = space.index[s][k][bi / 10][bi % 10];
                    if xb != usize::MAX {
                        triplets.push((xa, xb, gv));
                    }
                    let d = nets[s][k][bi / 10][bi % 10] - s9[s][k][bi / 10][bi % 10];
                    gb[0] += gv * d.x;
                    gb[1] += gv * d.y;
                    gb[2] += gv * d.z;
                }
                for axis in 0..3 {
                    rhs.write(xa, axis, -gb[axis]);
                }
            }
        }
    }
    for x in 0..count {
        triplets.push((x, x, ridge));
    }
    for (r, row) in space.rows.iter().enumerate() {
        for &(x, w) in row {
            triplets.push((count + r, x, w));
            triplets.push((x, count + r, w));
        }
        triplets.push((count + r, count + r, -eps));
    }
    let size = count + m;
    let a = SparseColMat::<usize, f64>::try_new_from_triplets(size, size, &triplets).ok()?;
    let symbolic = factorize_symbolic_lu(a.symbolic(), Default::default()).ok()?;
    let mut numeric = NumericLu::new();
    let req = symbolic
        .factorize_numeric_lu_req::<f64>(Parallelism::None)
        .ok()?;
    let lu = symbolic
        .factorize_numeric_lu(
            &mut numeric,
            a.as_ref(),
            Parallelism::None,
            PodStack::new(&mut GlobalPodBuffer::new(req)),
        )
        .ok()?;
    let solve_req = symbolic
        .solve_in_place_req::<f64>(3, Parallelism::None)
        .ok()?;
    lu.solve_in_place_with_conj(
        Conj::No,
        rhs.as_mut(),
        Parallelism::None,
        PodStack::new(&mut GlobalPodBuffer::new(solve_req)),
    );
    for s in 0..n {
        for k in 0..4 {
            for i in 0..10 {
                for j in 0..10 {
                    let x = space.index[s][k][i][j];
                    if x == usize::MAX {
                        continue;
                    }
                    let p = nets[s][k][i][j];
                    let (dx, dy, dz) = (rhs.read(x, 0), rhs.read(x, 1), rhs.read(x, 2));
                    if !(dx.is_finite() && dy.is_finite() && dz.is_finite()) {
                        return None;
                    }
                    nets[s][k][i][j] = crate::geometry::Point3::new(p.x + dx, p.y + dy, p.z + dz);
                }
            }
        }
    }
    Some(())
}

/// The fair-interior decision for one cap (PK-600).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20FairInterior {
    /// Whether the fair interior serves (it measured fairer).
    pub served: bool,
    /// Largest principal curvature of the resampled cap before.
    pub max_curvature_before: f64,
    /// Largest principal curvature of the fair candidate, `None` when
    /// its solve failed.
    pub max_curvature_after: Option<f64>,
}

/// Outcome of the fair-interior step: the candidate nets when served.
pub(super) struct FairStep {
    pub candidate: Option<Vec<[Net; 4]>>,
    pub record: Kp20FairInterior,
}

/// Largest principal curvature over the nets' interiors (24×24 samples
/// per piece; exact partials).
pub(super) fn max_principal_curvature(nets: &[[Net; 4]]) -> f64 {
    use crate::geometry::Vec3;
    let mut worst = 0.0f64;
    for sector in nets {
        for net in sector {
            for a in 1..24 {
                for c in 1..24 {
                    let (u, v) = (a as f64 / 24.0, c as f64 / 24.0);
                    let (bu, bv): (Vec<[f64; 3]>, Vec<[f64; 3]>) = (
                        (0..10).map(|i| b9(i, u)).collect(),
                        (0..10).map(|j| b9(j, v)).collect(),
                    );
                    let d = |p: usize, q: usize| {
                        let mut acc = Vec3::new(0.0, 0.0, 0.0);
                        for i in 0..10 {
                            for j in 0..10 {
                                let x = net[i][j];
                                acc = acc + Vec3::new(x.x, x.y, x.z) * (bu[i][p] * bv[j][q]);
                            }
                        }
                        acc
                    };
                    let (su, sv, suu, suv, svv) = (d(1, 0), d(0, 1), d(2, 0), d(1, 1), d(0, 2));
                    let Some(n) = su.cross(&sv).normalized() else {
                        continue;
                    };
                    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
                    let (l, m, o) = (n.dot(&suu), n.dot(&suv), n.dot(&svv));
                    let det = e * g - f * f;
                    if det.abs() < 1e-30 {
                        continue;
                    }
                    let (kg, h) = (
                        (l * o - m * m) / det,
                        (e * o + g * l - 2.0 * f * m) / (2.0 * det),
                    );
                    worst = worst.max(h.abs() + (h * h - kg).max(0.0).sqrt());
                }
            }
        }
    }
    worst
}

/// Solve the fair interior from `nets` and serve it only when its
/// largest principal curvature does not exceed the start's. Following
/// the accurate guide sampling is fairer wherever the guide is (cube,
/// saddle, pyramids); where the guide itself is poor (pentagon5
/// variable measured 4.71 → 7.60) the start is kept, recorded.
pub(super) fn fair_interior_step(
    n: usize,
    nets: &[[Net; 4]],
    s9: &[[Net; 4]],
    b: [[f64; 3]; 2],
    tau: &super::kp20_bi6_cap::Kp20MacroPatch,
) -> FairStep {
    let before = max_principal_curvature(nets);
    let mut candidate = nets.to_vec();
    if fair_interior(n, &mut candidate, s9, b, tau).is_none() {
        return FairStep {
            candidate: None,
            record: Kp20FairInterior {
                served: false,
                max_curvature_before: before,
                max_curvature_after: None,
            },
        };
    }
    let after = max_principal_curvature(&candidate);
    let served = after <= before;
    FairStep {
        candidate: served.then_some(candidate),
        record: Kp20FairInterior {
            served,
            max_curvature_before: before,
            max_curvature_after: Some(after),
        },
    }
}
