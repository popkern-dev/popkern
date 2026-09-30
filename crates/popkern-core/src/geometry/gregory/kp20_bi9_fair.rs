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
pub(super) fn metric_third_order_gram(tau: &super::kp20_bi6_cap::Kp20Bi6Patch) -> Vec<Vec<f64>> {
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

/// The energy a fair-interior candidate minimised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20FairEnergy {
    /// PK-601: the frozen-reference curvature variation (the shape
    /// operator's variation), with a tangential anchor and a small pull
    /// towards the guide sampling.
    CurvatureVariation,
    /// PK-600: the third-order departure from the guide sampling `S₉`.
    GuideThirdOrder,
}

/// One fair-interior candidate: its energy, its largest principal
/// curvature (`None` when the solve declined) and why it declined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20FairAttempt {
    /// The candidate's energy.
    pub energy: Kp20FairEnergy,
    /// Largest principal curvature of the candidate, `None` when its
    /// solve declined.
    pub max_curvature_after: Option<f64>,
    /// The curvature-variation solve's decline, when it declined.
    pub decline: Option<super::kp20_bi9_fair_curvature::Kp20CurvatureVariationDecline>,
    /// Refined residual of the coupled solve (curvature-variation
    /// candidate).
    pub refinement_residual: Option<f64>,
}

/// The fair-interior decision for one cap (PK-600, PK-601): both
/// candidates measured, the fairer served when it does not exceed the
/// resampled cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20FairInterior {
    /// The energy whose candidate serves, `None` when the resampled cap
    /// is kept.
    pub served: Option<Kp20FairEnergy>,
    /// Largest principal curvature of the resampled cap before.
    pub max_curvature_before: f64,
    /// The candidates: curvature variation, then guide third-order.
    pub attempts: [Option<Kp20FairAttempt>; 2],
}

impl Kp20FairInterior {
    /// The served candidate's largest principal curvature.
    pub fn max_curvature_served(&self) -> Option<f64> {
        let served = self.served?;
        self.attempts
            .iter()
            .flatten()
            .find(|a| a.energy == served)
            .and_then(|a| a.max_curvature_after)
    }
}

/// Outcome of the fair-interior step: the candidate nets when served.
pub(super) struct FairStep {
    pub candidate: Option<Vec<[Net; 4]>>,
    pub record: Kp20FairInterior,
}

/// Interior samples per direction of [`max_principal_curvature`].
const CURVATURE_SAMPLES: usize = 24;

/// The degree-9 basis at the curvature samples `k / 24`, `k = 1..24`,
/// evaluated once (PK-604: re-evaluating it per sample cost 0.03 s per
/// cap measurement at n = 8).
pub(super) fn sample_basis() -> &'static [[[f64; 3]; 10]] {
    static TABLE: std::sync::OnceLock<Vec<[[f64; 3]; 10]>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        (1..CURVATURE_SAMPLES)
            .map(|k| {
                let t = k as f64 / CURVATURE_SAMPLES as f64;
                std::array::from_fn(|i| b9(i, t))
            })
            .collect()
    })
}

/// Largest principal curvature over the nets' interiors (23×23 samples
/// per piece at `k / 24`; exact partials).
pub(super) fn max_principal_curvature(nets: &[[Net; 4]]) -> f64 {
    let table = sample_basis();
    let mut worst = 0.0f64;
    for sector in nets {
        for net in sector {
            for bu in table {
                for bv in table {
                    worst = worst.max(point_curvature(net, bu, bv));
                }
            }
        }
    }
    worst
}

/// Largest principal curvature of one bi-9 net at the parameter whose
/// basis is `bu × bv` (0 where the frame is degenerate).
fn point_curvature(net: &Net, bu: &[[f64; 3]; 10], bv: &[[f64; 3]; 10]) -> f64 {
    use crate::geometry::Vec3;
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
        return 0.0;
    };
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (l, m, o) = (n.dot(&suu), n.dot(&suv), n.dot(&svv));
    let det = e * g - f * f;
    if det.abs() < 1e-30 {
        return 0.0;
    }
    let (kg, h) = (
        (l * o - m * m) / det,
        (e * o + g * l - 2.0 * f * m) / (2.0 * det),
    );
    h.abs() + (h * h - kg).max(0.0).sqrt()
}

/// The fair-interior step (PK-601): both candidates, the curvature
/// variation and PK-600's guide third-order, are solved and measured;
/// the fairer serves when it does not exceed the resampled cap's largest
/// principal curvature, otherwise the resampled cap is kept. Both are
/// measured because neither dominates: curvature variation is fairer on
/// every fillet fixture (cube 2.39 → 1.79, blunt pentagon 5.19 → 1.71),
/// while on a strongly curved synthetic n = 3 border it measured 31
/// against the guide candidate's 10.8. Following the guide alone is
/// worse where the guide is poor (pentagon5 variable, PK-600: 4.71 →
/// 7.60), which is why the gate stays.
pub(super) fn fair_interior_step(
    n: usize,
    nets: &[[Net; 4]],
    s9: &[[Net; 4]],
    b: [[f64; 3]; 2],
    tau: &super::kp20_bi6_cap::Kp20MacroPatch,
) -> FairStep {
    let before = max_principal_curvature(nets);
    // The two candidates are independent pure solves: under the
    // `parallel` feature they run concurrently (identical results either
    // way; PK-604: the guide candidate's 0.13 s at n = 8 hides behind
    // the curvature candidate).
    let curvature_candidate = || {
        let mut variation = nets.to_vec();
        let attempt = match super::kp20_bi9_fair_curvature::fair_interior_curvature(
            n,
            &mut variation,
            s9,
            b,
            tau,
        ) {
            Ok(residual) => Kp20FairAttempt {
                energy: Kp20FairEnergy::CurvatureVariation,
                max_curvature_after: Some(max_principal_curvature(&variation)),
                decline: None,
                refinement_residual: Some(residual),
            },
            Err(decline) => Kp20FairAttempt {
                energy: Kp20FairEnergy::CurvatureVariation,
                max_curvature_after: None,
                decline: Some(decline),
                refinement_residual: None,
            },
        };
        (variation, attempt)
    };
    let guide_candidate = || {
        let mut guide = nets.to_vec();
        let solved = fair_interior(n, &mut guide, s9, b, tau).is_some();
        let attempt = Kp20FairAttempt {
            energy: Kp20FairEnergy::GuideThirdOrder,
            max_curvature_after: solved.then(|| max_principal_curvature(&guide)),
            decline: None,
            refinement_residual: None,
        };
        (guide, attempt)
    };
    #[cfg(feature = "parallel")]
    let ((variation, variation_attempt), (guide, guide_attempt)) =
        rayon::join(curvature_candidate, guide_candidate);
    #[cfg(not(feature = "parallel"))]
    let ((variation, variation_attempt), (guide, guide_attempt)) =
        (curvature_candidate(), guide_candidate());
    let mut record = Kp20FairInterior {
        served: None,
        max_curvature_before: before,
        attempts: [Some(variation_attempt), Some(guide_attempt)],
    };
    // The fairest candidate at or under the resampled cap; ties go to
    // the curvature variation (the ladder's first rung).
    let mut best: Option<(f64, Kp20FairEnergy)> = None;
    for attempt in [variation_attempt, guide_attempt] {
        if let Some(k) = attempt.max_curvature_after
            && k <= before
            && best.is_none_or(|(b, _)| k < b)
        {
            best = Some((k, attempt.energy));
        }
    }
    record.served = best.map(|(_, energy)| energy);
    let candidate = match record.served {
        Some(Kp20FairEnergy::CurvatureVariation) => Some(variation),
        Some(Kp20FairEnergy::GuideThirdOrder) => Some(guide),
        None => None,
    };
    FairStep { candidate, record }
}
