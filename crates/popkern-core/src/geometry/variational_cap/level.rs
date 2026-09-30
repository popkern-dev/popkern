//! One refinement level: the hard system, the fixed rim rows, the rim
//! samples with their neighbour jets, and the solves over the null
//! space.

use std::sync::Arc;

use faer::Mat;

use super::energy::{Anchor, Integrand, energy_gram};
use super::jets::{NeighbourJet, NeighbourSampler};
use super::linalg::{cholesky, solve_lower};
use super::measure::RimSample;
use super::nullspace::{HardSystem, hard_system};
use super::rim::side_rows;
use super::space::{SectorSpace, gauss_legendre};
use super::{RimHalf, VariationalCapError, VariationalCapSide};
use crate::geometry::{Point3, Vec3};

/// Model ↔ normalised frame: `p_model = origin + scale · p`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Frame {
    pub origin: Point3,
    pub scale: f64,
}

impl Frame {
    pub fn to_local(self, p: Point3) -> Vec3 {
        Vec3::new(
            (p.x - self.origin.x) / self.scale,
            (p.y - self.origin.y) / self.scale,
            (p.z - self.origin.z) / self.scale,
        )
    }
}

/// The per-level state.
pub(super) struct Level {
    pub space: SectorSpace,
    pub hs: Arc<HardSystem>,
    /// Fixed (rim) values per coordinate, in `hs.fixed` order.
    pub fixed: [Vec<f64>; 3],
    /// Particular free values per coordinate.
    pub xp: [Vec<f64>; 3],
    /// Rim samples with neighbour jets (normalised frame).
    pub samples: Vec<RimSample>,
}

/// Build a level: rim rows, consistency check, samples.
pub(super) fn build_level(
    space: SectorSpace,
    sides: &[VariationalCapSide],
    frame: &Frame,
) -> Result<Level, VariationalCapError> {
    let n = space.n;
    let m = space.m;
    let hs = hard_system(&space).ok_or(VariationalCapError::Singular {
        stage: "spoke null space",
    })?;
    let mut full: [Vec<f64>; 3] = std::array::from_fn(|_| vec![f64::NAN; space.total()]);
    fn put(full: &mut [Vec<f64>; 3], idx: usize, p: Vec3) {
        for (c, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            full[c][idx] = v;
        }
    }
    let mut params = Vec::with_capacity(n);
    for (s, side) in sides.iter().enumerate() {
        let rows = side_rows(&space, side, s)?;
        let prev = (s + n - 1) % n;
        for (j, &p) in rows.trailing.iter().enumerate() {
            put(&mut full, space.index(s, m - 1, j), frame.to_local(p));
        }
        // The leading row's last point is corner c_{s−1}, written by side
        // s−1's trailing row (closure was checked against the tolerance
        // before any level was built).
        for (i, &p) in rows.leading.iter().enumerate().take(m - 1) {
            put(&mut full, space.index(prev, i, m - 1), frame.to_local(p));
        }
        params.push(rows.params);
    }
    let fixed: [Vec<f64>; 3] =
        std::array::from_fn(|c| hs.fixed.iter().map(|&i| full[c][i]).collect());
    if fixed.iter().flatten().any(|v| !v.is_finite()) {
        return Err(VariationalCapError::Inconsistent {
            stage: "rim rows incomplete",
            residual: f64::INFINITY,
        });
    }
    let mut xp: [Vec<f64>; 3] = Default::default();
    for c in 0..3 {
        let (x, residual) = hs.particular(&fixed[c]);
        if residual > 1.0e-9 {
            return Err(VariationalCapError::Inconsistent {
                stage: "spoke rows against rim rows",
                residual,
            });
        }
        xp[c] = x;
    }
    let samples = rim_samples(&space, sides, frame, &params)?;
    Ok(Level {
        space,
        hs,
        fixed,
        xp,
        samples,
    })
}

/// Gauss samples on every rim half (`degree + 3` per span) plus both
/// ends, each with its neighbour jet in the normalised frame.
fn rim_samples(
    space: &SectorSpace,
    sides: &[VariationalCapSide],
    frame: &Frame,
    params: &[(f64, f64, f64)],
) -> Result<Vec<RimSample>, VariationalCapError> {
    let n = space.n;
    let nodes = gauss_legendre(space.degree + 3);
    let mut xs: Vec<f64> = vec![0.0];
    for (a, b) in space.spans() {
        xs.extend(nodes.iter().map(|&(x, _)| a + (b - a) * x));
    }
    xs.push(1.0);
    let mut out = Vec::new();
    for (s, side) in sides.iter().enumerate() {
        let (t0, tm, t1) = params[s];
        for half in [RimHalf::Leading, RimHalf::Trailing] {
            let mut sampler = NeighbourSampler::new(&side.neighbour, side.neighbour_reversed)
                .map_err(|reason| VariationalCapError::NeighbourJet {
                    side: s,
                    at: tm,
                    reason,
                })?;
            let (sector, sign, h) = match half {
                RimHalf::Trailing => (s, 1.0, t1 - tm),
                RimHalf::Leading => ((s + n - 1) % n, -1.0, tm - t0),
            };
            for &x in &xs {
                let t = tm + sign * x * h;
                let p = crate::geometry::curve::CurveEval::point_at(&side.curve, t);
                let model =
                    sampler
                        .jet_at(&p)
                        .map_err(|reason| VariationalCapError::NeighbourJet {
                            side: s,
                            at: t,
                            reason,
                        })?;
                out.push(RimSample {
                    side: s,
                    half,
                    sector,
                    x,
                    t,
                    jet: localise(&model, frame),
                });
            }
        }
    }
    Ok(out)
}

/// A model-frame jet in the normalised frame (curvature × scale).
fn localise(jet: &NeighbourJet, frame: &Frame) -> NeighbourJet {
    let mut form = jet.form;
    for row in &mut form {
        for v in row.iter_mut() {
            *v *= frame.scale;
        }
    }
    let foot = frame.to_local(jet.foot);
    NeighbourJet {
        foot: Point3::new(foot.x, foot.y, foot.z),
        normal: jet.normal,
        form,
    }
}

/// The energy reduced to the null space and Cholesky-compressed:
/// `E(y) ∝ ‖Lᵀ y + r‖²` over `y = [y_x; y_y; y_z]`.
pub(super) struct EnergyFactor {
    pub l: Mat<f64>,
    pub r: Vec<f64>,
    /// The normalised Gram `L Lᵀ` (for the normal-equation solve).
    pub gram: Mat<f64>,
}

/// One sector's reduced energy: `(c, d, Zₛᵀ E_cd Zₛ, constant-term
/// row)` per nonzero coordinate block, and its linear term.
type SectorReduction = (Vec<(usize, usize, Mat<f64>, Vec<f64>)>, Vec<f64>);

/// Reduce the energy onto `level`'s null space. The Gram is normalised
/// to unit mean diagonal, so the fairness weight is a relative weight
/// per unknown.
pub(super) fn energy_factor(
    level: &Level,
    integrand: &Integrand<'_>,
    anchor: Option<&Anchor<'_>>,
) -> Result<EnergyFactor, VariationalCapError> {
    use faer::linalg::matmul::matmul;
    let space = &level.space;
    let hs = &level.hs;
    let k = hs.dim();
    let gram = energy_gram(space, integrand, anchor);
    let per = space.m * space.m;
    let size = gram.size;
    // Each sector's reduction `Zₛᵀ E_cd Zₛ` is independent; sectors
    // reduce in parallel and their contributions accumulate below in
    // the serial loop's order (sector, then c, then d), so `g` and `h`
    // are bit-identical to it (PK-598: 0.05 s per refresh at n = 8).
    let reduce = |s: usize| -> SectorReduction {
        let block = &gram.blocks[s];
        let base = s * per;
        let free_local: Vec<(usize, usize)> = (0..per)
            .filter_map(|l| match hs.slots[base + l] {
                (false, r) => Some((l, r)),
                _ => None,
            })
            .collect();
        let fixed_local: Vec<(usize, usize)> = (0..per)
            .filter_map(|l| match hs.slots[base + l] {
                (true, r) => Some((l, r)),
                _ => None,
            })
            .collect();
        let fl = free_local.len();
        let zs = Mat::from_fn(fl, k, |a, j| hs.z.read(free_local[a].1, j));
        let mut pairs = Vec::new();
        for c in 0..3 {
            for d in 0..3 {
                let e_uu = Mat::from_fn(fl, fl, |a, b| {
                    block[(c * per + free_local[a].0) * size + d * per + free_local[b].0]
                });
                if (0..fl).all(|a| (0..fl).all(|b| e_uu.read(a, b) == 0.0)) {
                    continue;
                }
                let mut ez = Mat::<f64>::zeros(fl, k);
                matmul(
                    ez.as_mut(),
                    e_uu.as_ref(),
                    zs.as_ref(),
                    None,
                    1.0,
                    faer::Parallelism::None,
                );
                let mut zez = Mat::<f64>::zeros(k, k);
                matmul(
                    zez.as_mut(),
                    zs.transpose(),
                    ez.as_ref(),
                    None,
                    1.0,
                    faer::Parallelism::None,
                );
                // Constants: E_cd (x_p,d on free rows, f_d on fixed rows).
                let mut ex = vec![0.0; fl];
                for (a, &(la, _)) in free_local.iter().enumerate() {
                    let row = (c * per + la) * size + d * per;
                    let mut acc = 0.0;
                    for &(lb, rb) in &free_local {
                        acc += block[row + lb] * level.xp[d][rb];
                    }
                    for &(lb, rb) in &fixed_local {
                        acc += block[row + lb] * level.fixed[d][rb];
                    }
                    ex[a] = acc;
                }
                let hc: Vec<f64> = (0..k)
                    .map(|i| (0..fl).map(|a| zs.read(a, i) * ex[a]).sum::<f64>())
                    .collect();
                pairs.push((c, d, zez, hc));
            }
        }
        // The linear term (E = xᵀBx − 2bᵀx) on the free coefficients.
        let lin = &gram.linear[s];
        let linear: Vec<f64> = (0..3)
            .flat_map(|c| {
                let (zs, free_local) = (&zs, &free_local);
                (0..k).map(move |i| {
                    free_local
                        .iter()
                        .enumerate()
                        .map(|(a, &(la, _))| zs.read(a, i) * lin[c * per + la])
                        .sum::<f64>()
                })
            })
            .collect();
        (pairs, linear)
    };
    #[cfg(feature = "parallel")]
    let reductions: Vec<SectorReduction> = {
        use rayon::prelude::*;
        (0..gram.blocks.len()).into_par_iter().map(reduce).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let reductions: Vec<SectorReduction> = (0..gram.blocks.len()).map(reduce).collect();
    let mut g = Mat::<f64>::zeros(3 * k, 3 * k);
    let mut h = vec![0.0; 3 * k];
    for (pairs, linear) in reductions {
        for (c, d, zez, hc) in pairs {
            for i in 0..k {
                for j in 0..k {
                    let at = (c * k + i, d * k + j);
                    g.write(at.0, at.1, g.read(at.0, at.1) + zez.read(i, j));
                }
            }
            for (i, v) in hc.into_iter().enumerate() {
                h[c * k + i] += v;
            }
        }
        for (hi, v) in h.iter_mut().zip(linear) {
            *hi -= v;
        }
    }
    let dim = 3 * k;
    let mean_diag = (0..dim).map(|i| g.read(i, i)).sum::<f64>() / dim.max(1) as f64;
    if mean_diag.is_nan() || mean_diag <= 0.0 {
        return Err(VariationalCapError::Singular { stage: "energy" });
    }
    for i in 0..dim {
        for j in 0..dim {
            g.write(i, j, g.read(i, j) / mean_diag);
        }
        g.write(i, i, g.read(i, i) + 1e-12);
    }
    let gram = g.clone();
    let l = cholesky(g).ok_or(VariationalCapError::Singular { stage: "energy" })?;
    let hn: Vec<f64> = h.iter().map(|v| v / mean_diag).collect();
    let r = solve_lower(&l, &hn);
    Ok(EnergyFactor { l, r, gram })
}

/// Full coefficient vectors from null-space coordinates.
pub(super) fn assemble(level: &Level, y: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
    let hs = &level.hs;
    let k = hs.dim();
    std::array::from_fn(|c| {
        let mut x = vec![0.0; level.space.total()];
        for (r, &i) in hs.free.iter().enumerate() {
            let mut v = level.xp[c][r];
            for (j, &yj) in y[c].iter().enumerate().take(k) {
                v += hs.z.read(r, j) * yj;
            }
            x[i] = v;
        }
        for (r, &i) in hs.fixed.iter().enumerate() {
            x[i] = level.fixed[c][r];
        }
        x
    })
}
