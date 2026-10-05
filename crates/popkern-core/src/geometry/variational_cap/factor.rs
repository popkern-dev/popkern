//! The fairness energy over the null space as rows per sector (PK-620).
//!
//! Each sector's energy involves only its own coefficients:
//! `E = Σ_s (x_sᵀ E_s x_s + 2 c_sᵀ x_s) + const` over its free
//! coefficients `x_s = x_p,s + Z_s y` (the fixed values folded into
//! `c_s`). With `Ê_s = F_s F_sᵀ` the Cholesky factor of the normalised
//! block, sector s contributes the rows
//!
//! ```text
//! F_sᵀ Z_s y + (F_sᵀ x_p,s + F_s⁻¹ ĉ_s)
//! ```
//!
//! over the few null-space columns its coefficients touch (its interior,
//! the local modes of the spokes beside it, the global block). This is
//! the energy `‖Lᵀ y + r‖²` of the dense factor `L Lᵀ = Zᵀ E Z / mean`
//! up to a constant, which the constants' re-centring at the minimiser
//! removes (see [`energy_factor`]), without ever forming the dense
//! reduction: measured on the octagon h1.2
//! variable apex's 16-sector degree-9 fit (2754 unknowns), the dense
//! reduction and its Schur complement cost 1.0 s per refresh.
//!
//! The normalisation is the dense one's: the reduced Gram's mean
//! diagonal `tr(Zᵀ E Z) / 3k` (`Z` is orthonormal, so the trace is
//! `Σ_s tr(E_s Z_s Z_sᵀ)`), and the `1e-12` ridge on the reduced Gram is
//! the same ridge on each `Ê_s` (`x_p` is the minimum-norm particular
//! solution, orthogonal to the null space, so the ridge adds no linear
//! term).

use std::sync::Arc;

use faer::Mat;
use popkern_math::linalg;

use crate::geometry::dense::{dense, mat};

use super::VariationalCapError;
use super::energy::{Anchor, Integrand, energy_gram};
use super::level::Level;
use super::linalg::{cholesky, par_map};
use super::nullspace::HardSystem;

/// The ridge on the normalised Gram (the dense factor's).
const GRAM_RIDGE: f64 = 1e-12;

/// The energy as rows per sector: `E(y) ∝ Σ_s ‖a_s y_s + c_s‖²`, `y_s`
/// the sector's columns of `y`, with its minimum 0.
pub(super) struct EnergyFactor {
    /// The hard system the columns belong to.
    pub hs: Arc<HardSystem>,
    /// Null-space dimension per coordinate.
    pub k: usize,
    /// Each sector's rows.
    pub sectors: Vec<SectorRows>,
    /// The energy's minimiser.
    minimiser: [Vec<f64>; 3],
    /// The normalised linear term `h = Σ_s a_sᵀ c_s`, in the `[x; y; z]`
    /// order (the normal-equation solve's).
    pub h: Vec<f64>,
}

/// One sector's energy rows.
pub(super) struct SectorRows {
    /// The columns (in the `[x; y; z]` order: `c·k + j`) the rows touch,
    /// ascending within each coordinate.
    pub cols: Vec<usize>,
    /// `F_sᵀ Z_s` (rows × `cols`).
    pub a: Mat<f64>,
    /// The constant per row.
    pub c: Vec<f64>,
}

impl EnergyFactor {
    /// The energy at `y`, `Σ_s ‖a_s y_s + c_s‖²`: 0 at its minimiser.
    pub fn value(&self, y: &[Vec<f64>; 3]) -> f64 {
        let k = self.k;
        self.sectors
            .iter()
            .map(|sector| {
                let ys: Vec<f64> = sector.cols.iter().map(|&col| y[col / k][col % k]).collect();
                (0..sector.a.nrows())
                    .map(|i| {
                        let v = sector.c[i]
                            + ys.iter()
                                .enumerate()
                                .map(|(j, &yj)| sector.a.read(i, j) * yj)
                                .sum::<f64>();
                        v * v
                    })
                    .sum::<f64>()
            })
            .sum()
    }

    /// The energy's minimiser.
    pub fn minimiser(&self) -> Result<[Vec<f64>; 3], VariationalCapError> {
        Ok(self.minimiser.clone())
    }

    /// The normalised Gram, dense in the `[x; y; z]` order.
    pub fn dense_gram(&self) -> Mat<f64> {
        let mut g = Mat::<f64>::zeros(3 * self.k, 3 * self.k);
        for sector in &self.sectors {
            let a = dense(&sector.a);
            let ata = mat(&linalg::matmul_tn(&a, &a));
            for (a, &ca) in sector.cols.iter().enumerate() {
                for (b, &cb) in sector.cols.iter().enumerate() {
                    g.write(ca, cb, g.read(ca, cb) + ata.read(a, b));
                }
            }
        }
        g
    }

    /// Every sector's rows over the `[x; y; z]` columns, stacked, with
    /// their constants: the dense stacked QR's energy rows.
    pub fn dense_rows(&self) -> (Mat<f64>, Vec<f64>) {
        let total = self.sectors.iter().map(|s| s.a.nrows()).sum();
        let mut m = Mat::<f64>::zeros(total, 3 * self.k);
        let mut c = Vec::with_capacity(total);
        let mut row = 0;
        for sector in &self.sectors {
            for i in 0..sector.a.nrows() {
                for (j, &col) in sector.cols.iter().enumerate() {
                    m.write(row, col, sector.a.read(i, j));
                }
                c.push(sector.c[i]);
                row += 1;
            }
        }
        (m, c)
    }
}

/// One sector's raw energy on its free coefficients: `E_s` (three
/// coordinates, coordinate-major), the constant row `c_s` (fixed values
/// and the linear term), and `tr(E_s Z_s Z_sᵀ)`.
struct SectorEnergy {
    e: Mat<f64>,
    c: Vec<f64>,
    trace: f64,
}

/// Reduce the energy onto `level`'s null space as rows per sector.
pub(super) fn energy_factor(
    level: &Level,
    integrand: &Integrand<'_>,
    anchor: Option<&Anchor<'_>>,
) -> Result<EnergyFactor, VariationalCapError> {
    let space = &level.space;
    let hs = &level.hs;
    let k = hs.dim();
    let gram = energy_gram(space, integrand, anchor);
    let raw = par_map(space.n, |s| {
        sector_energy(level, &gram.blocks[s], &gram.linear[s], s, gram.size)
    });
    let trace: f64 = raw.iter().map(|r| r.trace).sum();
    let mean_diag = trace / (3 * k).max(1) as f64;
    if mean_diag.is_nan() || mean_diag <= 0.0 {
        return Err(VariationalCapError::Singular { stage: "energy" });
    }
    let rows = par_map(space.n, |s| sector_rows(level, &raw[s], s, mean_diag));
    let sectors: Vec<SectorRows> = rows
        .into_iter()
        .collect::<Option<_>>()
        .ok_or(VariationalCapError::Singular { stage: "energy" })?;
    let mut h = vec![0.0; 3 * k];
    for sector in &sectors {
        for (j, &col) in sector.cols.iter().enumerate() {
            h[col] += (0..sector.a.nrows())
                .map(|i| sector.a.read(i, j) * sector.c[i])
                .sum::<f64>();
        }
    }
    let mut factor = EnergyFactor {
        hs: hs.clone(),
        k,
        sectors,
        minimiser: Default::default(),
        h,
    };
    // The unconstrained minimum y*: the banded solve over the energy rows
    // alone. Each sector's constant is then re-centred to `−a_s y*_s`, so
    // the rows read `a (y − y*)`: the same quadratic (at y* the residual
    // `a y* + c` is orthogonal to every column, so the two differ by the
    // constant ‖a y* + c‖²) with its minimum 0, the dense factor's.
    // A least-squares solve's error grows with κ² times its residual, and
    // the irreducible residual was measured to cost accuracy: on the
    // octagon h0.4 variable apex's corner-cell fit the G1 step's
    // objective read 8% over the dense factor's, which one step of
    // iterative refinement recovered exactly (PK-620).
    let none = Mat::<f64>::zeros(0, 3 * k);
    let y = super::banded::solve(&factor, 1.0, &none, &[])
        .ok_or(VariationalCapError::Singular { stage: "energy" })?;
    for sector in &mut factor.sectors {
        let ys: Vec<f64> = sector.cols.iter().map(|&col| y[col / k][col % k]).collect();
        for (i, slot) in sector.c.iter_mut().enumerate() {
            *slot = -(0..ys.len())
                .map(|j| sector.a.read(i, j) * ys[j])
                .sum::<f64>();
        }
    }
    factor.minimiser = y;
    Ok(factor)
}

/// Sector `s`'s free-coefficient block of the energy Gram, its constant
/// row, and its share of the reduced trace.
fn sector_energy(
    level: &Level,
    block: &[f64],
    linear: &[f64],
    s: usize,
    size: usize,
) -> SectorEnergy {
    let hs = &level.hs;
    let per = level.space.per();
    let base = s * per;
    let rows = &hs.sector_rows[s];
    let nf = rows.len();
    let local: Vec<usize> = rows.iter().map(|&r| hs.free[r] - base).collect();
    let fixed: Vec<(usize, usize)> = (0..per)
        .filter_map(|l| match hs.slots[base + l] {
            (true, r) => Some((l, r)),
            _ => None,
        })
        .collect();
    let entry =
        |c: usize, la: usize, d: usize, lb: usize| block[(c * per + la) * size + d * per + lb];
    let e = Mat::from_fn(3 * nf, 3 * nf, |a, b| {
        entry(a / nf, local[a % nf], b / nf, local[b % nf])
    });
    // E_fF x_F − b on the free rows (the energy is xᵀBx − 2bᵀx).
    let c: Vec<f64> = (0..3 * nf)
        .map(|a| {
            let (cc, la) = (a / nf, local[a % nf]);
            let mut acc = -linear[cc * per + la];
            for d in 0..3 {
                for &(lb, r) in &fixed {
                    acc += entry(cc, la, d, lb) * level.fixed[d][r];
                }
            }
            acc
        })
        .collect();
    // tr(Z_sᵀ E_s Z_s) = Σ_c tr(E_cc P), P = Z_s Z_sᵀ.
    let z = &hs.sector_z[s];
    let zd = dense(z);
    let mut pd = linalg::DenseMat::zeros(nf, nf);
    linalg::matmul(&mut pd, None, &zd, &zd.transpose(), 1.0);
    let p = mat(&pd);
    let mut trace = 0.0;
    for cc in 0..3 {
        for a in 0..nf {
            for b in 0..nf {
                trace += e.read(cc * nf + a, cc * nf + b) * p.read(b, a);
            }
        }
    }
    SectorEnergy { e, c, trace }
}

/// Sector `s`'s rows `F_sᵀ Z_s y + (F_sᵀ x_p,s + F_s⁻¹ ĉ_s)`; `None`
/// when the normalised block is not positive definite.
fn sector_rows(level: &Level, raw: &SectorEnergy, s: usize, mean_diag: f64) -> Option<SectorRows> {
    let hs = &level.hs;
    let k = hs.dim();
    let rows = &hs.sector_rows[s];
    let nf = rows.len();
    let z = &hs.sector_z[s];
    let nc = z.ncols();
    let scale = 1.0 / mean_diag;
    let mut e = Mat::from_fn(3 * nf, 3 * nf, |a, b| raw.e.read(a, b) * scale);
    for a in 0..3 * nf {
        e.write(a, a, e.read(a, a) + GRAM_RIDGE);
    }
    let f = cholesky(e)?;
    // Ẑ = blockdiag(Z_s, Z_s, Z_s); a = Fᵀ Ẑ.
    let zhat = Mat::from_fn(3 * nf, 3 * nc, |a, b| {
        if a / nf == b / nc {
            z.read(a % nf, b % nc)
        } else {
            0.0
        }
    });
    let fd = dense(&f);
    let a = mat(&linalg::matmul_tn(&fd, &dense(&zhat)));
    // F⁻¹ ĉ, then + Fᵀ x_p.
    let mut fcd = linalg::DenseMat::from_fn(3 * nf, 1, |i, _| raw.c[i] * scale);
    linalg::solve_lower_in_place(&fd, &mut fcd);
    let fc = mat(&fcd);
    let xp: Vec<f64> = (0..3 * nf)
        .map(|i| level.xp[i / nf][rows[i % nf]])
        .collect();
    let c: Vec<f64> = (0..3 * nf)
        .map(|i| fc.read(i, 0) + (i..3 * nf).map(|j| f.read(j, i) * xp[j]).sum::<f64>())
        .collect();
    let cols = (0..3)
        .flat_map(|cc| hs.sector_cols[s].iter().map(move |&j| cc * k + j))
        .collect();
    Some(SectorRows { cols, a, c })
}
