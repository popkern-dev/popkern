//! The infill core's stress peak (PH-079, PH-110): the homogenised
//! gyroid's stress is the core cells' own field σ_h, read with the
//! core's material, and its strength check needs the worst of it.
//!
//! A homogenised core means nothing below its unit cell: Bean 2025's
//! strength is the yield of a whole periodic cell. So the core's stress
//! is averaged over a cube of the gyroid's cell size λ centred on each
//! core cell (core-volume weighted, over core cells only), and the peak
//! is the worst average. Unaveraged, the core's peak where a wall feeds
//! load into it grew from 2.48 to 2.92 times its parallel share from
//! 1 mm to 0.5 mm cells (PH-110).
//!
//! The core cells are full core cells and partial cells that are mostly
//! core (the cells that recover with the core's material), each read at
//! its core's centroid. Partial cells that are mostly wall recover with
//! the walls' material, so their core share is not read.

use serde::Serialize;

use crate::element::operators::Core;
use crate::grid::Grid;
use crate::moments::SolidMoments;

use super::cell::{StressPoly, centroid, stress_at};
use super::tensor::{Voigt, principal, von_mises};

/// The core's governing stress sample.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct CorePeak {
    /// The averaged stress, MPa (Voigt, tensor shears).
    pub stress: Voigt,
    /// The equivalent stress the core check reads, MPa:
    /// [`core_equivalent`] of `stress`.
    pub equivalent: f64,
    /// The centre of the averaging window, model coordinates.
    pub position: [f64; 3],
    /// Local relative ZZ error of the centre's cell, when estimated.
    pub error: Option<f64>,
    /// Core cells averaged.
    pub cells: usize,
}

/// The core's equivalent stress: the larger of von Mises and the
/// largest principal magnitude. Uniaxial stress reads its magnitude in
/// tension and compression alike (Bean 2025 models the gyroid's parent
/// as von Mises, symmetric in sign); the principal term keeps a
/// hydrostatic load, which von Mises cannot see but which collapses a
/// cellular core, from passing unchecked.
pub fn core_equivalent(s: &Voigt) -> f64 {
    let p = principal(s);
    von_mises(s).max(p[0].abs()).max(p[2].abs())
}

/// One core cell: grid indices, core volume, its stress and where.
struct Sample {
    ijk: [i64; 3],
    w: f64,
    s: Voigt,
    x: [f64; 3],
    error: Option<f64>,
}

/// A grid line of core cells along x: its (k, j), its cells' i and the
/// running sums of w and w·σ.
struct Line {
    key: (i64, i64),
    i: Vec<i64>,
    sums: Vec<[f64; 7]>,
}

impl Line {
    /// Σ (w, w·σ) over the line's cells with i in [lo, hi].
    fn range(&self, lo: i64, hi: i64) -> [f64; 7] {
        let a = self.i.partition_point(|&i| i < lo);
        let b = self.i.partition_point(|&i| i <= hi);
        let upper = self.sums[b];
        let lower = self.sums[a];
        std::array::from_fn(|q| upper[q] - lower[q])
    }
}

/// The worst `window`-averaged core stress of one solid, `None` without
/// a core cell.
pub fn core_peak(
    grid: &Grid,
    moments: &SolidMoments,
    cores: &[Core],
    polys: &[StressPoly],
    cell_error: &[Option<f64>],
    window: f64,
) -> Option<CorePeak> {
    let h = grid.cell_size();
    let mut samples = Vec::new();
    for (entry, ((cell, core), poly)) in moments.cells.iter().zip(cores).zip(polys).enumerate() {
        let (w, xi) = match core {
            Core::None => continue,
            Core::Full => (cell.m[0], centroid(&cell.m)),
            Core::Partial(m) if m[0] > 0.5 * cell.m[0] => (m[0], centroid(m)),
            Core::Partial(_) => continue,
        };
        let Some(ijk) = grid.cell_ijk(cell.cell) else {
            continue;
        };
        let lo = grid.cell_min(ijk);
        samples.push(Sample {
            ijk: ijk.map(|v| v as i64),
            w,
            s: stress_at(poly, xi),
            x: [0, 1, 2].map(|a| lo[a] + 0.5 * h * (xi[a] + 1.0)),
            error: cell_error.get(entry).copied().flatten(),
        });
    }
    samples.sort_by_key(|s| (s.ijk[2], s.ijk[1], s.ijk[0]));
    let mut lines: Vec<Line> = Vec::new();
    for s in &samples {
        let key = (s.ijk[2], s.ijk[1]);
        if lines.last().is_none_or(|l| l.key != key) {
            lines.push(Line {
                key,
                i: Vec::new(),
                sums: vec![[0.0; 7]],
            });
        }
        if let Some(line) = lines.last_mut() {
            let last = line.sums.last().copied().unwrap_or([0.0; 7]);
            let mut next = last;
            next[0] += s.w;
            for q in 0..6 {
                next[q + 1] += s.w * s.s[q];
            }
            line.i.push(s.ijk[0]);
            line.sums.push(next);
        }
    }
    // Cells either side of the centre: the window's half-width.
    let r = (0.5 * window / h).floor().max(0.0) as i64;
    let mut best: Option<CorePeak> = None;
    for s in &samples {
        let mut total = [0.0; 7];
        let mut cells = 0;
        for dk in -r..=r {
            for dj in -r..=r {
                let key = (s.ijk[2] + dk, s.ijk[1] + dj);
                let Ok(at) = lines.binary_search_by(|l| l.key.cmp(&key)) else {
                    continue;
                };
                let line = &lines[at];
                let part = line.range(s.ijk[0] - r, s.ijk[0] + r);
                for (t, p) in total.iter_mut().zip(part) {
                    *t += p;
                }
                cells += line.i.partition_point(|&i| i <= s.ijk[0] + r)
                    - line.i.partition_point(|&i| i < s.ijk[0] - r);
            }
        }
        if total[0] <= 0.0 {
            continue;
        }
        let stress: Voigt = std::array::from_fn(|q| total[q + 1] / total[0]);
        let equivalent = core_equivalent(&stress);
        if best.is_none_or(|b| equivalent > b.equivalent) {
            best = Some(CorePeak {
                stress,
                equivalent,
                position: s.x,
                error: s.error,
                cells,
            });
        }
    }
    best
}
