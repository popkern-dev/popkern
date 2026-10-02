//! The hard system: exact spoke rows over the free coefficients, with
//! the rim rows (every sector's `i = m−1` and `j = m−1` control rows)
//! fixed by the rim approximant.
//!
//! Per coordinate, `C_U x_U + C_F x_F = 0`. The null space `Z` of `C_U`
//! and the pseudo-inverse pieces for the particular solution depend
//! only on the space, so they are computed once and memoised.
//!
//! `Z` is built in the spokes' structure (PK-620; see
//! `nullspace_basis`): each sector's interior, each spoke's local
//! modes, and a small global block, orthonormal as a whole.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use faer::Mat;

use super::chart::CharacteristicMap;
use super::nullspace_basis::{BasisRow, ColumnLayout, ParticularSolve, nested_basis};
use super::space::{Functional, SectorSpace};
use super::spoke::spoke_rows;

/// Relative singular-value threshold for the numerical rank.
pub(super) const RANK_TOL: f64 = 1.0e-10;

/// The eliminated hard system for one space.
///
/// Columns of `Z` (per coordinate) are ordered sector interiors (sector
/// 0, 1, ...), then each spoke's local modes (spoke 0, 1, ...), then the
/// global block; see [`ColumnLayout`]. Sector s's free coefficients
/// touch only its interior, the modes of the spokes on either side of
/// it, and the global block, which is the structure the solves exploit
/// (see `banded`).
#[derive(Debug)]
pub(super) struct HardSystem {
    /// Global indices of the free coefficients.
    pub free: Vec<usize>,
    /// Global indices of the fixed (rim) coefficients.
    pub fixed: Vec<usize>,
    /// Per global index: `(fixed?, row in free or fixed)`.
    pub slots: Vec<(bool, usize)>,
    /// The column blocks.
    pub layout: ColumnLayout,
    /// Each free coefficient's basis row: `(first column, length,
    /// offset into `values`)` segments.
    rows: Vec<BasisRow>,
    values: Vec<f64>,
    /// Per sector: its free rows, the columns they touch (ascending),
    /// and `Z` restricted to them (free rows × columns).
    pub sector_rows: Vec<Vec<usize>>,
    pub sector_cols: Vec<Vec<usize>>,
    pub sector_z: Vec<Mat<f64>>,
    /// The spoke rows and the solve of `C_U x = b` in their structure.
    spoke_rows: Vec<Functional>,
    particular: ParticularSolve,
}

type Key = (usize, usize, bool, Vec<u64>);

static MEMO: Mutex<Option<HashMap<Key, Arc<HardSystem>>>> = Mutex::new(None);

/// The memoised hard system for `space`; `None` when an SVD workspace
/// cannot be sized.
pub(super) fn hard_system(space: &SectorSpace) -> Option<Arc<HardSystem>> {
    // The corner cells change the coefficient count, so they are part
    // of the key.
    let key: Key = (
        space.n,
        space.degree,
        space.corner_cells,
        space.knots.iter().map(|k| k.to_bits()).collect(),
    );
    if let Some(found) = MEMO
        .lock()
        .unwrap_or_else(|p| p.into_inner())
        .as_ref()
        .and_then(|m| m.get(&key).cloned())
    {
        return Some(found);
    }
    let built = Arc::new(build(space)?);
    let mut guard = MEMO.lock().unwrap_or_else(|p| p.into_inner());
    let map = guard.get_or_insert_with(HashMap::new);
    if map.len() > 64 {
        map.clear();
    }
    map.insert(key, built.clone());
    Some(built)
}

fn build(space: &SectorSpace) -> Option<HardSystem> {
    let m = space.m;
    let total = space.total();
    let per = space.per();
    let mut is_fixed = vec![false; total];
    for s in 0..space.n {
        for k in 0..m {
            is_fixed[space.index(s, m - 1, k)] = true;
            is_fixed[space.index(s, k, m - 1)] = true;
        }
    }
    let free: Vec<usize> = (0..total).filter(|&i| !is_fixed[i]).collect();
    let fixed: Vec<usize> = (0..total).filter(|&i| is_fixed[i]).collect();
    let mut slots = vec![(false, 0usize); total];
    for (c, &i) in free.iter().enumerate() {
        slots[i] = (false, c);
    }
    for (c, &i) in fixed.iter().enumerate() {
        slots[i] = (true, c);
    }
    let spoke_rows = spoke_rows(space, CharacteristicMap::new(space.n).c);
    let basis = nested_basis(space, &free, &slots, &spoke_rows)?;
    // Per sector: its free rows and the columns they touch.
    let mut sector_rows = vec![Vec::new(); space.n];
    for (r, &i) in free.iter().enumerate() {
        sector_rows[i / per].push(r);
    }
    let mut sector_cols = Vec::with_capacity(space.n);
    let mut sector_z = Vec::with_capacity(space.n);
    for rows in &sector_rows {
        let mut cols: Vec<usize> = rows
            .iter()
            .flat_map(|&r| {
                basis.rows[r]
                    .iter()
                    .flat_map(|&(start, len, _)| start..start + len)
            })
            .collect();
        cols.sort_unstable();
        cols.dedup();
        let mut z = Mat::<f64>::zeros(rows.len(), cols.len());
        for (a, &r) in rows.iter().enumerate() {
            for &(start, len, offset) in &basis.rows[r] {
                for t in 0..len {
                    if let Ok(c) = cols.binary_search(&(start + t)) {
                        z.write(a, c, basis.values[offset + t]);
                    }
                }
            }
        }
        sector_cols.push(cols);
        sector_z.push(z);
    }
    Some(HardSystem {
        free,
        fixed,
        slots,
        layout: basis.layout,
        rows: basis.rows,
        values: basis.values,
        sector_rows,
        sector_cols,
        sector_z,
        spoke_rows,
        particular: basis.particular,
    })
}

impl HardSystem {
    /// Null-space dimension (per coordinate).
    pub fn dim(&self) -> usize {
        self.layout.k
    }

    /// `a += w · Z[r, :]` for free row `r`.
    pub fn add_row(&self, r: usize, w: f64, a: &mut [f64]) {
        for &(start, len, offset) in &self.rows[r] {
            for (slot, &z) in a[start..start + len]
                .iter_mut()
                .zip(&self.values[offset..offset + len])
            {
                *slot += w * z;
            }
        }
    }

    /// `Z[r, :] · y` for free row `r`.
    pub fn dot_row(&self, r: usize, y: &[f64]) -> f64 {
        self.rows[r]
            .iter()
            .map(|&(start, len, offset)| {
                self.values[offset..offset + len]
                    .iter()
                    .zip(&y[start..start + len])
                    .map(|(z, v)| z * v)
                    .sum::<f64>()
            })
            .sum()
    }

    /// `Zᵀ d` for one coordinate's free-row values `d` (the basis has
    /// orthonormal columns, so this is the least-squares projection).
    pub fn project(&self, d: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.dim()];
        for (r, &v) in d.iter().enumerate() {
            self.add_row(r, v, &mut y);
        }
        y
    }

    /// Minimum-norm free coefficients satisfying `C_U x = −C_F f` for
    /// one coordinate's fixed values `f`, with the relative residual of
    /// that system (0 when consistent). A solution in the spokes'
    /// structure less its null-space part: the minimum-norm solution is
    /// the one orthogonal to the null space.
    pub fn particular(&self, fixed_values: &[f64]) -> (Vec<f64>, f64) {
        let nfree = self.free.len();
        let fixed_part = |row: &Functional| -> f64 {
            row.iter()
                .filter_map(|&(i, w)| match self.slots[i] {
                    (true, r) => Some(w * fixed_values[r]),
                    _ => None,
                })
                .sum()
        };
        let b: Vec<f64> = self.spoke_rows.iter().map(|row| -fixed_part(row)).collect();
        let x0 = self.particular.solve(&b, nfree);
        // The part of b the solution misses.
        let residual2: f64 = self
            .spoke_rows
            .iter()
            .map(|row| {
                let v: f64 = row
                    .iter()
                    .map(|&(i, w)| match self.slots[i] {
                        (false, r) => w * x0[r],
                        (true, r) => w * fixed_values[r],
                    })
                    .sum();
                v * v
            })
            .sum();
        let scale = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1.0);
        let y = self.project(&x0);
        let x = (0..nfree).map(|r| x0[r] - self.dot_row(r, &y)).collect();
        (x, residual2.sqrt() / scale)
    }
}
