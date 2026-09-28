//! The hard system: exact spoke rows over the free coefficients, with
//! the rim rows (every sector's `i = m−1` and `j = m−1` control rows)
//! fixed by the rim approximant.
//!
//! Per coordinate, `C_U x_U + C_F x_F = 0`. The null space `Z` of `C_U`
//! and the pseudo-inverse pieces for the particular solution depend
//! only on (n, degree, K), so they are computed once and memoised.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use faer::Mat;

use super::chart::CharacteristicMap;
use super::linalg::svd;
use super::space::SectorSpace;
use super::spoke::spoke_rows;

/// Relative singular-value threshold for the numerical rank.
const RANK_TOL: f64 = 1.0e-10;

/// The eliminated hard system for one space.
#[derive(Debug)]
pub(super) struct HardSystem {
    /// Global indices of the free coefficients (columns of `z`'s rows).
    pub free: Vec<usize>,
    /// Global indices of the fixed (rim) coefficients.
    pub fixed: Vec<usize>,
    /// Orthonormal null-space basis of `C_U`, `|free| × k`.
    pub z: Mat<f64>,
    /// Per global index: `(fixed?, row in free or fixed)`.
    pub slots: Vec<(bool, usize)>,
    c_fixed: Mat<f64>,
    u_r: Mat<f64>,
    s_r: Vec<f64>,
    v_r: Mat<f64>,
}

type Key = (usize, usize, Vec<u64>);

static MEMO: Mutex<Option<HashMap<Key, Arc<HardSystem>>>> = Mutex::new(None);

/// The memoised hard system for `space`; `None` when the SVD workspace
/// cannot be sized.
pub(super) fn hard_system(space: &SectorSpace) -> Option<Arc<HardSystem>> {
    let key: Key = (
        space.n,
        space.degree,
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
    let mut is_fixed = vec![false; total];
    for s in 0..space.n {
        for k in 0..m {
            is_fixed[space.index(s, m - 1, k)] = true;
            is_fixed[space.index(s, k, m - 1)] = true;
        }
    }
    let free: Vec<usize> = (0..total).filter(|&i| !is_fixed[i]).collect();
    let fixed: Vec<usize> = (0..total).filter(|&i| is_fixed[i]).collect();
    let mut col_of = vec![(false, 0usize); total];
    for (c, &i) in free.iter().enumerate() {
        col_of[i] = (false, c);
    }
    for (c, &i) in fixed.iter().enumerate() {
        col_of[i] = (true, c);
    }
    let rows = spoke_rows(space, CharacteristicMap::new(space.n).c);
    let mut c_free = Mat::<f64>::zeros(rows.len(), free.len());
    let mut c_fixed = Mat::<f64>::zeros(rows.len(), fixed.len());
    for (r, row) in rows.iter().enumerate() {
        for &(i, w) in row {
            match col_of[i] {
                (false, c) => c_free.write(r, c, c_free.read(r, c) + w),
                (true, c) => c_fixed.write(r, c, c_fixed.read(r, c) + w),
            }
        }
    }
    let d = svd(&c_free)?;
    let s_max = d.s.first().copied().unwrap_or(0.0);
    let rank = d.s.iter().filter(|&&s| s > RANK_TOL * s_max).count();
    let k = free.len() - rank;
    let z = Mat::from_fn(free.len(), k, |i, j| d.v.read(i, rank + j));
    let u_r = Mat::from_fn(rows.len(), rank, |i, j| d.u.read(i, j));
    let v_r = Mat::from_fn(free.len(), rank, |i, j| d.v.read(i, j));
    Some(HardSystem {
        free,
        fixed,
        z,
        slots: col_of,
        c_fixed,
        u_r,
        s_r: d.s[..rank].to_vec(),
        v_r,
    })
}

impl HardSystem {
    /// Null-space dimension.
    pub fn dim(&self) -> usize {
        self.z.ncols()
    }

    /// Minimum-norm free coefficients satisfying `C_U x = −C_F f` for
    /// one coordinate's fixed values `f`, with the relative residual of
    /// that system (0 when consistent).
    pub fn particular(&self, fixed_values: &[f64]) -> (Vec<f64>, f64) {
        let rows = self.c_fixed.nrows();
        let mut rhs = vec![0.0; rows];
        for (r, slot) in rhs.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (c, &f) in fixed_values.iter().enumerate() {
                acc -= self.c_fixed.read(r, c) * f;
            }
            *slot = acc;
        }
        let rank = self.s_r.len();
        let mut coef = vec![0.0; rank];
        for (j, slot) in coef.iter_mut().enumerate() {
            *slot = (0..rows).map(|r| self.u_r.read(r, j) * rhs[r]).sum();
        }
        // Residual: the part of rhs outside range(C_U).
        let mut residual2 = 0.0;
        for (r, &b) in rhs.iter().enumerate() {
            let proj: f64 = (0..rank).map(|j| self.u_r.read(r, j) * coef[j]).sum();
            residual2 += (b - proj) * (b - proj);
        }
        let scale = rhs.iter().map(|b| b * b).sum::<f64>().sqrt().max(1.0);
        let x = (0..self.free.len())
            .map(|i| {
                (0..rank)
                    .map(|j| self.v_r.read(i, j) * coef[j] / self.s_r[j])
                    .sum()
            })
            .collect();
        (x, residual2.sqrt() / scale)
    }
}
