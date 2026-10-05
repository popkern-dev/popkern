//! The sparse factorisation behind the KP20 bi-9 saddle solves
//! (PK-676): faer's approximate-minimum-degree ordering (integer work
//! only, the same on every machine) and `popkern_math`'s up-looking
//! `L D Lᵀ`, whose numeric phase runs in one fixed operation order.
//! faer's own supernodal numeric factorisation chose SIMD kernels by
//! runtime CPU detection, so the corner caps differed in their last
//! bits between machines.

use faer::dyn_stack::{GlobalPodBuffer, PodStack};
use faer::sparse::SymbolicSparseColMatRef;
use popkern_math::linalg::LdltSymbolic;

/// Order the `n × n` symmetric pattern given by its lower triangle
/// (compressed columns, rows sorted within each column) by AMD and
/// analyse it for `L D Lᵀ`. `None` for a malformed pattern.
pub(super) fn analyse(n: usize, col_ptr: &[usize], row_idx: &[usize]) -> Option<LdltSymbolic> {
    let pattern = SymbolicSparseColMatRef::new_checked(n, n, col_ptr, None, row_idx);
    let (mut perm, mut perm_inv) = (vec![0usize; n], vec![0usize; n]);
    let req = faer::sparse::linalg::amd::order_req::<usize>(n, row_idx.len()).ok()?;
    faer::sparse::linalg::amd::order(
        &mut perm,
        &mut perm_inv,
        pattern,
        faer::sparse::linalg::amd::Control::default(),
        PodStack::new(&mut GlobalPodBuffer::new(req)),
    )
    .ok()?;
    LdltSymbolic::new(n, col_ptr, row_idx, perm)
}

/// A symmetric system stored as its lower triangle in compressed
/// columns, with the positions of the diagonals regularisation shifts.
pub(super) struct LowerSystem {
    pub(super) size: usize,
    pub(super) col_ptr: Vec<usize>,
    pub(super) row_idx: Vec<usize>,
    pub(super) values: Vec<f64>,
}

impl LowerSystem {
    /// Collect `(row, col, value)` entries of a symmetric matrix (only
    /// those with `row >= col` are kept; repeated positions are summed
    /// in the order given). Every diagonal is present (zero if unset).
    pub(super) fn from_entries(size: usize, entries: &[(usize, usize, f64)]) -> LowerSystem {
        let mut lower: Vec<(usize, usize, f64)> = entries
            .iter()
            .copied()
            .filter(|&(r, c, _)| r >= c && r < size)
            .chain((0..size).map(|d| (d, d, 0.0)))
            .collect();
        lower.sort_by_key(|&(r, c, _)| (c, r));
        let mut col_ptr = vec![0usize; size + 1];
        let (mut row_idx, mut values): (Vec<usize>, Vec<f64>) = (Vec::new(), Vec::new());
        let mut last: Option<(usize, usize)> = None;
        for (r, c, v) in lower {
            if last == Some((r, c)) {
                if let Some(tail) = values.last_mut() {
                    *tail += v;
                }
                continue;
            }
            row_idx.push(r);
            values.push(v);
            col_ptr[c + 1] += 1;
            last = Some((r, c));
        }
        for c in 0..size {
            col_ptr[c + 1] += col_ptr[c];
        }
        LowerSystem {
            size,
            col_ptr,
            row_idx,
            values,
        }
    }

    /// The value index of diagonal `d` (each column starts with it:
    /// rows are sorted and the diagonal is the smallest row kept).
    pub(super) fn diagonal(&self, d: usize) -> usize {
        self.col_ptr[d]
    }

    /// `A x` for the symmetric matrix with these `values`.
    pub(super) fn product(&self, values: &[f64], x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.size];
        for col in 0..self.size {
            let span = self.col_ptr[col]..self.col_ptr[col + 1];
            for (&row, &w) in self.row_idx[span.clone()].iter().zip(&values[span]) {
                y[row] += w * x[col];
                if row != col {
                    y[col] += w * x[row];
                }
            }
        }
        y
    }
}
