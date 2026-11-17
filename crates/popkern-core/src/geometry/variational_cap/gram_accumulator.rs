//! A Gram matrix `Σ w · row rowᵀ` accumulated in batches (PK-620).
//!
//! The quadrature adds a handful of weighted rank-one terms per node,
//! each over every coefficient the node's basis supports. Adding them
//! one at a time is a scalar triple loop; stacking the scaled rows as
//! the columns of `W` and adding `W Wᵀ` is the same sum as one matrix
//! product. Measured on the octagon h1.2 variable apex's 16-sector
//! degree-9 fit: 0.7 s per energy refresh, 17% of the fit, with the
//! scalar loop.

/// Columns held before they are folded into the Gram.
const BATCH: usize = 256;

/// The accumulator for one `size × size` block.
pub(super) struct GramAccumulator {
    size: usize,
    columns: Vec<f64>,
    count: usize,
}

impl GramAccumulator {
    /// An empty accumulator for a `size × size` Gram.
    pub fn new(size: usize) -> Self {
        GramAccumulator {
            size,
            columns: vec![0.0; size * BATCH],
            count: 0,
        }
    }

    /// The next column, zeroed, for the caller to fill with `√w · row`;
    /// full batches are folded into `block` first.
    pub fn column(&mut self, block: &mut [f64]) -> &mut [f64] {
        if self.count == BATCH {
            self.flush(block);
        }
        let at = self.count * self.size;
        self.count += 1;
        &mut self.columns[at..at + self.size]
    }

    /// Add `√w · row` as a column, `row` sparse over the block's indices.
    pub fn push_sparse(&mut self, block: &mut [f64], row: &[(usize, f64)], w: f64) {
        let s = w.sqrt();
        let column = self.column(block);
        for &(i, a) in row {
            column[i] += s * a;
        }
    }

    /// Fold the held columns into `block` (row-major; the Gram is
    /// symmetric, so the layout does not matter).
    pub fn flush(&mut self, block: &mut [f64]) {
        if self.count == 0 {
            return;
        }
        // block += W Wᵀ, each entry summed over the columns in order
        // (one multiply and one add per term: the same bits on every
        // target, PK-676), a column whose `j` entry is zero adding
        // nothing to row `j`. Column-major on both sides. Each strip of
        // a block row stays in registers while every column passes
        // through it: sweeping the whole block once per column ran out
        // of cache (PK-1223); the terms and their order are unchanged.
        let size = self.size;
        let count = self.count;
        let w = &self.columns[..size * count];
        // Columns are taken four at a time where all four have a row's
        // entry, else in pairs: one pass over a block row adds every
        // column's term, each entry the columns' terms in their order
        // (as one column at a time), so the row is loaded and stored a
        // quarter as often (PK-1223). A zero entry adds nothing.
        let one = |target: &mut [f64], column: &[f64], wj: f64| {
            for (t, &wi) in target.iter_mut().zip(column) {
                *t += wi * wj;
            }
        };
        let two =
            |target: &mut [f64], a: &[f64], b: &[f64], j: usize| match (a[j] != 0.0, b[j] != 0.0) {
                (false, false) => {}
                (false, true) => one(target, b, b[j]),
                (true, false) => one(target, a, a[j]),
                (true, true) => {
                    let (wa, wb) = (a[j], b[j]);
                    for ((t, &x), &y) in target.iter_mut().zip(a).zip(b) {
                        *t = (*t + x * wa) + y * wb;
                    }
                }
            };
        let mut quads = w.chunks_exact(4 * size);
        for quad in &mut quads {
            let (a, rest) = quad.split_at(size);
            let (b, rest) = rest.split_at(size);
            let (c, d) = rest.split_at(size);
            for (j, target) in block.chunks_exact_mut(size).enumerate() {
                let (wa, wb, wc, wd) = (a[j], b[j], c[j], d[j]);
                if wa != 0.0 && wb != 0.0 && wc != 0.0 && wd != 0.0 {
                    for ((((t, &x), &y), &z), &u) in target.iter_mut().zip(a).zip(b).zip(c).zip(d) {
                        *t = (((*t + x * wa) + y * wb) + z * wc) + u * wd;
                    }
                } else {
                    two(target, a, b, j);
                    two(target, c, d, j);
                }
            }
        }
        for column in quads.remainder().chunks_exact(size) {
            for (j, target) in block.chunks_exact_mut(size).enumerate() {
                if column[j] != 0.0 {
                    one(target, column, column[j]);
                }
            }
        }
        self.columns[..size * self.count].fill(0.0);
        self.count = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::PortableMath;

    /// The per-column sweep `flush` ran before columns were grouped.
    fn flush_by_columns(columns: &[f64], size: usize, block: &mut [f64]) {
        for column in columns.chunks_exact(size) {
            for (j, &wj) in column.iter().enumerate() {
                if wj == 0.0 {
                    continue;
                }
                let target = &mut block[j * size..(j + 1) * size];
                for (t, &wi) in target.iter_mut().zip(column) {
                    *t += wi * wj;
                }
            }
        }
    }

    #[test]
    fn the_grouped_fold_is_the_column_sweep_bit_for_bit() {
        for (size, count) in [(81, 256), (13, 300), (121, 37), (8, 9), (3, 1)] {
            let value = |c: usize, i: usize| {
                if (c * 3 + i * 5).is_multiple_of(4) {
                    0.0
                } else {
                    ((c * 31 + i * 17) as f64 * 0.29).psin() * 10f64.powi((i % 5) as i32 - 2)
                }
            };
            let mut acc = GramAccumulator::new(size);
            let mut block = vec![0.0; size * size];
            let mut held = Vec::new();
            let mut reference = vec![0.0; size * size];
            for c in 0..count {
                // A full batch folds before the next column, as in use.
                if held.len() == size * BATCH {
                    flush_by_columns(&held, size, &mut reference);
                    held.clear();
                }
                let column = acc.column(&mut block);
                for (i, slot) in column.iter_mut().enumerate() {
                    *slot = value(c, i);
                }
                held.extend((0..size).map(|i| value(c, i)));
            }
            acc.flush(&mut block);
            flush_by_columns(&held, size, &mut reference);
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(
                bits(&block),
                bits(&reference),
                "size {size}, {count} columns"
            );
        }
    }
}
