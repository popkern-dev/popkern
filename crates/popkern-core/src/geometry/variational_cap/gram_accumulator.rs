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
        // target, PK-676). Column-major on both sides.
        let size = self.size;
        let w = &self.columns[..size * self.count];
        for c in 0..self.count {
            let column = &w[c * size..(c + 1) * size];
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
        self.columns[..size * self.count].fill(0.0);
        self.count = 0;
    }
}
