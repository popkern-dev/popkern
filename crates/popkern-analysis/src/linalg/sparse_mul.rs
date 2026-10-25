//! Sparse products A B, threaded by blocks of rows.

use super::{Csr, PARALLEL_MIN_ENTRIES, SparseError, available_workers};

impl Csr {
    /// A B, row by row: each output row accumulates A's entries in
    /// column order, each scaled row of B in column order. A product of
    /// at least [`PARALLEL_MIN_ENTRIES`] multiply-adds (estimated) splits
    /// A's rows over threads natively; each row sums in the same order
    /// either way, so the product is identical.
    pub fn mul(&self, b: &Csr) -> Result<Csr, SparseError> {
        let work = self.nnz() as f64 * b.nnz() as f64 / b.rows.max(1) as f64;
        let workers = if work >= PARALLEL_MIN_ENTRIES as f64 {
            available_workers()
        } else {
            1
        };
        self.mul_on(b, workers)
    }

    /// [`Csr::mul`] on `workers` threads, each taking one fixed,
    /// contiguous block of A's rows.
    pub fn mul_on(&self, b: &Csr, workers: usize) -> Result<Csr, SparseError> {
        if self.cols != b.rows {
            return Err(SparseError::Shape("mul"));
        }
        let parts = self.mul_blocks(b, workers);
        let mut row_start = Vec::with_capacity(self.rows + 1);
        row_start.push(0);
        let nnz = parts.iter().map(|p| p.1.len()).sum();
        let (mut col, mut val) = (Vec::with_capacity(nnz), Vec::with_capacity(nnz));
        for (lengths, c, v) in parts {
            for l in lengths {
                row_start.push(col.len() + l);
            }
            col.extend(c);
            val.extend(v);
        }
        let out = Csr {
            rows: self.rows,
            cols: b.cols,
            row_start,
            col,
            val,
        };
        out.check_columns("mul")?;
        Ok(out)
    }

    /// A B in `workers` contiguous blocks of A's rows, in order.
    #[allow(clippy::type_complexity)]
    fn mul_blocks(&self, b: &Csr, workers: usize) -> Vec<(Vec<usize>, Vec<usize>, Vec<f64>)> {
        #[cfg(not(target_arch = "wasm32"))]
        if workers > 1 && self.rows > workers {
            let block = self.rows.div_ceil(workers);
            return std::thread::scope(|scope| {
                let handles: Vec<_> = (0..self.rows)
                    .step_by(block)
                    .map(|first| {
                        let rows = first..(first + block).min(self.rows);
                        scope.spawn(move || self.mul_rows(b, rows))
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            });
        }
        let _ = workers;
        vec![self.mul_rows(b, 0..self.rows)]
    }

    /// Rows `rows` of A B: each row's end within the block, its columns
    /// and its values.
    fn mul_rows(
        &self,
        b: &Csr,
        rows: std::ops::Range<usize>,
    ) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let mut acc = vec![0.0; b.cols];
        let mut seen = vec![false; b.cols];
        let mut touched = Vec::new();
        let mut ends = Vec::with_capacity(rows.len());
        let mut col = Vec::new();
        let mut val = Vec::new();
        for r in rows {
            for (k, a) in self.row(r) {
                for (c, v) in b.row(k) {
                    if !seen[c] {
                        seen[c] = true;
                        touched.push(c);
                    }
                    acc[c] += a * v;
                }
            }
            touched.sort_unstable();
            for &c in &touched {
                col.push(c);
                val.push(acc[c]);
                acc[c] = 0.0;
                seen[c] = false;
            }
            touched.clear();
            ends.push(col.len());
        }
        (ends, col, val)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Csr, SparseError};

    /// PH-129: products and block-diagonal stacks pass the column
    /// check, and a matrix with a column out of range or out of order
    /// fails it, naming its constructor.
    #[test]
    fn the_column_invariant_is_checked() {
        let a = Csr::from_triplets(3, 4, &[(0, 1, 1.0), (0, 3, 2.0), (2, 0, 3.0)])
            .unwrap_or_else(|e| panic!("{e}"));
        let b = Csr::from_triplets(4, 2, &[(1, 0, 1.0), (3, 1, 1.0), (0, 1, 2.0)])
            .unwrap_or_else(|e| panic!("{e}"));
        let p = a.mul(&b).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(p.check_columns("mul"), Ok(()));
        let d = Csr::block_diagonal(vec![a.clone(), b]).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(d.check_columns("block_diagonal"), Ok(()));
        let narrow = Csr {
            cols: 3,
            ..a.clone()
        };
        assert_eq!(
            narrow.check_columns("test"),
            Err(SparseError::Invariant {
                constructor: "test",
                row: 0
            })
        );
        let mut swapped = a;
        swapped.col = vec![3, 1, 0];
        assert!(swapped.check_columns("test").is_err());
    }
}
