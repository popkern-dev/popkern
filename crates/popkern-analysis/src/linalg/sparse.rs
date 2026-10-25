//! Scalar sparse matrices in compressed sparse row form, with the
//! products multigrid needs (A·x, Aᵀ, A·B), all in a fixed order.

/// Why a sparse matrix could not be built or applied.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SparseError {
    /// A triplet lies outside the matrix.
    #[error("entry ({row}, {col}) is outside {rows}x{cols}")]
    OutOfRange {
        /// Row.
        row: usize,
        /// Column.
        col: usize,
        /// Rows of the matrix.
        rows: usize,
        /// Columns of the matrix.
        cols: usize,
    },
    /// A value is NaN or infinite.
    #[error("entry ({row}, {col}) is not finite")]
    NonFinite {
        /// Row.
        row: usize,
        /// Column.
        col: usize,
    },
    /// Operand shapes do not match.
    #[error("shape mismatch: {0}")]
    Shape(&'static str),
    /// A constructed matrix broke the CSR invariants (a column out of
    /// range or out of order): a defect in the constructor named.
    #[error("{constructor}: row {row} breaks the column invariant")]
    Invariant {
        /// The constructor.
        constructor: &'static str,
        /// The first row that breaks it.
        row: usize,
    },
}

/// Stored entries from which [`Csr::mul_vec`] uses threads: below it a
/// product takes well under a millisecond, about the cost of starting
/// them.
pub const PARALLEL_MIN_ENTRIES: usize = 1 << 20;

/// Threads for a large product: the cores available, at most 8 (the
/// product is memory-bound, and more threads gain nothing). 1 on wasm.
pub fn available_workers() -> usize {
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::thread::available_parallelism().map_or(1, |n| n.get().min(8))
    }
    #[cfg(target_arch = "wasm32")]
    {
        1
    }
}

/// A scalar CSR matrix. Columns within a row are ascending.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    rows: usize,
    cols: usize,
    row_start: Vec<usize>,
    col: Vec<usize>,
    val: Vec<f64>,
}

impl Csr {
    /// Builds a matrix from (row, col, value) triplets; duplicates are
    /// summed in input order.
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        triplets: &[(usize, usize, f64)],
    ) -> Result<Csr, SparseError> {
        for &(row, col, v) in triplets {
            if row >= rows || col >= cols {
                return Err(SparseError::OutOfRange {
                    row,
                    col,
                    rows,
                    cols,
                });
            }
            if !v.is_finite() {
                return Err(SparseError::NonFinite { row, col });
            }
        }
        // Bucket by row in input order (a counting sort), then sort each
        // row by column with a stable sort: the same order as a stable
        // sort by (row, col), so duplicates still sum in input order and
        // the result is bit-identical, in linear time for the rows.
        let mut start = vec![0usize; rows + 1];
        for &(row, _, _) in triplets {
            start[row + 1] += 1;
        }
        for r in 0..rows {
            start[r + 1] += start[r];
        }
        let mut next = start.clone();
        let mut order = vec![0usize; triplets.len()];
        for (i, &(row, _, _)) in triplets.iter().enumerate() {
            order[next[row]] = i;
            next[row] += 1;
        }
        for r in 0..rows {
            order[start[r]..start[r + 1]].sort_by_key(|&i| triplets[i].1);
        }
        let mut row_start = vec![0; rows + 1];
        let mut col = Vec::new();
        let mut val: Vec<f64> = Vec::new();
        let mut last = None;
        for i in order {
            let (r, c, v) = triplets[i];
            match val.last_mut() {
                Some(sum) if last == Some((r, c)) => *sum += v,
                _ => {
                    col.push(c);
                    val.push(v);
                    row_start[r + 1] += 1;
                    last = Some((r, c));
                }
            }
        }
        for r in 0..rows {
            row_start[r + 1] += row_start[r];
        }
        Ok(Csr {
            rows,
            cols,
            row_start,
            col,
            val,
        })
    }

    /// The block-diagonal matrix of `blocks`, in order: the same matrix
    /// as `from_triplets` of their entries offset by the blocks above
    /// (rows sorted by column, explicit zeros kept), without the sort.
    pub fn block_diagonal(blocks: Vec<Csr>) -> Result<Csr, SparseError> {
        let rows: usize = blocks.iter().map(|b| b.rows).sum();
        let cols: usize = blocks.iter().map(|b| b.cols).sum();
        let mut blocks = blocks;
        if blocks.len() == 1
            && let Some(only) = blocks.pop()
        {
            only.check_finite()?;
            return Ok(only);
        }
        let nnz: usize = blocks.iter().map(|b| b.val.len()).sum();
        let mut row_start = Vec::with_capacity(rows + 1);
        row_start.push(0);
        let (mut col, mut val) = (Vec::with_capacity(nnz), Vec::with_capacity(nnz));
        let (mut r0, mut c0) = (0, 0);
        for b in &blocks {
            b.check_finite().map_err(|e| match e {
                SparseError::NonFinite { row, col } => SparseError::NonFinite {
                    row: row + r0,
                    col: col + c0,
                },
                other => other,
            })?;
            let base = val.len();
            row_start.extend(b.row_start[1..].iter().map(|s| base + s));
            col.extend(b.col.iter().map(|c| c + c0));
            val.extend_from_slice(&b.val);
            r0 += b.rows;
            c0 += b.cols;
        }
        let out = Csr {
            rows,
            cols,
            row_start,
            col,
            val,
        };
        out.check_columns("block_diagonal")?;
        Ok(out)
    }

    /// `Invariant` at the first row whose columns are out of range or
    /// not strictly ascending, or whose row start runs past the entries
    /// (PH-129: a one-off column of 2⁶³ + 3.4e10 reached coarsening).
    pub(crate) fn check_columns(&self, constructor: &'static str) -> Result<(), SparseError> {
        let bad = |row| SparseError::Invariant { constructor, row };
        if self.row_start.len() != self.rows + 1 || self.row_start.last() != Some(&self.col.len()) {
            return Err(bad(0));
        }
        for r in 0..self.rows {
            let (lo, hi) = (self.row_start[r], self.row_start[r + 1]);
            let cols = self.col.get(lo..hi).ok_or_else(|| bad(r))?;
            if cols.iter().any(|&c| c >= self.cols) || cols.windows(2).any(|w| w[0] >= w[1]) {
                return Err(bad(r));
            }
        }
        Ok(())
    }

    /// `NonFinite` at the first stored entry that is not finite.
    fn check_finite(&self) -> Result<(), SparseError> {
        for r in 0..self.rows {
            for k in self.row_start[r]..self.row_start[r + 1] {
                if !self.val[k].is_finite() {
                    return Err(SparseError::NonFinite {
                        row: r,
                        col: self.col[k],
                    });
                }
            }
        }
        Ok(())
    }

    /// A zero matrix with the pattern of `elements` (each a list of
    /// node indices) on `nodes` nodes with `per` dofs each (dof =
    /// per·node + component): every pair of nodes sharing an element is
    /// coupled in every component.
    pub fn block_pattern(
        nodes: usize,
        per: usize,
        elements: &[[usize; 8]],
    ) -> Result<Csr, SparseError> {
        let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); nodes];
        for element in elements {
            for &a in element {
                let row = adjacent.get_mut(a).ok_or(SparseError::OutOfRange {
                    row: a,
                    col: 0,
                    rows: nodes,
                    cols: nodes,
                })?;
                row.extend_from_slice(element);
            }
        }
        let n = nodes * per;
        let mut row_start = Vec::with_capacity(n + 1);
        row_start.push(0);
        let mut col = Vec::new();
        for list in &mut adjacent {
            list.sort_unstable();
            list.dedup();
            for _ in 0..per {
                for &m in list.iter() {
                    col.extend((0..per).map(|j| per * m + j));
                }
                row_start.push(col.len());
            }
        }
        let val = vec![0.0; col.len()];
        Ok(Csr {
            rows: n,
            cols: n,
            row_start,
            col,
            val,
        })
    }

    /// Adds `v` to entry (`row`, `col`), which must be in the pattern.
    pub fn add_at(&mut self, row: usize, col: usize, v: f64) -> Result<(), SparseError> {
        if !v.is_finite() {
            return Err(SparseError::NonFinite { row, col });
        }
        let missing = SparseError::OutOfRange {
            row,
            col,
            rows: self.rows,
            cols: self.cols,
        };
        let (&lo, &hi) = (
            self.row_start.get(row).ok_or(missing.clone())?,
            self.row_start.get(row + 1).ok_or(missing.clone())?,
        );
        let k = self.col[lo..hi].binary_search(&col).map_err(|_| missing)?;
        self.val[lo + k] += v;
        Ok(())
    }

    /// Adds `block` at rows `row..row + 3` and columns `col..col + 3`,
    /// each of which must be stored (as `add_at` would add them one by
    /// one): one search per row instead of one per entry.
    pub fn add_block3(
        &mut self,
        row: usize,
        col: usize,
        block: &[[f64; 3]; 3],
    ) -> Result<(), SparseError> {
        for (i, values) in block.iter().enumerate() {
            let r = row + i;
            let (Some(&lo), Some(&hi)) = (self.row_start.get(r), self.row_start.get(r + 1)) else {
                return Err(SparseError::OutOfRange {
                    row: r,
                    col,
                    rows: self.rows,
                    cols: self.cols,
                });
            };
            let cols = &self.col[lo..hi];
            let k = match cols.binary_search(&col) {
                Ok(k) if cols.get(k + 2) == Some(&(col + 2)) => k,
                _ => {
                    for (j, v) in values.iter().enumerate() {
                        self.add_at(r, col + j, *v)?;
                    }
                    continue;
                }
            };
            for (j, v) in values.iter().enumerate() {
                if !v.is_finite() {
                    return Err(SparseError::NonFinite {
                        row: r,
                        col: col + j,
                    });
                }
                self.val[lo + k + j] += v;
            }
        }
        Ok(())
    }

    /// The index in the stored values of entry (`row`, `col`), when
    /// stored.
    pub fn position(&self, row: usize, col: usize) -> Option<usize> {
        let (&lo, &hi) = (self.row_start.get(row)?, self.row_start.get(row + 1)?);
        self.col[lo..hi].binary_search(&col).ok().map(|k| lo + k)
    }

    /// Whether `other` has the same shape, stored entries and zero
    /// entries (a graph read with its zeros, as coarsening reads it).
    pub fn same_graph(&self, other: &Csr) -> bool {
        self.rows == other.rows
            && self.cols == other.cols
            && self.row_start == other.row_start
            && self.col == other.col
            && self
                .val
                .iter()
                .zip(&other.val)
                .all(|(a, b)| (*a == 0.0) == (*b == 0.0))
    }

    /// The stored values, in row order (columns ascending), to refill
    /// in place on a fixed pattern.
    pub fn values_mut(&mut self) -> &mut [f64] {
        &mut self.val
    }

    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.val.len()
    }

    /// The stored entries of `row` as (column, value), ascending.
    pub fn row(&self, row: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let range = self.row_start[row]..self.row_start[row + 1];
        self.col[range.clone()]
            .iter()
            .copied()
            .zip(self.val[range].iter().copied())
    }

    /// y = A x. A matrix with at least [`PARALLEL_MIN_ENTRIES`] stored
    /// entries splits its rows over threads natively (never on wasm);
    /// each row sums in the same order either way, so y is identical.
    pub fn mul_vec(&self, x: &[f64], y: &mut [f64]) -> Result<(), SparseError> {
        let workers = if self.val.len() >= PARALLEL_MIN_ENTRIES {
            available_workers()
        } else {
            1
        };
        self.mul_vec_on(x, y, workers)
    }

    /// [`Csr::mul_vec`] on `workers` threads, each taking one fixed,
    /// contiguous block of rows.
    pub fn mul_vec_on(&self, x: &[f64], y: &mut [f64], workers: usize) -> Result<(), SparseError> {
        if x.len() != self.cols || y.len() != self.rows {
            return Err(SparseError::Shape("mul_vec"));
        }
        let rows_into = |first: usize, out: &mut [f64]| {
            for (k, slot) in out.iter_mut().enumerate() {
                *slot = self.row(first + k).map(|(c, v)| v * x[c]).sum();
            }
        };
        #[cfg(not(target_arch = "wasm32"))]
        if workers > 1 && self.rows > workers {
            let block = self.rows.div_ceil(workers);
            std::thread::scope(|scope| {
                for (i, out) in y.chunks_mut(block).enumerate() {
                    scope.spawn(move || rows_into(i * block, out));
                }
            });
            return Ok(());
        }
        let _ = workers;
        rows_into(0, y);
        Ok(())
    }

    /// The diagonal (zero where not stored).
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.rows.min(self.cols))
            .map(|r| self.row(r).find(|(c, _)| *c == r).map_or(0.0, |(_, v)| v))
            .collect()
    }

    /// Aᵀ.
    pub fn transpose(&self) -> Csr {
        // Counting sort by column; rows are visited in order, so each
        // transposed row comes out with ascending columns.
        let mut row_start = vec![0; self.cols + 1];
        for &c in &self.col {
            row_start[c + 1] += 1;
        }
        for c in 0..self.cols {
            row_start[c + 1] += row_start[c];
        }
        let mut next = row_start.clone();
        let mut col = vec![0; self.nnz()];
        let mut val = vec![0.0; self.nnz()];
        for r in 0..self.rows {
            for (c, v) in self.row(r) {
                col[next[c]] = r;
                val[next[c]] = v;
                next[c] += 1;
            }
        }
        Csr {
            rows: self.cols,
            cols: self.rows,
            row_start,
            col,
            val,
        }
    }

    /// The dense row-major copy (for small coarse levels).
    pub fn to_dense(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.rows * self.cols];
        for r in 0..self.rows {
            for (c, v) in self.row(r) {
                d[r * self.cols + c] = v;
            }
        }
        d
    }
}

#[path = "sparse_mul.rs"]
mod mul;
