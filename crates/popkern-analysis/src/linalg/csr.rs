//! Sparse matrices of 3×3 blocks in compressed sparse row form, for
//! coarse multigrid levels of the vector (displacement) problem.
//!
//! Every reduction has a fixed order: duplicate triplets are summed in
//! the order they were given, and the product accumulates each block
//! row over its blocks in ascending column, each block in row-major
//! order. The result's bits depend only on the input.

/// A 3×3 block, row-major.
pub type Block = [[f64; 3]; 3];

/// One block entry for assembly: block row, block column, value.
pub type BlockTriplet = (usize, usize, Block);

/// Why a block matrix could not be built or applied.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CsrError {
    /// A triplet lies outside the matrix.
    #[error("triplet {index} at block ({row}, {col}) is outside {rows}x{cols} blocks")]
    OutOfRange {
        /// Position of the triplet in the input.
        index: usize,
        /// Its block row.
        row: usize,
        /// Its block column.
        col: usize,
        /// Block rows of the matrix.
        rows: usize,
        /// Block columns of the matrix.
        cols: usize,
    },
    /// A triplet holds a NaN or infinite value.
    #[error("triplet {0} holds a non-finite value")]
    NonFinite(usize),
    /// A vector has the wrong length.
    #[error("vector has {len} entries, expected {expected}")]
    VectorLength {
        /// Entries given.
        len: usize,
        /// Entries required (3 per block).
        expected: usize,
    },
}

/// A sparse matrix of 3×3 blocks.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockCsr {
    rows: usize,
    cols: usize,
    row_start: Vec<usize>,
    col_index: Vec<usize>,
    blocks: Vec<Block>,
}

impl BlockCsr {
    /// Builds a `rows`×`cols` block matrix from triplets. Duplicates
    /// are summed in input order; block rows with no triplets are
    /// empty.
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        triplets: &[BlockTriplet],
    ) -> Result<BlockCsr, CsrError> {
        for (index, &(row, col, ref block)) in triplets.iter().enumerate() {
            if row >= rows || col >= cols {
                return Err(CsrError::OutOfRange {
                    index,
                    row,
                    col,
                    rows,
                    cols,
                });
            }
            if block.iter().flatten().any(|v| !v.is_finite()) {
                return Err(CsrError::NonFinite(index));
            }
        }
        // Stable: equal (row, col) keys keep their input order.
        let mut order: Vec<usize> = (0..triplets.len()).collect();
        order.sort_by_key(|&i| (triplets[i].0, triplets[i].1));

        let mut row_start = vec![0; rows + 1];
        let mut col_index: Vec<usize> = Vec::new();
        let mut blocks: Vec<Block> = Vec::new();
        let mut last: Option<(usize, usize)> = None;
        for &i in &order {
            let (row, col, ref block) = triplets[i];
            match blocks.last_mut() {
                Some(sum) if last == Some((row, col)) => {
                    for a in 0..3 {
                        for b in 0..3 {
                            sum[a][b] += block[a][b];
                        }
                    }
                }
                _ => {
                    col_index.push(col);
                    blocks.push(*block);
                    row_start[row + 1] += 1;
                    last = Some((row, col));
                }
            }
        }
        for r in 0..rows {
            row_start[r + 1] += row_start[r];
        }
        Ok(BlockCsr {
            rows,
            cols,
            row_start,
            col_index,
            blocks,
        })
    }

    /// Number of block rows.
    pub fn block_rows(&self) -> usize {
        self.rows
    }

    /// Number of block columns.
    pub fn block_cols(&self) -> usize {
        self.cols
    }

    /// Number of stored blocks.
    pub fn stored_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// The stored blocks of block row `row` as (column, block), in
    /// ascending column order. Empty for rows past the end.
    pub fn row(&self, row: usize) -> impl Iterator<Item = (usize, &Block)> {
        let range = match (self.row_start.get(row), self.row_start.get(row + 1)) {
            (Some(&s), Some(&e)) => s..e,
            _ => 0..0,
        };
        self.col_index[range.clone()]
            .iter()
            .copied()
            .zip(&self.blocks[range])
    }

    /// y = A·x, overwriting `y`.
    pub fn mul_vec(&self, x: &[f64], y: &mut [f64]) -> Result<(), CsrError> {
        if x.len() != 3 * self.cols {
            return Err(CsrError::VectorLength {
                len: x.len(),
                expected: 3 * self.cols,
            });
        }
        if y.len() != 3 * self.rows {
            return Err(CsrError::VectorLength {
                len: y.len(),
                expected: 3 * self.rows,
            });
        }
        for r in 0..self.rows {
            let mut acc = [0.0; 3];
            for (c, block) in self.row(r) {
                let xc = &x[3 * c..3 * c + 3];
                for a in 0..3 {
                    for b in 0..3 {
                        acc[a] += block[a][b] * xc[b];
                    }
                }
            }
            y[3 * r..3 * r + 3].copy_from_slice(&acc);
        }
        Ok(())
    }
}
