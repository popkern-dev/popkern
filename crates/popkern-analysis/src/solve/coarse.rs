//! The coarsest-level solve, factored per connected block.
//!
//! Solids that do not touch give a block-diagonal operator, and
//! coarsening stops at a few dofs per block, so a part made of hundreds
//! of separate bodies keeps thousands of coarsest dofs. Each connected
//! block (by the operator's sparsity, found in index order) is factored
//! densely on its own: the cost is the sum of the blocks' cubes and the
//! dense cap applies to the largest block, not the whole level.

use crate::linalg::sparse::Csr;
use crate::linalg::{Cholesky, CholeskyError};

/// Dense factors of the connected blocks of the coarsest operator.
#[derive(Debug, Clone)]
pub struct BlockCholesky {
    /// Each block's dofs (ascending) and its factor.
    blocks: Vec<(Vec<usize>, Cholesky)>,
    /// Size of the operator.
    order: usize,
}

/// The connected blocks of `a`'s sparsity, each ascending, in order of
/// their lowest dof.
pub fn blocks(a: &Csr) -> Vec<Vec<usize>> {
    let n = a.rows();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut block = vec![start];
        let mut k = 0;
        while k < block.len() {
            let r = block[k];
            for (c, v) in a.row(r) {
                if v != 0.0 && !seen[c] {
                    seen[c] = true;
                    block.push(c);
                }
            }
            k += 1;
        }
        block.sort_unstable();
        out.push(block);
    }
    out
}

impl BlockCholesky {
    /// Factors every block of `a`; `Err(size)` when a block exceeds
    /// `cap` dofs.
    pub fn factor(a: &Csr, cap: usize) -> Result<Result<BlockCholesky, CholeskyError>, usize> {
        let groups = blocks(a);
        if let Some(largest) = groups.iter().map(Vec::len).max().filter(|&m| m > cap) {
            return Err(largest);
        }
        let mut position = vec![0usize; a.rows()];
        let mut factored = Vec::with_capacity(groups.len());
        for block in groups {
            let m = block.len();
            for (i, &r) in block.iter().enumerate() {
                position[r] = i;
            }
            let mut dense = vec![0.0; m * m];
            for (i, &r) in block.iter().enumerate() {
                for (c, v) in a.row(r) {
                    dense[i * m + position[c]] = v;
                }
            }
            match Cholesky::factor(&dense, m, 0.0) {
                Ok(f) => factored.push((block, f)),
                Err(e) => return Ok(Err(e)),
            }
        }
        Ok(Ok(BlockCholesky {
            blocks: factored,
            order: a.rows(),
        }))
    }

    /// Solves A x = b block by block.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, CholeskyError> {
        let mut x = vec![0.0; self.order];
        for (block, factor) in &self.blocks {
            let rhs: Vec<f64> = block
                .iter()
                .map(|&r| b.get(r).copied().unwrap_or(0.0))
                .collect();
            let y = factor.solve(&rhs)?;
            for (&r, v) in block.iter().zip(y) {
                x[r] = v;
            }
        }
        Ok(x)
    }

    /// Number of blocks.
    pub fn block_count(&self) -> usize {
        self.blocks.len()
    }

    /// Dofs of the largest block.
    pub fn largest_block(&self) -> usize {
        self.blocks.iter().map(|(b, _)| b.len()).max().unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two SPD blocks {0, 2} and {1, 3}, interleaved.
    fn interleaved() -> Csr {
        let t = [
            (0, 0, 4.0),
            (0, 2, 1.0),
            (2, 0, 1.0),
            (2, 2, 3.0),
            (1, 1, 5.0),
            (1, 3, -2.0),
            (3, 1, -2.0),
            (3, 3, 2.0),
        ];
        Csr::from_triplets(4, 4, &t).unwrap_or_else(|e| panic!("{e}"))
    }

    #[test]
    fn blocks_are_solved_like_the_whole() {
        let a = interleaved();
        assert_eq!(blocks(&a), [vec![0, 2], vec![1, 3]]);
        let f = BlockCholesky::factor(&a, 2)
            .unwrap_or_else(|m| panic!("{m}"))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((f.block_count(), f.largest_block()), (2, 2));
        let b = [1.0, -2.0, 0.5, 3.0];
        let x = f.solve(&b).unwrap_or_else(|e| panic!("{e}"));
        let whole = Cholesky::factor(&a.to_dense(), 4, 0.0).unwrap_or_else(|e| panic!("{e}"));
        let y = whole.solve(&b).unwrap_or_else(|e| panic!("{e}"));
        for (p, q) in x.iter().zip(&y) {
            assert!((p - q).abs() <= 1e-14, "{x:?} vs {y:?}");
        }
    }

    #[test]
    fn the_cap_applies_to_the_largest_block() {
        let a = interleaved();
        assert!(BlockCholesky::factor(&a, 2).is_ok());
        // Negative control: a cap below one block's size declines.
        assert_eq!(BlockCholesky::factor(&a, 1).err(), Some(2));
    }
}
