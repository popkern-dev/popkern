//! A sparse `L D Lᵀ` factorisation (the up-looking algorithm of Davis'
//! LDL package) for symmetric quasi-definite systems, in one fixed
//! operation order (PK-676).
//!
//! A quasi-definite matrix `[[H, Aᵀ], [A, −G]]` (`H`, `G` positive
//! definite) has an `L D Lᵀ` factorisation for every symmetric
//! permutation, so no numeric pivoting is needed and the elimination
//! order is the caller's fill-reducing ordering alone. The symbolic
//! analysis (elimination tree, column counts) depends only on the
//! pattern and the ordering, so a caller factoring many systems with
//! one pattern analyses it once.

/// The pattern-dependent part: the permuted upper triangle's structure
/// and the elimination tree.
#[derive(Debug, Clone)]
pub struct LdltSymbolic {
    n: usize,
    perm: Vec<usize>,
    perm_inv: Vec<usize>,
    parent: Vec<usize>,
    lp: Vec<usize>,
    up_col_ptr: Vec<usize>,
    up_row: Vec<usize>,
    /// The input entry each upper entry takes its value from.
    up_source: Vec<usize>,
}

/// The numeric factor `P A Pᵀ = L D Lᵀ`.
#[derive(Debug, Clone)]
pub struct LdltFactor {
    li: Vec<usize>,
    lx: Vec<f64>,
    d: Vec<f64>,
}

const NONE: usize = usize::MAX;

impl LdltSymbolic {
    /// Analyse the `n × n` symmetric matrix given by its lower triangle
    /// in compressed columns (`row_idx[p] >= column`, each position
    /// once), eliminated in the order `perm` (`perm[k]` is the row of
    /// the input eliminated `k`-th). `None` for a malformed pattern or
    /// permutation.
    pub fn new(
        n: usize,
        col_ptr: &[usize],
        row_idx: &[usize],
        perm: Vec<usize>,
    ) -> Option<LdltSymbolic> {
        if col_ptr.len() != n + 1 || perm.len() != n {
            return None;
        }
        let mut perm_inv = vec![NONE; n];
        for (k, &i) in perm.iter().enumerate() {
            if i >= n || perm_inv[i] != NONE {
                return None;
            }
            perm_inv[i] = k;
        }
        // The permuted upper triangle, by column.
        let mut counts = vec![0usize; n + 1];
        let mut entries = Vec::with_capacity(row_idx.len());
        for c in 0..n {
            for p in col_ptr[c]..col_ptr[c + 1] {
                let r = *row_idx.get(p)?;
                if r < c || r >= n {
                    return None;
                }
                let (a, b) = (perm_inv[r], perm_inv[c]);
                let (row, col) = (a.min(b), a.max(b));
                counts[col + 1] += 1;
                entries.push((col, row, p));
            }
        }
        for k in 0..n {
            counts[k + 1] += counts[k];
        }
        let mut fill = counts.clone();
        let mut up_row = vec![0usize; entries.len()];
        let mut up_source = vec![0usize; entries.len()];
        for &(col, row, p) in &entries {
            up_row[fill[col]] = row;
            up_source[fill[col]] = p;
            fill[col] += 1;
        }
        // Elimination tree and column counts (LDL symbolic).
        let mut parent = vec![NONE; n];
        let mut flag = vec![NONE; n];
        let mut lnz = vec![0usize; n];
        for k in 0..n {
            flag[k] = k;
            for &row in &up_row[counts[k]..counts[k + 1]] {
                let mut i = row;
                while i < k && flag[i] != k {
                    if parent[i] == NONE {
                        parent[i] = k;
                    }
                    lnz[i] += 1;
                    flag[i] = k;
                    i = parent[i];
                }
            }
        }
        let mut lp = vec![0usize; n + 1];
        for k in 0..n {
            lp[k + 1] = lp[k] + lnz[k];
        }
        Some(LdltSymbolic {
            n,
            perm,
            perm_inv,
            parent,
            lp,
            up_col_ptr: counts,
            up_row,
            up_source,
        })
    }

    /// The order size.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Factor numerically with `values` aligned to the input pattern.
    /// `None` when a pivot is zero or not finite.
    pub fn factor(&self, values: &[f64]) -> Option<LdltFactor> {
        let n = self.n;
        let nnz = self.lp[n];
        let mut li = vec![0usize; nnz];
        let mut lx = vec![0.0f64; nnz];
        let mut d = vec![0.0f64; n];
        let mut y = vec![0.0f64; n];
        let mut pattern = vec![0usize; n];
        let mut flag = vec![NONE; n];
        let mut lnz = vec![0usize; n];
        for k in 0..n {
            let mut top = n;
            flag[k] = k;
            for p in self.up_col_ptr[k]..self.up_col_ptr[k + 1] {
                let mut i = self.up_row[p];
                y[i] += *values.get(self.up_source[p])?;
                let mut len = 0;
                while i < n && flag[i] != k {
                    pattern[len] = i;
                    len += 1;
                    flag[i] = k;
                    i = self.parent[i];
                }
                while len > 0 {
                    top -= 1;
                    len -= 1;
                    pattern[top] = pattern[len];
                }
            }
            d[k] = y[k];
            y[k] = 0.0;
            for &i in &pattern[top..n] {
                let yi = y[i];
                y[i] = 0.0;
                let end = self.lp[i] + lnz[i];
                for p in self.lp[i]..end {
                    y[li[p]] -= lx[p] * yi;
                }
                let l_ki = yi / d[i];
                d[k] -= l_ki * yi;
                li[end] = k;
                lx[end] = l_ki;
                lnz[i] += 1;
            }
            if d[k] == 0.0 || !d[k].is_finite() {
                return None;
            }
        }
        Some(LdltFactor { li, lx, d })
    }

    /// Solve `A x = b` in place with `factor` (one right-hand side).
    pub fn solve(&self, factor: &LdltFactor, b: &mut [f64]) {
        let n = self.n;
        let mut z: Vec<f64> = (0..n).map(|k| b[self.perm[k]]).collect();
        for j in 0..n {
            let zj = z[j];
            for p in self.lp[j]..self.lp[j + 1] {
                z[factor.li[p]] -= factor.lx[p] * zj;
            }
        }
        for (zj, dj) in z.iter_mut().zip(&factor.d) {
            *zj /= dj;
        }
        for j in (0..n).rev() {
            let mut s = z[j];
            for p in self.lp[j]..self.lp[j + 1] {
                s -= factor.lx[p] * z[factor.li[p]];
            }
            z[j] = s;
        }
        for (k, zk) in z.into_iter().enumerate() {
            b[self.perm[k]] = zk;
        }
        debug_assert!(self.perm_inv.len() == n);
    }
}
