//! A sparse supernodal `L D Lᵀ` factorisation for symmetric
//! quasi-definite systems, in one fixed operation order (PK-676,
//! PK-600).
//!
//! A quasi-definite matrix `[[H, Aᵀ], [A, −G]]` (`H`, `G` positive
//! definite) has an `L D Lᵀ` factorisation for every symmetric
//! permutation, so no numeric pivoting is needed and the elimination
//! order is the caller's fill-reducing ordering alone. The symbolic
//! analysis (elimination tree, column structures, supernodes) depends
//! only on the pattern and the ordering, so a caller factoring many
//! systems with one pattern analyses it once.
//!
//! Columns whose factor structure nests (a fundamental supernode: a
//! chain in the elimination tree whose column counts fall by one per
//! column) are factored together as one dense panel, left-looking:
//! each panel gathers the updates of the earlier panels that reach it
//! in ascending order, then factors its own diagonal block. The KP20
//! curvature system is per-piece dense Grams coupled by seams, so most
//! of its work lies in panels 8 to 300 columns wide, where dense loops
//! run several times faster than the column-at-a-time up-looking
//! elimination (PK-600: n = 8 cap 0.75 s back under its 0.5 s budget).
//! Every sum runs in an order fixed by the pattern, with no fused
//! multiply-add, so a factor depends only on its inputs.

/// The pattern-dependent part: the ordering, the supernode partition,
/// each panel's row structure and the panel-to-panel update lists.
#[derive(Debug, Clone)]
pub struct LdltSymbolic {
    pub(super) n: usize,
    pub(super) perm: Vec<usize>,
    /// The permuted lower triangle by column: `lower_row[p]` with the
    /// input entry `lower_source[p]` it takes its value from.
    pub(super) lower_col_ptr: Vec<usize>,
    pub(super) lower_row: Vec<usize>,
    pub(super) lower_source: Vec<usize>,
    /// Supernode `s` holds columns `sn_col[s]..sn_col[s + 1]`.
    pub(super) sn_col: Vec<usize>,
    /// Its panel rows (its own columns, then the rows below, ascending)
    /// are `sn_rows[sn_row_ptr[s]..sn_row_ptr[s + 1]]`; its dense
    /// column-major panel starts at `sn_val_ptr[s]`.
    pub(super) sn_row_ptr: Vec<usize>,
    pub(super) sn_rows: Vec<usize>,
    pub(super) sn_val_ptr: Vec<usize>,
    /// The earlier panels updating panel `s`:
    /// `updates[update_ptr[s]..update_ptr[s + 1]]`, each `(k, a, b)`
    /// meaning rows `a..b` of panel `k` fall in `s`'s columns (and rows
    /// `a..` of `k` all lie in `s`'s structure). Ascending in `k`.
    pub(super) update_ptr: Vec<usize>,
    pub(super) updates: Vec<(usize, usize, usize)>,
}

/// The numeric factor `P A Pᵀ = L D Lᵀ`: dense panels plus `D`.
#[derive(Debug, Clone)]
pub struct LdltFactor {
    pub(super) panels: Vec<f64>,
    pub(super) d: Vec<f64>,
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
        // The permuted entries as (column, row) of the lower triangle.
        let mut entries = Vec::with_capacity(row_idx.len());
        for c in 0..n {
            for p in col_ptr[c]..col_ptr[c + 1] {
                let r = *row_idx.get(p)?;
                if r < c || r >= n {
                    return None;
                }
                let (a, b) = (perm_inv[r], perm_inv[c]);
                entries.push((a.min(b), a.max(b), p));
            }
        }
        let (up_ptr, up_row, _) = compress(n, entries.iter().map(|&(c, r, p)| (r, c, p)));
        let (lower_col_ptr, lower_row, lower_source) = compress(n, entries.into_iter());
        // Elimination tree and each column's below-diagonal structure,
        // ascending (row `k` reaches column `i` through the row subtree).
        let mut parent = vec![NONE; n];
        let mut flag = vec![NONE; n];
        let mut col_rows: Vec<Vec<usize>> = vec![Vec::new(); n];
        for k in 0..n {
            flag[k] = k;
            for &row in &up_row[up_ptr[k]..up_ptr[k + 1]] {
                let mut i = row;
                while i < k && flag[i] != k {
                    if parent[i] == NONE {
                        parent[i] = k;
                    }
                    col_rows[i].push(k);
                    flag[i] = k;
                    i = parent[i];
                }
            }
        }
        // Fundamental supernodes: column k − 1 joins k when k is its
        // parent and only child and the structures nest.
        let mut children = vec![0usize; n];
        for &p in &parent {
            if p != NONE {
                children[p] += 1;
            }
        }
        let mut sn_col = vec![0usize];
        for k in 1..n {
            let nests = parent[k - 1] == k
                && children[k] == 1
                && col_rows[k - 1].len() == col_rows[k].len() + 1;
            if !nests {
                sn_col.push(k);
            }
        }
        if n > 0 {
            sn_col.push(n);
        }
        let supernodes = sn_col.len().saturating_sub(1);
        let mut sn_of = vec![0usize; n];
        let (mut sn_row_ptr, mut sn_rows, mut sn_val_ptr) =
            (vec![0usize], Vec::new(), vec![0usize]);
        for s in 0..supernodes {
            let (c0, c1) = (sn_col[s], sn_col[s + 1]);
            sn_of[c0..c1].fill(s);
            sn_rows.extend(c0..c1);
            sn_rows.extend_from_slice(&col_rows[c1 - 1]);
            sn_row_ptr.push(sn_rows.len());
            let m = sn_rows.len() - sn_row_ptr[s];
            sn_val_ptr.push(sn_val_ptr[s] + m * (c1 - c0));
        }
        // Update lists: the rows below each panel, grouped by the
        // supernode owning them.
        let mut lists: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); supernodes];
        for k in 0..supernodes {
            let rows = &sn_rows[sn_row_ptr[k]..sn_row_ptr[k + 1]];
            let mut a = sn_col[k + 1] - sn_col[k];
            while a < rows.len() {
                let target = sn_of[rows[a]];
                let mut b = a + 1;
                while b < rows.len() && sn_of[rows[b]] == target {
                    b += 1;
                }
                lists[target].push((k, a, b));
                a = b;
            }
        }
        let mut update_ptr = vec![0usize];
        let mut updates = Vec::new();
        for list in lists {
            updates.extend(list);
            update_ptr.push(updates.len());
        }
        Some(LdltSymbolic {
            n,
            perm,
            lower_col_ptr,
            lower_row,
            lower_source,
            sn_col,
            sn_row_ptr,
            sn_rows,
            sn_val_ptr,
            update_ptr,
            updates,
        })
    }

    /// The order size.
    pub fn n(&self) -> usize {
        self.n
    }

    /// The number of supernodal panels (diagnostic: one per column
    /// means no dense panel formed).
    pub fn supernodes(&self) -> usize {
        self.sn_col.len().saturating_sub(1)
    }
}

/// Compressed columns from `(column, row, source)` triples: rows
/// ascending within each column (sources in input order on ties).
fn compress(
    n: usize,
    triples: impl Iterator<Item = (usize, usize, usize)>,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut triples: Vec<(usize, usize, usize)> = triples.collect();
    triples.sort_by_key(|&(c, r, _)| (c, r));
    let mut ptr = vec![0usize; n + 1];
    for &(c, _, _) in &triples {
        ptr[c + 1] += 1;
    }
    for k in 0..n {
        ptr[k + 1] += ptr[k];
    }
    let rows = triples.iter().map(|&(_, r, _)| r).collect();
    let sources = triples.iter().map(|&(_, _, p)| p).collect();
    (ptr, rows, sources)
}
