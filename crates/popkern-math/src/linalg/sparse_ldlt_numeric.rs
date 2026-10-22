//! The numeric phase of the supernodal `L D Lᵀ` (PK-600): left-looking
//! panel factorisation and the triangular solves, each in an order
//! fixed by the pattern.

use super::sparse_ldlt::{LdltFactor, LdltSymbolic};

impl LdltSymbolic {
    /// Factor numerically with `values` aligned to the input pattern.
    /// `None` when a pivot is zero or not finite, or `values` is short.
    pub fn factor(&self, values: &[f64]) -> Option<LdltFactor> {
        let supernodes = self.supernodes();
        let mut panels = vec![0.0f64; self.sn_val_ptr[supernodes]];
        let mut d = vec![0.0f64; self.n];
        let mut position = vec![0usize; self.n];
        let mut block: Vec<f64> = Vec::new();
        let mut scaled: Vec<f64> = Vec::new();
        let mut relative: Vec<usize> = Vec::new();
        for s in 0..supernodes {
            let (c0, c1) = (self.sn_col[s], self.sn_col[s + 1]);
            let width = c1 - c0;
            let rows = &self.sn_rows[self.sn_row_ptr[s]..self.sn_row_ptr[s + 1]];
            let m = rows.len();
            for (t, &r) in rows.iter().enumerate() {
                position[r] = t;
            }
            let base = self.sn_val_ptr[s];
            // Assemble this panel's columns of P A Pᵀ.
            for j in 0..width {
                let c = c0 + j;
                for p in self.lower_col_ptr[c]..self.lower_col_ptr[c + 1] {
                    let value = *values.get(self.lower_source[p])?;
                    panels[base + j * m + position[self.lower_row[p]]] += value;
                }
            }
            // Updates from earlier panels, ascending.
            for &(k, a, b) in &self.updates[self.update_ptr[s]..self.update_ptr[s + 1]] {
                let k_rows = &self.sn_rows[self.sn_row_ptr[k]..self.sn_row_ptr[k + 1]];
                let (k_m, k_width) = (k_rows.len(), self.sn_col[k + 1] - self.sn_col[k]);
                let (k_base, k_c0) = (self.sn_val_ptr[k], self.sn_col[k]);
                let tall = k_m - a;
                // block[u, t] = Σ_q L_k[a + u, q] d_q L_k[a + t, q], for
                // t < b − a, accumulated over q ascending.
                block.clear();
                block.resize(tall * (b - a), 0.0);
                for q in 0..k_width {
                    let column = &panels[k_base + q * k_m..k_base + (q + 1) * k_m];
                    let dq = d[k_c0 + q];
                    scaled.clear();
                    scaled.extend(column[a..b].iter().map(|l| l * dq));
                    for (t, &f) in scaled.iter().enumerate() {
                        let out = &mut block[t * tall..(t + 1) * tall];
                        for (o, &l) in out[t..].iter_mut().zip(&column[a + t..]) {
                            *o += l * f;
                        }
                    }
                }
                // Scatter into this panel. Rows of `k` that are
                // consecutive in this panel's rows subtract as a slice.
                relative.clear();
                relative.extend(k_rows[a..].iter().map(|&r| position[r]));
                let consecutive = relative.windows(2).all(|w| w[1] == w[0] + 1);
                for t in 0..b - a {
                    let target = base + (k_rows[a + t] - c0) * m;
                    let source = &block[t * tall + t..(t + 1) * tall];
                    if consecutive {
                        let start = target + relative[t];
                        for (o, &v) in panels[start..start + source.len()].iter_mut().zip(source) {
                            *o -= v;
                        }
                    } else {
                        for (&r, &v) in relative[t..].iter().zip(source) {
                            panels[target + r] -= v;
                        }
                    }
                }
            }
            // Factor the panel's own columns, left-looking within it.
            for j in 0..width {
                for q in 0..j {
                    let f = d[c0 + q] * panels[base + q * m + j];
                    let (left, right) = panels[base..base + m * width].split_at_mut(j * m);
                    let source = &left[q * m..(q + 1) * m];
                    for (o, &l) in right[j..m].iter_mut().zip(&source[j..]) {
                        *o -= l * f;
                    }
                }
                let pivot = panels[base + j * m + j];
                if pivot == 0.0 || !pivot.is_finite() {
                    return None;
                }
                d[c0 + j] = pivot;
                panels[base + j * m + j] = 1.0;
                for o in &mut panels[base + j * m + j + 1..base + (j + 1) * m] {
                    *o /= pivot;
                }
            }
        }
        Some(LdltFactor { panels, d })
    }

    /// Solve `A x = b` in place with `factor` (one right-hand side).
    pub fn solve(&self, factor: &LdltFactor, b: &mut [f64]) {
        let supernodes = self.supernodes();
        let mut z: Vec<f64> = self.perm.iter().map(|&i| b[i]).collect();
        let panel = |s: usize| {
            let rows = &self.sn_rows[self.sn_row_ptr[s]..self.sn_row_ptr[s + 1]];
            (rows, self.sn_val_ptr[s], self.sn_col[s], self.sn_col[s + 1])
        };
        for s in 0..supernodes {
            let (rows, base, c0, c1) = panel(s);
            let m = rows.len();
            for j in 0..c1 - c0 {
                let zj = z[c0 + j];
                let column = &factor.panels[base + j * m..base + (j + 1) * m];
                for (&r, &l) in rows[j + 1..].iter().zip(&column[j + 1..]) {
                    z[r] -= l * zj;
                }
            }
        }
        for (zj, dj) in z.iter_mut().zip(&factor.d) {
            *zj /= dj;
        }
        for s in (0..supernodes).rev() {
            let (rows, base, c0, c1) = panel(s);
            let m = rows.len();
            for j in (0..c1 - c0).rev() {
                let column = &factor.panels[base + j * m..base + (j + 1) * m];
                let mut acc = z[c0 + j];
                for (&r, &l) in rows[j + 1..].iter().zip(&column[j + 1..]) {
                    acc -= l * z[r];
                }
                z[c0 + j] = acc;
            }
        }
        for (k, zk) in z.into_iter().enumerate() {
            b[self.perm[k]] = zk;
        }
    }
}
