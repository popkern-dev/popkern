//! Factorisations and triangular solves, each in one fixed operation
//! order.

use super::DenseMat;
use crate::PortableMath;

/// Lower Cholesky factor `L` of the symmetric positive definite `a`
/// (`a = L Lᵀ`; only the lower triangle of `a` is read). `None` when a
/// pivot is not positive or not finite.
pub fn cholesky(a: &DenseMat) -> Option<DenseMat> {
    let n = a.nrows();
    let mut l = DenseMat::zeros(n, n);
    for j in 0..n {
        let mut d = a.get(j, j);
        for k in 0..j {
            let ljk = l.get(j, k);
            d -= ljk * ljk;
        }
        if d.is_nan() || d <= 0.0 || d.is_infinite() {
            return None;
        }
        let pivot = d.sqrt();
        l.set(j, j, pivot);
        for i in j + 1..n {
            let mut s = a.get(i, j);
            for k in 0..j {
                s -= l.get(i, k) * l.get(j, k);
            }
            l.set(i, j, s / pivot);
        }
    }
    Some(l)
}

/// Solve `L X = B` in place (`L` lower triangular, every column of `B`).
pub fn solve_lower_in_place(l: &DenseMat, b: &mut DenseMat) {
    let n = l.nrows();
    for c in 0..b.ncols() {
        let x = b.col_mut(c);
        for i in 0..n {
            let mut s = x[i];
            for (k, xk) in x.iter().enumerate().take(i) {
                s -= l.get(i, k) * xk;
            }
            x[i] = s / l.get(i, i);
        }
    }
}

/// Solve `Lᵀ X = B` in place (`L` lower triangular).
pub fn solve_lower_transposed_in_place(l: &DenseMat, b: &mut DenseMat) {
    let n = l.nrows();
    for c in 0..b.ncols() {
        let x = b.col_mut(c);
        for i in (0..n).rev() {
            let mut s = x[i];
            for (k, xk) in x.iter().enumerate().skip(i + 1) {
                s -= l.get(k, i) * xk;
            }
            x[i] = s / l.get(i, i);
        }
    }
}

/// Solve `U X = B` in place (`U` upper triangular: the leading
/// `U.ncols()` rows of each column of `B`).
pub fn solve_upper_in_place(u: &DenseMat, b: &mut DenseMat) {
    let n = u.ncols();
    for c in 0..b.ncols() {
        let x = b.col_mut(c);
        for i in (0..n).rev() {
            let mut s = x[i];
            for (k, xk) in x.iter().enumerate().take(n).skip(i + 1) {
                s -= u.get(i, k) * xk;
            }
            x[i] = s / u.get(i, i);
        }
    }
}

/// Solve the square system `A X = B` by LU with partial pivoting (the
/// largest magnitude in the column, the lowest row on a tie). `None`
/// for a singular or non-finite pivot.
pub fn lu_solve(mut a: DenseMat, mut b: DenseMat) -> Option<DenseMat> {
    let n = a.nrows();
    if a.ncols() != n || b.nrows() != n {
        return None;
    }
    for k in 0..n {
        let mut pivot = k;
        let mut largest = a.get(k, k).abs();
        for i in k + 1..n {
            let v = a.get(i, k).abs();
            if v > largest {
                largest = v;
                pivot = i;
            }
        }
        if largest.is_nan() || largest <= 0.0 || largest.is_infinite() {
            return None;
        }
        if pivot != k {
            for j in 0..n {
                let (x, y) = (a.get(k, j), a.get(pivot, j));
                a.set(k, j, y);
                a.set(pivot, j, x);
            }
            for j in 0..b.ncols() {
                let (x, y) = (b.get(k, j), b.get(pivot, j));
                b.set(k, j, y);
                b.set(pivot, j, x);
            }
        }
        let d = a.get(k, k);
        for i in k + 1..n {
            let f = a.get(i, k) / d;
            if f == 0.0 {
                continue;
            }
            a.set(i, k, 0.0);
            for j in k + 1..n {
                a.set(i, j, a.get(i, j) - f * a.get(k, j));
            }
            for j in 0..b.ncols() {
                b.set(i, j, b.get(i, j) - f * b.get(k, j));
            }
        }
    }
    solve_upper_in_place(&a, &mut b);
    Some(b)
}

/// A Householder QR in the LAPACK layout: `R` on and above the
/// diagonal, each reflector's vector below it (its leading 1 implied),
/// `tau` per reflector; `H_k = I − τ_k v_k v_kᵀ`.
#[derive(Debug, Clone)]
pub struct HouseholderQr {
    /// `R` and the reflectors.
    pub qr: DenseMat,
    /// One scale per reflector.
    pub tau: Vec<f64>,
}

impl HouseholderQr {
    /// Factor the first `left` columns of `m` (`left <= rows`), applying
    /// every reflector to the remaining columns too: on return the
    /// columns past `left` hold `Qᵀ` times their input.
    pub fn factor_left(mut m: DenseMat, left: usize) -> HouseholderQr {
        let rows = m.nrows();
        let tau = (0..left.min(rows))
            .map(|k| HouseholderQr::reflect_column(&mut m, k))
            .collect();
        HouseholderQr { qr: m, tau }
    }

    /// The diagonal entry `R_kk`.
    pub fn r_diag(&self, k: usize) -> f64 {
        self.qr.get(k, k)
    }

    /// `Qᵀ B` in place.
    pub fn apply_qt(&self, b: &mut DenseMat) {
        let rows = self.qr.nrows();
        for (k, &t) in self.tau.iter().enumerate() {
            if t == 0.0 {
                continue;
            }
            for j in 0..b.ncols() {
                let mut w = b.get(k, j);
                for i in k + 1..rows {
                    w += self.qr.get(i, k) * b.get(i, j);
                }
                w *= t;
                b.set(k, j, b.get(k, j) - w);
                for i in k + 1..rows {
                    b.set(i, j, b.get(i, j) - w * self.qr.get(i, k));
                }
            }
        }
    }

    /// The least-squares solution `X` of `A X = B` (the factored `A`
    /// tall, full rank): `R X = (Qᵀ B)` on the leading rows.
    pub fn solve(&self, mut b: DenseMat) -> DenseMat {
        let cols = self.tau.len();
        self.apply_qt(&mut b);
        let r = DenseMat::from_fn(
            cols,
            cols,
            |i, j| if j >= i { self.qr.get(i, j) } else { 0.0 },
        );
        solve_upper_in_place(&r, &mut b);
        DenseMat::from_fn(cols, b.ncols(), |i, j| b.get(i, j))
    }
}

/// A Householder QR with column pivoting (`A P = Q R`): at each step the
/// remaining column of largest norm (the lowest index on a tie, norms
/// recomputed exactly) is moved forward.
#[derive(Debug, Clone)]
pub struct ColPivQr {
    qr: HouseholderQr,
    /// `perm[k]` is the column of `A` in position `k`.
    perm: Vec<usize>,
}

impl ColPivQr {
    /// Factor `a` (rows >= cols for a least-squares solve).
    pub fn new(mut a: DenseMat) -> ColPivQr {
        let (rows, cols) = (a.nrows(), a.ncols());
        let mut perm: Vec<usize> = (0..cols).collect();
        let steps = rows.min(cols);
        let mut tau = Vec::with_capacity(steps);
        for k in 0..steps {
            let mut best = k;
            let mut best_norm = -1.0;
            for j in k..cols {
                let mut s = 0.0;
                for i in k..rows {
                    let v = a.get(i, j);
                    s += v * v;
                }
                if s > best_norm {
                    best_norm = s;
                    best = j;
                }
            }
            if best != k {
                for i in 0..rows {
                    let (x, y) = (a.get(i, k), a.get(i, best));
                    a.set(i, k, y);
                    a.set(i, best, x);
                }
                perm.swap(k, best);
            }
            // One reflector on column k, applied to the columns after it.
            let single = HouseholderQr::reflect_column(&mut a, k);
            tau.push(single);
        }
        ColPivQr {
            qr: HouseholderQr { qr: a, tau },
            perm,
        }
    }

    /// The least-squares solution of `A X = B`: the full `R` solved
    /// (no rank truncation), then un-permuted.
    pub fn solve_lstsq(&self, b: DenseMat) -> DenseMat {
        let y = self.qr.solve(b);
        let mut x = DenseMat::zeros(y.nrows(), y.ncols());
        for (k, &column) in self.perm.iter().enumerate() {
            for j in 0..y.ncols() {
                x.set(column, j, y.get(k, j));
            }
        }
        x
    }
}

impl HouseholderQr {
    /// Reflect column `k` of `m` onto its diagonal and apply the
    /// reflector to every later column; the reflector's scale.
    fn reflect_column(m: &mut DenseMat, k: usize) -> f64 {
        let rows = m.nrows();
        let alpha = m.get(k, k);
        let mut tail = 0.0;
        for i in k + 1..rows {
            let v = m.get(i, k);
            tail += v * v;
        }
        if tail == 0.0 {
            return 0.0;
        }
        let norm = alpha.phypot(tail.sqrt());
        let beta = if alpha >= 0.0 { -norm } else { norm };
        let t = (beta - alpha) / beta;
        let scale = 1.0 / (alpha - beta);
        for i in k + 1..rows {
            m.set(i, k, m.get(i, k) * scale);
        }
        m.set(k, k, beta);
        for j in k + 1..m.ncols() {
            let mut w = m.get(k, j);
            for i in k + 1..rows {
                w += m.get(i, k) * m.get(i, j);
            }
            w *= t;
            m.set(k, j, m.get(k, j) - w);
            for i in k + 1..rows {
                m.set(i, j, m.get(i, j) - w * m.get(i, k));
            }
        }
        t
    }
}
