//! Eigenvalues of a symmetric matrix: Householder reduction to
//! tridiagonal form, then the implicit QL iteration (the EISPACK
//! `tred` / `tql` pair, eigenvalues only), in one fixed operation order
//! (PK-676). The singular values of a symmetric matrix are the absolute
//! values of its eigenvalues, so condition numbers of the KKT and
//! normal-equation matrices come from here: robust on badly conditioned
//! input, and O(n³) with a small constant where one-sided Jacobi needs
//! many sweeps.

use super::DenseMat;
use crate::PortableMath;

/// QL iterations per eigenvalue before giving up (a few suffice).
const MAX_ITERATIONS: usize = 60;

/// The eigenvalues of the symmetric `a` (only its lower triangle is
/// read), in no particular order. `None` when the QL iteration does not
/// converge.
// The reduction is written in the textbook's index form (rows and
// columns of the working triangle addressed together), which reads
// against the algorithm; iterator rewrites would obscure it.
#[allow(clippy::needless_range_loop)]
pub fn symmetric_eigenvalues(a: &DenseMat) -> Option<Vec<f64>> {
    let n = a.nrows();
    if n == 0 {
        return Some(Vec::new());
    }
    // Row-major working copy of the lower triangle (a[i][k], k <= i).
    let mut w: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..=i).map(|k| a.get(i, k)).collect())
        .collect();
    let mut d = vec![0.0f64; n];
    let mut e = vec![0.0f64; n];
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        if l > 0 {
            let scale: f64 = w[i][..=l].iter().map(|v| v.abs()).sum();
            if scale == 0.0 {
                e[i] = w[i][l];
            } else {
                for k in 0..=l {
                    w[i][k] /= scale;
                    h += w[i][k] * w[i][k];
                }
                let f = w[i][l];
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                w[i][l] = f - g;
                let mut f = 0.0;
                for j in 0..=l {
                    let mut g = 0.0;
                    for k in 0..=j {
                        g += w[j][k] * w[i][k];
                    }
                    for k in j + 1..=l {
                        g += w[k][j] * w[i][k];
                    }
                    e[j] = g / h;
                    f += e[j] * w[i][j];
                }
                let hh = f / (h + h);
                for j in 0..=l {
                    let f = w[i][j];
                    let g = e[j] - hh * f;
                    e[j] = g;
                    for k in 0..=j {
                        w[j][k] -= f * e[k] + g * w[i][k];
                    }
                }
            }
        } else {
            e[i] = w[i][l];
        }
        d[i] = h;
    }
    for (i, di) in d.iter_mut().enumerate() {
        *di = w[i][i];
    }
    tql(&mut d, &mut e)?;
    Some(d)
}

/// The implicit QL iteration on the tridiagonal `(d, e)` (`e[i]` the
/// subdiagonal entry below row `i`, `e[0]` unused): `d` becomes the
/// eigenvalues.
fn tql(d: &mut [f64], e: &mut [f64]) -> Option<()> {
    let n = d.len();
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iterations = 0;
        loop {
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iterations += 1;
            if iterations > MAX_ITERATIONS {
                return None;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.phypot(1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r.abs() } else { -r.abs() });
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut deflated = false;
            for i in (l..m).rev() {
                let f = s * e[i];
                let b = c * e[i];
                r = f.phypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    deflated = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
            }
            if deflated {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    Some(())
}
