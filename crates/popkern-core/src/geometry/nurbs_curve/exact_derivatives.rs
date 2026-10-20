//! Exact derivatives of a rational B-spline curve to any order (The
//! NURBS Book A2.3 basis-function derivatives, A4.2 rational quotient
//! rule). [`CurveEval::derivative_at`](crate::geometry::curve::CurveEval)
//! finite-differences orders above one; a third-order join test (G3,
//! `popkern_blend::curve_join`) needs them exact.

use super::NurbsCurve;
use crate::geometry::Vec3;

impl NurbsCurve {
    /// The point and derivatives `C(t), C'(t), …, C^(n)(t)` (`n + 1`
    /// vectors; the point as a vector from the origin), exact.
    /// Derivatives above the degree are zero on each span; at an
    /// interior knot the span on the right is used (a knot of
    /// multiplicity `m` leaves only `degree - m` derivatives continuous).
    pub fn derivatives(&self, t: f64, n: usize) -> Vec<Vec3> {
        let p = self.degree;
        let count = self.control_points.len();
        if count == 0 || self.knots.len() != count + p + 1 {
            return vec![Vec3::new(0.0, 0.0, 0.0); n + 1];
        }
        let (lo, hi) = self.parameter_range_inner();
        let t = t.clamp(lo, hi);
        let span = crate::geometry::nurbs_ops::find_span(p, &self.knots, count, t);
        let ders = basis_derivatives(p, &self.knots, span, t, n);
        // Homogeneous derivatives A^(k) = sum N^(k) w P and w^(k).
        let mut aw = vec![(Vec3::new(0.0, 0.0, 0.0), 0.0); n + 1];
        for (k, row) in ders.iter().enumerate() {
            for (j, &d) in row.iter().enumerate() {
                let i = span - p + j;
                let w = self.weights.get(i).copied().unwrap_or(1.0);
                let c = self.control_points[i];
                aw[k].0 = aw[k].0 + Vec3::new(c.x, c.y, c.z) * (d * w);
                aw[k].1 += d * w;
            }
        }
        // C^(k) = (A^(k) - sum_{i=1..k} binom(k, i) w^(i) C^(k-i)) / w.
        let w0 = aw[0].1;
        let mut out: Vec<Vec3> = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let mut v = aw[k].0;
            for i in 1..=k {
                v = v - out[k - i] * (binomial(k, i) * aw[i].1);
            }
            out.push(if w0.abs() > 1e-300 { v * (1.0 / w0) } else { v });
        }
        out
    }

    fn parameter_range_inner(&self) -> (f64, f64) {
        let p = self.degree;
        let n = self.control_points.len();
        (self.knots[p], self.knots[n])
    }
}

fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// Basis functions and their derivatives up to order `n` at `t` on
/// `span` (The NURBS Book A2.3): `out[k][j]` is the `k`-th derivative
/// of basis function `span - p + j`.
fn basis_derivatives(p: usize, knots: &[f64], span: usize, t: f64, n: usize) -> Vec<Vec<f64>> {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let (mut left, mut right) = (vec![0.0; p + 1], vec![0.0; p + 1]);
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            // Lower triangle: knot differences; upper: basis values.
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r] != 0.0 {
                ndu[r][j - 1] / ndu[j][r]
            } else {
                0.0
            };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut ders = vec![vec![0.0; p + 1]; n + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = vec![vec![0.0; p + 1]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=n.min(p) {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                let rk = rk as usize;
                a[s2][0] = if ndu[pk + 1][rk] != 0.0 {
                    a[s1][0] / ndu[pk + 1][rk]
                } else {
                    0.0
                };
                d = a[s2][0] * ndu[rk][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize {
                k - 1
            } else {
                p - r
            };
            for j in j1..=j2 {
                let idx = (rk + j as isize) as usize;
                a[s2][j] = if ndu[pk + 1][idx] != 0.0 {
                    (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx]
                } else {
                    0.0
                };
                d += a[s2][j] * ndu[idx][pk];
            }
            if r <= pk {
                a[s2][k] = if ndu[pk + 1][r] != 0.0 {
                    -a[s1][k - 1] / ndu[pk + 1][r]
                } else {
                    0.0
                };
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    // Multiply through by p! / (p - k)!.
    let mut factor = p as f64;
    for (k, row) in ders.iter_mut().enumerate().skip(1) {
        if k > p {
            row.iter_mut().for_each(|v| *v = 0.0);
            continue;
        }
        for v in row.iter_mut() {
            *v *= factor;
        }
        factor *= (p as f64) - k as f64;
    }
    ders
}

#[cfg(test)]
#[path = "exact_derivatives_tests.rs"]
mod tests;
