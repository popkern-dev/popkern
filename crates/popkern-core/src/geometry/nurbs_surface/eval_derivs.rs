//! Fused analytic evaluation of a NURBS surface point and both first
//! partial derivatives (PK-281).
//!
//! One basis computation per direction (values + first derivatives via
//! the Cox–de Boor triangle), one contraction over the local
//! `(degree_u + 1) × (degree_v + 1)` control window, then the rational
//! quotient rule. Replaces the 5-evaluation pattern (1 point + 2×2
//! central-difference partials) that dominated every Gauss–Newton
//! projection iteration, and removes the 1e-7 finite-difference
//! truncation error from the derivatives entirely.

use crate::geometry::nurbs_ops::{MAX_DE_BOOR_DEGREE, find_span};
use crate::geometry::{Point3, Vec3};

use super::NurbsSurface;

/// B-spline basis function values and first derivatives at `t`.
///
/// Returns `(values, derivs)` for the `degree + 1` basis functions
/// that are non-zero on `span` (local index `r` ↔ global index
/// `span - degree + r`). Values via the standard left/right
/// Cox–de Boor triangle; derivatives from the degree-(p−1) row:
/// `N′_i^p = p·(N_i^{p−1}/(k[i+p]−k[i]) − N_{i+1}^{p−1}/(k[i+p+1]−k[i+1]))`,
/// with zero-multiplicity denominators contributing zero.
pub(super) fn basis_funs_and_first_derivs(
    degree: usize,
    knots: &[f64],
    span: usize,
    t: f64,
) -> ([f64; MAX_DE_BOOR_DEGREE], [f64; MAX_DE_BOOR_DEGREE]) {
    debug_assert!(
        degree < MAX_DE_BOOR_DEGREE,
        "basis_funs_and_first_derivs: degree {degree} exceeds MAX_DE_BOOR_DEGREE"
    );
    let mut values = [0.0; MAX_DE_BOOR_DEGREE];
    let mut derivs = [0.0; MAX_DE_BOOR_DEGREE];
    let mut left = [0.0; MAX_DE_BOOR_DEGREE];
    let mut right = [0.0; MAX_DE_BOOR_DEGREE];
    // Degree-(p−1) row of the triangle, captured on the way up.
    let mut lower = [0.0; MAX_DE_BOOR_DEGREE];

    values[0] = 1.0;
    if degree == 0 {
        return (values, derivs);
    }
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() < 1e-30 {
                0.0
            } else {
                values[r] / denom
            };
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
        if j == degree - 1 {
            lower[..=j].copy_from_slice(&values[..=j]);
        }
    }
    if degree == 1 {
        // The captured row is the trivial degree-0 row.
        lower[0] = 1.0;
    }

    let p = degree as f64;
    for r in 0..=degree {
        let i = span - degree + r;
        let mut val = 0.0;
        if r >= 1 {
            let denom = knots[i + degree] - knots[i];
            if denom.abs() > 1e-30 {
                val += lower[r - 1] / denom;
            }
        }
        if r < degree {
            let denom = knots[i + degree + 1] - knots[i + 1];
            if denom.abs() > 1e-30 {
                val -= lower[r] / denom;
            }
        }
        derivs[r] = p * val;
    }

    (values, derivs)
}

/// Homogeneous 4-vector accumulator for the tensor contraction.
#[derive(Clone, Copy)]
struct H4 {
    x: f64,
    y: f64,
    z: f64,
    w: f64,
}

impl H4 {
    const ZERO: H4 = H4 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        w: 0.0,
    };

    #[inline]
    fn add_scaled(&mut self, other: &H4, s: f64) {
        self.x += s * other.x;
        self.y += s * other.y;
        self.z += s * other.z;
        self.w += s * other.w;
    }
}

impl NurbsSurface {
    /// Point and both first partial derivatives in one fused analytic
    /// evaluation.
    ///
    /// Exact rational derivatives via the quotient rule
    /// `S_u = (A_u − w_u·S)/w` on the homogeneous forms — no finite
    /// differencing. Cost is one basis+deriv computation per
    /// direction plus a single pass over the local control window,
    /// roughly the price of one `point_at` and ~5× cheaper than the
    /// point + 2×2-FD pattern it replaces in Gauss–Newton loops.
    pub fn point_and_first_partials(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        let span_u = find_span(self.degree_u, &self.knots_u, self.num_u, u);
        let span_v = find_span(self.degree_v, &self.knots_v, self.num_v, v);
        let (nu, du) = basis_funs_and_first_derivs(self.degree_u, &self.knots_u, span_u, u);
        let (nv, dv) = basis_funs_and_first_derivs(self.degree_v, &self.knots_v, span_v, v);

        let mut a = H4::ZERO;
        let mut a_u = H4::ZERO;
        let mut a_v = H4::ZERO;
        for j_u in 0..=self.degree_u {
            let i_u = span_u - self.degree_u + j_u;
            // Row sums over v first: Σ_jv N_jv·Pw and Σ_jv N′_jv·Pw.
            let mut row = H4::ZERO;
            let mut row_dv = H4::ZERO;
            for j_v in 0..=self.degree_v {
                let i_v = span_v - self.degree_v + j_v;
                let cp = self.cp(i_u, i_v);
                let wt = self.wt(i_u, i_v);
                let hw = H4 {
                    x: wt * cp.x,
                    y: wt * cp.y,
                    z: wt * cp.z,
                    w: wt,
                };
                row.add_scaled(&hw, nv[j_v]);
                row_dv.add_scaled(&hw, dv[j_v]);
            }
            a.add_scaled(&row, nu[j_u]);
            a_u.add_scaled(&row, du[j_u]);
            a_v.add_scaled(&row_dv, nu[j_u]);
        }

        if a.w.abs() < 1e-30 {
            return (
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
            );
        }
        let inv_w = 1.0 / a.w;
        let s = Point3::new(a.x * inv_w, a.y * inv_w, a.z * inv_w);
        let s_u = Vec3::new(
            (a_u.x - a_u.w * s.x) * inv_w,
            (a_u.y - a_u.w * s.y) * inv_w,
            (a_u.z - a_u.w * s.z) * inv_w,
        );
        let s_v = Vec3::new(
            (a_v.x - a_v.w * s.x) * inv_w,
            (a_v.y - a_v.w * s.y) * inv_w,
            (a_v.z - a_v.w * s.z) * inv_w,
        );
        (s, s_u, s_v)
    }
}

#[cfg(test)]
#[path = "eval_derivs_tests.rs"]
mod tests;
