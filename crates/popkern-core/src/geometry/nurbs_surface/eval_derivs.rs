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

use crate::geometry::nurbs_ops::basis::{basis_funs_and_first_derivs, basis_funs_and_two_derivs};
use crate::geometry::nurbs_ops::find_span;
use crate::geometry::{Point3, Vec3};

use super::NurbsSurface;

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

    /// Point, first, and second partial derivatives in one fused
    /// analytic evaluation: `(S, S_u, S_v, S_uu, S_uv, S_vv)`.
    ///
    /// Exact rational second derivatives from differentiating
    /// `A = w·S` twice (PK-285):
    /// `S_uu = (A_uu − 2·w_u·S_u − w_uu·S)/w`,
    /// `S_uv = (A_uv − w_u·S_v − w_v·S_u − w_uv·S)/w`.
    /// Replaces the 1e-5 finite differences (3–5 evaluations each)
    /// that fed curvature and G2 metrics.
    #[allow(clippy::type_complexity)]
    pub fn point_and_two_partials(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        let span_u = find_span(self.degree_u, &self.knots_u, self.num_u, u);
        let span_v = find_span(self.degree_v, &self.knots_v, self.num_v, v);
        let (nu, du, duu) = basis_funs_and_two_derivs(self.degree_u, &self.knots_u, span_u, u);
        let (nv, dv, dvv) = basis_funs_and_two_derivs(self.degree_v, &self.knots_v, span_v, v);

        let mut a = H4::ZERO;
        let mut a_u = H4::ZERO;
        let mut a_v = H4::ZERO;
        let mut a_uu = H4::ZERO;
        let mut a_uv = H4::ZERO;
        let mut a_vv = H4::ZERO;
        for j_u in 0..=self.degree_u {
            let i_u = span_u - self.degree_u + j_u;
            let mut row = H4::ZERO;
            let mut row_dv = H4::ZERO;
            let mut row_dvv = H4::ZERO;
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
                row_dvv.add_scaled(&hw, dvv[j_v]);
            }
            a.add_scaled(&row, nu[j_u]);
            a_u.add_scaled(&row, du[j_u]);
            a_uu.add_scaled(&row, duu[j_u]);
            a_v.add_scaled(&row_dv, nu[j_u]);
            a_uv.add_scaled(&row_dv, du[j_u]);
            a_vv.add_scaled(&row_dvv, nu[j_u]);
        }

        if a.w.abs() < 1e-30 {
            let z = Vec3::new(0.0, 0.0, 0.0);
            return (Point3::new(0.0, 0.0, 0.0), z, z, z, z, z);
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
        let s_uu = Vec3::new(
            (a_uu.x - 2.0 * a_u.w * s_u.x - a_uu.w * s.x) * inv_w,
            (a_uu.y - 2.0 * a_u.w * s_u.y - a_uu.w * s.y) * inv_w,
            (a_uu.z - 2.0 * a_u.w * s_u.z - a_uu.w * s.z) * inv_w,
        );
        let s_uv = Vec3::new(
            (a_uv.x - a_u.w * s_v.x - a_v.w * s_u.x - a_uv.w * s.x) * inv_w,
            (a_uv.y - a_u.w * s_v.y - a_v.w * s_u.y - a_uv.w * s.y) * inv_w,
            (a_uv.z - a_u.w * s_v.z - a_v.w * s_u.z - a_uv.w * s.z) * inv_w,
        );
        let s_vv = Vec3::new(
            (a_vv.x - 2.0 * a_v.w * s_v.x - a_vv.w * s.x) * inv_w,
            (a_vv.y - 2.0 * a_v.w * s_v.y - a_vv.w * s.y) * inv_w,
            (a_vv.z - 2.0 * a_v.w * s_v.z - a_vv.w * s.z) * inv_w,
        );
        (s, s_u, s_v, s_uu, s_uv, s_vv)
    }
}

#[cfg(test)]
#[path = "eval_derivs_tests.rs"]
mod tests;
