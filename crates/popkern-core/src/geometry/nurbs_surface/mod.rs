//! NURBS surface evaluation via tensor-product de Boor algorithm.

use crate::geometry::nurbs_ops::{Homogeneous, MAX_DE_BOOR_DEGREE, de_boor_1d_local, find_span};
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

mod degree_and_split;
mod eval_derivs;
mod iso_curve;
mod knot_ops;

/// A Non-Uniform Rational B-Spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NurbsSurface {
    /// Polynomial degree in the u direction.
    pub degree_u: usize,
    /// Polynomial degree in the v direction.
    pub degree_v: usize,
    /// Knot vector in the u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in the v direction.
    pub knots_v: Vec<f64>,
    /// Control point grid, stored row-major (u varies within each row).
    /// Length = num_u * num_v.
    pub control_points: Vec<Point3>,
    /// Number of control points in the u direction.
    pub num_u: usize,
    /// Number of control points in the v direction.
    pub num_v: usize,
    /// Weights, one per control point.
    pub weights: Vec<f64>,
}

impl NurbsSurface {
    /// Access control point at grid position (i_u, i_v) in the flat array.
    pub(super) fn cp(&self, i_u: usize, i_v: usize) -> &Point3 {
        &self.control_points[i_u * self.num_v + i_v]
    }

    /// Access weight at grid position (i_u, i_v) in the flat array.
    pub(super) fn wt(&self, i_u: usize, i_v: usize) -> f64 {
        self.weights[i_u * self.num_v + i_v]
    }

    /// Evaluate the surface at (u, v) using tensor-product de Boor in
    /// homogeneous coordinates.
    ///
    /// Span-local: a de Boor evaluation at (u, v) only ever reads the
    /// `(degree_u + 1) × (degree_v + 1)` control window around the two
    /// knot spans, so that is all this builds — one `find_span` per
    /// direction, `degree_u + 1` row recursions, one column recursion.
    /// Cost is O((p+1)·(q+1)) regardless of grid size, always
    /// allocation-free (the scratch is bounded by the degree cap, not
    /// by `num_u`/`num_v`), and bit-identical to evaluating the full
    /// grid because [`de_boor_1d_local`] runs the same lerps on the
    /// same values in the same order (PK-281; the earlier full-grid
    /// form evaluated every row and discarded all but `degree_u + 1`
    /// intermediates, which is what made `de_boor_1d` dominate the
    /// G2 cube profile at 74/100 samples).
    fn eval_homogeneous(&self, u: f64, v: f64) -> Homogeneous {
        let zero = Homogeneous {
            wx: 0.0,
            wy: 0.0,
            wz: 0.0,
            w: 0.0,
        };
        let span_u = find_span(self.degree_u, &self.knots_u, self.num_u, u);
        let span_v = find_span(self.degree_v, &self.knots_v, self.num_v, v);
        let mut intermediate = [zero; MAX_DE_BOOR_DEGREE];
        let mut row = [zero; MAX_DE_BOOR_DEGREE];
        for (j_u, intermediate_slot) in intermediate.iter_mut().enumerate().take(self.degree_u + 1)
        {
            let i_u = span_u - self.degree_u + j_u;
            for (j_v, row_slot) in row.iter_mut().enumerate().take(self.degree_v + 1) {
                let i_v = span_v - self.degree_v + j_v;
                *row_slot = Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v));
            }
            *intermediate_slot = de_boor_1d_local(
                self.degree_v,
                &self.knots_v,
                span_v,
                &row[..=self.degree_v],
                v,
            );
        }
        de_boor_1d_local(
            self.degree_u,
            &self.knots_u,
            span_u,
            &intermediate[..=self.degree_u],
            u,
        )
    }
}

impl SurfaceEval for NurbsSurface {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        self.eval_homogeneous(u, v).to_point3()
    }

    // First partials are exact rational derivatives via the fused
    // basis-contraction evaluation (PK-281) — previously 1e-7 central
    // differences, which cost two full evaluations each and carried
    // ~1e-9 truncation error into every Gauss–Newton projection.
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_first_partials(u, v).1
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_first_partials(u, v).2
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        let (_, du, dv) = self.point_and_first_partials(u, v);
        du.cross(&dv)
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        let p = self.degree_u;
        let n_u = self.num_u - 1;
        let q = self.degree_v;
        let n_v = self.num_v - 1;
        (
            (self.knots_u[p], self.knots_u[n_u + 1]),
            (self.knots_v[q], self.knots_v[n_v + 1]),
        )
    }

    // Second partials are exact rational derivatives via the fused
    // two-derivative basis contraction (PK-285) — previously 1e-5
    // central differences of point_at (3–5 evaluations each), whose
    // truncation error fed every curvature and G2 metric on NURBS.
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).3
    }

    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).4
    }

    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).5
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        use crate::geometry::surface::closest_point_gauss_newton;
        closest_point_gauss_newton(self, point)
    }

    fn bounding_box(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> crate::geometry::intersection::surface_subdivision::Aabb {
        use crate::geometry::surface::bounding_box_by_sampling;
        bounding_box_by_sampling(self, u_range, v_range)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod eval_span_local_tests;

#[cfg(test)]
mod iso_curve_tests;
