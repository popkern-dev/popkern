//! NURBS surface evaluation via tensor-product de Boor algorithm.

use crate::geometry::nurbs_ops::{Homogeneous, MAX_DE_BOOR_DEGREE, de_boor_1d};
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

mod degree_and_split;
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
    /// Hot path: when both `num_u` and `num_v` fit within
    /// [`MAX_DE_BOOR_DEGREE`], both scratch arrays (the `intermediate`
    /// column of partially-evaluated rows and the per-row buffer
    /// feeding `de_boor_1d`) live on the stack and the call is
    /// allocation-free. For larger grids (e.g. offset NURBS surfaces
    /// with 64+ control points per direction), the function falls
    /// back to heap-allocated scratch — those calls are not the
    /// allocation-profile bottleneck.
    ///
    /// Measured impact of the allocation-free path: a 12-edge cube
    /// fillet's `record_fallback_cap_support_continuity` path used to
    /// dominate the allocation profile because every NURBS surface
    /// evaluation built a fresh `Vec<Homogeneous>` per row plus one
    /// for the intermediate column. The cascade ran ~12,500
    /// allocations per cap/support seam edge measured (cube fillet:
    /// 269k → 612k after Stage C1+ pcurve authoring widened the
    /// "measurable seam" set). With the stack-scratch path the same
    /// cube fillet runs at 195k allocations — below the pre-Stage-
    /// C1+ baseline.
    fn eval_homogeneous(&self, u: f64, v: f64) -> Homogeneous {
        let zero = Homogeneous {
            wx: 0.0,
            wy: 0.0,
            wz: 0.0,
            w: 0.0,
        };
        if self.num_u <= MAX_DE_BOOR_DEGREE && self.num_v <= MAX_DE_BOOR_DEGREE {
            let mut intermediate = [zero; MAX_DE_BOOR_DEGREE];
            let mut row = [zero; MAX_DE_BOOR_DEGREE];
            for (i_u, intermediate_slot) in intermediate.iter_mut().enumerate().take(self.num_u) {
                for (i_v, row_slot) in row.iter_mut().enumerate().take(self.num_v) {
                    *row_slot =
                        Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v));
                }
                *intermediate_slot =
                    de_boor_1d(self.degree_v, &self.knots_v, &row[..self.num_v], v);
            }
            de_boor_1d(self.degree_u, &self.knots_u, &intermediate[..self.num_u], u)
        } else {
            // Heap fallback for large grids (offset surfaces, swept NURBS,
            // etc.). These are not the allocation hotspot — the cap-
            // support continuity sampler iterates analytic and Gregory
            // surfaces that fit the stack path.
            let mut intermediate: Vec<Homogeneous> = Vec::with_capacity(self.num_u);
            for i_u in 0..self.num_u {
                let row: Vec<Homogeneous> = (0..self.num_v)
                    .map(|i_v| Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v)))
                    .collect();
                intermediate.push(de_boor_1d(self.degree_v, &self.knots_v, &row, v));
            }
            de_boor_1d(self.degree_u, &self.knots_u, &intermediate, u)
        }
    }
}

impl SurfaceEval for NurbsSurface {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        self.eval_homogeneous(u, v).to_point3()
    }

    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let p_fwd = self.point_at(u + h, v);
        let p_bwd = self.point_at(u - h, v);
        let inv_2h = 1.0 / (2.0 * h);
        Vec3::new(
            (p_fwd.x - p_bwd.x) * inv_2h,
            (p_fwd.y - p_bwd.y) * inv_2h,
            (p_fwd.z - p_bwd.z) * inv_2h,
        )
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let p_fwd = self.point_at(u, v + h);
        let p_bwd = self.point_at(u, v - h);
        let inv_2h = 1.0 / (2.0 * h);
        Vec3::new(
            (p_fwd.x - p_bwd.x) * inv_2h,
            (p_fwd.y - p_bwd.y) * inv_2h,
            (p_fwd.z - p_bwd.z) * inv_2h,
        )
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        let du = self.partial_u(u, v);
        let dv = self.partial_v(u, v);
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

    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        // Central finite difference: f''(u) ≈ (f(u+h) - 2f(u) + f(u-h)) / h²
        let h = 1e-5;
        let origin = Point3::new(0.0, 0.0, 0.0);
        let fwd = self.point_at(u + h, v) - origin;
        let mid = self.point_at(u, v) - origin;
        let bwd = self.point_at(u - h, v) - origin;
        (fwd - mid * 2.0 + bwd) * (1.0 / (h * h))
    }

    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        // Central finite difference for mixed partial
        let h = 1e-5;
        let origin = Point3::new(0.0, 0.0, 0.0);
        let pp = self.point_at(u + h, v + h) - origin;
        let pm = self.point_at(u + h, v - h) - origin;
        let mp = self.point_at(u - h, v + h) - origin;
        let mm = self.point_at(u - h, v - h) - origin;
        (pp - pm - mp + mm) * (1.0 / (4.0 * h * h))
    }

    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-5;
        let origin = Point3::new(0.0, 0.0, 0.0);
        let fwd = self.point_at(u, v + h) - origin;
        let mid = self.point_at(u, v) - origin;
        let bwd = self.point_at(u, v - h) - origin;
        (fwd - mid * 2.0 + bwd) * (1.0 / (h * h))
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
mod iso_curve_tests;
