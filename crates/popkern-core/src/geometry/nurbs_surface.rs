//! NURBS surface evaluation via tensor-product de Boor algorithm.

use crate::geometry::nurbs_ops::{Homogeneous, de_boor_1d};
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

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
    fn cp(&self, i_u: usize, i_v: usize) -> &Point3 {
        &self.control_points[i_u * self.num_v + i_v]
    }

    /// Access weight at grid position (i_u, i_v) in the flat array.
    fn wt(&self, i_u: usize, i_v: usize) -> f64 {
        self.weights[i_u * self.num_v + i_v]
    }

    /// Evaluate the surface at (u, v) using tensor-product de Boor in homogeneous coordinates.
    fn eval_homogeneous(&self, u: f64, v: f64) -> Homogeneous {
        // Step 1: For each u-row, evaluate the 1D B-spline in v direction.
        let mut intermediate: Vec<Homogeneous> = Vec::with_capacity(self.num_u);
        for i_u in 0..self.num_u {
            let row: Vec<Homogeneous> = (0..self.num_v)
                .map(|i_v| Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v)))
                .collect();
            intermediate.push(de_boor_1d(self.degree_v, &self.knots_v, &row, v));
        }

        // Step 2: Evaluate the 1D B-spline in u direction on the intermediate points.
        de_boor_1d(self.degree_u, &self.knots_u, &intermediate, u)
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    /// Helper to build a bilinear patch for testing.
    fn bilinear_patch() -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                // Row u=0: (0,0,0), (0,1,0)
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                // Row u=1: (1,0,0), (1,1,0)
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        }
    }

    #[test]
    fn bilinear_patch_center() {
        let surface = bilinear_patch();
        let p = surface.point_at(0.5, 0.5);
        assert!((p.x - 0.5).abs() < 1e-10);
        assert!((p.y - 0.5).abs() < 1e-10);
        assert!((p.z).abs() < 1e-10);
    }

    #[test]
    fn bilinear_patch_corners() {
        let surface = bilinear_patch();
        let p00 = surface.point_at(0.0, 0.0);
        assert!((p00.x).abs() < 1e-10);
        assert!((p00.y).abs() < 1e-10);
        // Use slightly less than 1.0 to stay in valid knot span
        let p11 = surface.point_at(1.0 - 1e-12, 1.0 - 1e-12);
        assert!((p11.x - 1.0).abs() < 1e-6);
        assert!((p11.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bilinear_patch_normal_is_z() {
        let surface = bilinear_patch();
        let n = surface.normal_at(0.5, 0.5);
        // Normal should point in +z or -z direction for a flat XY patch
        let nlen = n.length();
        assert!(nlen > 1e-10, "Normal should be non-zero");
        let nz = n.z.abs() / nlen;
        assert!((nz - 1.0).abs() < 1e-6, "Normal should be ±z, got {:?}", n);
    }

    #[test]
    fn parameter_range() {
        let surface = bilinear_patch();
        let ((u0, u1), (v0, v1)) = surface.parameter_range();
        assert_eq!(u0, 0.0);
        assert_eq!(u1, 1.0);
        assert_eq!(v0, 0.0);
        assert_eq!(v1, 1.0);
    }

    #[test]
    fn bilinear_patch_partial_u() {
        let surface = bilinear_patch();
        let du = surface.partial_u(0.5, 0.5);
        // For a bilinear patch, d/du should be (1, 0, 0)
        assert!((du.x - 1.0).abs() < 1e-4);
        assert!((du.y).abs() < 1e-4);
    }

    #[test]
    fn bilinear_patch_partial_v() {
        let surface = bilinear_patch();
        let dv = surface.partial_v(0.5, 0.5);
        // For a bilinear patch, d/dv should be (0, 1, 0)
        assert!((dv.x).abs() < 1e-4);
        assert!((dv.y - 1.0).abs() < 1e-4);
    }

    #[test]
    fn bilinear_patch_at_end() {
        let surface = bilinear_patch();
        // Evaluate at exactly t=1.0 (clamped end)
        let p = surface.point_at(1.0, 1.0);
        assert!((p.x - 1.0).abs() < 1e-6);
        assert!((p.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bilinear_patch_at_before_start() {
        let surface = bilinear_patch();
        // Evaluate before start — clamped to degree span, extrapolates
        // Just verify it doesn't panic and returns finite values
        let p = surface.point_at(-0.1, -0.1);
        assert!(p.x.is_finite());
        assert!(p.y.is_finite());
    }
}
