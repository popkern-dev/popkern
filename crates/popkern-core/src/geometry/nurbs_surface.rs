//! NURBS surface evaluation (stub — full implementation in T17).

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

impl SurfaceEval for NurbsSurface {
    fn point_at(&self, _u: f64, _v: f64) -> Point3 {
        Point3::new(0.0, 0.0, 0.0) // stub
    }

    fn partial_u(&self, _u: f64, _v: f64) -> Vec3 {
        Vec3::zero() // stub
    }

    fn partial_v(&self, _u: f64, _v: f64) -> Vec3 {
        Vec3::zero() // stub
    }

    fn normal_at(&self, _u: f64, _v: f64) -> Vec3 {
        Vec3::new(0.0, 0.0, 1.0) // stub
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
