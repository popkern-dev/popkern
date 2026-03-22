//! Curve geometry types, evaluation trait, and topology binding.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Evaluation interface for parametric curves.
pub trait CurveEval {
    /// Point on the curve at parameter t.
    fn point_at(&self, t: f64) -> Point3;
    /// Derivative of order `order` at parameter t.
    fn derivative_at(&self, t: f64, order: usize) -> Vec3;
    /// Unit tangent vector at parameter t.
    fn tangent_at(&self, t: f64) -> Vec3;
    /// Valid parameter range [t_start, t_end].
    fn parameter_range(&self) -> (f64, f64);
    /// Arc length between two parameter values.
    fn length_between(&self, t0: f64, t1: f64) -> f64;
}

/// Geometric curve types.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CurveGeometry {
    /// A straight line.
    Line {
        /// A point on the line.
        origin: Point3,
        /// Direction vector (not necessarily unit length; length determines parameterization speed).
        direction: Vec3,
    },
    /// A full circle.
    Circle {
        /// Center point.
        center: Point3,
        /// Normal to the circle plane (unit vector).
        normal: Vec3,
        /// Radius.
        radius: f64,
    },
    /// A circular arc (portion of a circle).
    Arc {
        /// Center point.
        center: Point3,
        /// Normal to the arc plane (unit vector).
        normal: Vec3,
        /// Radius.
        radius: f64,
        /// Start angle in radians.
        start_angle: f64,
        /// End angle in radians.
        end_angle: f64,
    },
    /// A NURBS (Non-Uniform Rational B-Spline) curve.
    Nurbs(NurbsCurve),
}

/// Binding between an edge and its geometric curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurveOnEdge {
    /// The geometric curve.
    pub geometry: CurveGeometry,
    /// Parameter at the start vertex.
    pub t_start: f64,
    /// Parameter at the end vertex.
    pub t_end: f64,
}
