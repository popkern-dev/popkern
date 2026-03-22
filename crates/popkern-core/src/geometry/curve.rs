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

// ---------------------------------------------------------------------------
// Analytic-to-NURBS conversion
// ---------------------------------------------------------------------------

/// Build a full-circle NURBS (degree 2, 9 control points, 4 quarter arcs).
fn circle_to_nurbs(center: Point3, normal: Vec3, radius: f64) -> NurbsCurve {
    use crate::geometry::analytic_curves::circle_basis;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    let (u, v) = circle_basis(normal);
    let w = 1.0_f64 / 2.0_f64.sqrt(); // weight for midpoints

    // 9 control points: on-curve (weight 1) alternate with off-curve (weight 1/√2)
    let pts = [
        center + u * radius,        // angle 0
        center + (u + v) * radius,  // control pt (tangent intersection)
        center + v * radius,        // angle π/2
        center + (v - u) * radius,  // control pt
        center - u * radius,        // angle π
        center + (-u - v) * radius, // control pt
        center - v * radius,        // angle 3π/2
        center + (u - v) * radius,  // control pt
        center + u * radius,        // angle 2π
    ];

    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0,
        0.0,
        0.0,
        FRAC_PI_2,
        FRAC_PI_2,
        PI,
        PI,
        PI + FRAC_PI_2,
        PI + FRAC_PI_2,
        TAU,
        TAU,
        TAU,
    ];

    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights,
    }
}

/// Build a NURBS arc by constructing the full circle and splitting to the arc range.
fn arc_to_nurbs(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> NurbsCurve {
    let full = circle_to_nurbs(center, normal, radius);

    // Split the full circle to extract just the arc portion.
    // split_at returns Result; if the split parameter is at a boundary it may fail,
    // so we fall back to the full circle in degenerate cases.
    let right_of_start = if start_angle > 0.0 {
        match full.split_at(start_angle) {
            Ok((_left, right)) => right,
            Err(_) => full,
        }
    } else {
        full
    };

    if end_angle < std::f64::consts::TAU {
        match right_of_start.split_at(end_angle) {
            Ok((left, _right)) => left,
            Err(_) => right_of_start,
        }
    } else {
        right_of_start
    }
}

impl CurveGeometry {
    /// Convert to a [`NurbsCurve`] covering the given parameter range.
    /// The NURBS knot domain will match the range.
    pub fn to_nurbs_in_range(&self, t_range: (f64, f64)) -> NurbsCurve {
        match self {
            CurveGeometry::Line { origin, direction } => {
                let p0 = *origin + *direction * t_range.0;
                let p1 = *origin + *direction * t_range.1;
                NurbsCurve {
                    degree: 1,
                    knots: vec![t_range.0, t_range.0, t_range.1, t_range.1],
                    control_points: vec![p0, p1],
                    weights: vec![1.0, 1.0],
                }
            }
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => circle_to_nurbs(*center, *normal, *radius),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => arc_to_nurbs(*center, *normal, *radius, *start_angle, *end_angle),
            CurveGeometry::Nurbs(nurbs) => nurbs.clone(),
        }
    }

    /// Convert with the default range from [`CurveEval::parameter_range()`].
    pub fn to_nurbs(&self) -> NurbsCurve {
        let range = self.parameter_range();
        self.to_nurbs_in_range(range)
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_to_nurbs_in_range() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let nurbs = line.to_nurbs_in_range((-2.0, 4.0));
        assert_eq!(nurbs.degree, 1);
        for t in [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.99] {
            let p1 = line.point_at(t);
            let p2 = nurbs.point_at(t);
            assert!((p1.x - p2.x).abs() < 1e-10, "x differs at t={t}");
        }
    }

    #[test]
    fn circle_to_nurbs_shape() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let nurbs = circle.to_nurbs();

        // At knot boundaries the parameterizations agree exactly.
        use std::f64::consts::{FRAC_PI_2, PI, TAU};
        for t in [0.0, FRAC_PI_2, PI, PI + FRAC_PI_2] {
            let p1 = circle.point_at(t);
            let p2 = nurbs.point_at(t);
            assert!(
                (p1.x - p2.x).abs() < 1e-10,
                "x differs at t={t}: {} vs {}",
                p1.x,
                p2.x
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-10,
                "y differs at t={t}: {} vs {}",
                p1.y,
                p2.y
            );
        }

        // All NURBS-sampled points should lie on the unit circle (correct radius).
        for i in 0..20 {
            let t = TAU * (i as f64) / 20.0;
            let p = nurbs.point_at(t);
            let dist = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (dist - 1.0).abs() < 1e-10,
                "point not on circle at t={t}: dist={}",
                dist
            );
        }
    }

    #[test]
    fn nurbs_to_nurbs_is_clone() {
        let nurbs = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let geom = CurveGeometry::Nurbs(nurbs.clone());
        let promoted = geom.to_nurbs();
        assert_eq!(promoted.degree, 2);
        assert_eq!(promoted.control_points.len(), 3);
    }
}
