//! NURBS (Non-Uniform Rational B-Spline) curve evaluation via de Boor's algorithm.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_ops::{Homogeneous, de_boor_1d};
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A Non-Uniform Rational B-Spline curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NurbsCurve {
    /// Polynomial degree.
    pub degree: usize,
    /// Knot vector (non-decreasing, length = control_points.len() + degree + 1).
    pub knots: Vec<f64>,
    /// Control points.
    pub control_points: Vec<Point3>,
    /// Weights (one per control point, all positive).
    pub weights: Vec<f64>,
}

impl NurbsCurve {
    /// Evaluate the curve at parameter t using de Boor's algorithm in homogeneous coordinates.
    fn de_boor(&self, t: f64) -> Homogeneous {
        let homogeneous_points: Vec<Homogeneous> = self
            .control_points
            .iter()
            .zip(self.weights.iter())
            .map(|(p, &w)| Homogeneous::from_point_weight(p, w))
            .collect();
        de_boor_1d(self.degree, &self.knots, &homogeneous_points, t)
    }
}

impl CurveEval for NurbsCurve {
    fn point_at(&self, t: f64) -> Point3 {
        self.de_boor(t).to_point3()
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        if order == 0 {
            let p = self.point_at(t);
            return Vec3::new(p.x, p.y, p.z);
        }

        // Finite difference approximation.
        // For order 1: central difference (f(t+h) - f(t-h)) / (2h)
        // For higher orders: recurse on lower order derivatives.
        let h = 1e-7;
        let (t_min, t_max) = self.parameter_range();

        if order == 1 {
            let t_fwd = (t + h).min(t_max);
            let t_bwd = (t - h).max(t_min);
            let dt = t_fwd - t_bwd;
            if dt.abs() < 1e-30 {
                return Vec3::zero();
            }
            let p1 = self.point_at(t_fwd);
            let p0 = self.point_at(t_bwd);
            Vec3::new((p1.x - p0.x) / dt, (p1.y - p0.y) / dt, (p1.z - p0.z) / dt)
        } else {
            // Higher order: finite difference of (order-1) derivatives
            let t_fwd = (t + h).min(t_max);
            let t_bwd = (t - h).max(t_min);
            let dt = t_fwd - t_bwd;
            if dt.abs() < 1e-30 {
                return Vec3::zero();
            }
            let d1 = self.derivative_at(t_fwd, order - 1);
            let d0 = self.derivative_at(t_bwd, order - 1);
            Vec3::new((d1.x - d0.x) / dt, (d1.y - d0.y) / dt, (d1.z - d0.z) / dt)
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        let d = self.derivative_at(t, 1);
        d.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        let p = self.degree;
        let n = self.control_points.len() - 1;
        (self.knots[p], self.knots[n + 1])
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        // 5-point Gauss-Legendre quadrature
        let nodes: [f64; 5] = [
            -0.906_179_845_938_664,
            -0.538_469_310_105_683,
            0.0,
            0.538_469_310_105_683,
            0.906_179_845_938_664,
        ];
        let weights: [f64; 5] = [
            0.236_926_885_056_189_1,
            0.478_628_670_499_366_5,
            0.568_888_888_888_889,
            0.478_628_670_499_366_5,
            0.236_926_885_056_189_1,
        ];

        let mid = 0.5 * (t1 + t0);
        let half_len = 0.5 * (t1 - t0);

        let mut sum = 0.0;
        for i in 0..5 {
            let t = mid + half_len * nodes[i];
            let d = self.derivative_at(t, 1);
            sum += weights[i] * d.length();
        }
        sum * half_len
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    #[test]
    fn degree_1_nurbs_is_polyline() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let p = curve.point_at(0.5);
        assert!((p.x - 0.5).abs() < 1e-10);
        assert!((p.y).abs() < 1e-10);
    }

    #[test]
    fn degree_1_endpoints() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 3.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let p0 = curve.point_at(0.0);
        assert!((p0.x).abs() < 1e-10);
        let p1 = curve.point_at(1.0 - 1e-12);
        assert!((p1.x - 2.0).abs() < 1e-6);
    }

    #[test]
    fn degree_2_quarter_circle() {
        let w = 1.0_f64 / 2.0_f64.sqrt();
        let curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            weights: vec![1.0, w, 1.0],
        };
        // At t=0.5, point should be on the unit circle
        let p = curve.point_at(0.5);
        let radius = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (radius - 1.0).abs() < 1e-10,
            "Expected radius 1.0, got {}",
            radius
        );
    }

    #[test]
    fn degree_2_bezier_midpoint() {
        let curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let p = curve.point_at(0.5);
        assert!((p.x - 0.5).abs() < 1e-10);
        assert!((p.y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn derivative_degree_1_is_constant() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let d = curve.derivative_at(0.5, 1);
        assert!((d.x - 3.0).abs() < 1e-8);
        assert!((d.y - 4.0).abs() < 1e-8);
    }

    #[test]
    fn length_degree_1() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let len = curve.length_between(0.0, 1.0);
        assert!((len - 5.0).abs() < 1e-6, "Expected length 5.0, got {}", len);
    }

    #[test]
    fn parameter_range() {
        let curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let (t0, t1) = curve.parameter_range();
        assert_eq!(t0, 0.0);
        assert_eq!(t1, 1.0);
    }

    #[test]
    fn tangent_at_degree_1() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let t = curve.tangent_at(0.5);
        assert!((t.length() - 1.0).abs() < 1e-6, "Tangent should be unit");
    }

    #[test]
    fn derivative_order_0() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let d0 = curve.derivative_at(0.5, 0);
        // order 0 returns position as vec
        assert!((d0.x - 0.5).abs() < 1e-6);
    }

    #[test]
    fn derivative_order_2() {
        let curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        // Second derivative should be non-zero for a quadratic
        let d2 = curve.derivative_at(0.5, 2);
        assert!(d2.length() > 1e-3);
    }

    #[test]
    fn find_span_clamped_at_end() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        // Evaluate at exactly 1.0 (clamped)
        let p = curve.point_at(1.0);
        assert!((p.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn find_span_clamped_before_start() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        // Evaluate before start — clamped to span p, so result is at t=0 origin
        let p = curve.point_at(-0.5);
        // The find_span clamps to p=1, so de Boor uses span 1 with t=-0.5
        // which extrapolates. Just check it doesn't panic.
        assert!(p.x.is_finite());
    }

    #[test]
    fn multi_span_nurbs_curve() {
        // Degree 2 with 4 control points = 2 spans
        let curve = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.25, 1.0, 0.0),
                Point3::new(0.75, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };
        let p = curve.point_at(0.5);
        // Should be roughly at (0.5, something, 0)
        assert!((p.x - 0.5).abs() < 0.1);
        // Length should be computable
        let len = curve.length_between(0.0, 1.0);
        assert!(len > 1.0);
    }
}
