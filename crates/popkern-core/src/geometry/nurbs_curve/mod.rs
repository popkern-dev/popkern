//! NURBS (Non-Uniform Rational B-Spline) curve evaluation via de Boor's algorithm.

use crate::error::CoreError;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_ops::{self, Homogeneous, de_boor_1d};
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
    /// Insert a single knot at parameter `t`. The curve shape is mathematically unchanged.
    pub fn insert_knot(&self, t: f64) -> Result<NurbsCurve, CoreError> {
        let (t_min, t_max) = self.parameter_range();
        if t < t_min || t > t_max {
            return Err(CoreError::PreconditionFailed(format!(
                "knot {} outside parameter range [{}, {}]",
                t, t_min, t_max
            )));
        }
        let (knots, pts, wts) = nurbs_ops::knot_insertion::insert_knot_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        );
        Ok(NurbsCurve {
            degree: self.degree,
            knots,
            control_points: pts,
            weights: wts,
        })
    }

    /// Attempt to remove the knot at `knot_index` from the knot vector.
    ///
    /// Returns `Some(new_curve)` if removal is within `tolerance`, or `None` if
    /// the knot cannot be removed within the given tolerance.  Use `tolerance = 0.0`
    /// for exact removal only (e.g. when reversing a prior insertion).
    pub fn remove_knot(&self, knot_index: usize, tolerance: f64) -> Option<NurbsCurve> {
        let (k, p, w) = nurbs_ops::knot_removal::remove_knot_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            knot_index,
            tolerance,
        )?;
        Some(NurbsCurve {
            degree: self.degree,
            knots: k,
            control_points: p,
            weights: w,
        })
    }

    /// Remove the knot at `knot_index` up to `times` times.
    ///
    /// Returns `(new_curve, actual_count)` where `actual_count` is the number of
    /// times removal succeeded (may be less than `times` if tolerance is exceeded
    /// or the knot multiplicity is lower than requested).
    pub fn remove_knot_times(
        &self,
        knot_index: usize,
        times: usize,
        tolerance: f64,
    ) -> (NurbsCurve, usize) {
        let (k, p, w, count) = nurbs_ops::knot_removal::remove_knot_times_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            knot_index,
            times,
            tolerance,
        );
        (
            NurbsCurve {
                degree: self.degree,
                knots: k,
                control_points: p,
                weights: w,
            },
            count,
        )
    }

    /// Insert multiple sorted knots. The curve shape is mathematically unchanged.
    pub fn insert_knots(&self, knots_to_insert: &[f64]) -> Result<NurbsCurve, CoreError> {
        let (t_min, t_max) = self.parameter_range();
        for &t in knots_to_insert {
            if t < t_min || t > t_max {
                return Err(CoreError::PreconditionFailed(format!(
                    "knot {} outside parameter range [{}, {}]",
                    t, t_min, t_max
                )));
            }
        }
        let (knots, pts, wts) = nurbs_ops::knot_insertion::insert_knots_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            knots_to_insert,
        );
        Ok(NurbsCurve {
            degree: self.degree,
            knots,
            control_points: pts,
            weights: wts,
        })
    }

    /// Elevate degree by 1. Shape unchanged.
    pub fn elevate_degree(&self) -> NurbsCurve {
        let (deg, knots, pts, wts) = nurbs_ops::degree_elevation::elevate_degree_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
        );
        NurbsCurve {
            degree: deg,
            knots,
            control_points: pts,
            weights: wts,
        }
    }

    /// Elevate degree to `target`. Returns `Err` if `target <= current degree`.
    pub fn elevate_degree_to(&self, target: usize) -> Result<NurbsCurve, CoreError> {
        if target <= self.degree {
            return Err(CoreError::PreconditionFailed(format!(
                "target degree {} must be > current degree {}",
                target, self.degree
            )));
        }
        let mut curve = self.clone();
        for _ in 0..(target - self.degree) {
            curve = curve.elevate_degree();
        }
        Ok(curve)
    }

    /// Split the curve at parameter `t` into two independent curves.
    ///
    /// Returns `(left, right)` where `left` covers `[t_start, t]` and `right` covers `[t, t_end]`.
    /// The parameter `t` must be strictly inside the parameter range.
    pub fn split_at(&self, t: f64) -> Result<(NurbsCurve, NurbsCurve), CoreError> {
        let (lk, lp, lw, rk, rp, rw) = nurbs_ops::splitting::split_1d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )?;
        Ok((
            NurbsCurve {
                degree: self.degree,
                knots: lk,
                control_points: lp,
                weights: lw,
            },
            NurbsCurve {
                degree: self.degree,
                knots: rk,
                control_points: rp,
                weights: rw,
            },
        ))
    }

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
mod tests;
