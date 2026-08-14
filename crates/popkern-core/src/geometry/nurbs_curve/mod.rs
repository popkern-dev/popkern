//! NURBS (Non-Uniform Rational B-Spline) curve evaluation via de Boor's algorithm.

use crate::error::CoreError;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_ops::{self, Homogeneous, de_boor_1d_points};
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A Non-Uniform Rational B-Spline curve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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

    /// Concatenate two clamped curves of the same degree into one
    /// piecewise curve, C0 at the join.
    ///
    /// `self` covers `[a, b]`, `right` covers `[c, d]`; the result
    /// covers `[a, b + (d − c)]` with `right`'s knots SHIFTED so its
    /// start lands at `b` (the caller re-anchors wrapped windows of
    /// closed curves this way — a periodic parameterization's anchor
    /// is arbitrary). The shared control point is kept once; the join
    /// knot has multiplicity `degree` (positional continuity — the
    /// pieces meet with whatever tangent continuity they inherently
    /// have). `right`'s weights are uniformly rescaled so the shared
    /// point's weight matches (a global weight scale does not change
    /// a rational curve).
    ///
    /// # Errors
    ///
    /// [`CoreError::PreconditionFailed`] when the degrees differ or
    /// the join positions are farther apart than `tolerance`.
    pub fn concatenate_clamped(
        &self,
        right: &NurbsCurve,
        tolerance: f64,
    ) -> Result<NurbsCurve, CoreError> {
        if self.degree != right.degree {
            return Err(CoreError::PreconditionFailed(format!(
                "concatenate_clamped: degree mismatch {} vs {}",
                self.degree, right.degree
            )));
        }
        let (Some(l_end), Some(r_start)) =
            (self.control_points.last(), right.control_points.first())
        else {
            return Err(CoreError::PreconditionFailed(
                "concatenate_clamped: empty control polygon".into(),
            ));
        };
        if l_end.distance_to(r_start) > tolerance {
            return Err(CoreError::PreconditionFailed(format!(
                "concatenate_clamped: join gap {} exceeds tolerance {tolerance}",
                l_end.distance_to(r_start)
            )));
        }
        let (Some(&w_l_end), Some(&w_r_start)) = (self.weights.last(), right.weights.first())
        else {
            return Err(CoreError::PreconditionFailed(
                "concatenate_clamped: empty weights".into(),
            ));
        };
        if w_r_start.abs() < 1e-30 {
            return Err(CoreError::PreconditionFailed(
                "concatenate_clamped: zero weight at join".into(),
            ));
        }
        let deg = self.degree;
        let b = *self.knots.last().ok_or_else(|| {
            CoreError::PreconditionFailed("concatenate_clamped: empty knots".into())
        })?;
        let c = *right.knots.first().ok_or_else(|| {
            CoreError::PreconditionFailed("concatenate_clamped: empty knots".into())
        })?;
        let shift = b - c;
        let w_scale = w_l_end / w_r_start;

        // Left knots minus ONE trailing clamp knot (join multiplicity
        // becomes `degree`), then right's knots past its leading
        // clamp, shifted.
        let mut knots = self.knots[..self.knots.len() - 1].to_vec();
        knots.extend(right.knots[deg + 1..].iter().map(|k| k + shift));
        let mut control_points = self.control_points.clone();
        control_points.extend(right.control_points[1..].iter().copied());
        let mut weights = self.weights.clone();
        weights.extend(right.weights[1..].iter().map(|w| w * w_scale));

        Ok(NurbsCurve {
            degree: deg,
            knots,
            control_points,
            weights,
        })
    }

    /// Evaluate the curve at parameter t using de Boor's algorithm in homogeneous coordinates.
    ///
    /// Allocation-free hot path: [`de_boor_1d_points`] handles the
    /// homogeneous conversion inside a fixed-size stack scratch, so
    /// neither this function nor the downstream recursion touches
    /// the allocator. Measured impact on the §18 overlap-detect
    /// profile: every `point_at` call that used to build a temporary
    /// `Vec<Homogeneous>` now runs allocation-free.
    fn de_boor(&self, t: f64) -> Homogeneous {
        de_boor_1d_points(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }

    /// Exact rational first derivative (PK-285): basis values +
    /// derivatives at the span, one local contraction, then the
    /// quotient rule C′ = (A′ − w′·C)/w on the homogeneous forms.
    /// Replaces the 1e-7 central difference (2 evals + ~1e-9
    /// truncation).
    fn derivative_analytic_order1(&self, t: f64) -> Vec3 {
        let span = nurbs_ops::find_span(self.degree, &self.knots, self.control_points.len(), t);
        let (n, d) =
            nurbs_ops::basis::basis_funs_and_first_derivs(self.degree, &self.knots, span, t);
        let mut a = [0.0_f64; 4];
        let mut a_t = [0.0_f64; 4];
        for r in 0..=self.degree {
            let i = span - self.degree + r;
            let w = self.weights[i];
            let cp = &self.control_points[i];
            let hw = [w * cp.x, w * cp.y, w * cp.z, w];
            for c in 0..4 {
                a[c] += n[r] * hw[c];
                a_t[c] += d[r] * hw[c];
            }
        }
        if a[3].abs() < 1e-30 {
            return Vec3::zero();
        }
        let inv_w = 1.0 / a[3];
        let cx = a[0] * inv_w;
        let cy = a[1] * inv_w;
        let cz = a[2] * inv_w;
        Vec3::new(
            (a_t[0] - a_t[3] * cx) * inv_w,
            (a_t[1] - a_t[3] * cy) * inv_w,
            (a_t[2] - a_t[3] * cz) * inv_w,
        )
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

        if order == 1 {
            return self.derivative_analytic_order1(t);
        }

        // Higher orders: finite difference of (order-1) derivatives —
        // now a difference of EXACT first derivatives at order 2.
        let h = 1e-7;
        let (t_min, t_max) = self.parameter_range();

        {
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

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        use crate::geometry::curve::closest_point_on_curve;
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(
        &self,
        t_range: (f64, f64),
    ) -> crate::geometry::intersection::surface_subdivision::Aabb {
        use crate::geometry::curve::curve_bounding_box_by_sampling;
        curve_bounding_box_by_sampling(self, t_range)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod derivative_tests;
