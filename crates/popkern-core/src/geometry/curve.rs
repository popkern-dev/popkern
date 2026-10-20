//! Curve geometry types, evaluation trait, and topology binding.

use crate::geometry::intersection::surface_subdivision::Aabb;
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

    /// Find the parameter t of the closest point on the curve to `point`.
    ///
    /// Returns `None` if projection fails to converge.
    fn closest_point(&self, point: &Point3) -> Option<f64>;

    /// Compute an axis-aligned bounding box over the given parameter range.
    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb;
}

/// Newton-based closest-point projection for any `CurveEval`.
///
/// Multi-start from 5 evenly-spaced candidates within the parameter range.
pub fn closest_point_on_curve(curve: &dyn CurveEval, point: &Point3) -> Option<f64> {
    let (t0, t1) = curve.parameter_range();
    let t_lo = t0.max(-100.0);
    let t_hi = t1.min(100.0);

    let n_starts = 5;
    let max_iter = 20;
    let step_tol = 1e-12;

    let mut best: Option<(f64, f64)> = None; // (t, dist)

    for i in 0..n_starts {
        let mut t = t_lo + (t_hi - t_lo) * (i as f64) / ((n_starts - 1) as f64);

        for _ in 0..max_iter {
            let p = curve.point_at(t);
            let diff = p - *point;
            let d1 = curve.derivative_at(t, 1);

            // f(t) = diff · d1 = 0  =>  Newton: dt = -f / f'
            let f = diff.dot(&d1);
            let d2 = curve.derivative_at(t, 2);
            let fp = d1.dot(&d1) + diff.dot(&d2);

            if fp.abs() < 1e-30 {
                break;
            }

            let dt = -f / fp;
            t += dt;

            if dt.abs() < step_tol {
                break;
            }
        }

        let dist = curve.point_at(t).distance_to(point);
        match &best {
            Some((_, bd)) if dist >= *bd => {}
            _ => best = Some((t, dist)),
        }
    }

    // A closed curve's Newton can settle outside the curve's range: an
    // analytic circle on another period branch (the same point at
    // t ± period), a clamped NURBS circle past its end (evaluated at
    // the end). Callers clamp the result into the curve's window, and
    // an off-branch t clamped to a window END misreads the point —
    // measured, PK-462: a rim circle's closest t came back as −6.672 =
    // −0.389 − 2π, clamped to −π/2, and a point ON the rim read 3.3e-3
    // away. Return an in-range parameter that evaluates to the SAME
    // point — wrapped when the parameterization repeats, clamped when
    // evaluation clamps — and the raw one when neither does.
    best.map(|(t, _)| {
        if (t0..=t1).contains(&t) || !(t0.is_finite() && t1.is_finite() && t1 > t0) {
            return t;
        }
        let at = curve.point_at(t);
        let same = |c: f64| {
            curve.point_at(c).distance_to(&at)
                <= 1e-12 * (1.0 + at.distance_to(&curve.point_at(t0)))
        };
        let wrapped = t0 + (t - t0).rem_euclid(t1 - t0);
        if same(wrapped) {
            return wrapped;
        }
        let clamped = t.clamp(t0, t1);
        if same(clamped) { clamped } else { t }
    })
}

/// The closest point over the window `[t0, t1]` only (an edge's own
/// window, which may lie on another period than the curve's declared
/// range): the best of dense samples, polished by Newton steps clamped
/// to the window.
pub fn closest_point_in_window(curve: &dyn CurveEval, point: &Point3, t0: f64, t1: f64) -> f64 {
    let (t0, t1) = (t0.min(t1), t0.max(t1));
    let samples = 64;
    let dist = |t: f64| curve.point_at(t).distance_to(point);
    let polish = |mut t: f64| -> f64 {
        for _ in 0..30 {
            let diff = curve.point_at(t) - *point;
            let d1 = curve.derivative_at(t, 1);
            let fp = d1.dot(&d1) + diff.dot(&curve.derivative_at(t, 2));
            if fp.abs() < 1e-30 {
                break;
            }
            let next = (t - diff.dot(&d1) / fp).clamp(t0, t1);
            if dist(next) > dist(t) {
                break;
            }
            let done = (next - t).abs() < 1e-12;
            t = next;
            if done {
                break;
            }
        }
        t
    };
    // Seeds: the three best samples and both window ends. On a closed
    // edge a point just past the start ties with the end, whose Newton
    // step clamps there: the start's does not.
    let mut seeds: Vec<f64> = (0..=samples)
        .map(|k| t0 + (t1 - t0) * k as f64 / samples as f64)
        .collect();
    seeds.sort_by(|a, b| dist(*a).total_cmp(&dist(*b)));
    seeds.truncate(3);
    seeds.extend([t0, t1]);
    seeds
        .into_iter()
        .map(polish)
        .min_by(|a, b| dist(*a).total_cmp(&dist(*b)))
        .unwrap_or(t0)
}

/// Compute AABB by sampling a curve at evenly-spaced parameters.
pub fn curve_bounding_box_by_sampling(curve: &dyn CurveEval, t_range: (f64, f64)) -> Aabb {
    let n = 16;
    let mut points = Vec::with_capacity(n + 1);

    for i in 0..=n {
        let t = t_range.0 + (t_range.1 - t_range.0) * (i as f64) / (n as f64);
        points.push(curve.point_at(t));
    }

    Aabb::from_points(&points).unwrap_or(Aabb {
        min: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
        max: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
    })
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
    /// An elliptical arc. Evaluates as
    /// `center + cos(t) · u_basis + sin(t) · v_basis`. The two
    /// basis vectors carry both the major/minor axis directions
    /// and their lengths (the semi-radii). They must be coplanar
    /// (i.e. their cross product points along the ellipse's
    /// normal) but need not be orthogonal — the parameterisation
    /// still produces an ellipse for any non-degenerate basis pair.
    Ellipse {
        /// Centre point.
        center: Point3,
        /// First semi-axis vector (length = first semi-radius).
        u_basis: Vec3,
        /// Second semi-axis vector (length = second semi-radius).
        v_basis: Vec3,
        /// Sweep start angle (radians).
        start_angle: f64,
        /// Sweep end angle (radians).
        end_angle: f64,
    },
    /// A parabolic arc parameterised as
    /// `vertex + t · linear + t² · quadratic`.
    ///
    /// Canonical form for plane-cone intersections where the cutting
    /// plane is parallel to exactly one cone generator. The `linear`
    /// vector is the tangent at `t = 0`; the `quadratic` vector is
    /// half the second derivative and points along the parabola's
    /// axis of symmetry. Both vectors lie in the parabola's plane;
    /// they need not be orthogonal but must be linearly independent.
    Parabola {
        /// Vertex / base point: `point_at(0)`.
        vertex: Point3,
        /// Coefficient of `t`: tangent vector at `t = 0`. Length
        /// determines parameterisation speed.
        linear: Vec3,
        /// Coefficient of `t²`: defines the axis of symmetry direction
        /// and the curvature (related to focal parameter `p` by
        /// `|quadratic| = 1 / (4·p)` for the canonical form).
        quadratic: Vec3,
        /// Sweep start parameter.
        t_start: f64,
        /// Sweep end parameter.
        t_end: f64,
    },
    /// A hyperbolic arc (one branch) parameterised as
    /// `center + cosh(t) · semi_major + sinh(t) · semi_minor`.
    ///
    /// Canonical form for plane-cone intersections where the cutting
    /// plane is steeper than the cone's generator angle. The
    /// `semi_major` vector points from the centre toward the branch's
    /// vertex (nearest-approach point at `t = 0`); `semi_minor` is
    /// perpendicular within the hyperbola plane and scales the
    /// `sinh` arm. To select the opposite branch, negate `semi_major`.
    /// Both vectors must be linearly independent.
    Hyperbola {
        /// Centre point — lies on the line of symmetry between the
        /// two branches.
        center: Point3,
        /// Real semi-axis vector. `cosh(t)·semi_major` is the
        /// "horizontal" component. `|semi_major|` = real semi-axis
        /// `a`; the curve approaches `center + cosh(t)·semi_major`
        /// asymptotically as `|t|→∞`.
        semi_major: Vec3,
        /// Imaginary semi-axis vector. `sinh(t)·semi_minor` is the
        /// "vertical" component. `|semi_minor|` = imaginary semi-axis
        /// `b`.
        semi_minor: Vec3,
        /// Sweep start parameter.
        t_start: f64,
        /// Sweep end parameter.
        t_end: f64,
    },
    /// A NURBS (Non-Uniform Rational B-Spline) curve.
    Nurbs(NurbsCurve),
}

use crate::geometry::curve_to_nurbs::{
    arc_to_nurbs, circle_to_nurbs, conic_arc_to_nurbs, ellipse_to_nurbs, hyperbola_to_nurbs,
    parabola_to_nurbs,
};

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
            } => match angular_sub_range(t_range, (0.0, std::f64::consts::TAU)) {
                Some((a0, a1)) => arc_to_nurbs(*center, *normal, *radius, a0, a1),
                None => circle_to_nurbs(*center, *normal, *radius),
            },
            // The parameter of a circle, arc or ellipse is its angle, so
            // the requested range is the sub-arc an edge actually uses
            // (PK-820: a moved countersink's planar pcurves spanned the
            // whole declared arc, 4.9 off their edge).
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => {
                let (a0, a1) = angular_sub_range(t_range, (*start_angle, *end_angle))
                    .unwrap_or((*start_angle, *end_angle));
                arc_to_nurbs(*center, *normal, *radius, a0, a1)
            }
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => {
                // A sub-range takes the exact segmentation: the
                // shoulder construction splits at parameter = angle,
                // which a rational quadratic meets only at quarter turns.
                match angular_sub_range(t_range, (*start_angle, *end_angle)) {
                    Some((a0, a1)) => conic_arc_to_nurbs(*center, *u_basis, *v_basis, a0, a1),
                    None => ellipse_to_nurbs(*center, *u_basis, *v_basis, *start_angle, *end_angle),
                }
            }
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => parabola_to_nurbs(*vertex, *linear, *quadratic, *t_start, *t_end),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => hyperbola_to_nurbs(*center, *semi_major, *semi_minor, *t_start, *t_end),
            CurveGeometry::Nurbs(nurbs) => {
                nurbs_curve_in_range(nurbs, t_range).unwrap_or_else(|| nurbs.clone())
            }
        }
    }

    /// Convert with the default range from [`CurveEval::parameter_range()`].
    pub fn to_nurbs(&self) -> NurbsCurve {
        let range = self.parameter_range();
        self.to_nurbs_in_range(range)
    }
}

/// The ordered angular range `t_range` when it is a proper, non-empty
/// part of `declared` (a full turn or less), else `None`: the whole
/// declared curve is then converted as before.
fn angular_sub_range(t_range: (f64, f64), declared: (f64, f64)) -> Option<(f64, f64)> {
    let (a0, a1) = (t_range.0.min(t_range.1), t_range.0.max(t_range.1));
    let (d0, d1) = (declared.0.min(declared.1), declared.0.max(declared.1));
    let width = a1 - a0;
    let proper = width < (d1 - d0).min(std::f64::consts::TAU) - 1e-12;
    (a0.is_finite() && a1.is_finite() && width > 1e-12 && proper).then_some((a0, a1))
}

fn nurbs_curve_in_range(nurbs: &NurbsCurve, t_range: (f64, f64)) -> Option<NurbsCurve> {
    let (full_start, full_end) = nurbs.parameter_range();
    let (slice_start, slice_end, reversed) = if t_range.0 <= t_range.1 {
        (t_range.0, t_range.1, false)
    } else {
        (t_range.1, t_range.0, true)
    };
    let eps = 1.0e-12 * (full_end - full_start).abs().max(1.0);
    if slice_start < full_start - eps
        || slice_end > full_end + eps
        || slice_end <= slice_start + eps
    {
        return None;
    }

    let at_start = (slice_start - full_start).abs() <= eps;
    let at_end = (slice_end - full_end).abs() <= eps;
    let mut restricted = nurbs.clone();
    if !at_start {
        restricted = restricted.split_at(slice_start).ok()?.1;
    }
    if !at_end {
        let (current_start, _) = restricted.parameter_range();
        let split_at = if at_start {
            slice_end
        } else {
            current_start + (slice_end - slice_start)
        };
        restricted = restricted.split_at(split_at).ok()?.0;
    }

    Some(if reversed {
        reverse_nurbs_curve(&restricted)
    } else {
        restricted
    })
}

fn reverse_nurbs_curve(nurbs: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = nurbs.parameter_range();
    let mut knots = nurbs
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect::<Vec<_>>();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = nurbs.control_points.clone();
    control_points.reverse();
    let mut weights = nurbs.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: nurbs.degree,
        knots,
        control_points,
        weights,
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

impl CurveOnEdge {
    /// Construct an analytic bounded line segment from `start` to `end`.
    ///
    /// The stored curve is the infinite analytic line
    /// `start + (end - start) * t`, with the edge bounded by
    /// `t_start=0` and `t_end=1`.
    pub fn line_segment(start: Point3, end: Point3, tolerance: f64) -> Option<Self> {
        let direction = end - start;
        if direction.length() <= tolerance.max(1.0e-12) {
            return None;
        }
        Some(Self {
            geometry: CurveGeometry::Line {
                origin: start,
                direction,
            },
            t_start: 0.0,
            t_end: 1.0,
        })
    }
}

#[cfg(test)]
#[path = "curve_tests.rs"]
mod tests;
