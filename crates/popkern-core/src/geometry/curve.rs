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

    best.map(|(t, _)| t)
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
    arc_to_nurbs, circle_to_nurbs, ellipse_to_nurbs, hyperbola_to_nurbs, parabola_to_nurbs,
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
            } => circle_to_nurbs(*center, *normal, *radius),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => arc_to_nurbs(*center, *normal, *radius, *start_angle, *end_angle),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => ellipse_to_nurbs(*center, *u_basis, *v_basis, *start_angle, *end_angle),
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

    #[test]
    fn arc_to_nurbs_quarter_circle() {
        use std::f64::consts::FRAC_PI_2;
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let nurbs = arc.to_nurbs();
        // The NURBS should represent a quarter circle; sampled points should be on unit circle
        let p = nurbs.point_at(nurbs.parameter_range().0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (r - 1.0).abs() < 1e-6,
            "Start point not on unit circle: r={}",
            r
        );
    }

    #[test]
    fn arc_to_nurbs_half_circle() {
        use std::f64::consts::PI;
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: PI,
        };
        let nurbs = arc.to_nurbs();
        // Just verify it's a valid NURBS
        assert!(nurbs.degree >= 2);
        assert!(!nurbs.control_points.is_empty());
    }

    #[test]
    fn arc_to_nurbs_from_nonzero_start() {
        use std::f64::consts::{FRAC_PI_2, PI};
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: FRAC_PI_2,
            end_angle: PI,
        };
        let nurbs = arc.to_nurbs();
        assert!(nurbs.degree >= 2);
    }

    #[test]
    fn arc_to_nurbs_full_circle_range() {
        // Arc with end_angle = TAU (full circle)
        use std::f64::consts::TAU;
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: TAU,
        };
        let nurbs = arc.to_nurbs();
        assert!(nurbs.degree >= 2);
    }

    #[test]
    fn to_nurbs_default_range() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let nurbs = line.to_nurbs();
        assert_eq!(nurbs.degree, 1);
    }

    #[test]
    fn arc_to_nurbs_in_range() {
        use std::f64::consts::FRAC_PI_2;
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let nurbs = arc.to_nurbs_in_range((0.0, FRAC_PI_2));
        assert!(nurbs.degree >= 2);
    }

    #[test]
    fn arc_to_nurbs_boundary_split_fallback() {
        // Arc with start_angle at a knot boundary that might cause split_at to fail.
        // This is a degenerate case where the split parameter is at the boundary.
        use std::f64::consts::TAU;

        // Arc starting near TAU boundary
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: TAU - 0.001, // very close to TAU
            end_angle: TAU,
        };
        // Should not panic
        let nurbs = arc.to_nurbs();
        assert!(nurbs.degree >= 2);
    }
}
