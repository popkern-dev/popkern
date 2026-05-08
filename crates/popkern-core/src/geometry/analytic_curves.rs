//! CurveEval implementations for analytic curve types (Line, Circle, Arc).

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::{Point3, Vec3};

#[path = "analytic_curve_primitives.rs"]
mod analytic_curve_primitives;

#[cfg(test)]
pub(crate) use analytic_curve_primitives::arbitrary_perpendicular;
pub(crate) use analytic_curve_primitives::circle_basis;
use analytic_curve_primitives::{Arc, Circle, Ellipse, Hyperbola, Line, Parabola};

// ---------------------------------------------------------------------------
// CurveGeometry dispatch
// ---------------------------------------------------------------------------

impl CurveEval for CurveGeometry {
    fn point_at(&self, t: f64) -> Point3 {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .point_at(t),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).point_at(t),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle).point_at(t),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle).point_at(t),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end).point_at(t),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end).point_at(t),
            CurveGeometry::Nurbs(nurbs) => nurbs.point_at(t),
        }
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .derivative_at(t, order),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).derivative_at(t, order),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle)
                .derivative_at(t, order),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle)
                .derivative_at(t, order),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end)
                .derivative_at(t, order),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end)
                .derivative_at(t, order),
            CurveGeometry::Nurbs(nurbs) => nurbs.derivative_at(t, order),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .tangent_at(t),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).tangent_at(t),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle).tangent_at(t),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle).tangent_at(t),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end).tangent_at(t),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end).tangent_at(t),
            CurveGeometry::Nurbs(nurbs) => nurbs.tangent_at(t),
        }
    }

    fn parameter_range(&self) -> (f64, f64) {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .parameter_range(),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).parameter_range(),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle).parameter_range(),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle)
                .parameter_range(),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end).parameter_range(),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end)
                .parameter_range(),
            CurveGeometry::Nurbs(nurbs) => nurbs.parameter_range(),
        }
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .length_between(t0, t1),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).length_between(t0, t1),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => {
                Arc::new(*center, *normal, *radius, *start_angle, *end_angle).length_between(t0, t1)
            }
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle)
                .length_between(t0, t1),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => {
                Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end).length_between(t0, t1)
            }
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end)
                .length_between(t0, t1),
            CurveGeometry::Nurbs(nurbs) => nurbs.length_between(t0, t1),
        }
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .closest_point(point),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).closest_point(point),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle).closest_point(point),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle)
                .closest_point(point),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end).closest_point(point),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end)
                .closest_point(point),
            CurveGeometry::Nurbs(nurbs) => nurbs.closest_point(point),
        }
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        match self {
            CurveGeometry::Line { origin, direction } => Line {
                origin: *origin,
                direction: *direction,
            }
            .curve_bounding_box(t_range),
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => Circle::new(*center, *normal, *radius).curve_bounding_box(t_range),
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => Arc::new(*center, *normal, *radius, *start_angle, *end_angle)
                .curve_bounding_box(t_range),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                start_angle,
                end_angle,
            } => Ellipse::new(*center, *u_basis, *v_basis, *start_angle, *end_angle)
                .curve_bounding_box(t_range),
            CurveGeometry::Parabola {
                vertex,
                linear,
                quadratic,
                t_start,
                t_end,
            } => Parabola::new(*vertex, *linear, *quadratic, *t_start, *t_end)
                .curve_bounding_box(t_range),
            CurveGeometry::Hyperbola {
                center,
                semi_major,
                semi_minor,
                t_start,
                t_end,
            } => Hyperbola::new(*center, *semi_major, *semi_minor, *t_start, *t_end)
                .curve_bounding_box(t_range),
            CurveGeometry::Nurbs(nurbs) => nurbs.curve_bounding_box(t_range),
        }
    }
}

#[cfg(test)]
#[path = "analytic_curves_tests.rs"]
mod tests;
