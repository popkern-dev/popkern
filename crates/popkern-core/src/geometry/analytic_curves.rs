//! CurveEval implementations for analytic curve types (Line, Circle, Arc).

use crate::geometry::curve::{
    CurveEval, CurveGeometry, closest_point_on_curve, curve_bounding_box_by_sampling,
};
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::{Point3, Vec3};

/// Returns an arbitrary vector perpendicular to `normal` (assumed non-zero).
pub(crate) fn arbitrary_perpendicular(normal: Vec3) -> Vec3 {
    // Pick an axis not parallel to normal, then cross to get a perpendicular.
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    normal.cross(&axis)
}

/// Compute the two orthonormal basis vectors (u, v) for a circle/arc in the plane
/// defined by `normal`. `u` and `v` satisfy:
///   - u and v are unit vectors
///   - u, v are perpendicular to each other and to normal
///   - normal ≈ u × v (right-hand rule)
pub(crate) fn circle_basis(normal: Vec3) -> (Vec3, Vec3) {
    let perp = arbitrary_perpendicular(normal);
    // u = normalize(perp)
    // SAFETY: arbitrary_perpendicular guarantees a non-zero vector,
    // and normal × u is non-zero when both are non-zero and not parallel.
    // unwrap_or fallback to unit X/Y is defensive but should never trigger.
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

// ---------------------------------------------------------------------------
// Line
// ---------------------------------------------------------------------------

struct Line {
    origin: Point3,
    direction: Vec3,
}

impl CurveEval for Line {
    fn point_at(&self, t: f64) -> Point3 {
        self.origin + self.direction * t
    }

    fn derivative_at(&self, _t: f64, order: usize) -> Vec3 {
        if order == 1 {
            self.direction
        } else {
            Vec3::zero()
        }
    }

    fn tangent_at(&self, _t: f64) -> Vec3 {
        self.direction.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.direction.length() * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        // Analytic: project point onto line
        let d = self.direction;
        let len_sq = d.dot(&d);
        if len_sq < 1e-30 {
            return None;
        }
        let diff = *point - self.origin;
        Some(diff.dot(&d) / len_sq)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Circle
// ---------------------------------------------------------------------------

struct Circle {
    center: Point3,
    radius: f64,
    u: Vec3,
    v: Vec3,
}

impl Circle {
    fn new(center: Point3, normal: Vec3, radius: f64) -> Self {
        let (u, v) = circle_basis(normal);
        Self {
            center,
            radius,
            u,
            v,
        }
    }
}

impl CurveEval for Circle {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u * (self.radius * t.cos()) + self.v * (self.radius * t.sin())
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u * (self.radius * t.cos()) + self.v * (self.radius * t.sin()),
            1 => self.u * (-self.radius * t.sin()) + self.v * (self.radius * t.cos()),
            2 => self.u * (-self.radius * t.cos()) + self.v * (-self.radius * t.sin()),
            3 => self.u * (self.radius * t.sin()) + self.v * (-self.radius * t.cos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        let d = self.derivative_at(t, 1);
        d.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (0.0, std::f64::consts::TAU)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.radius * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Arc
// ---------------------------------------------------------------------------

struct Arc {
    center: Point3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    u: Vec3,
    v: Vec3,
}

impl Arc {
    fn new(center: Point3, normal: Vec3, radius: f64, start_angle: f64, end_angle: f64) -> Self {
        let (u, v) = circle_basis(normal);
        Self {
            center,
            radius,
            start_angle,
            end_angle,
            u,
            v,
        }
    }
}

impl CurveEval for Arc {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u * (self.radius * t.cos()) + self.v * (self.radius * t.sin())
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u * (self.radius * t.cos()) + self.v * (self.radius * t.sin()),
            1 => self.u * (-self.radius * t.sin()) + self.v * (self.radius * t.cos()),
            2 => self.u * (-self.radius * t.cos()) + self.v * (-self.radius * t.sin()),
            3 => self.u * (self.radius * t.sin()) + self.v * (-self.radius * t.cos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        let d = self.derivative_at(t, 1);
        d.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.start_angle, self.end_angle)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.radius * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Ellipse — `center + cos(t) · u_basis + sin(t) · v_basis`.
// Used directly by `CurveGeometry::Ellipse` dispatch; the (u_basis,
// v_basis) form keeps both the axis directions and the semi-radii
// in one struct, so non-orthogonal bases work too (the trace is
// still an ellipse for any non-degenerate basis pair).
// ---------------------------------------------------------------------------

struct Ellipse {
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    start_angle: f64,
    end_angle: f64,
}

impl Ellipse {
    fn new(center: Point3, u_basis: Vec3, v_basis: Vec3, start_angle: f64, end_angle: f64) -> Self {
        Self {
            center,
            u_basis,
            v_basis,
            start_angle,
            end_angle,
        }
    }
}

impl CurveEval for Ellipse {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u_basis * t.cos() + self.v_basis * t.sin()
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u_basis * t.cos() + self.v_basis * t.sin(),
            1 => self.u_basis * (-t.sin()) + self.v_basis * t.cos(),
            2 => self.u_basis * (-t.cos()) + self.v_basis * (-t.sin()),
            3 => self.u_basis * t.sin() + self.v_basis * (-t.cos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        self.derivative_at(t, 1)
            .normalized()
            .unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.start_angle, self.end_angle)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        // No closed form — sample-and-sum. 32 samples are enough
        // for tessellation-grade length to ~1e-9 over a half-arc.
        let n = 32;
        let mut sum = 0.0;
        let mut prev = self.point_at(t0);
        for i in 1..=n {
            let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
            let p = self.point_at(t);
            sum += prev.distance_to(&p);
            prev = p;
        }
        sum
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

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
            CurveGeometry::Nurbs(nurbs) => nurbs.curve_bounding_box(t_range),
        }
    }
}

#[cfg(test)]
#[path = "analytic_curves_tests.rs"]
mod tests;
