//! CurveEval implementations for analytic curve types (Line, Circle, Arc).

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};

/// Returns an arbitrary vector perpendicular to `normal` (assumed non-zero).
fn arbitrary_perpendicular(normal: Vec3) -> Vec3 {
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
fn circle_basis(normal: Vec3) -> (Vec3, Vec3) {
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
        // nth derivative of (r*cos(t), r*sin(t)) follows a 4-cycle of sin/cos.
        match order % 4 {
            0 => {
                // 0th derivative = position offset from center
                self.u * (self.radius * t.cos()) + self.v * (self.radius * t.sin())
            }
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
            CurveGeometry::Nurbs(nurbs) => nurbs.length_between(t0, t1),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line_point_at() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let p = line.point_at(0.5);
        assert!((p.x - 0.5).abs() < 1e-10);
        assert!((p.y).abs() < 1e-10);
    }

    #[test]
    fn line_derivative() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(2.0, 0.0, 0.0),
        };
        let d = line.derivative_at(0.5, 1);
        assert!((d.x - 2.0).abs() < 1e-10);
        let d2 = line.derivative_at(0.5, 2);
        assert!((d2.length()).abs() < 1e-10);
    }

    #[test]
    fn line_length() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        assert!((line.length_between(0.0, 1.0) - 1.0).abs() < 1e-10);
        assert!((line.length_between(0.0, 0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn circle_point_at_start() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let p = circle.point_at(0.0);
        // At t=0, should be at distance 1 from center in the XY plane
        let dist = ((p.x - 0.0).powi(2) + (p.y - 0.0).powi(2)).sqrt();
        assert!((dist - 1.0).abs() < 1e-10);
    }

    #[test]
    fn circle_point_at_quarter() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let p = circle.point_at(FRAC_PI_2);
        let dist = ((p.x).powi(2) + (p.y).powi(2)).sqrt();
        assert!((dist - 1.0).abs() < 1e-10);
    }

    #[test]
    fn circle_length_full() {
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let len = circle.length_between(0.0, 2.0 * PI);
        assert!((len - 4.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn arc_length() {
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let len = arc.length_between(0.0, FRAC_PI_2);
        assert!((len - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn line_tangent_is_normalized() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(3.0, 4.0, 0.0),
        };
        let t = line.tangent_at(0.0);
        assert!((t.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn line_parameter_range() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let (t0, t1) = line.parameter_range();
        assert!(t0.is_infinite() && t0 < 0.0);
        assert!(t1.is_infinite() && t1 > 0.0);
    }

    #[test]
    fn circle_all_curve_eval_methods() {
        let c = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        // derivative_at all 4 branches (order % 4)
        assert!((c.derivative_at(0.0, 0).length() - 2.0).abs() < 1e-10); // order 0
        assert!((c.derivative_at(0.0, 1).length() - 2.0).abs() < 1e-10); // order 1
        assert!((c.derivative_at(0.0, 2).length() - 2.0).abs() < 1e-10); // order 2
        assert!((c.derivative_at(FRAC_PI_2, 3).length() - 2.0).abs() < 1e-10); // order 3
        // tangent_at
        assert!((c.tangent_at(0.0).length() - 1.0).abs() < 1e-10);
        // parameter_range
        let (t0, t1) = c.parameter_range();
        assert!((t0).abs() < 1e-10);
        assert!((t1 - std::f64::consts::TAU).abs() < 1e-10);
        // length_between
        assert!((c.length_between(0.0, FRAC_PI_2) - FRAC_PI_2 * 2.0).abs() < 1e-10);
    }

    #[test]
    fn arc_all_curve_eval_methods() {
        let a = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.5,
            end_angle: 2.0,
        };
        // point_at
        let p = a.point_at(0.0);
        assert!(((p.x).powi(2) + (p.y).powi(2)).sqrt() - 1.0 < 1e-10);
        // derivative_at all 4 branches
        assert!((a.derivative_at(0.0, 0).length() - 1.0).abs() < 1e-10);
        assert!((a.derivative_at(0.0, 1).length() - 1.0).abs() < 1e-10);
        assert!((a.derivative_at(0.0, 2).length() - 1.0).abs() < 1e-10);
        assert!((a.derivative_at(FRAC_PI_2, 3).length() - 1.0).abs() < 1e-10);
        // tangent_at
        assert!((a.tangent_at(FRAC_PI_2).length() - 1.0).abs() < 1e-10);
        // parameter_range
        let (t0, t1) = a.parameter_range();
        assert!((t0 - 0.5).abs() < 1e-10);
        assert!((t1 - 2.0).abs() < 1e-10);
        // length_between
        assert!((a.length_between(0.0, PI) - PI).abs() < 1e-10);
    }

    #[test]
    fn nurbs_dispatch_all_methods() {
        use crate::geometry::nurbs_curve::NurbsCurve;
        let nurbs = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        assert!((nurbs.point_at(0.5).x - 0.5).abs() < 1e-6);
        assert!(nurbs.derivative_at(0.5, 1).length() > 0.5);
        assert!((nurbs.tangent_at(0.5).length() - 1.0).abs() < 1e-6);
        let (t0, t1) = nurbs.parameter_range();
        assert!((t0).abs() < 1e-10 && (t1 - 1.0).abs() < 1e-10);
        assert!(nurbs.length_between(0.0, 1.0) > 0.9);
    }

    #[test]
    fn arbitrary_perpendicular_both_branches() {
        // x < 0.9 => picks X axis
        let perp1 = arbitrary_perpendicular(Vec3::new(0.0, 0.0, 1.0));
        assert!(perp1.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() < 1e-10);
        // x >= 0.9 => picks Y axis
        let perp2 = arbitrary_perpendicular(Vec3::new(1.0, 0.0, 0.0));
        assert!(perp2.dot(&Vec3::new(1.0, 0.0, 0.0)).abs() < 1e-10);
    }
}
