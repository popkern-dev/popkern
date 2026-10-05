//! Curve splitting for sub-edges.
//!
//! When a graph edge is split at a point, the geometric curve must also be
//! subdivided. This module provides [`split_curve`] to divide a
//! [`CurveOnEdge`] into two sub-curves at a given parameter.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::geometry::{
    Point3, Vec3,
    curve::{CurveGeometry, CurveOnEdge},
};
use popkern_math::PortableMath;

/// Split a [`CurveOnEdge`] at parameter `t_split` into two sub-curves.
///
/// Returns `(first_half, second_half)` where the first half covers the
/// portion before the split and the second half covers the portion after.
pub(crate) fn split_curve(
    curve: &CurveOnEdge,
    t_split: f64,
    split_point: Point3,
) -> (CurveOnEdge, CurveOnEdge) {
    match &curve.geometry {
        CurveGeometry::Line { origin, direction } => split_line(
            *origin,
            *direction,
            curve.t_start,
            curve.t_end,
            t_split,
            split_point,
        ),
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => split_circle(
            *center,
            *normal,
            *radius,
            curve.t_start,
            curve.t_end,
            split_point,
        ),
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => split_arc(
            *center,
            *normal,
            *radius,
            *start_angle,
            *end_angle,
            curve.t_start,
            curve.t_end,
            split_point,
        ),
        CurveGeometry::Ellipse { .. }
        | CurveGeometry::Parabola { .. }
        | CurveGeometry::Hyperbola { .. }
        | CurveGeometry::Nurbs(_) => split_nurbs(curve, t_split),
    }
}

/// Split a line segment at t_split.
///
/// Produces two line segments: origin→split_point and split_point→end,
/// each with parameter range [0, 1].
fn split_line(
    origin: Point3,
    direction: Vec3,
    t_start: f64,
    t_end: f64,
    _t_split: f64,
    split_point: Point3,
) -> (CurveOnEdge, CurveOnEdge) {
    // End point of the original line segment.
    let end_point = origin + direction * t_end;
    let start_point = origin + direction * t_start;

    let first = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: start_point,
            direction: split_point - start_point,
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let second = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: split_point,
            direction: end_point - split_point,
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    (first, second)
}

/// Compute the angle of a point on a circle in the circle's local frame.
fn angle_in_circle_frame(point: Point3, center: Point3, normal: Vec3) -> f64 {
    let (u, v) = circle_basis(normal);
    let d = point - center;
    d.dot(&v).patan2(d.dot(&u))
}

/// Compute an orthonormal basis (u, v) for a plane with the given normal.
///
/// Mirrors the logic in `popkern_core::geometry::analytic_curves::circle_basis`
/// so that angles are consistent.
fn circle_basis(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = normal.cross(&axis);
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

/// Split a full circle at a point, producing two arcs.
fn split_circle(
    center: Point3,
    normal: Vec3,
    radius: f64,
    t_start: f64,
    t_end: f64,
    split_point: Point3,
) -> (CurveOnEdge, CurveOnEdge) {
    let split_angle = angle_in_circle_frame(split_point, center, normal);
    // Normalize into [0, 2pi).
    let split_angle = split_angle.rem_euclid(std::f64::consts::TAU);

    let first = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle: t_start,
            end_angle: split_angle,
        },
        t_start,
        t_end: split_angle,
    };
    let second = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle: split_angle,
            end_angle: t_end,
        },
        t_start: split_angle,
        t_end,
    };
    (first, second)
}

/// Split an arc at a point, producing two smaller arcs.
#[allow(clippy::too_many_arguments)]
fn split_arc(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    t_start: f64,
    t_end: f64,
    split_point: Point3,
) -> (CurveOnEdge, CurveOnEdge) {
    let split_angle = angle_in_circle_frame(split_point, center, normal);
    let split_angle = split_angle.rem_euclid(std::f64::consts::TAU);

    let first = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle: split_angle,
        },
        t_start,
        t_end: split_angle,
    };
    let second = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle: split_angle,
            end_angle,
        },
        t_start: split_angle,
        t_end,
    };
    (first, second)
}

/// Split a NURBS curve by clamping parameter ranges.
fn split_nurbs(curve: &CurveOnEdge, t_split: f64) -> (CurveOnEdge, CurveOnEdge) {
    let first = CurveOnEdge {
        geometry: curve.geometry.clone(),
        t_start: curve.t_start,
        t_end: t_split,
    };
    let second = CurveOnEdge {
        geometry: curve.geometry.clone(),
        t_start: t_split,
        t_end: curve.t_end,
    };
    (first, second)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

    #[test]
    fn split_line_at_midpoint() {
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(2.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        };
        let split_pt = Point3::new(1.0, 0.0, 0.0);
        let (first, second) = split_curve(&curve, 0.5, split_pt);

        // First half: (0,0,0) → (1,0,0)
        match &first.geometry {
            CurveGeometry::Line { origin, direction } => {
                assert!(origin.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-10);
                assert!((direction.x - 1.0).abs() < 1e-10);
            }
            _ => panic!("Expected Line"),
        }

        // Second half: (1,0,0) → (2,0,0)
        match &second.geometry {
            CurveGeometry::Line { origin, direction } => {
                assert!(origin.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-10);
                assert!((direction.x - 1.0).abs() < 1e-10);
            }
            _ => panic!("Expected Line"),
        }
    }

    #[test]
    fn split_circle_produces_two_arcs() {
        use std::f64::consts::TAU;

        let curve = CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
            },
            t_start: 0.0,
            t_end: TAU,
        };
        // Split at (-1, 0, 0). The exact angle depends on the local frame
        // but the two arcs must partition the original range.
        let split_pt = Point3::new(-1.0, 0.0, 0.0);
        let (first, second) = split_curve(&curve, std::f64::consts::PI, split_pt);

        match (&first.geometry, &second.geometry) {
            (
                CurveGeometry::Arc {
                    start_angle: s1,
                    end_angle: e1,
                    ..
                },
                CurveGeometry::Arc {
                    start_angle: s2,
                    end_angle: e2,
                    ..
                },
            ) => {
                // First arc starts at 0.
                assert!(s1.abs() < 1e-10, "first start should be 0, got {s1}");
                // Arcs must be contiguous: first.end == second.start.
                assert!((e1 - s2).abs() < 1e-10, "arcs not contiguous: {e1} vs {s2}");
                // Second arc ends at TAU.
                assert!(
                    (e2 - TAU).abs() < 1e-10,
                    "second end should be TAU, got {e2}"
                );
                // Split angle must be positive and less than TAU.
                assert!(*e1 > 0.0 && *e1 < TAU, "split angle out of range: {e1}");
            }
            _ => panic!("Expected two Arcs"),
        }
    }

    #[test]
    fn split_nurbs_at_midpoint() {
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
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start: 0.0,
            t_end: 1.0,
        };
        let split_pt = Point3::new(0.5, 0.5, 0.0);
        let (first, second) = split_curve(&curve, 0.5, split_pt);

        assert!((first.t_start - 0.0).abs() < 1e-10);
        assert!((first.t_end - 0.5).abs() < 1e-10);
        assert!((second.t_start - 0.5).abs() < 1e-10);
        assert!((second.t_end - 1.0).abs() < 1e-10);
    }

    #[test]
    fn split_arc_preserves_bounds() {
        use std::f64::consts::PI;

        let curve = CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
                start_angle: 0.0,
                end_angle: PI,
            },
            t_start: 0.0,
            t_end: PI,
        };
        // Split at (-1, 0, 0). The exact split angle depends on the basis,
        // but bounds and continuity must hold.
        let split_pt = Point3::new(-1.0, 0.0, 0.0);
        let (first, second) = split_curve(&curve, PI / 2.0, split_pt);

        match (&first.geometry, &second.geometry) {
            (
                CurveGeometry::Arc {
                    start_angle: s1,
                    end_angle: e1,
                    ..
                },
                CurveGeometry::Arc {
                    start_angle: s2,
                    end_angle: e2,
                    ..
                },
            ) => {
                assert!(s1.abs() < 1e-10, "first start should be 0, got {s1}");
                assert!((e1 - s2).abs() < 1e-10, "arcs not contiguous: {e1} vs {s2}");
                assert!((e2 - PI).abs() < 1e-10, "second end should be PI, got {e2}");
            }
            _ => panic!("Expected two Arcs"),
        }
    }
}
