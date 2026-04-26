//! Analytic-to-NURBS conversion helpers for [`CurveGeometry`].
//!
//! Extracted from `curve.rs` to keep that file under the 500-line
//! structural cap. Builds degree-2 rational NURBS approximations of
//! analytic primitives (full circles, arcs, full ellipses, elliptical
//! arcs) using the standard 9-control-point construction with weight
//! `1/√2` on the alternating tangent-intersection points.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

/// Build a full-circle NURBS (degree 2, 9 control points, 4 quarter arcs).
pub(super) fn circle_to_nurbs(center: Point3, normal: Vec3, radius: f64) -> NurbsCurve {
    use crate::geometry::analytic_curves::circle_basis;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    let (u, v) = circle_basis(normal);
    let w = 1.0_f64 / 2.0_f64.sqrt(); // weight for midpoints

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

/// Build a degree-2 rational NURBS for a full ellipse using
/// `(u_basis, v_basis)` as the two semi-axis vectors. Same 9-CP /
/// w=1/√2 structure as the circle, just with `radius·u_axis` and
/// `radius·v_axis` replaced by `u_basis` and `v_basis` directly.
pub(super) fn ellipse_full_to_nurbs(center: Point3, u_basis: Vec3, v_basis: Vec3) -> NurbsCurve {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let w = 1.0_f64 / 2.0_f64.sqrt();
    let pts = [
        center + u_basis,           // angle 0
        center + u_basis + v_basis, // tangent intersection at 0
        center + v_basis,           // angle π/2
        center + v_basis - u_basis, // tangent intersection
        center - u_basis,           // angle π
        center - u_basis - v_basis,
        center - v_basis, // angle 3π/2
        center + u_basis - v_basis,
        center + u_basis, // angle 2π
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
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
        ],
        control_points: pts.to_vec(),
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build a NURBS elliptical arc by constructing the full ellipse
/// and splitting to the arc range.
pub(super) fn ellipse_to_nurbs(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    start_angle: f64,
    end_angle: f64,
) -> NurbsCurve {
    let full = ellipse_full_to_nurbs(center, u_basis, v_basis);
    let right = if start_angle > 0.0 {
        match full.split_at(start_angle) {
            Ok((_l, r)) => r,
            Err(_) => full,
        }
    } else {
        full
    };
    if end_angle < std::f64::consts::TAU {
        match right.split_at(end_angle) {
            Ok((l, _r)) => l,
            Err(_) => right,
        }
    } else {
        right
    }
}

/// Build an exact degree-2 Bézier-equivalent NURBS for a parabolic
/// arc `vertex + t·linear + t²·quadratic` over `[t_start, t_end]`.
///
/// A parabola is exactly representable as a non-rational degree-2
/// Bézier curve. Endpoints are evaluations at `t_start` and `t_end`;
/// the middle control point sits at the intersection of the two
/// endpoint tangent lines. Equivalently, `b1 = b0 + tangent(t_start)
/// · (t_end − t_start) / 2`.
pub(super) fn parabola_to_nurbs(
    vertex: Point3,
    linear: Vec3,
    quadratic: Vec3,
    t_start: f64,
    t_end: f64,
) -> NurbsCurve {
    let span = t_end - t_start;
    let b0 = vertex + linear * t_start + quadratic * (t_start * t_start);
    let b2 = vertex + linear * t_end + quadratic * (t_end * t_end);
    let tangent_at_start = linear + quadratic * (2.0 * t_start);
    let b1 = b0 + tangent_at_start * (span * 0.5);
    NurbsCurve {
        degree: 2,
        knots: vec![t_start, t_start, t_start, t_end, t_end, t_end],
        control_points: vec![b0, b1, b2],
        weights: vec![1.0, 1.0, 1.0],
    }
}

/// Build a NURBS arc by constructing the full circle and splitting to the arc range.
pub(super) fn arc_to_nurbs(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> NurbsCurve {
    let full = circle_to_nurbs(center, normal, radius);

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
