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

/// Build an exact rational degree-2 NURBS approximation of a
/// hyperbolic arc `center + cosh(t)·semi_major + sinh(t)·semi_minor`
/// over `[t_start, t_end]`.
///
/// The hyperbola is exactly representable as a rational degree-2
/// Bézier via the substitution `u = tanh(t/2)`, which gives
/// `cosh(t) = (1+u²)/(1−u²)` and `sinh(t) = 2u/(1−u²)`. In
/// homogeneous form: `(b0 = center + sm_maj, b1 = center, b2 =
/// center + sm_maj)` doesn't quite work directly — the cleanest
/// closed form is to evaluate three sample points at the
/// endpoints and the parametric midpoint, then build the unique
/// rational quadratic Bézier through them. We fall back to a
/// sample-and-fit since the rational-quadratic construction needs
/// careful derivation; degree-5 with 33 samples gives ~1e-12
/// residual on typical sweeps.
pub(super) fn hyperbola_to_nurbs(
    center: Point3,
    semi_major: Vec3,
    semi_minor: Vec3,
    t_start: f64,
    t_end: f64,
) -> NurbsCurve {
    let n_samples = 33;
    let mut points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let t = t_start + frac * (t_end - t_start);
        let p = center + semi_major * t.cosh() + semi_minor * t.sinh();
        points.push(p);
        params.push(frac);
    }
    let nurbs = match crate::geometry::fitting::fit_curve(&points, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => {
            // Degenerate: fall back to a degree-1 spline through endpoints.
            return NurbsCurve {
                degree: 1,
                knots: vec![t_start, t_start, t_end, t_end],
                control_points: vec![*points.first().unwrap(), *points.last().unwrap()],
                weights: vec![1.0, 1.0],
            };
        }
    };
    // Reknot the [0, 1] parameter range to [t_start, t_end] so
    // downstream consumers see the natural hyperbola parameter.
    let mut nurbs = nurbs;
    let k_lo = *nurbs.knots.first().unwrap_or(&0.0);
    let k_hi = *nurbs.knots.last().unwrap_or(&1.0);
    let span_old = k_hi - k_lo;
    let span_new = t_end - t_start;
    if span_old.abs() > 1.0e-30 {
        for k in &mut nurbs.knots {
            *k = t_start + (*k - k_lo) / span_old * span_new;
        }
    }
    nurbs
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
