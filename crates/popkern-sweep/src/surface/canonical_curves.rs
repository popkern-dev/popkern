//! Exact edge curves for sweep boundaries: a straight NURBS becomes a
//! `Line`, a full circle an `Arc` anchored at the curve's start point,
//! everything else stays NURBS (`analytic_or_nurbs_curve`); the
//! revolution builders keep circles as the surface's own iso-NURBS so
//! edge and surface stay co-parameterised (`line_or_nurbs_curve`).

use crate::surface::canonical::recognize_circle;

/// Recognize a straight NURBS curve: every control point collinear
/// with the endpoint chord (relative 1e-9). Returns the origin, unit
/// direction and length.
pub(crate) fn recognize_line(
    curve: &popkern_core::geometry::NurbsCurve,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Vec3,
    f64,
)> {
    let first = *curve.control_points.first()?;
    let last = *curve.control_points.last()?;
    let chord = last - first;
    let len = chord.length();
    if len < 1e-12 {
        return None;
    }
    let dir = chord * (1.0 / len);
    for p in &curve.control_points {
        let d = *p - first;
        let off = d - dir * d.dot(&dir);
        if off.length() > 1e-9 * len.max(1.0) {
            return None;
        }
    }
    Some((first, dir, len))
}

/// A line when the curve is straight (its length parameter is a linear
/// reparameterisation of the NURBS parameter), else the NURBS itself:
/// the edge stays co-parameterised with the surface iso-curve it came
/// from, which a circle's angle parameter would not be (PK-408: the
/// mesher maps edge-parameter fractions onto pcurves linearly).
pub(crate) fn line_or_nurbs_curve(
    curve: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::CurveOnEdge {
    use popkern_core::geometry::curve::CurveEval;
    if let Some((origin, direction, length)) = recognize_line(curve) {
        return popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Line { origin, direction },
            t_start: 0.0,
            t_end: length,
        };
    }
    let (t_start, t_end) = curve.parameter_range();
    popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    }
}

pub(crate) fn analytic_or_nurbs_curve(
    curve: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::CurveOnEdge {
    use popkern_core::geometry::curve::CurveEval;
    if let Some((origin, direction, length)) = recognize_line(curve) {
        return popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Line { origin, direction },
            t_start: 0.0,
            t_end: length,
        };
    }
    if let Some((center, normal, radius)) = recognize_circle(curve) {
        let (t0, _) = curve.parameter_range();
        let start = curve.point_at(t0);
        // Angles MUST be measured in the Arc evaluator's own frame
        // (circle_basis) — orthonormal_basis is a different
        // convention and a mismatched frame rotates the arc's
        // parameterization off its seam vertex (measured: closed cap
        // rims whose trim start landed 90° from the vertex).
        let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(normal);
        let d = start - center;
        let theta0 = d.dot(&e2).atan2(d.dot(&e1));
        return popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle: theta0,
                end_angle: theta0 + std::f64::consts::TAU,
            },
            t_start: theta0,
            t_end: theta0 + std::f64::consts::TAU,
        };
    }
    let (t_start, t_end) = curve.parameter_range();
    popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    }
}

#[cfg(test)]
mod recognition_phase_tests {
    use popkern_core::geometry::curve::CurveEval;

    /// PK-479 pin: circle recognition preserves the wire's start
    /// phase for BOTH windings — the recognized curve evaluated at
    /// its own t_start must land on the wire's start point (the
    /// winding-derived normal changes the evaluator frame, and a
    /// frame-inconsistent theta rotates every trim off its seam
    /// vertex).
    #[test]
    fn circle_recognition_preserves_start_phase_both_windings() {
        let circle = popkern_core::geometry::curve::CurveGeometry::Circle {
            center: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let fwd = circle.to_nurbs();
        let rev = {
            let mut c = fwd.clone();
            let (t0, t1) = c.parameter_range();
            let mut knots: Vec<f64> = c.knots.iter().rev().map(|&k| t0 + t1 - k).collect();
            knots.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
            c.knots = knots;
            c.control_points.reverse();
            c.weights.reverse();
            c
        };
        for (name, wire) in [("fwd", fwd), ("rev", rev)] {
            let (t0, _) = wire.parameter_range();
            let s = wire.point_at(t0);
            let rec = super::analytic_or_nurbs_curve(&wire);
            let rs = rec.geometry.point_at(rec.t_start);
            assert!(
                rs.distance_to(&s) < 1e-9,
                "{name}: recognized start {rs:?} != wire start {s:?}"
            );
        }
    }
}
