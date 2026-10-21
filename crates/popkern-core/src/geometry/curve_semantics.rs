//! Semantic curve classification helpers.
//!
//! These helpers answer geometric questions without requiring callers to
//! care about the storage variant.  For example, topology imported from
//! construction stages may store a straight segment as a degree-1 NURBS;
//! apex construction should be able to treat that as a line without
//! rewriting the BRep's curve representation.

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};

/// A bounded straight-line curve span.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineSegmentParams {
    /// Start point of the detected span.
    pub start: Point3,
    /// End point of the detected span.
    pub end: Point3,
    /// Unit direction from `start` to `end`.
    pub direction: Vec3,
    /// Span length.
    pub length: f64,
    /// Parameter at `start`.
    pub t_start: f64,
    /// Parameter at `end`.
    pub t_end: f64,
}

impl CurveGeometry {
    /// The declared angular span `(start_angle, end_angle)` of an Arc
    /// or Ellipse; `None` for curves without one.
    pub fn angular_span(&self) -> Option<(f64, f64)> {
        match self {
            CurveGeometry::Arc {
                start_angle,
                end_angle,
                ..
            }
            | CurveGeometry::Ellipse {
                start_angle,
                end_angle,
                ..
            } => Some((*start_angle, *end_angle)),
            _ => None,
        }
    }

    /// The same curve declaring `window` as its angular span. An Arc
    /// or Ellipse edge moved onto the complement (or another turn)
    /// still evaluates there, but projection and parameter-range
    /// queries answer for the declared span (PK-906). Curves without a
    /// declared span are returned unchanged.
    pub fn spanning(&self, window: (f64, f64)) -> CurveGeometry {
        let (lo, hi) = (window.0.min(window.1), window.0.max(window.1));
        let mut out = self.clone();
        if let CurveGeometry::Arc {
            start_angle,
            end_angle,
            ..
        }
        | CurveGeometry::Ellipse {
            start_angle,
            end_angle,
            ..
        } = &mut out
        {
            (*start_angle, *end_angle) = (lo, hi);
        }
        out
    }

    /// Detect whether this curve is geometrically a straight segment over its
    /// native parameter range.
    ///
    /// This recognizes analytic [`CurveGeometry::Line`] curves and NURBS curves
    /// whose control polygon is collinear over the full curve range.  It does
    /// not rewrite the curve; callers can use the returned parameters for
    /// predicates, diagnostics, or local canonicalization.
    pub fn as_line_segment(&self, tolerance: f64) -> Option<LineSegmentParams> {
        self.as_line_segment_in_range(self.parameter_range(), tolerance)
    }

    /// Detect whether this curve is geometrically a straight segment over
    /// `t_range`.
    pub fn as_line_segment_in_range(
        &self,
        t_range: (f64, f64),
        tolerance: f64,
    ) -> Option<LineSegmentParams> {
        match self {
            CurveGeometry::Line { origin, direction } => {
                let start = *origin + *direction * t_range.0;
                let end = *origin + *direction * t_range.1;
                line_segment_from_endpoints(start, end, t_range, tolerance)
            }
            CurveGeometry::Nurbs(nurbs) => {
                if nurbs.control_points.len() < 2
                    || nurbs.weights.len() != nurbs.control_points.len()
                {
                    return None;
                }
                if nurbs.weights.iter().any(|w| !w.is_finite() || *w <= 0.0) {
                    return None;
                }
                let params = line_segment_from_endpoints(
                    nurbs.point_at(t_range.0),
                    nurbs.point_at(t_range.1),
                    t_range,
                    tolerance,
                )?;
                for point in &nurbs.control_points {
                    if point_to_segment_line_distance(*point, params.start, params.direction)
                        > tolerance
                    {
                        return None;
                    }
                }
                Some(params)
            }
            CurveGeometry::Circle { .. }
            | CurveGeometry::Arc { .. }
            | CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. }
            | CurveGeometry::Hyperbola { .. } => None,
        }
    }

    /// Return a canonical analytic line for a semantically straight span.
    pub fn canonical_line_in_range(
        &self,
        t_range: (f64, f64),
        tolerance: f64,
    ) -> Option<CurveGeometry> {
        let line = self.as_line_segment_in_range(t_range, tolerance)?;
        let direction = (line.end - line.start) * (1.0 / (line.t_end - line.t_start));
        Some(CurveGeometry::Line {
            origin: line.start - direction * line.t_start,
            direction,
        })
    }
}

fn line_segment_from_endpoints(
    start: Point3,
    end: Point3,
    t_range: (f64, f64),
    tolerance: f64,
) -> Option<LineSegmentParams> {
    if !t_range.0.is_finite() || !t_range.1.is_finite() || (t_range.1 - t_range.0).abs() < 1e-30 {
        return None;
    }
    let chord = end - start;
    let length = chord.length();
    if length <= tolerance.max(1e-12) {
        return None;
    }
    Some(LineSegmentParams {
        start,
        end,
        direction: chord * (1.0 / length),
        length,
        t_start: t_range.0,
        t_end: t_range.1,
    })
}

fn point_to_segment_line_distance(point: Point3, origin: Point3, direction: Vec3) -> f64 {
    let offset = point - origin;
    let projected = origin + direction * offset.dot(&direction);
    point.distance_to(&projected)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::NurbsCurve;

    #[test]
    fn analytic_line_reports_bounded_segment() {
        let curve = CurveGeometry::Line {
            origin: Point3::new(1.0, 2.0, 3.0),
            direction: Vec3::new(2.0, 0.0, 0.0),
        };
        let line = curve
            .as_line_segment_in_range((0.25, 1.25), 1e-9)
            .expect("line segment");
        assert_eq!(line.start, Point3::new(1.5, 2.0, 3.0));
        assert_eq!(line.end, Point3::new(3.5, 2.0, 3.0));
        assert!((line.length - 2.0).abs() < 1e-12);
    }

    #[test]
    fn degree_one_nurbs_line_reports_segment() {
        let curve = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let line = curve.as_line_segment(1e-9).expect("line segment");
        assert_eq!(line.start, Point3::new(0.0, 0.0, 0.0));
        assert_eq!(line.end, Point3::new(2.0, 0.0, 0.0));
        assert!((line.direction.x - 1.0).abs() < 1e-12);
    }

    #[test]
    fn elevated_nurbs_line_reports_segment() {
        let curve = CurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        let line = curve.as_line_segment(1e-9).expect("line segment");
        assert_eq!(line.start, Point3::new(0.0, 1.0, 0.0));
        assert_eq!(line.end, Point3::new(2.0, 1.0, 0.0));
    }

    #[test]
    fn canonical_line_in_nonzero_range_preserves_endpoint_parameters() {
        let curve = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![2.0, 2.0, 6.0, 6.0],
            control_points: vec![Point3::new(10.0, 0.0, 0.0), Point3::new(18.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let line = curve
            .canonical_line_in_range((3.0, 5.0), 1e-9)
            .expect("canonical line");
        assert_eq!(line.point_at(3.0), curve.point_at(3.0));
        assert_eq!(line.point_at(5.0), curve.point_at(5.0));
    }

    #[test]
    fn bent_nurbs_polyline_is_not_line() {
        let curve = CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.2, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        assert!(curve.as_line_segment(1e-3).is_none());
    }

    #[test]
    fn circular_arc_is_not_line() {
        let curve = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: std::f64::consts::FRAC_PI_2,
        };
        assert!(curve.as_line_segment(1e-3).is_none());
    }
}
