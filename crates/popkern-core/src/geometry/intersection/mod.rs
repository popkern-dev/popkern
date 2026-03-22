//! Curve-curve intersection algorithms.

mod bezier_clipping;
mod circle_circle;
mod line_circle;
mod line_line;
mod newton;
mod nurbs_nurbs;
mod overlap;

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveGeometry;
use crate::tolerance::Tolerance;
use serde::{Deserialize, Serialize};

/// A single intersection result between two curves.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CurveIntersection {
    /// Transversal or tangential crossing at a point.
    Point {
        /// The intersection point in 3D space.
        point: Point3,
        /// Parameter value on the first curve.
        t1: f64,
        /// Parameter value on the second curve.
        t2: f64,
    },
    /// Curves coincide over a parameter interval.
    Overlap {
        /// Parameter range [start, end] on the first curve.
        t1_range: (f64, f64),
        /// Parameter range [start, end] on the second curve.
        t2_range: (f64, f64),
    },
}

/// Result of a curve-curve intersection operation.
#[derive(Debug, Clone)]
pub struct CurveIntersectionResult {
    /// All intersections found, sorted by t1.
    pub intersections: Vec<CurveIntersection>,
}

/// Intersect two curves within the given tolerance.
pub fn intersect_curves(
    curve1: &CurveGeometry,
    t1_range: (f64, f64),
    curve2: &CurveGeometry,
    t2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<CurveIntersectionResult, CoreError> {
    if t1_range.0 >= t1_range.1 {
        return Err(CoreError::PreconditionFailed(format!(
            "t1_range invalid: {} >= {}",
            t1_range.0, t1_range.1
        )));
    }
    if t2_range.0 >= t2_range.1 {
        return Err(CoreError::PreconditionFailed(format!(
            "t2_range invalid: {} >= {}",
            t2_range.0, t2_range.1
        )));
    }
    // Suppress unused parameter warnings until dispatch is implemented in Task 10.
    let _ = (curve1, curve2, tolerance);
    // Dispatch stub — will be filled in Task 10
    Ok(CurveIntersectionResult {
        intersections: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_t1_range_fails() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        };
        let result = intersect_curves(&line, (1.0, 0.0), &line, (0.0, 1.0), &Tolerance::default());
        assert!(result.is_err());
    }

    #[test]
    fn valid_range_returns_empty_stub() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        };
        let result =
            intersect_curves(&line, (0.0, 1.0), &line, (0.0, 1.0), &Tolerance::default()).unwrap();
        assert!(result.intersections.is_empty()); // stub returns empty
    }
}
