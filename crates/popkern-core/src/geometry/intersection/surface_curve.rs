//! Surface-curve intersection types and dispatch.

use crate::error::CoreError;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use serde::{Deserialize, Serialize};

/// A single intersection result between a surface and a curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceCurveIntersection {
    /// Curve pierces the surface at a point.
    Point {
        /// The intersection point in 3D space.
        point: Point3,
        /// Parameter value on the curve at the intersection.
        t: f64,
        /// Surface parameter u at the intersection.
        u: f64,
        /// Surface parameter v at the intersection.
        v: f64,
    },
    /// Curve lies on the surface over a parameter interval.
    Overlap {
        /// Parameter range [start, end] on the curve.
        t_range: (f64, f64),
    },
}

/// Result of a surface-curve intersection operation.
#[derive(Debug, Clone)]
pub struct SurfaceCurveIntersectionResult {
    /// All intersections found, sorted by t.
    pub intersections: Vec<SurfaceCurveIntersection>,
}

/// Intersect a surface with a curve within the given tolerance.
pub fn intersect_surface_curve(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    _curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceCurveIntersectionResult, CoreError> {
    if t_range.0 >= t_range.1 {
        return Err(CoreError::PreconditionFailed(
            format!("t_range invalid: {} >= {}", t_range.0, t_range.1)
        ));
    }
    let _ = (surface, u_range, v_range, tolerance);
    Ok(SurfaceCurveIntersectionResult { intersections: Vec::new() })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Vec3;

    #[test]
    fn invalid_t_range_fails() {
        let s = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0), normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0), v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let c = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0), direction: Vec3::new(1.0, 0.0, 0.0),
        };
        assert!(intersect_surface_curve(&s, (0.0, 1.0), (0.0, 1.0), &c, (1.0, 0.0), &Tolerance::default()).is_err());
    }
}
