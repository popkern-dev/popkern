//! Curve-curve intersection algorithms.

pub(crate) mod adaptive_step;
mod bezier_clipping;
mod circle_circle;
pub mod coincident;
pub mod coincident_analytic;
pub mod degeneracy;
mod line_circle;
mod line_cone;
mod line_cylinder;
mod line_line;
mod line_plane;
mod line_sphere;
mod line_torus;
pub mod marching;
mod newton;
mod nurbs_nurbs;
mod overlap;
pub(crate) mod plane_cone;
pub(crate) mod plane_cylinder;
pub(crate) mod plane_plane;
pub(crate) mod plane_sphere;
pub mod starting_points;
pub mod surface_curve;
mod surface_curve_general;
mod surface_curve_newton;
mod surface_curve_overlap;
pub mod surface_subdivision;
pub mod surface_surface;

pub use surface_curve::{
    SurfaceCurveIntersection, SurfaceCurveIntersectionResult, intersect_surface_curve,
};
pub use surface_surface::{
    SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult, intersect_surfaces,
    intersect_surfaces_with_handler,
};

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

/// Swap t1/t2 (and t1_range/t2_range) in a list of intersection results.
/// Used when the caller reverses the curve order for an asymmetric solver.
fn swap_parameters(results: Vec<CurveIntersection>) -> Vec<CurveIntersection> {
    results
        .into_iter()
        .map(|r| match r {
            CurveIntersection::Point { point, t1, t2 } => CurveIntersection::Point {
                point,
                t1: t2,
                t2: t1,
            },
            CurveIntersection::Overlap { t1_range, t2_range } => CurveIntersection::Overlap {
                t1_range: t2_range,
                t2_range: t1_range,
            },
        })
        .collect()
}

/// Returns true if `angle` lies within the arc range `[start, end]`.
/// Handles angle wrapping: if start > end, accepts [start, 2pi) union [0, end].
fn angle_in_arc_range(angle: f64, start: f64, end: f64) -> bool {
    let tau = std::f64::consts::TAU;
    let a = ((angle % tau) + tau) % tau;
    let s = ((start % tau) + tau) % tau;
    let e = ((end % tau) + tau) % tau;

    if s <= e {
        a >= s && a <= e
    } else {
        a >= s || a <= e // wraps around 0/2pi
    }
}

/// Filter intersection results to only include those where a curve's angle
/// parameter falls within the given arc range.
///
/// When `filter_curve1` is true, the arc is curve 1 (filter on t1);
/// when false, the arc is curve 2 (filter on t2).
fn filter_arc_range(
    results: Vec<CurveIntersection>,
    filter_curve1: bool,
    start_angle: f64,
    end_angle: f64,
) -> Vec<CurveIntersection> {
    results
        .into_iter()
        .filter(|r| match r {
            CurveIntersection::Point { t1, t2, .. } => {
                let angle = if filter_curve1 { *t1 } else { *t2 };
                angle_in_arc_range(angle, start_angle, end_angle)
            }
            CurveIntersection::Overlap { .. } => {
                // Keep overlaps — they are already constrained by the solver
                true
            }
        })
        .collect()
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

    let intersections = match (curve1, curve2) {
        // Line-Line
        (
            CurveGeometry::Line {
                origin: o1,
                direction: d1,
            },
            CurveGeometry::Line {
                origin: o2,
                direction: d2,
            },
        ) => line_line::intersect_line_line(o1, d1, t1_range, o2, d2, t2_range, tolerance),

        // Line-Circle
        (
            CurveGeometry::Line { origin, direction },
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            },
        ) => line_circle::intersect_line_circle(
            origin, direction, t1_range, center, normal, *radius, t2_range, tolerance,
        ),

        // Circle-Line (swap)
        (
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            },
            CurveGeometry::Line { origin, direction },
        ) => {
            let results = line_circle::intersect_line_circle(
                origin, direction, t2_range, center, normal, *radius, t1_range, tolerance,
            );
            swap_parameters(results)
        }

        // Line-Arc: treat arc as circle, filter to arc range
        (
            CurveGeometry::Line { origin, direction },
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            },
        ) => {
            let results = line_circle::intersect_line_circle(
                origin,
                direction,
                t1_range,
                center,
                normal,
                *radius,
                (*start_angle, *end_angle),
                tolerance,
            );
            filter_arc_range(results, false, *start_angle, *end_angle)
        }

        // Arc-Line (swap)
        (CurveGeometry::Arc { .. }, CurveGeometry::Line { .. }) => {
            let swapped = intersect_curves(curve2, t2_range, curve1, t1_range, tolerance)?;
            swap_parameters(swapped.intersections)
        }

        // Circle-Circle
        (
            CurveGeometry::Circle {
                center: c1,
                normal: n1,
                radius: r1,
            },
            CurveGeometry::Circle {
                center: c2,
                normal: n2,
                radius: r2,
            },
        ) => circle_circle::intersect_circle_circle(
            c1, n1, *r1, t1_range, c2, n2, *r2, t2_range, tolerance,
        ),

        // Circle-Arc: run as circle-circle, filter curve 2 arc range
        (
            CurveGeometry::Circle {
                center: c1,
                normal: n1,
                radius: r1,
            },
            CurveGeometry::Arc {
                center: c2,
                normal: n2,
                radius: r2,
                start_angle,
                end_angle,
            },
        ) => {
            let results = circle_circle::intersect_circle_circle(
                c1,
                n1,
                *r1,
                t1_range,
                c2,
                n2,
                *r2,
                (*start_angle, *end_angle),
                tolerance,
            );
            filter_arc_range(results, false, *start_angle, *end_angle)
        }

        // Arc-Circle (swap)
        (CurveGeometry::Arc { .. }, CurveGeometry::Circle { .. }) => {
            let swapped = intersect_curves(curve2, t2_range, curve1, t1_range, tolerance)?;
            swap_parameters(swapped.intersections)
        }

        // Arc-Arc: run as circle-circle, filter both arc ranges
        (
            CurveGeometry::Arc {
                center: c1,
                normal: n1,
                radius: r1,
                start_angle: s1,
                end_angle: e1,
            },
            CurveGeometry::Arc {
                center: c2,
                normal: n2,
                radius: r2,
                start_angle: s2,
                end_angle: e2,
            },
        ) => {
            let results = circle_circle::intersect_circle_circle(
                c1,
                n1,
                *r1,
                (*s1, *e1),
                c2,
                n2,
                *r2,
                (*s2, *e2),
                tolerance,
            );
            let results = filter_arc_range(results, true, *s1, *e1);
            filter_arc_range(results, false, *s2, *e2)
        }

        // Anything involving NURBS, or fallback: promote to NURBS
        _ => {
            let n1 = curve1.to_nurbs_in_range(t1_range);
            let n2 = curve2.to_nurbs_in_range(t2_range);
            nurbs_nurbs::intersect_nurbs_nurbs(&n1, t1_range, &n2, t2_range, tolerance)
        }
    };

    Ok(CurveIntersectionResult { intersections })
}

#[cfg(test)]
mod tests;
