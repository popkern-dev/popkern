//! Analytic plane-plane intersection solver.
//!
//! When both surfaces are planes the intersection is either empty (parallel,
//! separated), a coincident region (parallel, within tolerance), or a line
//! (non-parallel).  Computing this analytically is far cheaper than the
//! general NURBS marching pipeline.

use crate::error::CoreError;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect two planes.
///
/// Returns `Some(Vec<SurfaceSurfaceIntersection>)` when both surfaces are
/// planes (the fast path), or `None` when the caller should fall back to the
/// general NURBS marching pipeline.
///
/// The parameter ranges `(u1_range, v1_range)` and `(u2_range, v2_range)`
/// are used to clip the result line and to compute UV coordinates.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_plane(
    origin1: Point3,
    normal1: Vec3,
    u_axis1: Vec3,
    v_axis1: Vec3,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    origin2: Point3,
    normal2: Vec3,
    u_axis2: Vec3,
    v_axis2: Vec3,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<Vec<SurfaceSurfaceIntersection>, CoreError> {
    // d = n1 × n2
    let d = normal1.cross(&normal2);
    let d_len_sq = d.length_squared();

    // Check for parallel planes: |d|² < sin²(angular_tol) ≈ angular_tol² for small angles.
    let sin_tol = tolerance.angular.sin();
    if d_len_sq < sin_tol * sin_tol {
        // Planes are parallel (or anti-parallel). Measure the signed distance
        // from origin2 to plane 1: dist = (origin2 - origin1) · n1.
        let diff = origin2 - origin1;
        let dist = diff.dot(&normal1).abs();
        if dist < tolerance.linear {
            // Coincident planes.
            return Ok(vec![SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range,
                v1_range,
                u2_range,
                v2_range,
            }]);
        }
        // Parallel but separated — no intersection.
        return Ok(vec![]);
    }

    // Non-parallel planes: intersection is a line.
    let direction = d
        .normalized()
        .ok_or_else(|| CoreError::PreconditionFailed("degenerate intersection direction".into()))?;

    // Find a point on the intersection line by solving the 2×2 system:
    //   n1 · p = d1   where d1 = n1 · origin1
    //   n2 · p = d2   where d2 = n2 · origin2
    // We fix the coordinate with the largest |direction| component to 0 and
    // solve for the other two.
    let d1 = normal1.dot(&Vec3::new(origin1.x, origin1.y, origin1.z));
    let d2 = normal2.dot(&Vec3::new(origin2.x, origin2.y, origin2.z));

    let point_on_line = solve_point_on_line(normal1, d1, normal2, d2, direction)?;

    // Clip the line to the parameter ranges of both planes.
    // We express the line parameter t in terms of the 3D point:
    //   P(t) = point_on_line + t * direction
    // Then restrict t so that P(t) lies within the UV boxes of both planes.
    let (t_min, t_max) = clip_line_to_uv_boxes(
        point_on_line,
        direction,
        origin1,
        u_axis1,
        v_axis1,
        u1_range,
        v1_range,
        origin2,
        u_axis2,
        v_axis2,
        u2_range,
        v2_range,
    );

    if t_max - t_min < tolerance.linear {
        // Intersection line does not pass through the parameter regions.
        return Ok(vec![]);
    }

    // Build the 3D NURBS line segment.
    let line_3d = CurveGeometry::Line {
        origin: point_on_line,
        direction,
    };
    let curve_3d = line_3d.to_nurbs_in_range((t_min, t_max));

    // Build the UV curve on surface 1.
    // A 3D point P lies on surface 1 at:
    //   u1 = (P - origin1) · u_axis1,   v1 = (P - origin1) · v_axis1
    // Substituting P = point_on_line + t * direction:
    //   u1(t) = u1_0 + t * (direction · u_axis1)
    //   v1(t) = v1_0 + t * (direction · v_axis1)
    let offset1 = point_on_line - origin1;
    let u1_0 = offset1.dot(&u_axis1);
    let v1_0 = offset1.dot(&v_axis1);
    let du1 = direction.dot(&u_axis1);
    let dv1 = direction.dot(&v_axis1);
    let curve_uv1 = linear_nurbs_uv(u1_0, v1_0, du1, dv1, t_min, t_max);

    // Build the UV curve on surface 2.
    let offset2 = point_on_line - origin2;
    let u2_0 = offset2.dot(&u_axis2);
    let v2_0 = offset2.dot(&v_axis2);
    let du2 = direction.dot(&u_axis2);
    let dv2 = direction.dot(&v_axis2);
    let curve_uv2 = linear_nurbs_uv(u2_0, v2_0, du2, dv2, t_min, t_max);

    let cos = CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: false,
    };

    Ok(vec![SurfaceSurfaceIntersection::Curve(cos)])
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Solve for a point on the intersection line given:
///   n1 · p = d1
///   n2 · p = d2
/// by zeroing the coordinate corresponding to the largest component of
/// `direction` and solving the resulting 2×2 linear system.
fn solve_point_on_line(
    n1: Vec3,
    d1: f64,
    n2: Vec3,
    d2: f64,
    direction: Vec3,
) -> Result<Point3, CoreError> {
    // Identify the axis with the largest absolute direction component.
    let ax = direction.x.abs();
    let ay = direction.y.abs();
    let az = direction.z.abs();

    let (x, y, z) = if ax >= ay && ax >= az {
        // Fix x = 0; solve for y and z.
        // n1.y * y + n1.z * z = d1
        // n2.y * y + n2.z * z = d2
        let (y, z) = solve_2x2(n1.y, n1.z, d1, n2.y, n2.z, d2)?;
        (0.0, y, z)
    } else if ay >= ax && ay >= az {
        // Fix y = 0; solve for x and z.
        let (x, z) = solve_2x2(n1.x, n1.z, d1, n2.x, n2.z, d2)?;
        (x, 0.0, z)
    } else {
        // Fix z = 0; solve for x and y.
        let (x, y) = solve_2x2(n1.x, n1.y, d1, n2.x, n2.y, d2)?;
        (x, y, 0.0)
    };

    Ok(Point3::new(x, y, z))
}

/// Solve the 2×2 system:
///   a * x + b * y = c1
///   d * x + e * y = c2
pub(crate) fn solve_2x2(
    a: f64,
    b: f64,
    c1: f64,
    d: f64,
    e: f64,
    c2: f64,
) -> Result<(f64, f64), CoreError> {
    let det = a * e - b * d;
    if det.abs() < f64::EPSILON * 1e6 {
        return Err(CoreError::PreconditionFailed(
            "degenerate 2x2 system in plane-plane intersection".into(),
        ));
    }
    let x = (c1 * e - b * c2) / det;
    let y = (a * c2 - c1 * d) / det;
    Ok((x, y))
}

/// Build a degree-1 NURBS curve in UV space for the line:
///   u(t) = u0 + t * du
///   v(t) = v0 + t * dv
/// encoded as `Point3(u, v, 0.0)` control points on `[t_min, t_max]`.
fn linear_nurbs_uv(u0: f64, v0: f64, du: f64, dv: f64, t_min: f64, t_max: f64) -> NurbsCurve {
    let p0 = Point3::new(u0 + du * t_min, v0 + dv * t_min, 0.0);
    let p1 = Point3::new(u0 + du * t_max, v0 + dv * t_max, 0.0);
    NurbsCurve {
        degree: 1,
        knots: vec![t_min, t_min, t_max, t_max],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Clip the parametric line `P(t) = point_on_line + t * direction` to the
/// intersection of the UV boxes of both planes.
///
/// Returns `(t_min, t_max)`.  If the clipped range is empty, `t_max < t_min`.
#[allow(clippy::too_many_arguments)]
fn clip_line_to_uv_boxes(
    point: Point3,
    direction: Vec3,
    origin1: Point3,
    u_axis1: Vec3,
    v_axis1: Vec3,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    origin2: Point3,
    u_axis2: Vec3,
    v_axis2: Vec3,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> (f64, f64) {
    let mut t_lo = f64::NEG_INFINITY;
    let mut t_hi = f64::INFINITY;

    // Clip against surface 1 UV box.
    let offset1 = point - origin1;
    clip_axis(
        offset1.dot(&u_axis1),
        direction.dot(&u_axis1),
        u1_range,
        &mut t_lo,
        &mut t_hi,
    );
    clip_axis(
        offset1.dot(&v_axis1),
        direction.dot(&v_axis1),
        v1_range,
        &mut t_lo,
        &mut t_hi,
    );

    // Clip against surface 2 UV box.
    let offset2 = point - origin2;
    clip_axis(
        offset2.dot(&u_axis2),
        direction.dot(&u_axis2),
        u2_range,
        &mut t_lo,
        &mut t_hi,
    );
    clip_axis(
        offset2.dot(&v_axis2),
        direction.dot(&v_axis2),
        v2_range,
        &mut t_lo,
        &mut t_hi,
    );

    (t_lo, t_hi)
}

/// Clip a 1-D interval `[t_lo, t_hi]` using the constraint:
///   range.0 <= param0 + t * dparam <= range.1
fn clip_axis(param0: f64, dparam: f64, range: (f64, f64), t_lo: &mut f64, t_hi: &mut f64) {
    if dparam.abs() < f64::EPSILON * 1e6 {
        // Direction is (nearly) perpendicular to this axis — check the fixed value.
        if param0 < range.0 - f64::EPSILON || param0 > range.1 + f64::EPSILON {
            // The line is outside this slab entirely.
            *t_hi = f64::NEG_INFINITY; // force empty interval
        }
        return;
    }
    let t_enter = (range.0 - param0) / dparam;
    let t_exit = (range.1 - param0) / dparam;
    let (t_a, t_b) = if t_enter < t_exit {
        (t_enter, t_exit)
    } else {
        (t_exit, t_enter)
    };
    if t_a > *t_lo {
        *t_lo = t_a;
    }
    if t_b < *t_hi {
        *t_hi = t_b;
    }
}

#[cfg(test)]
#[path = "plane_plane_tests.rs"]
mod tests;
