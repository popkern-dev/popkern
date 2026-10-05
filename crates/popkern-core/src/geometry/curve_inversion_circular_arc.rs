use super::{CircularArcParams, circle_basis_for};
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

// ---------------------------------------------------------------------------
// CurveGeometry::as_circular_arc — sample-based detector
// ---------------------------------------------------------------------------

impl CurveGeometry {
    /// Detect whether this curve is geometrically a circular arc.
    ///
    /// Returns `Some(params)` if 5 evenly-spaced sample points all lie on
    /// the circle fitted through 3 of them (within `tolerance.linear()`);
    /// otherwise `None`.
    ///
    /// Works for `Arc` / `Circle` (returns the analytic representation
    /// directly) and for `Nurbs` curves that geometrically trace a
    /// circular arc (returns the fitted representation).  Returns `None`
    /// for `Line` / `Ellipse` / `Parabola` / `Hyperbola`.
    pub fn as_circular_arc(&self, tolerance: &Tolerance) -> Option<CircularArcParams> {
        match self {
            CurveGeometry::Line { .. } => None,
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => {
                let n = normal.normalized()?;
                let (u, _v) = circle_basis_for(n);
                let start = *center + u * *radius;
                Some(CircularArcParams {
                    center: *center,
                    radius: *radius,
                    normal: n,
                    start_angle: 0.0,
                    sweep_angle: std::f64::consts::TAU,
                    start,
                    end: start,
                })
            }
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => {
                let n = normal.normalized()?;
                let (u, v) = circle_basis_for(n);
                let start_dir = u * start_angle.pcos() + v * start_angle.psin();
                let end_dir = u * end_angle.pcos() + v * end_angle.psin();
                Some(CircularArcParams {
                    center: *center,
                    radius: *radius,
                    normal: n,
                    start_angle: *start_angle,
                    sweep_angle: end_angle - start_angle,
                    start: *center + start_dir * *radius,
                    end: *center + end_dir * *radius,
                })
            }
            CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. }
            | CurveGeometry::Hyperbola { .. } => None,
            CurveGeometry::Nurbs(_) => detect_circular_arc_from_samples(self, tolerance),
        }
    }
}

/// Sample-based circular-arc detector.  See `as_circular_arc`.
fn detect_circular_arc_from_samples(
    curve: &CurveGeometry,
    tolerance: &Tolerance,
) -> Option<CircularArcParams> {
    let (t_lo, t_hi) = curve.parameter_range();
    if !t_lo.is_finite() || !t_hi.is_finite() || t_hi <= t_lo {
        return None;
    }
    let frac = |x: f64| t_lo + (t_hi - t_lo) * x;
    let samples = [
        curve.point_at(frac(0.0)),
        curve.point_at(frac(0.25)),
        curve.point_at(frac(0.5)),
        curve.point_at(frac(0.75)),
        curve.point_at(frac(1.0)),
    ];

    let (center, radius, normal) = circle_through_three_points(samples[0], samples[2], samples[4])?;

    // Verify all 5 samples lie on the fitted circle within tolerance.
    let lin_tol = tolerance.linear();
    for (i, p) in samples.iter().enumerate() {
        let radial = (*p - center).length() - radius;
        let along_normal = (*p - center).dot(&normal).abs();
        let total = (radial * radial + along_normal * along_normal).sqrt();
        if total > lin_tol {
            // Sample i is off the candidate circle. Not an arc within tolerance.
            let _ = i;
            return None;
        }
    }

    // Compute start/end angles in the (u, v) basis.
    let (u, v) = circle_basis_for(normal);
    let start_vec = samples[0] - center;
    let end_vec = samples[4] - center;
    let start_angle = start_vec.dot(&v).patan2(start_vec.dot(&u));
    let end_angle = end_vec.dot(&v).patan2(end_vec.dot(&u));

    // Determine sweep direction by checking which sense puts the midpoint
    // of the parameter range close to the midpoint of the arc.
    let mid_vec = samples[2] - center;
    let mid_angle = mid_vec.dot(&v).patan2(mid_vec.dot(&u));

    let sweep_ccw = signed_arc_sweep(start_angle, end_angle, mid_angle);

    Some(CircularArcParams {
        center,
        radius,
        normal,
        start_angle,
        sweep_angle: sweep_ccw,
        start: samples[0],
        end: samples[4],
    })
}

/// Compute the signed sweep from `start` to `end` (around the chosen
/// `normal`, CCW positive) such that `mid_angle` lies between them.
///
/// Returns the sweep in radians; negate the caller's normal to flip
/// direction if needed.
fn signed_arc_sweep(start: f64, end: f64, mid: f64) -> f64 {
    use std::f64::consts::TAU;
    // Shift everything so start is 0; CCW sweep is in [0, 2π).
    let end_shifted = (end - start).rem_euclid(TAU);
    let mid_shifted = (mid - start).rem_euclid(TAU);
    if mid_shifted <= end_shifted {
        // Mid lies in CCW arc from start to end; sweep is +end_shifted.
        end_shifted
    } else {
        // Mid lies in the other arc; sweep is the CW direction.
        end_shifted - TAU
    }
}

/// Fit a circle through three points in 3-D.
///
/// Returns `(center, radius, unit_normal)` where `unit_normal` follows
/// the right-hand rule for traversal `p0 → p1 → p2`.  Returns `None`
/// if the points are collinear.
pub(super) fn circle_through_three_points(
    p0: Point3,
    p1: Point3,
    p2: Point3,
) -> Option<(Point3, f64, Vec3)> {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let n = v1.cross(&v2);
    let n_sq = n.dot(&n);
    if n_sq < 1e-30 {
        return None; // collinear or coincident
    }
    let v1_sq = v1.dot(&v1);
    let v2_sq = v2.dot(&v2);
    // Circumcenter formula in 3-D:
    //   centre = p0 + ( (|v2|² · (n × v1) + |v1|² · (v2 × n)) / (2 · |n|²) )
    let to_center = (n.cross(&v1) * v2_sq + v2.cross(&n) * v1_sq) * (1.0 / (2.0 * n_sq));
    let center = p0 + to_center;
    let radius = (p0 - center).length();
    let normal = n.normalized()?;
    Some((center, radius, normal))
}
