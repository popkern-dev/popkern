//! Revolution geometry helpers — point rotation, revolution curves, revolution surfaces.

use crate::error::CoreError;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{FRAC_PI_2, TAU};

/// Rotate a point around an axis using Rodrigues' rotation formula.
pub fn rotate_point(point: &Point3, axis_origin: &Point3, axis_dir: &Vec3, angle: f64) -> Point3 {
    let v = *point - *axis_origin;
    let k = axis_dir.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let rotated = v * angle.cos() + k.cross(&v) * angle.sin() + k * k.dot(&v) * (1.0 - angle.cos());
    *axis_origin + rotated
}

/// Create a revolution curve for a point revolving around an axis.
/// Returns Circle for full 360°, Arc for partial.
pub fn revolution_curve(
    point: &Point3,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    angle: f64,
) -> CurveGeometry {
    let axis = axis_dir.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let v = *point - *axis_origin;
    let proj = axis * v.dot(&axis);
    let center = *axis_origin + proj;
    let radial = v - proj;
    let radius = radial.length();

    if (angle - TAU).abs() < 1e-10 {
        CurveGeometry::Circle { center, normal: axis, radius }
    } else {
        CurveGeometry::Arc {
            center,
            normal: axis,
            radius,
            start_angle: 0.0,
            end_angle: angle,
        }
    }
}

/// Build the arc control points and weights for revolving a single point around
/// an axis by the given angle. The arc is degree-2 rational.
///
/// Returns `(arc_points, arc_weights, knots_v)`.
fn arc_control_points_for_revolution(
    point: &Point3,
    center: &Point3,
    axis_dir: &Vec3, // must be unit
    angle: f64,
) -> (Vec<Point3>, Vec<f64>, Vec<f64>) {
    let n_arcs = if angle <= FRAC_PI_2 + 1e-10 {
        1
    } else if angle <= std::f64::consts::PI + 1e-10 {
        2
    } else if angle <= 3.0 * FRAC_PI_2 + 1e-10 {
        3
    } else {
        4
    };

    let d_angle = angle / n_arcs as f64;
    let w_mid = (d_angle / 2.0).cos(); // weight for the mid control point of each arc

    let mut points = Vec::with_capacity(2 * n_arcs + 1);
    let mut weights = Vec::with_capacity(2 * n_arcs + 1);

    // Start point
    let mut p0 = *point;
    points.push(p0);
    weights.push(1.0);

    for i in 0..n_arcs {
        let angle_end = d_angle * (i as f64 + 1.0);

        let p2 = rotate_point(point, center, axis_dir, angle_end);

        // Mid control point: intersection of tangent lines from p0 and p2
        // Tangent at p0: T0 = axis_dir × (p0 - center)
        // Tangent at p2: T2 = axis_dir × (p2 - center)
        // P1 = p0 + s * T0 = p2 + t * T2
        // For a symmetric arc: P1 = (p0 + p2) / (2 * cos(d_angle/2)) projected
        // Actually, the standard formula: P1 lies at the intersection of the tangents.
        // For a circular arc with half-angle α = d_angle/2:
        //   P1 = center + (r / cos(α)) * direction_at_midangle
        // which equals: (P0 + P2) / (2 * cos(α)) relative to center

        let radial = p0 - *center;
        let r = radial.length();

        if r < 1e-15 {
            // Point is on the axis — degenerate case, all arc points collapse to the same point
            points.push(p0);
            weights.push(w_mid);
            points.push(p0);
            weights.push(1.0);
            p0 = p2;
            continue;
        }

        // Use tangent-line intersection approach
        let t0 = axis_dir.cross(&(p0 - *center));
        let t2 = axis_dir.cross(&(p2 - *center));

        // Solve p0 + s*t0 = p2 + t*t2 using least-squares on the 3 equations
        // Dot both sides with t0: (p0-p2)·t0 + s*(t0·t0) = t*(t2·t0)
        // Dot both sides with t2: (p0-p2)·t2 + s*(t0·t2) = t*(t2·t2)
        let dp = p0 - p2;
        let t0_dot_t0 = t0.dot(&t0);
        let t2_dot_t2 = t2.dot(&t2);
        let t0_dot_t2 = t0.dot(&t2);
        let dp_dot_t0 = dp.dot(&t0);

        let s = if t0_dot_t0.abs() > 1e-30 {
            // From first equation: s = (t*(t0·t2) - dp·t0) / (t0·t0)
            // From second: (dp·t2) + s*(t0·t2) = t*(t2·t2)
            // Substituting: det = t0·t0 * t2·t2 - (t0·t2)^2
            let det = t0_dot_t0 * t2_dot_t2 - t0_dot_t2 * t0_dot_t2;
            if det.abs() > 1e-30 {
                let dp_dot_t2 = dp.dot(&t2);
                (-dp_dot_t0 * t2_dot_t2 + dp_dot_t2 * t0_dot_t2) / det
            } else {
                // Fallback: use midpoint rotation approach
                1.0 / (2.0 * w_mid)
            }
        } else {
            0.0
        };

        let p1 = p0 + t0 * s;

        points.push(p1);
        weights.push(w_mid);
        points.push(p2);
        weights.push(1.0);

        p0 = p2;
    }

    // Build knot vector for the v direction
    let knots_v = match n_arcs {
        1 => vec![0.0, 0.0, 0.0, angle, angle, angle],
        2 => {
            let mid = angle / 2.0;
            vec![0.0, 0.0, 0.0, mid, mid, angle, angle, angle]
        }
        3 => {
            let a1 = angle / 3.0;
            let a2 = 2.0 * angle / 3.0;
            vec![0.0, 0.0, 0.0, a1, a1, a2, a2, angle, angle, angle]
        }
        4 => {
            let a1 = angle / 4.0;
            let a2 = angle / 2.0;
            let a3 = 3.0 * angle / 4.0;
            vec![
                0.0, 0.0, 0.0, a1, a1, a2, a2, a3, a3, angle, angle, angle,
            ]
        }
        _ => unreachable!(),
    };

    (points, weights, knots_v)
}

/// Construct a NURBS surface of revolution.
///
/// The profile curve defines the u-direction.
/// The revolution angle defines the v-direction (degree 2, rational).
///
/// # Arguments
/// * `profile` - The NURBS profile curve to revolve.
/// * `axis_origin` - A point on the revolution axis.
/// * `axis_dir` - Direction of the revolution axis (need not be unit length).
/// * `angle` - Revolution angle in radians, must be in (0, 2pi].
///
/// # Errors
/// Returns `CoreError::InvalidNurbs` if the angle is out of range or the axis
/// direction is zero-length.
pub fn revolution_surface(
    profile: &NurbsCurve,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    angle: f64,
) -> Result<NurbsSurface, CoreError> {
    if angle <= 0.0 || angle > TAU + 1e-10 {
        return Err(CoreError::InvalidNurbs(format!(
            "revolution angle must be in (0, 2*pi], got {}",
            angle
        )));
    }

    let axis_unit = axis_dir
        .normalized()
        .ok_or_else(|| CoreError::InvalidNurbs("axis direction is zero-length".into()))?;

    let angle = angle.min(TAU); // clamp to exactly TAU

    let n_profile = profile.control_points.len();

    // Build arc control points for the first profile control point to determine v structure
    let first_cp = &profile.control_points[0];
    let v_proj = (*first_cp - *axis_origin).dot(&axis_unit);
    let first_center = *axis_origin + axis_unit * v_proj;
    let (first_arc_pts, first_arc_wts, knots_v) =
        arc_control_points_for_revolution(first_cp, &first_center, &axis_unit, angle);
    let n_arc = first_arc_pts.len();

    // Build the full control point grid: n_profile (u) x n_arc (v)
    let mut control_points = Vec::with_capacity(n_profile * n_arc);
    let mut weights = Vec::with_capacity(n_profile * n_arc);

    // First row already computed
    let w_profile_0 = profile.weights[0];
    for j in 0..n_arc {
        control_points.push(first_arc_pts[j]);
        weights.push(w_profile_0 * first_arc_wts[j]);
    }

    // Remaining profile control points
    for i in 1..n_profile {
        let cp = &profile.control_points[i];
        let w_i = profile.weights[i];
        let proj = (*cp - *axis_origin).dot(&axis_unit);
        let center = *axis_origin + axis_unit * proj;
        let (arc_pts, arc_wts, _) =
            arc_control_points_for_revolution(cp, &center, &axis_unit, angle);

        for j in 0..n_arc {
            control_points.push(arc_pts[j]);
            weights.push(w_i * arc_wts[j]);
        }
    }

    Ok(NurbsSurface {
        degree_u: profile.degree,
        degree_v: 2,
        knots_u: profile.knots.clone(),
        knots_v,
        control_points,
        num_u: n_profile,
        num_v: n_arc,
        weights,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    #[test]
    fn rotate_90_degrees() {
        let p = Point3::new(1.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(0.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), FRAC_PI_2);
        assert!((r.x).abs() < 1e-10);
        assert!((r.y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rotate_360_returns_to_start() {
        let p = Point3::new(1.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(0.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), TAU);
        assert!((r.x - 1.0).abs() < 1e-10);
        assert!((r.y).abs() < 1e-10);
    }

    #[test]
    fn rotate_around_offset_axis() {
        let p = Point3::new(2.0, 0.0, 0.0);
        let r = rotate_point(&p, &Point3::new(1.0, 0.0, 0.0), &Vec3::new(0.0, 0.0, 1.0), PI);
        assert!((r.x).abs() < 1e-10);
        assert!((r.y).abs() < 1e-10);
    }

    #[test]
    fn revolution_curve_full_is_circle() {
        let c = revolution_curve(
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        assert!(matches!(c, CurveGeometry::Circle { .. }));
    }

    #[test]
    fn revolution_curve_partial_is_arc() {
        let c = revolution_curve(
            &Point3::new(1.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        );
        assert!(matches!(c, CurveGeometry::Arc { .. }));
    }

    #[test]
    fn revolution_curve_radius_correct() {
        let c = revolution_curve(
            &Point3::new(3.0, 0.0, 0.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        if let CurveGeometry::Circle { radius, .. } = c {
            assert!((radius - 3.0).abs() < 1e-10);
        } else {
            panic!("Expected Circle");
        }
    }

    #[test]
    fn revolution_curve_offset_point() {
        // Point at (2, 0, 3), axis is Z through origin. Should orbit at radius 2, center at (0,0,3)
        let c = revolution_curve(
            &Point3::new(2.0, 0.0, 3.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        );
        if let CurveGeometry::Circle { center, radius, .. } = c {
            assert!((center.z - 3.0).abs() < 1e-10);
            assert!((radius - 2.0).abs() < 1e-10);
        } else {
            panic!("Expected Circle");
        }
    }

    #[test]
    fn revolution_surface_cylinder() {
        use crate::geometry::surface::SurfaceEval;
        // Revolve a vertical line segment around Z axis → cylinder of radius 1
        let profile = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        };
        let surface = revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        )
        .expect("should construct cylinder");

        // Sample: all points should be at radius 1 from Z axis
        for u in [0.0, 0.5, 0.99] {
            for v_frac in [0.0, 0.25, 0.5, 0.75] {
                let (v_min, v_max) = (
                    surface.knots_v[surface.degree_v],
                    surface.knots_v[surface.knots_v.len() - surface.degree_v - 1],
                );
                let v = v_min + v_frac * (v_max - v_min);
                let p = surface.point_at(u, v);
                let r = (p.x * p.x + p.y * p.y).sqrt();
                assert!(
                    (r - 1.0).abs() < 0.15,
                    "Not on cylinder: r={r} at u={u}, v={v}"
                );
            }
        }
    }

    #[test]
    fn revolution_surface_partial_90() {
        use crate::geometry::surface::SurfaceEval;
        let profile = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        };
        let surface = revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            FRAC_PI_2,
        )
        .expect("should construct quarter-cylinder");

        // Should produce a quarter-cylinder
        let p_start = surface.point_at(0.0, surface.knots_v[surface.degree_v]);
        assert!(
            (p_start.x - 1.0).abs() < 0.15,
            "Start should be at x=1, got x={}",
            p_start.x
        );
    }

    #[test]
    fn revolution_surface_invalid_angle() {
        let profile = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
            weights: vec![1.0, 1.0],
        };
        assert!(revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            0.0
        )
        .is_err());
    }
}
