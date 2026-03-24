//! Helper functions for center curve computation.
//!
//! Provides surface projection, footprint clipping, NURBS fitting,
//! and parameter range clamping used by the [`super::center_curve`] module.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};

/// Clip footprint points to the face's UV parameter bounds.
///
/// Projects each 3D footprint point into the surface UV domain and retains
/// only those points whose UV coordinates fall within (or very close to)
/// `u_range` and `v_range`. This prevents trim curves from extending
/// beyond face boundaries, which would cause `split_face` to receive
/// mismatched anchor vertices.
///
/// If clipping would remove all points, the original footprint is returned
/// unchanged to avoid degenerate results.
pub fn clip_footprint_to_bounds(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    footprint: &[Point3],
) -> Vec<Point3> {
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    // Small tolerance for boundary inclusion.
    let eps = 1e-6;

    let clipped: Vec<Point3> = footprint
        .iter()
        .filter(|pt| {
            let (u, v) = match project_point_to_uv(surface, pt) {
                Some(uv) => uv,
                None => return true, // Keep points we can't project (conservative).
            };
            u >= u_range.0 - eps
                && u <= u_range.1 + eps
                && v >= v_range.0 - eps
                && v <= v_range.1 + eps
        })
        .copied()
        .collect();

    // Never return an empty footprint — fall back to unclipped.
    if clipped.len() < 2 {
        footprint.to_vec()
    } else {
        clipped
    }
}

/// Clamp a potentially infinite parameter range to a fallback finite range.
pub fn clamp_finite_range(range: (f64, f64), fallback: (f64, f64)) -> (f64, f64) {
    let lo = if range.0.is_finite() {
        range.0
    } else {
        fallback.0
    };
    let hi = if range.1.is_finite() {
        range.1
    } else {
        fallback.1
    };
    (lo, hi)
}

/// Build a degree-1 NURBS (polyline) through the given points.
///
/// Uses uniform chord-length parameterization. This is the simplest
/// possible interpolation; higher-degree fitting can be added later.
pub fn fit_nurbs_through_points(points: &[Point3]) -> Result<NurbsCurve, String> {
    let n = points.len();
    if n < 2 {
        return Err("need at least 2 points to fit a NURBS curve".into());
    }

    // Chord-length parameterization
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        params[i] = params[i - 1] + points[i].distance_to(&points[i - 1]);
    }
    let total = params[n - 1];
    if total < 1e-30 {
        return Err("degenerate: all points coincide".into());
    }
    for p in &mut params {
        *p /= total;
    }

    // Degree-1 NURBS: knots = clamped uniform, control_points = the points themselves.
    let degree = 1;
    let mut knots = Vec::with_capacity(n + degree + 1);
    knots.push(0.0);
    knots.push(0.0);
    for param in params.iter().take(n - 1).skip(1) {
        knots.push(*param);
    }
    knots.push(1.0);
    knots.push(1.0);

    Ok(NurbsCurve {
        degree,
        knots,
        control_points: points.to_vec(),
        weights: vec![1.0; n],
    })
}

/// Project a 3D point onto a surface by finding the closest surface point.
///
/// For planes, uses the analytic closed-form projection. For other surface
/// types, samples a coarse grid and picks the nearest point.
pub fn project_onto_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    point: &Point3,
) -> Result<Point3, String> {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => {
            // Project onto the plane: p_proj = p - ((p - origin) . n) * n
            let d = *point - *origin;
            let dist = d.dot(normal);
            let projected = *point - *normal * dist;

            // Clamp to parameter range
            let dp = projected - *origin;
            let u = dp.dot(u_axis);
            let v = dp.dot(v_axis);
            let u_clamped = u.clamp(u_range.0, u_range.1);
            let v_clamped = v.clamp(v_range.0, v_range.1);
            Ok(*origin + *u_axis * u_clamped + *v_axis * v_clamped)
        }
        _ => {
            // Coarse grid search for general surfaces.
            let n_u = 20;
            let n_v = 20;
            let mut best_pt =
                surface.point_at((u_range.0 + u_range.1) * 0.5, (v_range.0 + v_range.1) * 0.5);
            let mut best_dist = best_pt.distance_to(point);

            for iu in 0..n_u {
                let u =
                    u_range.0 + (u_range.1 - u_range.0) * (iu as f64) / ((n_u - 1).max(1) as f64);
                for iv in 0..n_v {
                    let v = v_range.0
                        + (v_range.1 - v_range.0) * (iv as f64) / ((n_v - 1).max(1) as f64);
                    let sp = surface.point_at(u, v);
                    let d = sp.distance_to(point);
                    if d < best_dist {
                        best_dist = d;
                        best_pt = sp;
                    }
                }
            }

            Ok(best_pt)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::curve::CurveEval;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn fit_nurbs_through_points_simple() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ];
        let curve = fit_nurbs_through_points(&pts).unwrap();
        assert_eq!(curve.degree, 1);
        assert_eq!(curve.control_points.len(), 3);

        // Endpoints should match
        let (t0, t1) = curve.parameter_range();
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        assert!(p0.distance_to(&pts[0]) < 1e-10);
        assert!(p1.distance_to(&pts[2]) < 1e-10);
    }

    #[test]
    fn fit_nurbs_too_few_points() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0)];
        assert!(fit_nurbs_through_points(&pts).is_err());
    }

    #[test]
    fn project_onto_plane() {
        let plane = xy_plane();
        let pt = Point3::new(1.0, 2.0, 5.0);
        let proj = project_onto_surface(&plane, (-10.0, 10.0), (-10.0, 10.0), &pt).unwrap();
        assert!((proj.z).abs() < 1e-10, "should project to z=0");
        assert!((proj.x - 1.0).abs() < 1e-10);
        assert!((proj.y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let pt = Point3::new(3.0, 0.0, 0.0);
        let proj = project_onto_surface(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (0.0, std::f64::consts::PI),
            &pt,
        )
        .unwrap();
        // Should be on the sphere, closest to (3,0,0) -> (1,0,0)
        let r = (proj.x * proj.x + proj.y * proj.y + proj.z * proj.z).sqrt();
        assert!(
            (r - 1.0).abs() < 0.2,
            "projected point should be on sphere, r={}",
            r
        );
    }
}
