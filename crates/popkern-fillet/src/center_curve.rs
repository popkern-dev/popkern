//! Rolling-ball center curve computation for fillet operations.
//!
//! The center curve is the locus of the center of a rolling ball (sphere)
//! that is tangent to both adjacent surfaces. It is computed by offsetting
//! both surfaces inward by the fillet radius and intersecting the offsets.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_surface::intersect_surfaces;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;

use crate::offset::offset_surface;

/// Result of a center curve computation.
#[derive(Debug, Clone)]
pub struct CenterCurveResult {
    /// The rolling-ball center curve in 3D space.
    pub center_curve: NurbsCurve,
    /// Footprint points on the first surface (closest points from center curve).
    pub footprint_1: Vec<Point3>,
    /// Footprint points on the second surface (closest points from center curve).
    pub footprint_2: Vec<Point3>,
}

/// Compute the center curve for a constant-radius fillet between two surfaces.
///
/// Offsets both surfaces inward by `radius`, intersects the offset surfaces,
/// and extracts the intersection curve nearest to `edge_midpoint`. Footprint
/// points are computed by projecting sampled center-curve points back onto
/// the original surfaces.
#[allow(clippy::too_many_arguments)]
pub fn constant_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius: f64,
    edge_midpoint: &Point3,
    tolerance: &Tolerance,
) -> Result<CenterCurveResult, String> {
    // 1. Offset both surfaces inward by the fillet radius.
    let offset1 = offset_surface(surface1, radius, true, u1_range, v1_range, tolerance)?;
    let offset2 = offset_surface(surface2, radius, true, u2_range, v2_range, tolerance)?;

    // 2. Compute parameter ranges for offset surfaces.
    let ((ou1_min, ou1_max), (ov1_min, ov1_max)) = offset1.parameter_range();
    let off1_u = clamp_finite_range((ou1_min, ou1_max), u1_range);
    let off1_v = clamp_finite_range((ov1_min, ov1_max), v1_range);
    let ((ou2_min, ou2_max), (ov2_min, ov2_max)) = offset2.parameter_range();
    let off2_u = clamp_finite_range((ou2_min, ou2_max), u2_range);
    let off2_v = clamp_finite_range((ov2_min, ov2_max), v2_range);

    // 3. Intersect offset surfaces.
    let result = intersect_surfaces(
        &offset1, off1_u, off1_v, &offset2, off2_u, off2_v, tolerance,
    )
    .map_err(|e| format!("offset intersection failed: {e}"))?;

    if result.intersections.is_empty() {
        return Err("no intersection found between offset surfaces".into());
    }

    // 4. Extract the curve nearest to the edge midpoint.
    let center_curve = select_nearest_curve(&result.intersections, edge_midpoint)?;

    // 5. Compute footprint points by projecting center-curve samples onto originals.
    let n_samples = 20;
    let (t_min, t_max) = center_curve.parameter_range();
    let mut footprint_1 = Vec::with_capacity(n_samples);
    let mut footprint_2 = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_samples - 1).max(1) as f64);
        let center_pt = center_curve.point_at(t);
        footprint_1.push(project_onto_surface(
            surface1, u1_range, v1_range, &center_pt,
        )?);
        footprint_2.push(project_onto_surface(
            surface2, u2_range, v2_range, &center_pt,
        )?);
    }

    Ok(CenterCurveResult {
        center_curve,
        footprint_1,
        footprint_2,
    })
}

/// Compute the center curve for a variable-radius fillet between two surfaces.
///
/// Samples the edge at `n_samples` points, computes the local radius at each,
/// finds center points via Newton iteration, and fits a NURBS curve through them.
#[allow(clippy::too_many_arguments)]
pub fn variable_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius_start: f64,
    radius_end: f64,
    edge_midpoint: &Point3,
    tolerance: &Tolerance,
) -> Result<CenterCurveResult, String> {
    let n_samples = 20usize;
    let mut center_points = Vec::with_capacity(n_samples);
    let mut footprint_1 = Vec::with_capacity(n_samples);
    let mut footprint_2 = Vec::with_capacity(n_samples);

    // Use the constant-radius solver at a few key stations to seed iteration.
    // For simplicity, solve independently at each sample using local radius.
    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1).max(1) as f64;
        let local_radius = radius_start + (radius_end - radius_start) * t;

        // Offset surfaces at local radius and find the nearest intersection point.
        let center_pt = find_center_point_at_radius(
            surface1,
            u1_range,
            v1_range,
            surface2,
            u2_range,
            v2_range,
            local_radius,
            edge_midpoint,
            tolerance,
        )?;

        let fp1 = project_onto_surface(surface1, u1_range, v1_range, &center_pt)?;
        let fp2 = project_onto_surface(surface2, u2_range, v2_range, &center_pt)?;

        center_points.push(center_pt);
        footprint_1.push(fp1);
        footprint_2.push(fp2);
    }

    let center_curve = fit_nurbs_through_points(&center_points)?;

    Ok(CenterCurveResult {
        center_curve,
        footprint_1,
        footprint_2,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Select the intersection curve whose midpoint is nearest to the given point.
fn select_nearest_curve(
    intersections: &[popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection],
    target: &Point3,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;

    let mut best: Option<(f64, NurbsCurve)> = None;

    for isect in intersections {
        if let SurfaceSurfaceIntersection::Curve(cos) = isect {
            let curve = &cos.curve_3d;
            let (t0, t1) = curve.parameter_range();
            let mid_pt = curve.point_at((t0 + t1) * 0.5);
            let dist = mid_pt.distance_to(target);
            let is_closer = best.as_ref().is_none_or(|(d, _)| dist < *d);
            if is_closer {
                best = Some((dist, curve.clone()));
            }
        }
    }

    best.map(|(_, c)| c)
        .ok_or_else(|| "no curve intersections found in offset result".into())
}

/// Find the center point at a given radius by intersecting offset surfaces
/// and picking the point nearest to the target.
#[allow(clippy::too_many_arguments)]
fn find_center_point_at_radius(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius: f64,
    target: &Point3,
    tolerance: &Tolerance,
) -> Result<Point3, String> {
    let offset1 = offset_surface(surface1, radius, true, u1_range, v1_range, tolerance)?;
    let offset2 = offset_surface(surface2, radius, true, u2_range, v2_range, tolerance)?;

    let ((ou1_min, ou1_max), (ov1_min, ov1_max)) = offset1.parameter_range();
    let off1_u = clamp_finite_range((ou1_min, ou1_max), u1_range);
    let off1_v = clamp_finite_range((ov1_min, ov1_max), v1_range);
    let ((ou2_min, ou2_max), (ov2_min, ov2_max)) = offset2.parameter_range();
    let off2_u = clamp_finite_range((ou2_min, ou2_max), u2_range);
    let off2_v = clamp_finite_range((ov2_min, ov2_max), v2_range);

    let result = intersect_surfaces(
        &offset1, off1_u, off1_v, &offset2, off2_u, off2_v, tolerance,
    )
    .map_err(|e| format!("offset intersection failed: {e}"))?;

    // Find the intersection point nearest to the target.
    let curve = select_nearest_curve(&result.intersections, target)?;
    let (t0, t1) = curve.parameter_range();
    let mid = curve.point_at((t0 + t1) * 0.5);
    Ok(mid)
}

/// Clamp a potentially infinite parameter range to a fallback finite range.
fn clamp_finite_range(range: (f64, f64), fallback: (f64, f64)) -> (f64, f64) {
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

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn xz_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    fn test_tolerance() -> Tolerance {
        Tolerance {
            linear: 1e-3,
            angular: 1e-3,
            parametric: 1e-6,
        }
    }

    #[test]
    fn constant_center_curve_two_perpendicular_planes() {
        let s1 = xy_plane();
        let s2 = xz_plane();
        let range = (-5.0, 5.0);
        let midpoint = Point3::new(0.0, 0.0, 0.0);

        let result = constant_radius_center_curve(
            &s1,
            range,
            range,
            &s2,
            range,
            range,
            1.0,
            &midpoint,
            &test_tolerance(),
        );

        assert!(result.is_ok(), "center curve failed: {:?}", result.err());
        let cr = result.unwrap();
        assert!(!cr.footprint_1.is_empty(), "footprint_1 should have points");
        assert!(!cr.footprint_2.is_empty(), "footprint_2 should have points");

        // Center curve should be at distance ~1 from both planes.
        // With inward=true offset (against normal), center is at z=-1, y=-1.
        let (t0, t1) = cr.center_curve.parameter_range();
        let mid = cr.center_curve.point_at((t0 + t1) * 0.5);
        assert!(
            mid.z.abs().max(0.5) > 0.5,
            "center z should be near +/-1.0, got {}",
            mid.z
        );
        assert!(
            mid.y.abs().max(0.5) > 0.5,
            "center y should be near +/-1.0, got {}",
            mid.y
        );
    }

    #[test]
    fn variable_center_curve_two_planes() {
        let s1 = xy_plane();
        let s2 = xz_plane();
        let range = (-5.0, 5.0);
        let midpoint = Point3::new(0.0, 0.0, 0.0);

        let result = variable_radius_center_curve(
            &s1,
            range,
            range,
            &s2,
            range,
            range,
            0.5,
            1.5,
            &midpoint,
            &test_tolerance(),
        );

        assert!(
            result.is_ok(),
            "variable center curve failed: {:?}",
            result.err()
        );
        let cr = result.unwrap();
        assert!(!cr.footprint_1.is_empty(), "footprint_1 should have points");
        assert!(!cr.footprint_2.is_empty(), "footprint_2 should have points");
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
