//! Fillet surface construction.
//!
//! Given a center curve and footprint points, builds the rolling-ball fillet
//! surface. For two-plane cases an analytic torus/cylinder section is returned;
//! otherwise a NURBS swept surface is constructed by skinning cross-section arcs.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use crate::center_curve::CenterCurveResult;

/// Result of a fillet surface construction.
#[derive(Debug, Clone)]
pub struct FilletSurfaceResult {
    /// The constructed fillet surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on the first surface (3D curve along footprint_1).
    pub trim_curve_1: NurbsCurve,
    /// Trim curve on the second surface (3D curve along footprint_2).
    pub trim_curve_2: NurbsCurve,
}

/// Build the fillet surface from center-curve data and the two adjacent surfaces.
///
/// Attempts an analytic construction for the plane-plane case, falling back
/// to a NURBS swept surface for general surface pairs.
pub fn build_fillet_surface(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
) -> Result<FilletSurfaceResult, String> {
    // Try analytic torus for plane-plane case.
    if let (SurfaceGeometry::Plane { normal: n1, .. }, SurfaceGeometry::Plane { normal: n2, .. }) =
        (surface1, surface2)
        && let Some(result) = try_plane_plane_torus(n1, n2, center_result, radius)
    {
        return Ok(result);
    }

    // General case: NURBS swept surface.
    build_nurbs_sweep(center_result, radius, arc_samples)
}

/// Try to build an analytic torus section for a plane-plane fillet.
///
/// Returns `None` if the geometry is degenerate (e.g., parallel planes).
fn try_plane_plane_torus(
    n1: &Vec3,
    n2: &Vec3,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    // Compute the dihedral angle between the two planes.
    let cos_angle = n1.dot(n2).clamp(-1.0, 1.0);
    let angle = cos_angle.acos(); // angle between normals

    // The half-angle of the dihedral measured from the inside.
    let half_angle = angle * 0.5;
    let sin_half = half_angle.sin();

    // If planes are nearly parallel, can't build a torus fillet.
    if sin_half.abs() < 1e-10 {
        return None;
    }

    // Torus axis = cross product of the two normals (edge direction).
    let axis = n1.cross(n2);
    let axis = axis.normalized()?;

    // Center of the torus: the fillet center curve midpoint.
    let fp1 = center_result.footprint_1.first()?;
    let fp2 = center_result.footprint_2.first()?;
    let center = Point3::new(
        (fp1.x + fp2.x) * 0.5,
        (fp1.y + fp2.y) * 0.5,
        (fp1.z + fp2.z) * 0.5,
    );

    // For a fillet between two planes at angle `angle`:
    // Major radius = r / sin(half_angle) for a general dihedral.
    // For 90-degree planes (angle = PI/2), major_radius = r / sin(PI/4) = r * sqrt(2).
    // However, for a simple rolling-ball fillet the surface is always a
    // cylinder section when angle = PI/2, or a torus section otherwise.
    let major_radius = radius / sin_half;

    let surface = SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius: radius,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}

/// Build a NURBS swept surface by skinning cross-section arcs along the center curve.
///
/// Each cross-section is a polyline arc from footprint_1[i] to footprint_2[i],
/// centered at the corresponding center-curve sample. The arcs are sampled
/// at `arc_samples` points and skinned into a tensor-product surface.
fn build_nurbs_sweep(
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
) -> Result<FilletSurfaceResult, String> {
    let n_sections = center_result.footprint_1.len();
    if n_sections < 2 {
        return Err("need at least 2 cross-sections for sweep".into());
    }
    let n_arc = arc_samples.max(3); // at least 3 points per arc

    // Sample center curve at n_sections points.
    let (t_min, t_max) =
        popkern_core::geometry::curve::CurveEval::parameter_range(&center_result.center_curve);

    let mut control_points = Vec::with_capacity(n_sections * n_arc);

    for i in 0..n_sections {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_sections - 1).max(1) as f64);
        let center =
            popkern_core::geometry::curve::CurveEval::point_at(&center_result.center_curve, t);
        let fp1 = &center_result.footprint_1[i];
        let fp2 = &center_result.footprint_2[i];

        // Build the arc from fp1 to fp2 centered at `center`.
        let arc_pts = sample_arc(&center, fp1, fp2, radius, n_arc);
        control_points.extend_from_slice(&arc_pts);
    }

    // Build the NurbsSurface: U = along edge (n_sections), V = around arc (n_arc).
    let degree_u = 1.min(n_sections - 1);
    let degree_v = 1.min(n_arc - 1);

    let knots_u = build_clamped_uniform_knots(n_sections, degree_u);
    let knots_v = build_clamped_uniform_knots(n_arc, degree_v);
    let weights = vec![1.0; n_sections * n_arc];

    let nurbs = NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n_sections,
        num_v: n_arc,
        weights,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1)?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2)?;

    Ok(FilletSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}

/// Sample a circular arc from `start` to `end` centered at `center`.
///
/// The arc lies in the plane defined by the three points. If `start` and
/// `end` are equidistant from `center` at approximately `radius`, the arc
/// sweeps the shorter angle between them.
fn sample_arc(
    center: &Point3,
    start: &Point3,
    end: &Point3,
    _radius: f64,
    n: usize,
) -> Vec<Point3> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();

    // Handle degenerate cases.
    if r1 < 1e-15 || r2 < 1e-15 || n < 2 {
        return vec![*start; n];
    }

    // Normalize to unit vectors.
    let u = Vec3::new(v1.x / r1, v1.y / r1, v1.z / r1);
    let w_raw = v2 - u * v2.dot(&u);
    let w_len = w_raw.length();

    if w_len < 1e-15 {
        // start and end are nearly coincident or antipodal — just linearly interpolate.
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1).max(1) as f64;
            pts.push(Point3::new(
                start.x + (end.x - start.x) * t,
                start.y + (end.y - start.y) * t,
                start.z + (end.z - start.z) * t,
            ));
        }
        return pts;
    }

    let v = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);

    // Compute sweep angle.
    let cos_theta = (v1.dot(&v2) / (r1 * r2)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    // Use average radius for the arc.
    let r = (r1 + r2) * 0.5;

    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1).max(1) as f64;
        let angle = theta * t;
        let x = r * angle.cos();
        let y = r * angle.sin();
        pts.push(Point3::new(
            center.x + u.x * x + v.x * y,
            center.y + u.y * x + v.y * y,
            center.z + u.z * x + v.z * y,
        ));
    }

    pts
}

/// Build clamped uniform knot vector for `n` control points at given degree.
fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_knots = n + degree + 1;
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::center_curve::CenterCurveResult;
    use popkern_core::geometry::{Point3, Vec3};

    fn simple_center_result() -> CenterCurveResult {
        let fp1 = (0..5)
            .map(|i| Point3::new(i as f64, 0.0, 0.0))
            .collect::<Vec<_>>();
        let fp2 = (0..5)
            .map(|i| Point3::new(i as f64, 0.0, 1.0))
            .collect::<Vec<_>>();
        let center_pts: Vec<Point3> = (0..5).map(|i| Point3::new(i as f64, 0.0, 0.5)).collect();
        let center_curve = crate::center_curve::fit_nurbs_through_points(&center_pts).unwrap();
        CenterCurveResult {
            center_curve,
            footprint_1: fp1,
            footprint_2: fp2,
        }
    }

    #[test]
    fn fillet_surface_plane_plane() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let cr = simple_center_result();

        let result = build_fillet_surface(&s1, &s2, &cr, 1.0, 5);
        assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
        let fs = result.unwrap();
        // Plane-plane should produce a torus.
        assert!(
            matches!(fs.surface, SurfaceGeometry::Torus { .. }),
            "expected torus for plane-plane, got {:?}",
            std::mem::discriminant(&fs.surface)
        );
    }

    #[test]
    fn fillet_surface_general() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 5.0,
        };
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cr = simple_center_result();

        let result = build_fillet_surface(&sphere, &plane, &cr, 1.0, 5);
        assert!(result.is_ok(), "fillet surface failed: {:?}", result.err());
        let fs = result.unwrap();
        // Non-plane-plane should produce NURBS.
        assert!(
            matches!(fs.surface, SurfaceGeometry::Nurbs(_)),
            "expected NURBS for general case"
        );
    }

    #[test]
    fn sample_arc_basic() {
        let center = Point3::new(0.0, 0.0, 0.0);
        let start = Point3::new(1.0, 0.0, 0.0);
        let end = Point3::new(0.0, 1.0, 0.0);
        let pts = sample_arc(&center, &start, &end, 1.0, 5);
        assert_eq!(pts.len(), 5);

        // First and last should be close to start and end.
        assert!(pts[0].distance_to(&start) < 1e-10);
        let end_dist = pts[4].distance_to(&end);
        assert!(
            end_dist < 0.2,
            "last arc point should be near end, dist={end_dist}"
        );
    }

    #[test]
    fn build_clamped_knots_basic() {
        let knots = build_clamped_uniform_knots(4, 1);
        assert_eq!(knots.len(), 4 + 1 + 1); // n + degree + 1
        assert!((knots[0] - 0.0).abs() < 1e-15);
        assert!((knots[knots.len() - 1] - 1.0).abs() < 1e-15);
    }
}
