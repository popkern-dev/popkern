//! Trim curve computation for fillet operations.
//!
//! Projects 3D footprint points into surface UV parameter space, then fits
//! both 3D and 2D (parametric) NURBS curves through the projected points.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3};

use crate::center_curve::fit_nurbs_through_points;

/// Result of trim curve computation for both sides of a fillet.
#[derive(Debug, Clone)]
pub struct TrimCurves {
    /// 3D trim curve on surface 1 (through footprint_1 points).
    pub curve_3d_1: NurbsCurve,
    /// Parameter-space (UV) trim curve on surface 1.
    /// Uses `Point3(u, v, 0.0)` convention.
    pub pcurve_1: NurbsCurve,
    /// 3D trim curve on surface 2 (through footprint_2 points).
    pub curve_3d_2: NurbsCurve,
    /// Parameter-space (UV) trim curve on surface 2.
    /// Uses `Point3(u, v, 0.0)` convention.
    pub pcurve_2: NurbsCurve,
    /// Exact circle through footprint_1 `(center, normal, radius)`
    /// when the contacts are circular (closed rims) — the NURBS fit
    /// above carries ~3e-3 interpolation error; consumers needing
    /// watertight-grade trim geometry (closed-fillet edges) use
    /// this.
    pub circle_1: Option<(Point3, popkern_core::geometry::Vec3, f64)>,
    /// Exact circle through footprint_2; see `circle_1`.
    pub circle_2: Option<(Point3, popkern_core::geometry::Vec3, f64)>,
}

/// Compute trim curves from footprint points projected onto two surfaces.
///
/// For analytic surface types (plane, cylinder, sphere, cone), UV projection
/// uses the closed-form inverse parametrisation. For NURBS and torus surfaces,
/// a coarse-grid Newton search is used as a fallback.
pub fn compute_trim_curves(
    surface1: &SurfaceGeometry,
    footprint_1: &[Point3],
    surface2: &SurfaceGeometry,
    footprint_2: &[Point3],
) -> Result<TrimCurves, String> {
    if footprint_1.len() < 2 || footprint_2.len() < 2 {
        return Err("need at least 2 footprint points per surface".into());
    }

    // 1. Project footprint points into UV space.
    let uv1 = project_footprint(surface1, footprint_1)?;
    let uv2 = project_footprint(surface2, footprint_2)?;

    // 2. Fit 3D NURBS curves through the footprint points.
    let curve_3d_1 = fit_nurbs_through_points(footprint_1)?;
    let curve_3d_2 = fit_nurbs_through_points(footprint_2)?;

    // 3. Fit UV pcurves (using Point3(u, v, 0) convention).
    let pcurve_1 = fit_nurbs_through_points(&uv1)?;
    let pcurve_2 = fit_nurbs_through_points(&uv2)?;

    let circle_of = |fp: &[Point3]| -> Option<(Point3, popkern_core::geometry::Vec3, f64)> {
        let (center, radius) = crate::fillet_surface_analytic::fit_exact_circle_center(fp)?;
        // A circle through (near-)collinear contacts is no circle: a
        // straight edge's footprint fitted an 8e15-radius "arc" that
        // reached the trims once exact projections stopped failing on
        // its far points (PK-737).
        let extent = fp
            .iter()
            .map(|p| p.distance_to(&fp[0]))
            .fold(0.0_f64, f64::max);
        if radius > 1e3 * extent {
            return None;
        }
        let mut nx = 0.0;
        let mut ny = 0.0;
        let mut nz = 0.0;
        for i in 0..fp.len() {
            let a = fp[i];
            let b = fp[(i + 1) % fp.len()];
            nx += (a.y - b.y) * (a.z + b.z);
            ny += (a.z - b.z) * (a.x + b.x);
            nz += (a.x - b.x) * (a.y + b.y);
        }
        let normal = popkern_core::geometry::Vec3::new(nx, ny, nz).normalized()?;
        Some((center, normal, radius))
    };
    let circle_1 = circle_of(footprint_1);
    let circle_2 = circle_of(footprint_2);

    Ok(TrimCurves {
        curve_3d_1,
        pcurve_1,
        curve_3d_2,
        pcurve_2,
        circle_1,
        circle_2,
    })
}

/// Project a sequence of 3D points onto a surface's UV parameter space.
///
/// Uses the analytic projection from `popkern_core::geometry::uv_projection`
/// when available, falling back to a grid-search approach for surface types
/// that lack an analytic inverse (torus, NURBS).
fn project_footprint(surface: &SurfaceGeometry, points: &[Point3]) -> Result<Vec<Point3>, String> {
    let mut uv_points = Vec::with_capacity(points.len());

    for pt in points {
        let (u, v) = match project_point_to_uv(surface, pt) {
            Some(uv) => uv,
            None => project_uv_grid_search(surface, pt)?,
        };
        uv_points.push(Point3::new(u, v, 0.0));
    }

    Ok(uv_points)
}

/// Fallback UV projection via grid search for surfaces without analytic inverse.
///
/// Samples the surface on a coarse grid and picks the (u, v) whose surface
/// point is closest to the query point.
fn project_uv_grid_search(surface: &SurfaceGeometry, point: &Point3) -> Result<(f64, f64), String> {
    use popkern_core::geometry::surface::SurfaceEval;

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // Use finite fallback ranges if the surface has infinite extent.
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    let n = 30;
    let mut best_u = (u_lo + u_hi) * 0.5;
    let mut best_v = (v_lo + v_hi) * 0.5;
    let mut best_dist = f64::INFINITY;

    for iu in 0..n {
        let u = u_lo + (u_hi - u_lo) * (iu as f64) / ((n - 1).max(1) as f64);
        for iv in 0..n {
            let v = v_lo + (v_hi - v_lo) * (iv as f64) / ((n - 1).max(1) as f64);
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    Ok((best_u, best_v))
}

/// A torus strip's contact circles, exactly (PK-738): each is a
/// v-iso circle of the torus, and where the support is a coaxial
/// cylinder or a plane normal to the axis it is known in closed form —
/// a coaxial cylinder of radius `Rc` touches the tube on the circle
/// centred at the torus centre, radius `Rc`; a plane at axial height
/// `h` touches it on the circle centred `h` along the axis, radius `R`.
/// The fitted circle (`circle_of`) carried roundoff in its centre
/// (z = 1.7999999999999987 for an exact 1.8): a boundary lying on the
/// torus's v = 0 seam then sat a hair either side of it by chance, and
/// a later boolean read the strip's chart as wound (a loop left open by
/// 0.87). A fitted circle is replaced only when the exact one agrees
/// with it within 1e-6 of the torus size.
pub(crate) fn exact_torus_contact_circles(
    trims: &mut TrimCurves,
    strip: &SurfaceGeometry,
    analysis: &crate::analysis::EdgeAnalysis,
) {
    let SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = *strip
    else {
        return;
    };
    let Some(a) = axis.normalized() else {
        return;
    };
    let scale = major_radius + minor_radius;
    let exact = |support: &SurfaceGeometry| -> Option<(Point3, popkern_core::geometry::Vec3, f64)> {
        match *support {
            SurfaceGeometry::Plane { origin, normal, .. } => {
                let n = normal.normalized()?;
                if n.cross(&a).length() > 1e-9 {
                    return None;
                }
                let h = (origin - center).dot(&a);
                ((h.abs() - minor_radius).abs() <= 1e-6 * scale)
                    .then(|| (center + a * h, axis, major_radius))
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis: c_axis,
                radius,
            } => {
                let ca = c_axis.normalized()?;
                let off = center - origin;
                let off_axis = (off - ca * off.dot(&ca)).length();
                (ca.cross(&a).length() <= 1e-9 && off_axis <= 1e-9 * scale)
                    .then_some((center, axis, radius))
            }
            _ => None,
        }
    };
    for (circle, support) in [
        (&mut trims.circle_1, &analysis.surface_1),
        (&mut trims.circle_2, &analysis.surface_2),
    ] {
        let (Some(fitted), Some(e)) = (*circle, exact(support)) else {
            continue;
        };
        if fitted.0.distance_to(&e.0) <= 1e-6 * scale && (fitted.2 - e.2).abs() <= 1e-6 * scale {
            *circle = Some(e);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::{Point3, Vec3};

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn trim_curves_on_plane() {
        let plane = xy_plane();
        let footprint = vec![
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(3.0, 4.0, 0.0),
            Point3::new(5.0, 6.0, 0.0),
        ];

        let uv_pts = project_footprint(&plane, &footprint).unwrap();
        assert_eq!(uv_pts.len(), 3);

        // UV should match (x, y) for z=0 plane.
        assert!((uv_pts[0].x - 1.0).abs() < 1e-10, "u[0] should be 1.0");
        assert!((uv_pts[0].y - 2.0).abs() < 1e-10, "v[0] should be 2.0");
        assert!((uv_pts[1].x - 3.0).abs() < 1e-10, "u[1] should be 3.0");
        assert!((uv_pts[1].y - 4.0).abs() < 1e-10, "v[1] should be 4.0");
        assert!((uv_pts[2].x - 5.0).abs() < 1e-10, "u[2] should be 5.0");
        assert!((uv_pts[2].y - 6.0).abs() < 1e-10, "v[2] should be 6.0");
    }

    #[test]
    fn trim_curves_fit() {
        let plane = xy_plane();
        let fp1 = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let fp2 = vec![
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ];

        let result = compute_trim_curves(&plane, &fp1, &plane, &fp2);
        assert!(
            result.is_ok(),
            "compute_trim_curves failed: {:?}",
            result.err()
        );
        let tc = result.unwrap();

        // 3D trim curves should pass through the footprint points.
        let (t0, t1) = tc.curve_3d_1.parameter_range();
        let p_start = tc.curve_3d_1.point_at(t0);
        let p_end = tc.curve_3d_1.point_at(t1);
        assert!(
            p_start.distance_to(&fp1[0]) < 1e-6,
            "3D curve should start at fp1[0]"
        );
        assert!(
            p_end.distance_to(&fp1[2]) < 1e-6,
            "3D curve should end at fp1[2]"
        );

        // Pcurves should have z=0 everywhere.
        let (pt0, _pt1) = tc.pcurve_1.parameter_range();
        let uv_start = tc.pcurve_1.point_at(pt0);
        assert!(
            uv_start.z.abs() < 1e-10,
            "pcurve z should be 0, got {}",
            uv_start.z
        );

        // For points on the xy-plane, pcurve u should match x, v should match y.
        assert!(
            (uv_start.x - fp1[0].x).abs() < 1e-6,
            "pcurve u should match x"
        );
        assert!(
            (uv_start.y - fp1[0].y).abs() < 1e-6,
            "pcurve v should match y"
        );
    }

    #[test]
    fn trim_curves_too_few_points() {
        let plane = xy_plane();
        let fp = vec![Point3::new(0.0, 0.0, 0.0)];
        let fp2 = vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
        assert!(compute_trim_curves(&plane, &fp, &plane, &fp2).is_err());
    }

    #[test]
    fn grid_search_fallback_for_torus() {
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 3.0,
            minor_radius: 1.0,
        };
        // A point on the torus at u=0, v=0 should be at (major+minor, 0, 0) = (4, 0, 0).
        let pt = Point3::new(4.0, 0.0, 0.0);
        let (u, v) = project_uv_grid_search(&torus, &pt).unwrap();
        // Should find something reasonable (not crash).
        assert!(u.is_finite());
        assert!(v.is_finite());
    }
}
