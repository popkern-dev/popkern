//! Analytic fast-path surface producers for `build_fillet_surface`.
//!
//! Each `try_*` function matches a specific surface pair and returns
//! `Some(FilletSurfaceResult)` when it can produce an **exact**
//! analytic surface (torus, cylinder, …), or `None` otherwise so the
//! caller falls through to the generic NURBS sweep.
//!
//! Split out of `fillet_surface.rs` to keep both files under the
//! 500-line structural cap and to give the §11.1 analytic-path
//! catalogue a dedicated home as more cases land.
//!
//! ## §11.1 catalogue
//!
//! | Pair | Surface | Function |
//! |---|---|---|
//! | plane-plane (dihedral) | torus | [`try_plane_plane_torus`] |
//! | plane-sphere | torus | [`try_circular_center_torus`] |
//! | plane-cylinder (⊥ axis) | torus | [`try_circular_center_torus`] |
//! | plane-cone (⊥ axis) | torus | [`try_circular_center_torus`] |
//! | cylinder-cylinder (∥ axes) | cylinder | [`try_parallel_cylinder_cylinder`] |
//!
//! **Plane-cone reasoning.** A plane perpendicular to a cone's axis
//! cuts the cone in a circle; by rotational symmetry the rolling-ball
//! center curve is a 3D circle too, so the fillet surface is an
//! exact torus. Oblique plane-cone intersections (ellipse, parabola,
//! hyperbola) produce canal surfaces that popkern has no primitive
//! for, so the detector gates on `|plane_normal · cone_axis| > 0.999`
//! up front and leaves non-perpendicular cases to the generic NURBS
//! sweep.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface::FilletSurfaceResult;

/// Try to build an analytic torus section for a plane-plane fillet.
///
/// Returns `None` if the geometry is degenerate (e.g., parallel planes).
pub(crate) fn try_plane_plane_torus(
    n1: &Vec3,
    n2: &Vec3,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    // Compute the dihedral angle between the two planes.
    let cos_angle = n1.dot(n2).clamp(-1.0, 1.0);
    let angle = cos_angle.acos();
    let half_angle = angle * 0.5;
    let sin_half = half_angle.sin();
    if sin_half.abs() < 1e-10 {
        return None;
    }

    // Torus axis = cross product of the two normals (edge direction).
    let axis = n1.cross(n2);
    let axis = axis.normalized()?;

    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let spine_mid = center_result.center_curve.point_at((t0 + t1) * 0.5);

    let major_radius = radius / sin_half;

    // The torus centre must be on the axis, not at the spine
    // midpoint. The spine is on the major circle at distance
    // `major_radius` from the axis. Project back to the axis
    // via the bisector: for inward-pointing normals, (n1+n2)
    // points from the axis toward the spine.
    let bisector = (*n1 + *n2).normalized().unwrap_or(*n1);
    let center = spine_mid - bisector * major_radius;

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

/// Try to build an analytic torus for a fillet whose center curve is
/// a 3D circle (§11.1 `plane-sphere`, `plane-cylinder`/`plane-cone`
/// with the plane perpendicular to the axis, and any other surface
/// pair with a circular center curve).
///
/// **Geometry.** Whenever the rolling-ball center locus is a 3D
/// circle, the resulting fillet surface is an exact torus:
///
/// - **Axis** = normal to the plane containing the circle
/// - **Center** = centroid of the circular center curve
/// - **Major radius** = radius of the circular center curve
/// - **Minor radius** = the fillet radius
///
/// For a **plane-sphere** fillet the center curve is always a circle
/// (intersection of a parallel plane with a concentric sphere).
///
/// For a **plane-cylinder** or **plane-cone** fillet the center curve
/// is a circle iff the plane is perpendicular to the curved surface's
/// axis. Oblique planes cut a cylinder in an ellipse and a cone in
/// an ellipse / parabola / hyperbola — none of which produce a
/// circular center curve — so the detector gates on
/// `|plane_normal · surface_axis| > 0.999` up front.
///
/// Rather than derive the torus center + major radius from raw
/// sphere/cylinder/cone/plane math (which needs sign-juggling over
/// which side of the plane the ball sits on), we read them directly
/// from the already-computed `center_result.center_curve`. The
/// approach mirrors [`try_plane_plane_torus`] above.
pub(crate) fn try_circular_center_torus(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let plane_normal = match (surface1, surface2) {
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Sphere { .. })
        | (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Plane { normal, .. }) => *normal,
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cylinder { axis, .. })
        | (SurfaceGeometry::Cylinder { axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            let n_unit = normal.normalized()?;
            let a_unit = axis.normalized()?;
            if n_unit.dot(&a_unit).abs() < 0.999 {
                return None;
            }
            *normal
        }
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cone { axis, .. })
        | (SurfaceGeometry::Cone { axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            let n_unit = normal.normalized()?;
            let a_unit = axis.normalized()?;
            // Oblique plane-cone intersections are conic sections
            // other than circles (ellipse / parabola / hyperbola).
            // Their fillet surfaces are canal surfaces with non-
            // circular spines — popkern has no canal primitive in
            // `SurfaceGeometry`, so we leave those cases to the
            // generic NURBS sweep. The deferral is documented in
            // `docs/decisions/2026-04-14-oblique-plane-cone-canal-
            // deferred.md` with the trigger criteria for revisiting.
            if n_unit.dot(&a_unit).abs() < 0.999 {
                return None;
            }
            *normal
        }
        _ => return None,
    };

    use popkern_core::geometry::curve::CurveEval;
    let n = 16usize;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let mut samples: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let t = t0 + (t1 - t0) * (i as f64) / ((n - 1).max(1) as f64);
        samples.push(center_result.center_curve.point_at(t));
    }
    if samples.len() < 3 {
        return None;
    }
    let cx = samples.iter().map(|p| p.x).sum::<f64>() / samples.len() as f64;
    let cy = samples.iter().map(|p| p.y).sum::<f64>() / samples.len() as f64;
    let cz = samples.iter().map(|p| p.z).sum::<f64>() / samples.len() as f64;
    let center = Point3::new(cx, cy, cz);

    let radii: Vec<f64> = samples.iter().map(|p| (*p - center).length()).collect();
    let major_radius = radii.iter().sum::<f64>() / radii.len() as f64;
    if !major_radius.is_finite() || major_radius <= 0.0 {
        return None;
    }

    // The center curve must be a CIRCLE — not just any 3D curve
    // that happens to average to a point. Enforce two conditions:
    //   1. planarity — max out-of-plane residual < 0.1 % of the
    //      major radius.
    //   2. uniform radius — max − min radial distance < 10 % of the
    //      major radius. The second check is what catches a linear
    //      center curve, which has min = 0 (centroid lies on the
    //      line) and max > 0.
    let plane_normal_n = plane_normal.normalized()?;
    let max_planar_residual = samples
        .iter()
        .map(|p| (*p - center).dot(&plane_normal_n).abs())
        .fold(0.0_f64, f64::max);
    if max_planar_residual > major_radius * 1e-3 {
        return None;
    }
    let radial_min = radii.iter().cloned().fold(f64::INFINITY, f64::min);
    let radial_max = radii.iter().cloned().fold(0.0_f64, f64::max);
    if (radial_max - radial_min) > major_radius * 0.1 {
        return None;
    }

    let surface = SurfaceGeometry::Torus {
        center,
        axis: plane_normal_n,
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

/// Try to build an analytic cylinder for a fillet between two cylinders
/// with **parallel axes** (§11.1).
///
/// **Geometry.** Two parallel infinite cylinders intersect in 0, 1, or
/// 2 straight lines, parallel to the common axis direction. Along one
/// such intersection line the rolling ball of radius `r` sits at a
/// fixed 2D position (x_c, y_c) in the cross-section plane, determined
/// by the circle-circle tangency problem `|center − axis_A| = R_A ± r`
/// and similarly for cylinder B. That position is constant in the
/// axial direction, so the center curve is a straight line parallel to
/// both cylinder axes and the fillet surface swept by the ball is an
/// exact cylinder:
///
/// - **Axis** = cylinder_axis_direction (shared by both supports)
/// - **Origin** = any point on the center curve
/// - **Radius** = the fillet radius
///
/// The detector reads the center-curve direction from sampling the
/// already-computed `center_result.center_curve` and bails if:
///
/// - the cylinders are not parallel (dot < 0.999);
/// - the center curve is not a straight line (max lateral residual
///   > 0.1 % of the curve span); or
/// - the center-curve direction is not parallel to the cylinder axes
///   (dot < 0.999) — which also catches cylinder-cylinder cases where
///   the axes happen to be parallel but the center curve isn't
///   (degenerate setups).
pub(crate) fn try_parallel_cylinder_cylinder(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let (axis1, axis2) = match (surface1, surface2) {
        (
            SurfaceGeometry::Cylinder { axis: a1, .. },
            SurfaceGeometry::Cylinder { axis: a2, .. },
        ) => (*a1, *a2),
        _ => return None,
    };
    let a1_n = axis1.normalized()?;
    let a2_n = axis2.normalized()?;
    if a1_n.dot(&a2_n).abs() < 0.999 {
        return None;
    }

    use popkern_core::geometry::curve::CurveEval;
    let n = 16usize;
    let (t0, t1) = center_result.center_curve.parameter_range();
    let samples: Vec<Point3> = (0..n)
        .map(|i| {
            let t = t0 + (t1 - t0) * (i as f64) / ((n - 1).max(1) as f64);
            center_result.center_curve.point_at(t)
        })
        .collect();
    if samples.len() < 3 {
        return None;
    }

    let p0 = samples[0];
    let p_last = samples[samples.len() - 1];
    let span = p_last - p0;
    let span_len = span.length();
    if span_len < 1e-12 {
        return None;
    }
    let dir = span.normalized()?;

    // The center curve direction must be parallel to the cylinder
    // axes. This is the load-bearing check — without it, a non-
    // straight center curve could still have a matching span
    // direction if the endpoints happen to line up.
    if dir.dot(&a1_n).abs() < 0.999 {
        return None;
    }

    // All intermediate samples must lie on the line through p0 in
    // direction `dir`. Max perpendicular residual < 0.1 % of the
    // span length. A straight line has residual ≈ 0 (floating-point
    // noise only); a curved center curve (e.g. from two cylinders
    // with non-parallel axes that accidentally passed the axis
    // check) will have larger residual and bail.
    let max_residual = samples
        .iter()
        .map(|p| {
            let v = *p - p0;
            let along = v.dot(&dir);
            let perp = v - dir * along;
            perp.length()
        })
        .fold(0.0_f64, f64::max);
    if max_residual > span_len * 1e-3 {
        return None;
    }

    let surface = SurfaceGeometry::Cylinder {
        origin: p0,
        axis: a1_n,
        radius,
    };

    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;

    Some(FilletSurfaceResult {
        surface,
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
    })
}
