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

use crate::center_curve_helpers::{
    clamp_finite_range, clip_footprint_to_bounds, project_onto_surface,
};
use crate::offset::offset_surface;

// Re-export commonly used items from helpers so external callers
// (e.g. blend::surface tests) keep working with crate::center_curve::.
pub use crate::center_curve_helpers::fit_nurbs_through_points;

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
    orientation_1: bool,
    orientation_2: bool,
) -> Result<CenterCurveResult, String> {
    // 1. Offset both surfaces inward by the fillet radius.
    // "Inward" means toward the solid interior.
    // In popkern, face normals point INTO the solid. When orientation=true,
    // surface normal = face normal (inward). offset_surface's inward flag
    // moves AGAINST the normal (sign=-1), so we need inward=false to move
    // WITH the inward-pointing normal (toward the interior).
    let inward_1 = !orientation_1;
    let inward_2 = !orientation_2;
    let offset1 = offset_surface(surface1, radius, inward_1, u1_range, v1_range, tolerance)?;
    let offset2 = offset_surface(surface2, radius, inward_2, u2_range, v2_range, tolerance)?;

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

    // 6. Clip footprint points to face UV bounds so that trim curves
    //    do not extend beyond the face domain.
    footprint_1 = clip_footprint_to_bounds(surface1, u1_range, v1_range, &footprint_1);
    footprint_2 = clip_footprint_to_bounds(surface2, u2_range, v2_range, &footprint_2);

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
    orientation_1: bool,
    orientation_2: bool,
) -> Result<CenterCurveResult, String> {
    parametric_radius_center_curve(
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        |t| radius_start + (radius_end - radius_start) * t,
        edge_midpoint,
        tolerance,
        orientation_1,
        orientation_2,
    )
}

/// Build a center curve for an arbitrary local-radius profile
/// `radius_at(t)` where `t` ranges over `[0, 1]` along the sharp
/// edge. Round 5l: ByPoints uses this with a piecewise-linear
/// closure so the full profile (not just the endpoint radii)
/// drives the rolling-ball center curve and footprints.
#[allow(clippy::too_many_arguments)]
pub fn parametric_radius_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius_at: impl Fn(f64) -> f64,
    edge_midpoint: &Point3,
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
) -> Result<CenterCurveResult, String> {
    let n_samples = 20usize;
    let mut center_points = Vec::with_capacity(n_samples);
    let mut footprint_1 = Vec::with_capacity(n_samples);
    let mut footprint_2 = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1).max(1) as f64;
        let local_radius = radius_at(t);

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
            orientation_1,
            orientation_2,
        )?;

        let fp1 = project_onto_surface(surface1, u1_range, v1_range, &center_pt)?;
        let fp2 = project_onto_surface(surface2, u2_range, v2_range, &center_pt)?;

        center_points.push(center_pt);
        footprint_1.push(fp1);
        footprint_2.push(fp2);
    }

    let center_curve = fit_nurbs_through_points(&center_points)?;

    let footprint_1 = clip_footprint_to_bounds(surface1, u1_range, v1_range, &footprint_1);
    let footprint_2 = clip_footprint_to_bounds(surface2, u2_range, v2_range, &footprint_2);

    Ok(CenterCurveResult {
        center_curve,
        footprint_1,
        footprint_2,
    })
}

/// Piecewise-linear interpolation of a `ByPoints` radius profile.
/// `pts` must be sorted by parameter ascending. Out-of-range `t`
/// is clamped to the endpoint radii.
pub fn bypoints_local_radius(pts: &[(f64, f64)], t: f64) -> f64 {
    if pts.is_empty() {
        return 0.0;
    }
    if t <= pts[0].0 {
        return pts[0].1;
    }
    if t >= pts[pts.len() - 1].0 {
        return pts[pts.len() - 1].1;
    }
    for window in pts.windows(2) {
        let (t0, r0) = window[0];
        let (t1, r1) = window[1];
        if t >= t0 && t <= t1 {
            let span = (t1 - t0).max(1e-12);
            let frac = (t - t0) / span;
            return r0 + (r1 - r0) * frac;
        }
    }
    pts[pts.len() - 1].1
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
    orientation_1: bool,
    orientation_2: bool,
) -> Result<Point3, String> {
    let offset1 = offset_surface(
        surface1,
        radius,
        !orientation_1,
        u1_range,
        v1_range,
        tolerance,
    )?;
    let offset2 = offset_surface(
        surface2,
        radius,
        !orientation_2,
        u2_range,
        v2_range,
        tolerance,
    )?;

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

// Helper functions (clip, clamp, fit, project) are in center_curve_helpers.rs.

#[cfg(test)]
#[path = "center_curve_tests.rs"]
mod tests;
