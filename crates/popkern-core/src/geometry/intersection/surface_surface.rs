//! Surface-surface intersection types and API.

use crate::error::CoreError;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::fitting::fit_curve;
use crate::geometry::intersection::coincident::detect_coincident_regions;
use crate::geometry::intersection::degeneracy::{DefaultDegeneracyHandler, DegeneracyHandler};
use crate::geometry::intersection::marching::{march, MarchSample};
use crate::geometry::intersection::starting_points::find_starting_points;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use serde::{Deserialize, Serialize};

/// A single intersection result between two surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceSurfaceIntersection {
    /// The two surfaces intersect along a curve.
    Curve(CurveOnSurface),
    /// The two surfaces are coincident over a rectangular parameter region.
    CoincidentRegion {
        /// Parameter range in u for surface 1.
        u1_range: (f64, f64),
        /// Parameter range in v for surface 1.
        v1_range: (f64, f64),
        /// Parameter range in u for surface 2.
        u2_range: (f64, f64),
        /// Parameter range in v for surface 2.
        v2_range: (f64, f64),
    },
}

/// Result of a surface-surface intersection operation.
#[derive(Debug, Clone)]
pub struct SurfaceSurfaceIntersectionResult {
    /// All intersections found between the two surfaces.
    pub intersections: Vec<SurfaceSurfaceIntersection>,
}

/// Intersect two surfaces with default degeneracy handling.
pub fn intersect_surfaces(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    intersect_surfaces_with_handler(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        &DefaultDegeneracyHandler,
    )
}

/// Intersect two surfaces with a custom degeneracy handler.
#[allow(clippy::too_many_arguments)]
pub fn intersect_surfaces_with_handler(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    let mut results = Vec::new();

    // 1. Find starting points
    let starting_points = find_starting_points(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );

    // 2. Handle tangent starting points -> coincident region detection
    let has_tangent = starting_points.iter().any(|sp| sp.is_tangent);
    if has_tangent {
        let coincident = detect_coincident_regions(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        );
        results.extend(coincident);
    }

    // 3. March from non-tangent, unconsumed starting points
    let mut consumed: Vec<Point3> = Vec::new();
    for sp in &starting_points {
        if sp.is_tangent {
            continue;
        }
        // Skip if this starting point is already consumed by a previous march
        if consumed
            .iter()
            .any(|c| c.distance_to(&sp.sample.point) < tolerance.linear * 10.0)
        {
            continue;
        }

        let (samples, is_closed) = march(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, &sp.sample,
            tolerance, handler,
        );

        if samples.len() < 2 {
            continue; // too few samples to fit
        }

        // Mark all starting points near this curve as consumed
        for s in &samples {
            consumed.push(s.point);
        }

        // 4. Fit NURBS curves
        match fit_march_samples(&samples, is_closed) {
            Ok(curve_on_surface) => {
                results.push(SurfaceSurfaceIntersection::Curve(curve_on_surface));
            }
            Err(_) => {
                // Fitting failed -- skip this curve (don't crash)
            }
        }
    }

    Ok(SurfaceSurfaceIntersectionResult {
        intersections: results,
    })
}

/// Fit NURBS curves (3D + UV1 + UV2) to a sequence of march samples.
fn fit_march_samples(
    samples: &[MarchSample],
    is_closed: bool,
) -> Result<CurveOnSurface, CoreError> {
    if samples.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "Too few samples to fit".into(),
        ));
    }

    // Chord-length parameterization
    let mut params = vec![0.0];
    let mut total_length = 0.0;
    for i in 1..samples.len() {
        total_length += samples[i].point.distance_to(&samples[i - 1].point);
        params.push(total_length);
    }
    if total_length < 1e-30 {
        return Err(CoreError::PreconditionFailed(
            "Degenerate curve (zero length)".into(),
        ));
    }
    // Normalize to [0, 1]
    for p in &mut params {
        *p /= total_length;
    }

    let degree = 3.min(samples.len() - 1);
    let n_cp = (samples.len() / 5).clamp(degree + 1, 200);

    // 3D curve
    let pts_3d: Vec<Point3> = samples.iter().map(|s| s.point).collect();
    let curve_3d = fit_curve(&pts_3d, &params, degree, n_cp)?;

    // UV1 curve (using Point3(u1, v1, 0.0) convention)
    let pts_uv1: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u1, s.v1, 0.0))
        .collect();
    let curve_uv1 = fit_curve(&pts_uv1, &params, degree, n_cp)?;

    // UV2 curve
    let pts_uv2: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u2, s.v2, 0.0))
        .collect();
    let curve_uv2 = fit_curve(&pts_uv2, &params, degree, n_cp)?;

    Ok(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Point3, Vec3};
    use std::f64::consts::{PI, TAU};

    /// A coarser tolerance so marching takes fewer, larger steps in tests.
    fn test_tolerance() -> Tolerance {
        Tolerance {
            linear: 1e-3,
            angular: 1e-3,
            parametric: 1e-6,
        }
    }

    #[test]
    fn two_perpendicular_planes_one_curve() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let result = intersect_surfaces(
            &s1,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &s2,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &test_tolerance(),
        )
        .unwrap();
        assert!(
            !result.intersections.is_empty(),
            "Should find intersection"
        );
        assert!(matches!(
            &result.intersections[0],
            SurfaceSurfaceIntersection::Curve(_)
        ));
    }

    #[test]
    fn coincident_planes_region() {
        let s = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = intersect_surfaces(
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &s,
            (0.0, 1.0),
            (0.0, 1.0),
            &test_tolerance(),
        )
        .unwrap();
        assert!(result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. })));
    }

    #[test]
    fn non_intersecting_empty() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 10.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = intersect_surfaces(
            &s1,
            (0.0, 1.0),
            (0.0, 1.0),
            &s2,
            (0.0, 1.0),
            (0.0, 1.0),
            &test_tolerance(),
        )
        .unwrap();
        assert!(result.intersections.is_empty());
    }

    #[test]
    fn plane_and_sphere() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let result = intersect_surfaces(
            &plane,
            (-2.0, 2.0),
            (-2.0, 2.0),
            &sphere,
            (0.0, TAU),
            (0.0, PI),
            &test_tolerance(),
        )
        .unwrap();
        // Should find an intersection curve (circle of radius 1 in XY plane)
        assert!(
            !result.intersections.is_empty(),
            "Should find plane-sphere intersection"
        );
    }
}
