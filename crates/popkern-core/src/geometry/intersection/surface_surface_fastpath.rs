//! The pre-marching coincidence fast path — split from
//! `surface_surface.rs` (500-line cap rule).

use super::{SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult};
use crate::geometry::intersection::coincident::{
    closest_point_distance_pub, detect_coincident_regions,
};
use crate::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Attempt to detect coincident surfaces before the expensive marching solver.
///
/// **Tier A** — analytic O(1) comparison for same-type analytic pairs.
/// **Tier B** — sample-based: 4-corner quick rejection followed by the
/// existing 8×8 grid check for NURBS / mixed-type pairs.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_coincident_fast(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<SurfaceSurfaceIntersectionResult> {
    // Tier A: analytic comparison for same-type analytic surfaces.
    if are_surfaces_coincident_analytic(surface1, surface2, tolerance) {
        return Some(SurfaceSurfaceIntersectionResult {
            intersections: vec![SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range,
                v1_range,
                u2_range,
                v2_range,
            }],
        });
    }

    // Tier B: sample-based detection for NURBS / mixed types.
    // Quick rejection: sample the 4 corners of surface1's parameter domain and
    // project each onto surface2. If any corner is far away, bail out early.
    let corners = [
        (u1_range.0, v1_range.0),
        (u1_range.0, v1_range.1),
        (u1_range.1, v1_range.0),
        (u1_range.1, v1_range.1),
    ];
    let rejection_threshold = tolerance.linear() * 10.0;
    for (u, v) in &corners {
        let p = surface1.point_at(*u, *v);
        let dist = closest_point_distance_pub(p, surface2, u2_range, v2_range);
        if dist > rejection_threshold {
            return None;
        }
    }

    // All corners are close — run the full 8×8 grid check.
    let regions = detect_coincident_regions(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );
    if regions.is_empty() {
        None
    } else {
        Some(SurfaceSurfaceIntersectionResult {
            intersections: regions,
        })
    }
}
