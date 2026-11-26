//! A closed curve round a touch is the touch: the face-pair pass drops
//! it (`popkern_core::geometry::intersection::touch_ring`), so the pair
//! meets in no curve and the contact pass finds the touch
//! (`ff_parallel_contacts`).

use popkern_core::geometry::intersection::touch_ring::{is_touch_band, window_extent};
use popkern_core::geometry::intersection::{
    SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult,
};
use popkern_core::tolerance::Tolerance;

use super::PairInputs;

/// `result` without the closed curves that are a touch's tolerance band.
pub(super) fn drop_touch_rings(
    mut result: SurfaceSurfaceIntersectionResult,
    inp: &PairInputs<'_>,
    tolerance: &Tolerance,
) -> SurfaceSurfaceIntersectionResult {
    let extent = window_extent(inp.surf_a, inp.r1).max(window_extent(inp.surf_b, inp.r2));
    let scale = (tolerance.linear(), tolerance.linear_reach(), extent);
    result.intersections.retain(|i| match i {
        SurfaceSurfaceIntersection::Curve(c) => {
            !is_touch_band(c, (inp.surf_a, inp.r1), (inp.surf_b, inp.r2), scale)
        }
        _ => true,
    });
    result
}
