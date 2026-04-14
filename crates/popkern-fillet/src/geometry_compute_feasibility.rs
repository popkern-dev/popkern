//! Feasibility pre-flight helpers extracted from
//! `geometry_compute.rs` to keep the parent under the 500-line
//! structural cap.
//!
//! Contains the self-intersection curvature check that the fillet
//! engine runs before attempting any real geometry work — if the
//! requested radius would make the offset of either support surface
//! self-intersect, we short-circuit with a structured
//! `BlendFailure::RadiusTooLarge` carrying the tighter bound so the
//! caller can reduce the radius and retry.

use popkern_core::geometry::surface::SurfaceGeometry;

use crate::offset::self_intersect::check_offset_self_intersection;
use crate::types::FilletRadius;

/// Largest requested radius in a radius spec. Used as the "requested"
/// value in `BlendFailure::RadiusTooLarge`, so the user sees the
/// worst-case instance of their spec, not an arbitrary midpoint.
pub(crate) fn largest_requested_radius(radius: &FilletRadius) -> f64 {
    match radius {
        FilletRadius::Constant(r) => *r,
        FilletRadius::Variable { start, end } => start.max(*end),
        FilletRadius::ByPoints(pts) => pts
            .iter()
            .map(|(_, r)| *r)
            .fold(0.0_f64, |acc, r| acc.max(r)),
    }
}

/// Pre-flight self-intersection check for a rolling-ball fillet at the
/// requested radius. Returns `Some(max_safe_radius)` when the
/// requested radius would exceed the feasibility envelope on either
/// support surface, or `None` when the radius is safe (i.e. the
/// caller should proceed with the normal geometry pipeline).
///
/// For [`FilletRadius::Constant`] we test the single value. For
/// variable/by-points specs we test the peak requested radius —
/// if the peak is unsafe, the spec is unconditionally unsafe.
#[allow(clippy::too_many_arguments)]
pub(crate) fn self_intersection_safe_radius(
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    radius: &FilletRadius,
    orientation_1: bool,
    orientation_2: bool,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> Option<f64> {
    let peak = largest_requested_radius(radius);
    if peak <= 0.0 {
        return None;
    }

    // "inward" matches the logic in center_curve.rs — offset_surface's
    // inward flag moves AGAINST the normal, and face normals point
    // into the solid, so we need `inward = !orientation`.
    let inward_1 = !orientation_1;
    let inward_2 = !orientation_2;

    let max_1 = check_offset_self_intersection(surface_1, peak, inward_1, u1_range, v1_range);
    let max_2 = check_offset_self_intersection(surface_2, peak, inward_2, u2_range, v2_range);

    // Both surfaces reported "safe" — nothing to flag.
    match (max_1, max_2) {
        (None, None) => None,
        // One or both reported a safety bound — take the tighter
        // (smaller) bound since both surfaces have to accommodate
        // the rolling ball simultaneously.
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
    }
}
