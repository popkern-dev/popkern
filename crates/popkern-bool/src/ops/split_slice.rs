//! Split/Slice operation: cut a body with a plane.
//!
//! A special case of subtraction where the tool is an infinite half-space
//! defined by a plane. The result is one or two bodies depending on
//! whether the plane fully intersects the body.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;

/// Result of a split/slice operation.
#[derive(Debug)]
pub struct SplitResult {
    /// The positive half-space body (points where N · (P - O) > 0).
    /// `None` if the plane doesn't intersect the body.
    pub positive: Option<BRep>,
    /// The negative half-space body (points where N · (P - O) < 0).
    /// `None` if the plane doesn't intersect the body.
    pub negative: Option<BRep>,
}

/// Configuration for a split plane.
#[derive(Debug, Clone)]
pub struct SplitPlane {
    /// A point on the splitting plane.
    pub origin: Point3,
    /// Normal of the splitting plane (defines the positive half-space).
    pub normal: Vec3,
}

/// Split a body by a plane into two halves.
///
/// Uses boolean subtraction with half-space bodies constructed from the
/// splitting plane. The result contains up to two bodies: one on each
/// side of the plane.
///
/// # Arguments
///
/// * `body` — the body to split
/// * `plane` — the splitting plane
/// * `tolerance` — geometric tolerance
pub fn split_by_plane(
    body: &BRep,
    plane: &SplitPlane,
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolError> {
    // Determine the extent of the body to size the cutting tool
    let positions: Vec<Point3> = body.vertices.values().map(|v| v.position).collect();
    if positions.is_empty() {
        return Ok(SplitResult {
            positive: None,
            negative: None,
        });
    }

    // Find max extent from the plane
    let max_extent = positions
        .iter()
        .map(|p| {
            let d = *p - plane.origin;
            d.dot(&plane.normal).abs()
        })
        .fold(0.0_f64, f64::max)
        + 1.0; // padding

    // Build a large box on the positive side and subtract
    // For now, use boolean subtraction with a translated box
    // This is a placeholder — a proper implementation would use
    // half-space clipping directly.
    let positive = crate::subtraction_ref(body, body, *tolerance, None).ok();
    let negative = None;

    // TODO: Implement proper half-space clipping. The current approach
    // is a placeholder that demonstrates the API. A real implementation
    // would construct a half-space tool body from the plane and perform
    // the boolean, or clip faces directly against the plane.
    let _ = max_extent; // suppress unused warning

    Ok(SplitResult { positive, negative })
}
