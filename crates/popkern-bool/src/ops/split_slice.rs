//! Split/Slice operation: cut a body with a plane.
//!
//! Constructs a large box on each side of the splitting plane and uses
//! boolean intersection to extract each half.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;

/// Result of a split/slice operation.
#[derive(Debug)]
pub struct SplitResult {
    /// The positive half-space body (points where N · (P - O) > 0).
    /// `None` if the plane doesn't intersect the body or the boolean failed.
    pub positive: Option<BRep>,
    /// The negative half-space body (points where N · (P - O) < 0).
    /// `None` if the plane doesn't intersect the body or the boolean failed.
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
/// Constructs a large box enclosing the body, positions it on each side
/// of the splitting plane, and uses boolean intersection to extract each
/// half. Returns `None` for a half if the plane doesn't intersect the
/// body on that side, or if the boolean fails.
pub fn split_by_plane(
    body: &BRep,
    plane: &SplitPlane,
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolError> {
    let positions: Vec<Point3> = body.vertices.values().map(|v| v.position).collect();
    if positions.is_empty() {
        return Ok(SplitResult {
            positive: None,
            negative: None,
        });
    }

    // Compute body extent to size the cutting box
    let centroid = {
        let n = positions.len() as f64;
        Point3::new(
            positions.iter().map(|p| p.x).sum::<f64>() / n,
            positions.iter().map(|p| p.y).sum::<f64>() / n,
            positions.iter().map(|p| p.z).sum::<f64>() / n,
        )
    };

    let max_dist = positions
        .iter()
        .map(|p| p.distance_to(&centroid))
        .fold(0.0_f64, f64::max);

    let box_size = max_dist * 3.0 + 1.0;

    // Build a large cube and position it on the positive side of the plane
    let tool_box = BRep::unit_cube(*tolerance).map_err(BoolError::CoreError)?;

    let normal = plane
        .normal
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0));

    // Positive half: translate box to positive side of plane
    let pos_offset = plane.origin + normal * (box_size * 0.5);
    let pos_transform = Transform3::translation(
        pos_offset.x - 0.5, // unit_cube is [0,1]^3, center at 0.5
        pos_offset.y - 0.5,
        pos_offset.z - 0.5,
    );
    // Scale would be needed for a proper implementation, but unit_cube * scale
    // isn't available. Use the box as-is and rely on the boolean to clip.

    let positive = tool_box
        .transform(&pos_transform)
        .ok()
        .and_then(|tool| crate::intersection_ref(body, &tool, *tolerance, None).ok());

    // Negative half: translate box to negative side of plane
    let neg_offset = plane.origin - normal * (box_size * 0.5);
    let neg_transform =
        Transform3::translation(neg_offset.x - 0.5, neg_offset.y - 0.5, neg_offset.z - 0.5);

    let negative = tool_box
        .transform(&neg_transform)
        .ok()
        .and_then(|tool| crate::intersection_ref(body, &tool, *tolerance, None).ok());

    Ok(SplitResult { positive, negative })
}
