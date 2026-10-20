//! Split/Slice operation: cut a body with a plane.
//!
//! Constructs a large box on each side of the splitting plane and uses
//! boolean intersection to extract each half.

use popkern_core::geometry::{Point3, Vec3};
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
/// Each half is the body intersected with the plane's half-space on
/// that side ([`super::half_space`]: a prism whose cap lies exactly on
/// the plane, oriented to it, sized from the body's conservative face
/// boxes). A half is `None` when no part of the body lies on its side,
/// or when its boolean declines. The cutter used to be a unit cube
/// translated half a body-size off the plane, unscaled and axis-aligned:
/// for any body larger than it, the cube missed the plane entirely.
pub fn split_by_plane(
    body: &BRep,
    plane: &SplitPlane,
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolError> {
    let half = |normal: Vec3| -> Result<Option<BRep>, BoolError> {
        Ok(
            super::half_space::half_space(body, plane.origin, normal, *tolerance)?
                .and_then(|cutter| crate::intersection_ref(body, &cutter, *tolerance, None).ok())
                .filter(|b| !b.faces.is_empty()),
        )
    };
    Ok(SplitResult {
        positive: half(plane.normal * -1.0)?,
        negative: half(plane.normal)?,
    })
}
