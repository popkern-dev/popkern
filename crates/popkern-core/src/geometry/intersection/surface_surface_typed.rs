//! Surface-surface intersection with its contacts (PK-1126).
//!
//! [`intersect_surfaces`] returns the curves where two surfaces cross
//! and the regions where they coincide; where they only touch it
//! returns nothing, and its callers rely on that (the setback corner
//! dispatch, the Gregory self-intersection check). This entry point
//! returns the same result beside the contacts
//! ([`find_surface_contacts`]), so a caller that needs to know where
//! the surfaces touch asks for it, and no other caller changes.

use crate::error::CoreError;
use crate::geometry::intersection::contact_classify::ContactScale;
use crate::geometry::intersection::contact_discovery::{ContactSearch, find_surface_contacts};
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersectionResult, intersect_surfaces,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// The curves and coincident regions of a surface pair, and its
/// contacts.
#[derive(Debug, Clone)]
pub struct TypedSurfaceIntersection {
    /// What [`intersect_surfaces`] returns.
    pub curves: SurfaceSurfaceIntersectionResult,
    /// Where the surfaces touch, typed, and their near misses.
    pub contacts: ContactSearch,
}

/// [`intersect_surfaces`] with the contacts, judged at `tolerance`'s
/// linear tolerance over `extent` (the size of the faces meeting).
#[allow(clippy::too_many_arguments)]
pub fn intersect_surfaces_typed(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    extent: f64,
) -> Result<TypedSurfaceIntersection, CoreError> {
    let curves = intersect_surfaces(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    )?;
    let contacts = find_surface_contacts(
        surface1,
        (u1_range, v1_range),
        surface2,
        (u2_range, v2_range),
        ContactScale {
            linear: tolerance.linear(),
            extent,
        },
    );
    Ok(TypedSurfaceIntersection { curves, contacts })
}
