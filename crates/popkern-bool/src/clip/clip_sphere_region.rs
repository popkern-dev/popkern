//! The zero-crossing vote on a sphere face holding its chart's poles.
//!
//! The clip's periodic containment counts boundary crossings along the
//! axial line above a sample, the top pole counted only where the
//! face's contour runs to it. A ball face whose one loop is the small
//! circle where a neighbouring ball meets it holds both poles and every
//! point but that cap's: no crossing lies above a sample on the far side,
//! and the count read it outside. The circle where a third ball of a
//! chain meets that face then dropped whole, the operands were left
//! with no intersection, and the union declined undecided. The core
//! trim classifier reads a contractible loop on the closed chart by
//! the side it bounds, so on a sphere its verdict decides the sample.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// Whether the chart sample `q` lies on the sphere face `face`, as the
/// core trim classifier reads it; `None` for a face on another surface,
/// a sample within its boundary band, or no answer (the clip's own probe
/// then decides).
pub(super) fn sphere_sample_inside(brep: &BRep, face: FaceHandle, q: Point3) -> Option<bool> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    if !matches!(surface, SurfaceGeometry::Sphere { .. }) {
        return None;
    }
    match brep.face_uv_trim_containment(face, surface, q.x, q.y)? {
        FaceUvContainment::Interior => Some(true),
        FaceUvContainment::Outside => Some(false),
        FaceUvContainment::Boundary => None,
    }
}

/// The vote's verdict for sample `q`: [`sphere_sample_inside`] where it
/// answers, else the clip's own region probe.
pub(super) fn sample_inside(
    brep: &BRep,
    face: FaceHandle,
    q: Point3,
    probe: &crate::face_region_probe::FaceRegionProbe<'_>,
) -> bool {
    sphere_sample_inside(brep, face, q).unwrap_or_else(|| probe.contains(q))
}
