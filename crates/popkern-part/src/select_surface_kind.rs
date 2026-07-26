//! [`SurfaceKind`] — the surface family a selector matches on.
//!
//! Split out of `select.rs`, which sits on the structural line-count
//! ratchet: the remedy for a file at its frozen size is a submodule,
//! not a smaller edit.

use popkern_core::geometry::surface::SurfaceGeometry;

/// The analytic family of a face's surface, for adjacency selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    /// Analytic plane.
    Plane,
    /// Analytic cylinder.
    Cylinder,
    /// Analytic sphere.
    Sphere,
    /// Analytic cone.
    Cone,
    /// Analytic torus.
    Torus,
    /// Freeform NURBS surface.
    Nurbs,
    /// Gregory blend patch.
    Gregory,
}

impl SurfaceKind {
    pub(crate) fn of(geometry: &SurfaceGeometry) -> SurfaceKind {
        match geometry {
            SurfaceGeometry::Plane { .. } => SurfaceKind::Plane,
            SurfaceGeometry::Cylinder { .. } => SurfaceKind::Cylinder,
            SurfaceGeometry::Sphere { .. } => SurfaceKind::Sphere,
            SurfaceGeometry::Cone { .. } => SurfaceKind::Cone,
            SurfaceGeometry::Torus { .. } => SurfaceKind::Torus,
            SurfaceGeometry::Nurbs(_) => SurfaceKind::Nurbs,
            // Both free-form cap families report as `Gregory` here.
            // This selector answers "what kind of surface is this for
            // picking purposes", and a user selecting cap faces wants
            // both; the construction that produced a patch is carried
            // by the geometry itself, not by this enum.
            SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
                SurfaceKind::Gregory
            }
        }
    }
}
