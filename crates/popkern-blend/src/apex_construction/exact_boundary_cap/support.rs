//! Support-surface readiness classification for exact-boundary caps.

use popkern_core::geometry::surface::SurfaceGeometry;

/// Analytic or trusted support-surface family used by support ribbons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundarySupportSurfaceFamily {
    /// Infinite plane support.
    Plane,
    /// Periodic cylindrical support.
    Cylinder,
    /// Conical support.
    Cone,
    /// Periodic spherical support.
    Sphere,
    /// Trusted NURBS support.
    Nurbs,
    /// Currently unsupported support family.
    Unsupported,
}

/// Classify a support surface for exact-boundary ribbon preparation.
pub(crate) fn support_surface_family(
    surface: &SurfaceGeometry,
) -> ExactBoundarySupportSurfaceFamily {
    match surface {
        SurfaceGeometry::Plane { .. } => ExactBoundarySupportSurfaceFamily::Plane,
        SurfaceGeometry::Cylinder { .. } => ExactBoundarySupportSurfaceFamily::Cylinder,
        SurfaceGeometry::Cone { .. } => ExactBoundarySupportSurfaceFamily::Cone,
        SurfaceGeometry::Sphere { .. } => ExactBoundarySupportSurfaceFamily::Sphere,
        SurfaceGeometry::Nurbs(_) => ExactBoundarySupportSurfaceFamily::Nurbs,
        SurfaceGeometry::Torus { .. }
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => ExactBoundarySupportSurfaceFamily::Unsupported,
    }
}
