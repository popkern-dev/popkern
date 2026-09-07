//! The `ssi` trace topic — split from `surface_surface.rs`
//! (500-line cap rule): one line per intersect call (surface types,
//! windows, result counts).

use super::{SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult};
use crate::error::CoreError;
use crate::geometry::surface::SurfaceGeometry;

#[allow(clippy::too_many_arguments)]
pub(super) fn trace_ssi(
    result: &Result<SurfaceSurfaceIntersectionResult, CoreError>,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) {
    if !crate::trace::trace_enabled("ssi") {
        return;
    }

    let kind = |s: &SurfaceGeometry| -> &'static str {
        match s {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cylinder",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            _ => "Other",
        }
    };
    match &result {
        Ok(r) => {
            let curves = r
                .intersections
                .iter()
                .filter(|i| matches!(i, SurfaceSurfaceIntersection::Curve(_)))
                .count();
            let regions = r.intersections.len() - curves;
            eprintln!(
                "  [ssi] {}x{} u1[{:.3},{:.3}] v1[{:.3},{:.3}] u2[{:.3},{:.3}] v2[{:.3},{:.3}] -> {curves} curve(s), {regions} region(s)",
                kind(surface1),
                kind(surface2),
                u1_range.0,
                u1_range.1,
                v1_range.0,
                v1_range.1,
                u2_range.0,
                u2_range.1,
                v2_range.0,
                v2_range.1,
            );
        }
        Err(e) => eprintln!("  [ssi] {}x{} -> ERROR {e}", kind(surface1), kind(surface2)),
    }
}
