//! Result types for surface-surface intersection (extracted from
//! `surface_surface.rs` under the file-size ratchet).

use crate::geometry::curve_on_surface::CurveOnSurface;
use serde::{Deserialize, Serialize};

/// A single intersection result between two surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceSurfaceIntersection {
    /// The two surfaces intersect along a curve.
    Curve(CurveOnSurface),
    /// The two surfaces are coincident over a rectangular parameter region.
    CoincidentRegion {
        /// Parameter range in u for surface 1.
        u1_range: (f64, f64),
        /// Parameter range in v for surface 1.
        v1_range: (f64, f64),
        /// Parameter range in u for surface 2.
        u2_range: (f64, f64),
        /// Parameter range in v for surface 2.
        v2_range: (f64, f64),
    },
}

/// Result of a surface-surface intersection operation.
#[derive(Debug, Clone)]
pub struct SurfaceSurfaceIntersectionResult {
    /// All intersections found between the two surfaces.
    pub intersections: Vec<SurfaceSurfaceIntersection>,
}
