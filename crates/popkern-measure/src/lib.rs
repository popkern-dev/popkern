#![warn(missing_docs)]

//! popkern-measure: Geometric measurement queries for popkern BRep solids.
//!
//! Pure computation functions that take typed topology handles. No knowledge
//! of entity names, feature trees, or expressions.

mod error;

pub use error::MeasureError;

mod face_helpers;
mod planar_area;

mod edge;

pub use edge::{edge_length, edge_tangent};

mod vertex;

pub use vertex::vertex_position;

mod face;

pub use face::{
    face_area, face_centroid, face_curvature, face_normal, face_perimeter, face_radius,
};

mod angle;

pub use angle::{dihedral_angle, edge_angle};

mod body;

pub use body::{
    Axis, BoundingBox, body_bounding_box, body_centroid, body_moment_of_inertia, body_surface_area,
    body_volume,
};

mod distance;

/// Per-vertex surface-curvature fields over a tessellation (PK-520).
mod surface_curvature;

pub use distance::{TopoHandle, min_distance};
pub use surface_curvature::{CurvatureField, CurvatureKind, curvature_at, curvature_field};
