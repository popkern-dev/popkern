#![warn(missing_docs)]

//! popkern-measure: Geometric measurement queries for popkern BRep solids.
//!
//! Pure computation functions that take typed topology handles. No knowledge
//! of entity names, feature trees, or expressions.

mod error;

pub use error::MeasureError;

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
