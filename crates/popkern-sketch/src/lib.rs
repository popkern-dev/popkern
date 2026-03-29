#![warn(missing_docs)]

//! popkern-sketch: Sketch-to-BRep bridge for the popkern CAD kernel.
//!
//! Converts solved 2D parametric sketches (from `popkern-constraint`) on 3D
//! planes into BRep faces (in `popkern-core`). The pipeline is:
//!
//! 1. Define a [`SketchPlane`] in 3D space
//! 2. Build and solve a 2D sketch using `popkern-constraint`
//! 3. Convert the sketch into a BRep face via edge graph construction,
//!    loop tracing, and topology building

pub(crate) mod edge_graph;
pub mod error;
pub(crate) mod loop_tracer;
pub mod plane;

pub use error::{SketchError, SketchWarning};
pub use plane::SketchPlane;

// Re-export frequently used types from dependencies for convenience.
pub use popkern_constraint::{Constraint, Sketch};
