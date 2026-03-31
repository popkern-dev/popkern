#![warn(missing_docs)]

//! popkern-measure: Geometric measurement queries for popkern BRep solids.
//!
//! Pure computation functions that take typed topology handles. No knowledge
//! of entity names, feature trees, or expressions.

mod error;

pub use error::MeasureError;

mod edge;

pub use edge::{edge_length, edge_tangent};
