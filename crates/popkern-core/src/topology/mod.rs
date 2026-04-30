//! Topological entities for boundary representation.

/// Edge topology entity.
pub mod edge;
/// Face topology entity.
pub mod face;
/// Half-edge topology entity.
pub mod half_edge;
pub mod handles;
/// Loop topology entity.
pub mod loop_;
/// Shell topology entity.
pub mod shell;
/// Solid topology entity.
pub mod solid;
/// Vertex topology entity.
pub mod vertex;

pub use edge::{CanonicalForm, Edge};
pub use face::Face;
pub use half_edge::HalfEdge;
pub use handles::*;
pub use loop_::Loop;
pub use shell::Shell;
pub use solid::Solid;
pub use vertex::Vertex;

pub mod brep;
pub use brep::BRep;

/// Non-manifold radial-cycle helpers for [`BRep`].
mod brep_radial;

/// Pcurve recomputation for BRep faces.
mod pcurve_recompute;
mod pcurve_recompute_helpers;

/// Identity queries and recovery operations for precision tracking.
mod brep_precision_queries;

/// Geometric deviation measurement for edges and vertices.
mod brep_measurement;

/// Edge and face splitting operations.
pub mod split;
mod split_boundary;
pub(crate) mod split_helpers;
pub use split::{SplitEdgeResult, SplitFaceResult};

/// Inner loop (hole) creation for face trimming.
pub mod trim;
pub use trim::AddInnerLoopResult;

/// Vertex-welding primitives for [`BRep`].
mod brep_weld;
