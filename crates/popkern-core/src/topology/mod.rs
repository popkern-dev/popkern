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

pub use edge::Edge;
pub use face::Face;
pub use half_edge::HalfEdge;
pub use handles::*;
pub use loop_::Loop;
pub use shell::Shell;
pub use solid::Solid;
pub use vertex::Vertex;

pub mod brep;
pub use brep::BRep;

/// Pcurve recomputation for BRep faces.
mod pcurve_recompute;

/// Edge and face splitting operations.
pub mod split;
pub use split::{SplitEdgeResult, SplitFaceResult};

/// Inner loop (hole) creation for face trimming.
pub mod trim;
pub use trim::AddInnerLoopResult;
