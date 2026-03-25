//! Edge and vertex recomputation after draft surface changes.

pub mod edge;
pub mod vertex;

pub use edge::recompute_draft_edges;
pub use vertex::recompute_draft_vertices;
