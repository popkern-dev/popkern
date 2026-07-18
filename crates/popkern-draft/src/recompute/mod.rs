//! Edge and vertex recomputation after draft surface changes.

pub mod curve_snap;
pub mod edge;
pub mod vertex;

pub use curve_snap::snap_vertices_to_recomputed_curves;
pub use edge::recompute_draft_edges;
pub use vertex::recompute_draft_vertices;
