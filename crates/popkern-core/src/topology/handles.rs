//! Handle types for all topological entities.
//!
//! Handles are opaque IDs with generational indices (via slotmap).
//! They are Copy, Eq, Hash, and detect use-after-delete.

use slotmap::new_key_type;

new_key_type! {
    /// Handle to a vertex in a BRep.
    pub struct VertexHandle;
    /// Handle to a half-edge in a BRep.
    pub struct HalfEdgeHandle;
    /// Handle to an edge in a BRep.
    pub struct EdgeHandle;
    /// Handle to a loop in a BRep.
    pub struct LoopHandle;
    /// Handle to a face in a BRep.
    pub struct FaceHandle;
    /// Handle to a shell in a BRep.
    pub struct ShellHandle;
    /// Handle to a solid in a BRep.
    pub struct SolidHandle;
}
