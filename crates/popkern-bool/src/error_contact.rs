//! Non-manifold contacts and body-kind labels for boolean errors
//! (split from `error.rs` under the 500-line rule).

/// "solids", "solids or sheets": the kinds an operation accepts.
pub(super) fn accepted_label(kinds: &[popkern_core::topology::BodyKind]) -> String {
    popkern_core::topology::BodyKind::list_label(kinds)
}

/// The kind of a zero-thickness contact in a non-manifold result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactKind {
    /// More than two result faces would meet at an edge.
    Edge,
    /// The result faces at a vertex would form more than one fan.
    Vertex,
}

/// One zero-thickness contact of a non-manifold result.
#[derive(Debug, Clone, PartialEq)]
pub struct NonManifoldContact {
    /// Edge or vertex contact.
    pub kind: ContactKind,
    /// A point on the contact (the vertex, or the edge's midpoint).
    pub location: popkern_core::geometry::Point3,
}
