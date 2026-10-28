//! Boolean options (PK-964): how a boolean treats a result that would
//! be non-manifold.
//!
//! Operands touching along an edge, at a point, or along a line make a
//! result that is non-manifold where they touch. No manifold solid
//! represents it, so the boolean either refuses — the default, and
//! what CAD systems do with a zero-thickness result — or, when asked,
//! serves the result as manifold pieces meeting at recorded contacts
//! (`BRep::topology_contacts`), the form a slicer or a mesher can take.

/// What a boolean does when its result would be non-manifold.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum NonManifoldPolicy {
    /// Decline with [`crate::BoolError::NonManifoldResult`], naming
    /// every contact.
    #[default]
    Decline,
    /// Serve manifold pieces that meet at the contacts: a contact edge
    /// or vertex becomes one copy per piece, and every contact is
    /// recorded on the result (`BRep::topology_contacts`), so
    /// `validate_all` accepts the deliberate touch and still rejects
    /// any other. A union whose operands share area at the touch is
    /// not split (they merge there); it declines if the merge fails.
    SplitToManifold,
}

/// Options for [`crate::boolean_with_options`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BoolOptions {
    /// Non-manifold result handling.
    pub non_manifold: NonManifoldPolicy,
}
