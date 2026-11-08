//! Boolean options (PK-964): how a boolean treats a result that would
//! be non-manifold.
//!
//! Operands touching along an edge, at a point, or along a line make a
//! result that is non-manifold where they touch. No manifold solid
//! represents it, so the boolean either refuses — the default, and
//! what CAD systems do with a zero-thickness result — or, when asked,
//! serves it as one general body holding the contacts in its topology,
//! or as separate manifold bodies that touch there (PK-1115: no single
//! body holds touching lumps, and nothing is recorded).

/// What a boolean does when its result would be non-manifold.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum NonManifoldPolicy {
    /// Decline with [`crate::BoolError::NonManifoldResult`], naming
    /// every contact.
    #[default]
    Decline,
    /// Serve separate manifold bodies that touch at the contacts, each
    /// with its own copy of a contact edge or vertex and valid on its
    /// own: the `_bodies` entry points ([`crate::boolean_bodies`],
    /// [`crate::nary_boolean_bodies`], ...) return them, the single-body
    /// ones decline with [`crate::BoolError::SplitIntoBodies`]. A result
    /// whose pieces touch nowhere is one body. A union whose operands
    /// share area at the touch is not split (they merge there); it
    /// declines if the merge fails.
    SplitToManifold,
    /// Serve one general body: the touching pieces share the contact
    /// vertex, the contact edge's ring holds every piece's faces in
    /// radial order, a contact curve across a face is imprinted on it
    /// (`BRep::glue_touching` over the split pieces, from their
    /// geometry). Operands may be general bodies (taken apart into
    /// their manifold pieces first, `BRep::separate_at_contacts`).
    /// Where the glue declines, [`crate::BoolError::GeneralBodyDeclined`]
    /// with a `GeneralDeclined` shadow line.
    General,
}

impl NonManifoldPolicy {
    /// Whether this is the default, [`NonManifoldPolicy::Decline`] (a
    /// recorded policy is left out of a saved document when it is).
    pub fn is_decline(&self) -> bool {
        *self == NonManifoldPolicy::Decline
    }

    /// Whether a non-manifold result serves (as pieces, or as the
    /// general body glued from them) rather than declining.
    pub fn serves(&self) -> bool {
        *self != NonManifoldPolicy::Decline
    }
}

/// Options for [`crate::boolean_with_options`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BoolOptions {
    /// Non-manifold result handling.
    pub non_manifold: NonManifoldPolicy,
}
