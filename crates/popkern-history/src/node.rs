//! Core node types for the feature tree.
//!
//! A [`FeatureNode`] is a single entry in the timeline — it pairs a
//! [`Feature`] (the operation's parameters) with all the bookkeeping needed
//! by the tree: dependencies, suppression state, cached result, etc.

use crate::error::FeatureResult;
use crate::feature::Feature;
use popkern_core::naming::OperationId;
use serde::{Deserialize, Serialize};
use slotmap::new_key_type;

new_key_type! {
    /// Stable handle for a feature node in the [`crate::FeatureTree`].
    pub struct FeatureId;
}

/// Identifies a body within the feature tree.
///
/// Multi-body workflows (e.g. Boolean tool bodies) assign each discrete solid
/// its own `BodyId`.  The main solid starts at `BodyId(0)`.
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub struct BodyId(pub u32);

/// A single node in the parametric feature timeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureNode {
    /// The stable handle for this node.
    pub id: FeatureId,
    /// The operation and its parameters.
    pub feature: Feature,
    /// The [`OperationId`] assigned when this feature last evaluated
    /// successfully, used to resolve stable names.
    pub operation_id: OperationId,
    /// Human-readable name shown in the feature manager UI.
    pub label: String,
    /// Features this node depends on (must be evaluated before this one).
    pub depends_on: Vec<FeatureId>,
    /// When `true` the feature is skipped during replay.
    pub suppressed: bool,
    /// Which body this feature belongs to / operates on.
    pub body: BodyId,
    /// Cached outcome of the most recent evaluation, or `None` if the node
    /// has never been evaluated.
    pub last_result: Option<FeatureResult>,
}
