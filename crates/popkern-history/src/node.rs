//! Core node types for the feature tree.
//!
//! A [`FeatureNode`] is a single entry in the timeline — it pairs a
//! [`Feature`] (the operation's parameters) with all the bookkeeping needed
//! by the tree: dependencies, suppression state, cached result, etc.

use crate::error::FeatureResult;
use crate::expression::Expression;
use crate::feature::Feature;
use popkern_core::naming::OperationId;
use serde::{Deserialize, Serialize};
use slotmap::new_key_type;
use std::collections::HashMap;

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
///
/// The node's [`FeatureId`] is the SlotMap key — callers already hold it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureNode {
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
    /// Auto-suppress condition. When this expression evaluates to non-zero,
    /// the feature is skipped during replay.
    #[serde(default)]
    pub suppress_if: Option<Expression>,
    /// Persisted OperationIds for pattern instances.
    ///
    /// Key: `(instance_index, seed_feature_index_within_chain)`.
    /// Only populated for pattern features (LinearPattern, CircularPattern, etc.).
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub instance_op_ids: HashMap<(u32, u32), OperationId>,
    /// Cached outcome of the most recent evaluation, or `None` if the node
    /// has never been evaluated.
    ///
    /// Not serialized — rebuilt by [`crate::FeatureTree::replay`] after loading.
    #[serde(skip)]
    pub last_result: Option<FeatureResult>,
}
