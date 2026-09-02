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
    #[serde(
        default,
        skip_serializing_if = "HashMap::is_empty",
        serialize_with = "serialize_op_ids",
        deserialize_with = "deserialize_op_ids"
    )]
    pub instance_op_ids: HashMap<(u32, u32), OperationId>,
    /// Cached outcome of the most recent evaluation, or `None` if the node
    /// has never been evaluated.
    ///
    /// Not serialized — rebuilt by [`crate::FeatureTree::replay`] after loading.
    #[serde(skip)]
    pub last_result: Option<FeatureResult>,
    /// How each selection argument was authored (scripting design D3):
    /// the query's canonical spelling and the labels it cited. Pure
    /// diagnostics — the feature stores resolved names and never
    /// re-evaluates the text.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub selections: Vec<SelectionRecord>,
}

/// The authored form of one selection argument of a feature.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SelectionRecord {
    /// Which argument (`edges`, `faces_to_remove`, `face`, …).
    pub role: String,
    /// Canonical spelling of the query at record time.
    pub describe: String,
    /// Entity labels the query cited.
    #[serde(default)]
    pub labels_used: Vec<String>,
    /// Seam edges the query matched and the record dropped.
    #[serde(default)]
    pub seams_skipped: usize,
}

/// Serializes `instance_op_ids` with tuple keys encoded as `"inst_seed"` strings.
///
/// JSON object keys must be strings; this converts `(u32, u32)` to `"<inst>_<seed>"`.
fn serialize_op_ids<S: serde::Serializer>(
    map: &HashMap<(u32, u32), OperationId>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    use serde::ser::SerializeMap;
    let mut s = serializer.serialize_map(Some(map.len()))?;
    // Sorted: the document bytes must not depend on hash order.
    let mut entries: Vec<_> = map.iter().collect();
    entries.sort_by_key(|entry| *entry.0);
    for (&(inst, seed), op_id) in entries {
        let key = format!("{inst}_{seed}");
        s.serialize_entry(&key, op_id)?;
    }
    s.end()
}

/// Deserializes `instance_op_ids` from `"inst_seed"` string keys back to `(u32, u32)`.
fn deserialize_op_ids<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<HashMap<(u32, u32), OperationId>, D::Error> {
    use serde::de::Error;
    let string_map: HashMap<String, OperationId> = HashMap::deserialize(deserializer)?;
    let mut result = HashMap::new();
    for (key, op_id) in string_map {
        let parts: Vec<&str> = key.split('_').collect();
        if parts.len() != 2 {
            return Err(D::Error::custom(format!("invalid op_id key: {key}")));
        }
        let inst: u32 = parts[0].parse().map_err(D::Error::custom)?;
        let seed: u32 = parts[1].parse().map_err(D::Error::custom)?;
        result.insert((inst, seed), op_id);
    }
    Ok(result)
}
