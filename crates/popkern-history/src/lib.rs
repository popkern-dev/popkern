#![warn(missing_docs)]

//! popkern-history: parametric feature tree for the popkern CAD kernel.
//!
//! This crate stores modeling operations as a timeline of [`Feature`] nodes
//! and can replay them in order to rebuild a [`popkern_core::topology::BRep`]
//! from scratch.  The tree supports rollback to any point in the timeline,
//! suppression of individual features, and stable-name resolution across
//! replays.
//!
//! # Quick overview
//!
//! ```text
//! FeatureTree
//!   ├── features: SlotMap<FeatureId, FeatureNode>
//!   ├── timeline: Vec<FeatureId>          // ordered evaluation sequence
//!   └── rollback_index: Option<usize>     // Some(n) → replay only [0..n]
//! ```

pub mod error;
pub mod feature;
pub mod node;
pub mod timeline;

pub use error::{FeatureResult, HistoryError, ReplayError};
pub use feature::{BoolOpKind, Feature, ShellDirectionKind, SketchPlaneDesc, Vec3Desc};
pub use node::{BodyId, FeatureId, FeatureNode};
pub use timeline::{feature_entity_names, infer_dependencies};

use slotmap::SlotMap;

/// The parametric feature tree.
///
/// Stores all features in a [`SlotMap`] keyed by [`FeatureId`], and maintains
/// a separate `timeline` [`Vec`] that specifies evaluation order.  Methods for
/// inserting, suppressing, reordering, and replaying features will be added in
/// subsequent tasks.
#[derive(Debug, Clone)]
#[allow(dead_code)] // fields used in subsequent tasks (timeline, replay, cache)
pub struct FeatureTree {
    /// All feature nodes, keyed by their stable [`FeatureId`].
    pub(crate) features: SlotMap<FeatureId, FeatureNode>,
    /// Ordered sequence of feature handles defining the evaluation timeline.
    pub(crate) timeline: Vec<FeatureId>,
    /// When `Some(n)`, replay stops after the feature at index `n − 1` in the
    /// timeline (i.e. the tree is "rolled back" to that point).
    pub(crate) rollback_index: Option<usize>,
    /// Monotonically increasing counter for [`popkern_core::naming::OperationId`]
    /// assignment.  Starts at 1 (0 is reserved for "unassigned").
    pub(crate) next_op_id: u64,
}

impl FeatureTree {
    /// Create an empty feature tree.
    pub fn new() -> Self {
        Self {
            features: SlotMap::with_key(),
            timeline: Vec::new(),
            rollback_index: None,
            next_op_id: 1,
        }
    }
}

impl Default for FeatureTree {
    fn default() -> Self {
        Self::new()
    }
}
