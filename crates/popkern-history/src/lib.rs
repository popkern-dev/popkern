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
pub(crate) mod replay;
pub(crate) mod resolve;
pub mod timeline;

pub use error::{FeatureResult, HistoryError, ReplayError};
pub use feature::{BoolOpKind, Feature, ShellDirectionKind, SketchPlaneDesc, Vec3Desc};
pub use node::{BodyId, FeatureId, FeatureNode};

use std::collections::HashMap;
use std::time::Instant;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use slotmap::SlotMap;

/// The result of replaying a feature tree.
///
/// Contains per-feature outcomes and the total elapsed time.
#[derive(Debug)]
pub struct ReplayResult {
    /// Per-feature outcomes, keyed by [`FeatureId`].
    pub feature_results: HashMap<FeatureId, FeatureResult>,
    /// Wall-clock time spent replaying.
    pub elapsed: std::time::Duration,
}

/// The parametric feature tree.
///
/// Stores all features in a [`SlotMap`] keyed by [`FeatureId`], and maintains
/// a separate `timeline` [`Vec`] that specifies evaluation order.  The tree
/// supports replay, rollback, and suppression.
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
    /// The BRep produced by the most recent replay, or `None` if the tree
    /// has not been replayed yet.
    pub(crate) current_brep: Option<BRep>,
}

impl FeatureTree {
    /// Create an empty feature tree.
    pub fn new() -> Self {
        Self {
            features: SlotMap::with_key(),
            timeline: Vec::new(),
            rollback_index: None,
            next_op_id: 1,
            current_brep: None,
        }
    }

    /// Replay the entire feature tree from scratch.
    ///
    /// Evaluates features in timeline order up to the rollback index (or the
    /// end of the timeline). Suppressed features are skipped. Failed features
    /// record an error but do not abort the replay — the BRep continues from
    /// the last good state.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] only for fatal infrastructure failures.
    /// Individual feature failures are recorded in
    /// [`ReplayResult::feature_results`].
    pub fn replay(&mut self, tolerance: Tolerance) -> Result<ReplayResult, ReplayError> {
        let start = Instant::now();
        let mut brep = BRep::new(tolerance);
        let mut results = HashMap::new();

        let end = self.rollback_index.unwrap_or(self.timeline.len());

        for i in 0..end {
            let id = self.timeline[i];
            let node = &self.features[id];

            if node.suppressed {
                results.insert(id, FeatureResult::Suppressed);
                continue;
            }

            let mut ctx = TrackingContext::new(node.operation_id);

            // Clone the feature so we don't hold a borrow on self.features
            // while mutating brep.
            let feature = node.feature.clone();

            match replay::execute_feature(&mut brep, &feature, &mut ctx, &tolerance) {
                Ok(()) => {
                    let history = ctx.finish();
                    results.insert(
                        id,
                        FeatureResult::Success {
                            history: history.clone(),
                        },
                    );
                    self.features[id].last_result = Some(FeatureResult::Success { history });
                }
                Err(e) => {
                    let error_msg = e.to_string();
                    results.insert(
                        id,
                        FeatureResult::Failed {
                            error: error_msg.clone(),
                        },
                    );
                    self.features[id].last_result =
                        Some(FeatureResult::Failed { error: error_msg });
                    // Continue with the BRep as-is (last good state).
                }
            }
        }

        self.current_brep = Some(brep);

        Ok(ReplayResult {
            feature_results: results,
            elapsed: start.elapsed(),
        })
    }

    /// Replay starting from a specific feature.
    ///
    /// Without caching (Task 4 adds that), this is equivalent to a full
    /// replay. Features before `id` are still re-executed to build up the
    /// BRep state that `id` depends on.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] if `id` is not in the timeline.
    pub fn replay_from(
        &mut self,
        id: FeatureId,
        tolerance: Tolerance,
    ) -> Result<ReplayResult, ReplayError> {
        // Validate that id exists in the timeline.
        let _pos =
            self.timeline
                .iter()
                .position(|&fid| fid == id)
                .ok_or(ReplayError::FeatureFailed {
                    label: "replay_from".into(),
                    error: "feature not in timeline".into(),
                })?;

        // Without caching we must replay from scratch.
        self.replay(tolerance)
    }

    /// Return a reference to the BRep produced by the most recent replay.
    ///
    /// Returns `None` if the tree has not been replayed yet.
    pub fn current_brep(&self) -> Option<&BRep> {
        self.current_brep.as_ref()
    }
}

impl Default for FeatureTree {
    fn default() -> Self {
        Self::new()
    }
}
