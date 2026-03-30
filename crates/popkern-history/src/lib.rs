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
//!   ├── rollback_index: Option<usize>     // Some(n) → replay only [0..n]
//!   └── cache: BRepCache                  // intermediate BRep snapshots
//! ```

pub(crate) mod body;
pub mod cache;
pub mod delete;
pub mod error;
pub mod feature;
pub mod node;
pub mod queries;
pub(crate) mod replay;
pub(crate) mod resolve;
pub mod serialize;
pub mod timeline;
pub(crate) mod undo;

pub use cache::CacheStrategy;
pub use delete::{DeleteBehavior, DeleteResult};
pub use error::{FeatureResult, HistoryError, ReplayError};
pub use feature::{BoolOpKind, Feature, ShellDirectionKind, SketchPlaneDesc, Vec3Desc};
pub use node::{BodyId, FeatureId, FeatureNode};
pub use queries::ImpactReport;

use std::collections::HashMap;
use std::time::Instant;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;

use crate::body::BodyManager;
use crate::cache::BRepCache;
use crate::undo::TreeCommand;

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
/// supports replay, rollback, suppression, and intermediate-state caching.
///
/// # Serialization
///
/// Only the structural state (features, timeline, rollback index, next op id)
/// is persisted.  Session-local fields (cached BReps, undo log) are skipped —
/// call [`replay`](Self::replay) after loading to rebuild runtime state.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    ///
    /// Not serialized — call [`replay`](Self::replay) after loading.
    #[serde(skip)]
    pub(crate) current_brep: Option<BRep>,
    /// Cache of intermediate BRep snapshots.
    ///
    /// Not serialized — rebuilt on demand after loading.
    #[serde(skip)]
    pub(crate) cache: BRepCache,
    /// Tracks the mapping from [`BodyId`]s to features and BRep snapshots.
    ///
    /// Not serialized — rebuilt from features on replay after loading.
    #[serde(skip)]
    pub(crate) body_manager: BodyManager,
    /// Log of reversible tree mutations, oldest first.
    ///
    /// Not serialized — session-local undo history is not persisted.
    #[serde(skip)]
    pub(crate) undo_log: Vec<TreeCommand>,
    /// Index into `undo_log` pointing to the next command that would be
    /// redone.  Commands at `[0..undo_cursor]` have been applied; commands
    /// at `[undo_cursor..]` are "future" (available for redo).
    ///
    /// Not serialized — reset to 0 on load.
    #[serde(skip)]
    pub(crate) undo_cursor: usize,
}

impl FeatureTree {
    /// Create an empty feature tree with the default [`CacheStrategy::EveryStep`] cache.
    pub fn new() -> Self {
        Self {
            features: SlotMap::with_key(),
            timeline: Vec::new(),
            rollback_index: None,
            next_op_id: 1,
            current_brep: None,
            cache: BRepCache::new(CacheStrategy::EveryStep),
            body_manager: BodyManager::new(),
            undo_log: Vec::new(),
            undo_cursor: 0,
        }
    }

    // -----------------------------------------------------------------------
    // Replay
    // -----------------------------------------------------------------------

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
        let end = self.rollback_index.unwrap_or(self.timeline.len());

        let (brep, results) = self.replay_range(BRep::new(tolerance), 0, end, &tolerance)?;

        self.current_brep = Some(brep);

        Ok(ReplayResult {
            feature_results: results,
            elapsed: start.elapsed(),
        })
    }

    /// Replay starting from a specific feature.
    ///
    /// If a cached BRep is available for any ancestor of `id`, replay resumes
    /// from the nearest such ancestor rather than from scratch.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] if `id` is not in the timeline.
    pub fn replay_from(
        &mut self,
        id: FeatureId,
        tolerance: Tolerance,
    ) -> Result<ReplayResult, ReplayError> {
        let pos = self
            .timeline
            .iter()
            .position(|&fid| fid == id)
            .ok_or(HistoryError::FeatureNotFound)?;

        // Find the nearest cached BRep at or before `pos`.
        let mut start_pos = 0;
        let mut start_brep = BRep::new(tolerance);
        for i in (0..pos).rev() {
            let fid = self.timeline[i];
            if let Some(cached) = self.cache.get(fid) {
                start_pos = i + 1;
                start_brep = cached.clone();
                break;
            }
        }

        let timer = Instant::now();
        let end = self.rollback_index.unwrap_or(self.timeline.len());

        let (brep, results) = self.replay_range(start_brep, start_pos, end, &tolerance)?;

        self.current_brep = Some(brep);

        Ok(ReplayResult {
            feature_results: results,
            elapsed: timer.elapsed(),
        })
    }

    // -----------------------------------------------------------------------
    // Rollback / scrubbing
    // -----------------------------------------------------------------------

    /// Move the rollback bar to `index`.
    ///
    /// After this call the timeline is logically "rolled back" to just before
    /// the feature at `index`.  The next [`replay`](Self::replay) will stop at
    /// that point.
    ///
    /// Returns the cached BRep for the feature immediately *before* the bar
    /// (i.e. the BRep state at `index − 1`) if one exists in the cache.
    pub fn rollback_to(&mut self, index: usize) -> Option<&BRep> {
        let old = self.rollback_index;
        self.rollback_index = Some(index);
        self.record_undo(TreeCommand::MoveRollback {
            old,
            new: Some(index),
        });
        if index > 0 {
            let id = self.timeline[index - 1];
            self.cache.get(id)
        } else {
            None
        }
    }

    /// Remove the rollback bar, making all features in the timeline active.
    pub fn rollback_to_end(&mut self) {
        self.rollback_index = None;
    }

    /// Return the current rollback position, or `None` if no bar is set.
    pub fn rollback_position(&self) -> Option<usize> {
        self.rollback_index
    }

    /// Return the BRep state after the feature at `index` in the timeline.
    ///
    /// If the snapshot is cached it is returned instantly; otherwise a
    /// partial replay is executed from the nearest cached ancestor and the
    /// new snapshot is stored.
    ///
    /// `index` is one-based (i.e. `brep_at(1)` returns the BRep after the
    /// first feature).  Pass `0` to get an empty BRep.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] if any feature in the required range fails
    /// fatally, or if `index` exceeds the timeline length.
    pub fn brep_at(&mut self, index: usize, tolerance: Tolerance) -> Result<&BRep, ReplayError> {
        if index > self.timeline.len() {
            return Err(ReplayError::HistoryError(HistoryError::InvalidPosition(
                index,
            )));
        }

        // Fast path: already cached.
        if index > 0 {
            let id = self.timeline[index - 1];
            if self.cache.get(id).is_some() {
                // SAFETY: we just checked `is_some()`. The two-step pattern
                // (check then re-borrow) is required to avoid holding an
                // immutable borrow on `self.cache` across the mutable
                // `self.replay()` call below.
                return self
                    .cache
                    .get(id)
                    .ok_or(HistoryError::FeatureNotFound.into());
            }
        }

        // Slow path: replay up to `index` (preserving the caller's rollback
        // setting afterwards).
        let saved_rollback = self.rollback_index;
        self.rollback_index = Some(index);
        self.replay(tolerance)?;
        self.rollback_index = saved_rollback;

        // replay() sets current_brep; return it.
        self.current_brep
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "brep_at".into(),
                error: "replay did not produce a BRep".into(),
            })
    }

    // -----------------------------------------------------------------------
    // Cache management
    // -----------------------------------------------------------------------

    /// Replace the caching strategy.  All existing cached snapshots are evicted.
    pub fn set_cache_strategy(&mut self, strategy: CacheStrategy) {
        self.cache.set_strategy(strategy);
    }

    /// Evict all cached BRep snapshots.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    // -----------------------------------------------------------------------
    // Undo / redo
    // -----------------------------------------------------------------------

    /// Undo the last tree mutation. Returns `true` if something was undone.
    pub fn undo(&mut self) -> bool {
        if self.undo_cursor == 0 {
            return false;
        }
        self.undo_cursor -= 1;
        let cmd = self.undo_log[self.undo_cursor].clone();
        undo::reverse_command(self, &cmd);
        true
    }

    /// Redo the last undone mutation. Returns `true` if something was redone.
    pub fn redo(&mut self) -> bool {
        if self.undo_cursor >= self.undo_log.len() {
            return false;
        }
        let cmd = self.undo_log[self.undo_cursor].clone();
        undo::apply_command(self, &cmd);
        self.undo_cursor += 1;
        true
    }

    /// Whether undo is available.
    pub fn can_undo(&self) -> bool {
        self.undo_cursor > 0
    }

    /// Whether redo is available.
    pub fn can_redo(&self) -> bool {
        self.undo_cursor < self.undo_log.len()
    }

    /// Record a command in the undo log. Truncates any redo history.
    pub(crate) fn record_undo(&mut self, cmd: TreeCommand) {
        self.undo_log.truncate(self.undo_cursor);
        self.undo_log.push(cmd);
        self.undo_cursor = self.undo_log.len();
    }

    // -----------------------------------------------------------------------
    // BRep accessor
    // -----------------------------------------------------------------------

    /// Return a reference to the BRep produced by the most recent replay.
    ///
    /// Returns `None` if the tree has not been replayed yet.
    pub fn current_brep(&self) -> Option<&BRep> {
        self.current_brep.as_ref()
    }

    // -----------------------------------------------------------------------
    // Body accessors
    // -----------------------------------------------------------------------

    /// Return all body IDs currently tracked by this tree.
    pub fn bodies(&self) -> Vec<BodyId> {
        self.body_manager.bodies()
    }

    /// Return the current BRep snapshot for a body.
    ///
    /// Returns `None` if no snapshot has been stored yet (i.e. the tree has
    /// not been replayed, or the body does not exist).
    pub fn body_brep(&self, body: BodyId) -> Option<&BRep> {
        self.body_manager.get_brep(body)
    }

    /// Return the features assigned to a body, in insertion order.
    ///
    /// Returns an empty [`Vec`] if the body does not exist.
    pub fn body_features(&self, body: BodyId) -> Vec<FeatureId> {
        self.body_manager.body_features(body)
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Execute features `from..to` against `start_brep`, caching results and
    /// returning the final BRep together with per-feature outcomes.
    ///
    /// This is the shared core used by both [`replay`](Self::replay) and
    /// [`replay_from`](Self::replay_from).
    fn replay_range(
        &mut self,
        start_brep: BRep,
        from: usize,
        to: usize,
        tolerance: &Tolerance,
    ) -> Result<(BRep, HashMap<FeatureId, FeatureResult>), ReplayError> {
        let mut brep = start_brep;
        let mut results = HashMap::new();

        for i in from..to {
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

            match replay::execute_feature(&mut brep, &feature, &mut ctx, tolerance) {
                Ok(()) => {
                    let history = ctx.finish();
                    results.insert(
                        id,
                        FeatureResult::Success {
                            history: history.clone(),
                        },
                    );
                    self.features[id].last_result = Some(FeatureResult::Success { history });
                    // Store snapshot in cache.
                    self.cache.store(id, brep.clone(), i);
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

        Ok((brep, results))
    }
}

impl Default for FeatureTree {
    fn default() -> Self {
        Self::new()
    }
}
