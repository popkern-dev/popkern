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

pub(crate) mod auto_merge;
pub(crate) mod body;
mod body_query;
pub mod cache;
pub mod config;
pub mod config_api;
pub(crate) mod config_eval;
pub mod csv;
pub mod deferred_certification;
pub mod delete;
mod document;
pub mod error;
pub mod expression;
pub mod failure;
pub mod failure_blend;
pub mod feature;
pub mod feature_blend;
pub mod feature_draft;
pub mod feature_extrude;
pub mod feature_fastener;
mod feature_params;
pub mod feature_sketch;
pub mod feature_snapshot;
pub mod feature_transform;
pub mod format;
pub mod health;
pub mod labels;
pub(crate) mod measure_cache;
pub mod measure_query;
pub(crate) mod measure_resolve;
pub(crate) mod modify_inplace_route;
pub mod node;
pub mod param_api;
pub mod parameter;
pub mod parse;
pub mod pattern;
pub mod pattern_api;
pub mod pattern_api_extras;
pub(crate) mod pattern_execute;
pub(crate) mod pattern_replay;
pub(crate) mod pattern_transform;
pub(crate) mod pattern_validate;
pub mod queries;
pub(crate) mod replay;
mod replay_exact_extrude;
mod replay_exact_extrude_names;
pub mod replay_report;
mod replay_side_effects;
mod tree_replay_config;
mod tree_state;
pub(crate) use replay_report::ReplayRangeOutput;
pub(crate) mod resolve;
pub mod serialize;
pub mod timeline;
mod tree_replay;
pub(crate) mod undo;
pub(crate) mod undo_command;
pub(crate) mod undo_config;
pub mod unit;

pub use config::{ConfigId, Configuration, DeleteConfigResult, ParamOverride};
pub use delete::{DeleteBehavior, DeleteResult};
pub use error::{FeatureResult, HistoryError, ReplayError};
pub use expression::{BinOp, CmpOp, Expression, UnaryFunc};
pub use failure::{
    BlendDeclineKind, BlendDeclined, BlendFailureReason, BlendFailureRecord, FeatureFailure,
    RemediationRecord,
};
pub use feature::{
    AxisDesc, BasePlane, BoolOpKind, ExtrudeEnd, ExtrudeSecondDir, Feature, JoinStyleKind,
    OperationMode, RevolveEnd, ShellDirectionKind, SketchPlaneDesc, TrimPolicyKind, Vec3Desc,
};
pub use feature_blend::{
    BlendProfileKind, CapKindDesc, ChamferKindDesc, ContinuityKind, CornerModeDesc,
    FilletRadiusLaw, FilletRecipeDesc, InterferencePolicyKind, RadiusStation,
};
pub use feature_snapshot::SnapshotOrigin;
pub use format::{format_expression, format_expression_with_unit};
pub use health::{HealthReport, ReferenceIssue};
pub use labels::LabelTable;
pub use measure_query::{BoundExtent, Component, MeasureContext, MeasureQuery};
pub use node::{BodyId, FeatureId, FeatureNode, SelectionRecord};
pub use parameter::{EvaluationWarning, ParamScope, Parameter, ParameterId, ParameterTable};
pub use parse::{EntityResolver, ParseError, parse_expression, parse_expression_with_resolver};
pub use queries::ImpactReport;
pub use replay_report::{
    BlendEvidenceRecord, FallbackRecord, FeatureReport, KernelWarning, RecordSource, ReplayOutcome,
    ReplayReport, ReplayResult, SkipRecord, SuppressionReason,
};
pub use serialize::{CURRENT_VERSION, MIN_SUPPORTED_VERSION};
pub use unit::{DimensionKind, Unit};
pub use {cache::CacheStrategy, csv::CsvImportResult};
pub use {feature_draft::NeutralPlaneDesc, feature_transform::RigidMotionDesc};
pub use {
    feature_sketch::AnalyticSketchDesc, feature_sketch::ScalarBinding,
    feature_sketch::SketchPayload,
};

use std::collections::HashMap;

use popkern_core::naming::OperationId;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;

use crate::body::BodyManager;
use crate::cache::BRepCache;
use crate::undo::TreeCommand;

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
    /// Timeline index below which features are "frozen" — never rebuilt.
    /// `None` means nothing is frozen.
    #[serde(default)]
    pub(crate) frozen_index: Option<usize>,
    /// Base values of literal parameters currently pinned by the
    /// active configuration's overrides. While a configuration is
    /// active the parameter table shows the EFFECTIVE values (what
    /// replayed); the next base replay restores these. Persisted so a
    /// document saved under a configuration keeps its base values.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub(crate) config_pinned_base: HashMap<parameter::ParameterId, f64>,
    /// The BRep produced by the most recent replay, or `None` if the tree
    /// has not been replayed yet.
    ///
    /// Not serialized — call [`replay`](Self::replay) after loading.
    #[serde(skip)]
    pub(crate) current_brep: Option<BRep>,
    /// The most recent report for every feature in the timeline
    /// (session-local; rebuilt by replay).
    #[serde(skip)]
    pub(crate) last_report: ReplayReport,
    /// The body each Boss/Cut feature merged into during the last
    /// replay (transient; rebuilt by every replay).
    #[serde(skip)]
    pub(crate) merged_into: HashMap<FeatureId, BodyId>,
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
    /// Session-local deferred seam certification (never serialized).
    #[serde(skip)]
    pub(crate) deferred: deferred_certification::DeferredCertification,
    /// Index from [`OperationId`] to [`FeatureId`] for O(1) dependency lookup.
    ///
    /// Not serialized — rebuilt from features on load via
    /// [`rebuild_indices`](Self::rebuild_indices).
    #[serde(skip)]
    pub(crate) op_index: HashMap<OperationId, FeatureId>,
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
    /// All named parameters used by features in this tree.
    ///
    /// Serialized as part of the tree's persistent state.
    pub(crate) parameter_table: ParameterTable,
    /// Named configurations (design table variants).
    #[serde(default)]
    pub(crate) configs: SlotMap<ConfigId, Configuration>,
    /// Currently active configuration. None = default (no overrides).
    /// Serialized — restores the last-active config when the file is opened.
    #[serde(default)]
    pub(crate) active_config: Option<ConfigId>,
    /// Entity labels and auto-label counters (see [`labels`]).
    #[serde(default)]
    pub(crate) labels: LabelTable,
    /// The tolerance the document is built with (see [`document`]
    /// attributes); `None` means [`Tolerance::default`].
    #[serde(default)]
    pub(crate) build_tolerance: Option<Tolerance>,
}

impl FeatureTree {
    /// Create an empty feature tree with the default [`CacheStrategy::EveryStep`] cache.
    pub fn new() -> Self {
        Self {
            features: SlotMap::with_key(),
            timeline: Vec::new(),
            rollback_index: None,
            next_op_id: 1,
            frozen_index: None,
            current_brep: None,
            last_report: ReplayReport::default(),

            merged_into: HashMap::new(),
            cache: BRepCache::new(CacheStrategy::EveryStep),
            body_manager: BodyManager::new(),
            deferred: Default::default(),
            op_index: HashMap::new(),
            undo_log: Vec::new(),
            undo_cursor: 0,
            parameter_table: ParameterTable::new(),
            configs: SlotMap::with_key(),
            active_config: None,
            labels: LabelTable::default(),
            build_tolerance: None,
            config_pinned_base: HashMap::new(),
        }
    }

    /// Return a shared reference to the parameter table.
    pub fn parameters(&self) -> &ParameterTable {
        &self.parameter_table
    }

    /// Return a mutable reference to the parameter table.
    ///
    /// Prefer the undo-safe methods (`set_parameter_value`, `add_global_parameter`,
    /// etc.) for mutations that should be reversible.
    pub fn parameters_mut(&mut self) -> &mut ParameterTable {
        &mut self.parameter_table
    }

    // -- Replay ---------------------------------------------------------------

    /// Replay the entire feature tree from scratch.
    ///
    /// Evaluates features in timeline order up to the rollback index (or the
    /// end of the timeline). Suppressed features are skipped. Failed features
    /// record an error but do not abort the replay -- the BRep continues from
    /// the last good state.
    ///
    /// # Errors
    ///
    /// Returns [`ReplayError`] only for fatal infrastructure failures.
    /// Individual feature failures are recorded in
    /// [`ReplayResult::feature_results`].
    pub fn replay(&mut self, tolerance: Tolerance) -> Result<ReplayResult, ReplayError> {
        let start = popkern_core::time::Stopwatch::start();

        // Evaluate parameters — config-aware when a config is active.
        let eval_warnings = if self.active_config.is_some() {
            self.apply_config_and_evaluate()
        } else {
            self.restore_config_pinned_base();
            expression::evaluate_all(
                &mut self.parameter_table,
                &measure_query::MeasureContext::empty(),
            )
        };

        let end = self.rollback_index.unwrap_or(self.timeline.len());

        // If frozen, try to start from the cached BRep at the freeze point.
        let (start_bodies, from) = if let Some(frozen) = self.frozen_index {
            let frozen = frozen.min(end);
            if frozen > 0 {
                let freeze_id = self.timeline[frozen - 1];
                if let Some(cached) = self.cache.get(self.active_config, freeze_id) {
                    let body = self.features[freeze_id].body;
                    let mut bodies = HashMap::new();
                    bodies.insert(body, cached.clone());
                    (bodies, frozen)
                } else {
                    // No cache at freeze point — replay from scratch.
                    (HashMap::new(), 0)
                }
            } else {
                (HashMap::new(), 0)
            }
        } else {
            (HashMap::new(), 0)
        };

        let out = self.replay_range(start_bodies, from, end, &tolerance)?;
        self.store_body_breps(&out.bodies);
        self.last_report.absorb(out.report.clone(), &self.timeline);
        Ok(ReplayResult {
            feature_results: out.results,
            elapsed: start.elapsed(),
            eval_warnings,
            report: out.report,
            tolerance,
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
        let mut start_bodies: HashMap<BodyId, BRep> = HashMap::new();
        for i in (0..pos).rev() {
            let fid = self.timeline[i];
            if let Some(cached) = self.cache.get(self.active_config, fid) {
                start_pos = i + 1;
                let body = self.features[fid].body;
                start_bodies.insert(body, cached.clone());
                break;
            }
        }

        let timer = popkern_core::time::Stopwatch::start();
        if self.active_config.is_none() {
            self.restore_config_pinned_base();
        }
        let eval_warnings = expression::evaluate_all(
            &mut self.parameter_table,
            &measure_query::MeasureContext::empty(),
        );
        let end = self.rollback_index.unwrap_or(self.timeline.len());
        let out = self.replay_range(start_bodies, start_pos, end, &tolerance)?;
        self.store_body_breps(&out.bodies);
        self.last_report.absorb(out.report.clone(), &self.timeline);
        Ok(ReplayResult {
            feature_results: out.results,
            elapsed: timer.elapsed(),
            eval_warnings,
            report: out.report,
            tolerance,
        })
    }

    // -- Rollback / scrubbing -------------------------------------------------

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
            self.cache.get(self.active_config, id)
        } else {
            None
        }
    }

    /// Remove the rollback bar, making all features in the timeline active.
    pub fn rollback_to_end(&mut self) {
        let old = self.rollback_index;
        self.rollback_index = None;
        self.record_undo(TreeCommand::MoveRollback { old, new: None });
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
            if self.cache.get(self.active_config, id).is_some() {
                // SAFETY: we just checked `is_some()`. The two-step pattern
                // (check then re-borrow) is required to avoid holding an
                // immutable borrow on `self.cache` across the mutable
                // `self.replay()` call below.
                return self
                    .cache
                    .get(self.active_config, id)
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

    // -- BRep accessor --------------------------------------------------------

    /// Return a reference to the BRep produced by the most recent replay.
    ///
    /// Returns `None` if the tree has not been replayed yet.
    pub fn last_report(&self) -> &ReplayReport {
        &self.last_report
    }

    /// The current BRep after the last replay.
    pub fn current_brep(&self) -> Option<&BRep> {
        self.current_brep.as_ref()
    }
}

impl Default for FeatureTree {
    fn default() -> Self {
        Self::new()
    }
}
