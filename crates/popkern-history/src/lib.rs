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
//!   └── cache: StepCache                  // per-step document states
//! ```

pub mod appearance;
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
pub mod document_state;
pub mod error;
pub mod expression;
pub mod failure;
pub mod failure_blend;
pub mod feature;
pub mod feature_blend;
pub mod feature_draft;
pub mod feature_extrude;
pub mod feature_fastener;
pub mod feature_hinge;
pub(crate) mod feature_names;
mod feature_params;
mod feature_shell_kinds;
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
pub mod reference_rebind;
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
mod step_key;
pub mod timeline;
mod tree_replay;
mod tree_scrub;
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
    BlendProfileKind, CapKindDesc, ChamferKindDesc, ContinuityKind, CornerCapDesc,
    CornerConstructionDesc, CornerModeDesc, EdgeRadiusGroup, FilletRadiusLaw, FilletRecipeDesc,
    InterferencePolicyKind, RadiusStation,
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
    BlendEvidenceRecord, CornerEvidence, FallbackRecord, FeatureReport, KernelWarning,
    RecordSource, ReplayOutcome, ReplayReport, ReplayResult, SkipRecord, SuppressionReason,
};
pub use serialize::{CURRENT_VERSION, MIN_SUPPORTED_VERSION};
pub use unit::{DimensionKind, Unit};
pub use {cache::CacheStrategy, csv::CsvImportResult, document_state::DocumentState};
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
use crate::cache::StepCache;
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
    /// The part's colour and finish, when set (display and exchange
    /// only; replay ignores it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) appearance: Option<appearance::Appearance>,
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
    /// The step cache: every replayed step's document state and
    /// published outputs, content-addressed (see [`cache`]).
    ///
    /// Not serialized — rebuilt on demand after loading.
    #[serde(skip)]
    pub(crate) cache: StepCache,
    /// Tracks the mapping from [`BodyId`]s to features and BRep snapshots.
    ///
    /// Not serialized — rebuilt from the feature nodes' bodies on load.
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
    /// The naming scheme this tree's references were minted under, and
    /// the rules its replays mint with (PK-677). A new tree mints under
    /// [`popkern_core::naming::NAMING_SCHEME_VERSION`]; a tree saved
    /// before the field existed was scheme 1.
    #[serde(default = "scheme_one")]
    pub(crate) naming_scheme: u32,
}

fn scheme_one() -> u32 {
    1
}

impl FeatureTree {
    /// The part's appearance, when set.
    pub fn appearance(&self) -> Option<appearance::Appearance> {
        self.appearance
    }

    /// Set (or clear) the part's appearance. Geometry is untouched, so
    /// nothing replays.
    pub fn set_appearance(&mut self, appearance: Option<appearance::Appearance>) {
        self.appearance = appearance;
    }

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
            cache: StepCache::new(CacheStrategy::EveryStep),
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
            appearance: None,
            naming_scheme: popkern_core::naming::NAMING_SCHEME_VERSION,
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

    // -- BRep accessor --------------------------------------------------------

    /// Return a reference to the BRep produced by the most recent replay.
    ///
    /// Returns `None` if the tree has not been replayed yet.
    pub fn last_report(&self) -> &ReplayReport {
        &self.last_report
    }

    /// The naming scheme this tree mints under (PK-677).
    pub fn naming_scheme(&self) -> u32 {
        self.naming_scheme
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

#[cfg(test)]
#[path = "quantity_literal_tests.rs"]
mod quantity_literal_tests;
