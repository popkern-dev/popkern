//! [`TreeCommand`] enum definition — all reversible tree mutation variants.
//!
//! This module is separated from [`crate::undo`] purely to keep individual
//! source files below 500 lines.  All variants are consumed by
//! [`crate::undo::reverse_command`] and [`crate::undo::apply_command`].

use crate::config::{ConfigId, Configuration, ParamOverride};
use crate::expression::Expression;
use crate::feature::Feature;
use crate::node::{FeatureId, FeatureNode};
use crate::parameter::{Parameter, ParameterId};
use crate::unit::Unit;

/// A single reversible tree mutation recorded in the undo log.
#[derive(Debug, Clone)]
#[allow(dead_code)] // DeleteFeature wired in Task 7; snapshot kept for spec completeness
pub(crate) enum TreeCommand {
    /// A feature was added at a specific timeline index.
    AddFeature {
        /// The ID that was assigned to the new feature.
        id: FeatureId,
        /// Full snapshot of the node at the time of insertion.
        snapshot: FeatureNode,
        /// The timeline index at which the feature was inserted.
        index: usize,
    },
    /// One or more features were deleted.
    DeleteFeature {
        /// Each deleted feature's ID, node snapshot, and former timeline index.
        snapshots: Vec<(FeatureId, FeatureNode, usize)>,
    },
    /// A feature's parameters were modified.
    ModifyFeature {
        /// The ID of the modified feature.
        id: FeatureId,
        /// The previous feature parameters.
        old: Feature,
        /// The new feature parameters.
        new: Feature,
    },
    /// A feature was moved in the timeline.
    ReorderFeature {
        /// The ID of the moved feature.
        id: FeatureId,
        /// The timeline index before the move.
        old_index: usize,
        /// The timeline index after the move.
        new_index: usize,
    },
    /// A feature was suppressed.
    Suppress {
        /// The ID of the suppressed feature.
        id: FeatureId,
    },
    /// A feature was unsuppressed.
    Unsuppress {
        /// The ID of the unsuppressed feature.
        id: FeatureId,
    },
    /// The rollback bar was moved.
    MoveRollback {
        /// The previous rollback position.
        old: Option<usize>,
        /// The new rollback position.
        new: Option<usize>,
    },
    /// A feature was renamed.
    RenameFeature {
        /// The ID of the renamed feature.
        id: FeatureId,
        /// The label before the rename.
        old_label: String,
        /// The label after the rename.
        new_label: String,
    },
    /// An expression was set or cleared on a parameter.
    SetExpression {
        /// The parameter whose expression changed.
        param: ParameterId,
        /// The expression before the change (`None` = was literal).
        old: Option<Expression>,
        /// The expression after the change (`None` = cleared).
        new: Option<Expression>,
    },
    /// A parameter's value and/or unit was changed.
    SetParamValue {
        /// The parameter whose value changed.
        param: ParameterId,
        /// The value before the change (SI units).
        old_value: f64,
        /// The unit before the change.
        old_unit: Unit,
        /// The value after the change (SI units).
        new_value: f64,
        /// The unit after the change.
        new_unit: Unit,
    },
    /// Parameter bounds were changed or cleared.
    SetBounds {
        /// The parameter whose bounds changed.
        param: ParameterId,
        /// The bounds before the change.
        old: Option<(f64, f64)>,
        /// The bounds after the change.
        new: Option<(f64, f64)>,
    },
    /// A global parameter was added.
    AddGlobalParam {
        /// The stable ID assigned to the new parameter.
        id: ParameterId,
        /// Full snapshot of the parameter at insertion time.
        snapshot: Parameter,
    },
    /// A global parameter was removed.
    RemoveGlobalParam {
        /// The ID of the removed parameter.
        id: ParameterId,
        /// Full snapshot of the parameter before removal.
        snapshot: Parameter,
    },
    /// A parameter was renamed.
    RenameParam {
        /// The parameter that was renamed.
        id: ParameterId,
        /// The name before the rename.
        old_name: String,
        /// The name after the rename.
        new_name: String,
    },
    /// A suppress_if expression was set or cleared on a feature.
    SetSuppressIf {
        /// The feature whose suppress_if changed.
        feature: FeatureId,
        /// The expression before the change.
        old: Option<Expression>,
        /// The expression after the change.
        new: Option<Expression>,
    },
    /// A configuration was added (via add_config, add_derived_config).
    AddConfig {
        /// The ID assigned to the new configuration.
        id: ConfigId,
        /// Full snapshot of the configuration at insertion time.
        snapshot: Configuration,
    },
    /// A configuration was deleted.
    DeleteConfig {
        /// The ID of the deleted configuration.
        id: ConfigId,
        /// Full snapshot of the configuration before deletion.
        snapshot: Configuration,
        /// Children that were reparented: `(child_id, old_parent)`.
        reparented: Vec<(ConfigId, Option<ConfigId>)>,
    },
    /// A configuration was renamed.
    RenameConfig {
        /// The ID of the renamed configuration.
        id: ConfigId,
        /// The name before the rename.
        old_name: String,
        /// The name after the rename.
        new_name: String,
    },
    /// A configuration's description was changed.
    SetConfigDescription {
        /// The ID of the configuration.
        id: ConfigId,
        /// The description before the change.
        old: String,
        /// The description after the change.
        new: String,
    },
    /// A configuration was duplicated.
    DuplicateConfig {
        /// The ID assigned to the duplicate configuration.
        id: ConfigId,
        /// Full snapshot of the duplicate at insertion time.
        snapshot: Configuration,
    },
    /// A config parameter override was set or cleared.
    SetConfigOverride {
        /// The configuration being modified.
        config: ConfigId,
        /// The parameter whose override changed.
        param: ParameterId,
        /// The override before the change (`None` = no override existed).
        old: Option<ParamOverride>,
        /// The override after the change (`None` = override was cleared).
        new: Option<ParamOverride>,
    },
    /// A config feature suppression override was set or cleared.
    SetConfigSuppression {
        /// The configuration being modified.
        config: ConfigId,
        /// The feature whose suppression override changed.
        feature: FeatureId,
        /// The suppression state before the change (`None` = no override existed).
        old: Option<bool>,
        /// The suppression state after the change (`None` = override was cleared).
        new: Option<bool>,
    },
}
