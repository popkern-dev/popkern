//! The transactional record loop (spec §1.6): append one feature,
//! replay, and either serve or remove the feature again.

use popkern_history::error::FeatureResult;
use popkern_history::feature::Feature;
use popkern_history::node::FeatureId;
use popkern_history::{DeleteBehavior, FeatureFailure, FeatureTree, ReplayError};

use crate::error::ScriptError;

/// Append `feature`, replay the document, and roll the append back on
/// failure. Returns the new feature's id.
///
/// `geometry_current` is the recorder's freshness bit for the append
/// fast path ([`FeatureTree::replay_appended`]): true means the
/// tree's stored body snapshots reflect a completed replay of the
/// whole timeline with no parameter/config/timeline mutation since.
/// Recording appends preserve it; the parameter and configuration
/// setters clear it; every full replay a recorder verb performs
/// restores it. When it is false (or the record inserts at a
/// rollback bar) the record replays the full document exactly as
/// before.
pub(crate) fn record_feature(
    tree: &mut FeatureTree,
    tolerance: popkern_core::tolerance::Tolerance,
    feature: Feature,
    label: &str,
    body: Option<popkern_history::node::BodyId>,
    selections: Vec<popkern_history::SelectionRecord>,
    geometry_current: &mut bool,
) -> Result<FeatureId, ScriptError> {
    let operation = feature.operation_label().to_string();
    // Recording at a rollback bar inserts there (migration `insert_*`).
    let at_bar = tree.rollback_position().is_some();
    let id = if at_bar {
        tree.insert_feature_at_rollback(feature, label)?
    } else {
        tree.add_feature_on(feature, label, body)
    };
    let _ = tree.set_selection_records(id, selections);
    let minted = tree
        .feature(id)
        .map(|n| n.label.clone())
        .unwrap_or_default();
    let outcome = if *geometry_current && !at_bar {
        tree.replay_appended(tolerance)
    } else {
        tree.replay(tolerance)
    };
    *geometry_current = tree.rollback_position().is_none();
    let failure = match outcome {
        Ok(result) => match result.feature_results.get(&id) {
            Some(FeatureResult::Failed { error }) => Some(Box::new(error.clone())),
            _ => None,
        },
        Err(ReplayError::FeatureFailed {
            label: failed,
            error,
        }) if failed == minted => Some(error),
        Err(other) => Some(Box::new(FeatureFailure::from(other))),
    };
    if let Some(error) = failure {
        let report = Box::new(tree.last_report().clone());
        // Transactional: the document is as it was before the call.
        let _ = tree.delete_feature(id, DeleteBehavior::Cascade);
        let _ = tree.replay(tolerance);
        return Err(ScriptError::FeatureFailed {
            label: minted,
            operation,
            error,
            report,
        });
    }
    Ok(id)
}
