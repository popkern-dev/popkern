//! The transactional record loop (spec §1.6): append one feature,
//! replay, and either serve or remove the feature again.

use popkern_history::error::FeatureResult;
use popkern_history::feature::Feature;
use popkern_history::node::FeatureId;
use popkern_history::{DeleteBehavior, FeatureTree, ReplayError};

use crate::error::ScriptError;

/// Append `feature`, replay the document, and roll the append back on
/// failure. Returns the new feature's id.
pub(crate) fn record_feature(
    tree: &mut FeatureTree,
    tolerance: popkern_core::tolerance::Tolerance,
    feature: Feature,
    label: &str,
    body: Option<popkern_history::node::BodyId>,
    selections: Vec<popkern_history::SelectionRecord>,
) -> Result<FeatureId, ScriptError> {
    let operation = feature.operation_label().to_string();
    // Recording at a rollback bar inserts there (migration `insert_*`).
    let id = if tree.rollback_position().is_some() {
        tree.insert_feature_at_rollback(feature, label)?
    } else {
        tree.add_feature_on(feature, label, body)
    };
    let _ = tree.set_selection_records(id, selections);
    let minted = tree
        .feature(id)
        .map(|n| n.label.clone())
        .unwrap_or_default();
    let outcome = tree.replay(tolerance);
    let failure = match outcome {
        Ok(result) => match result.feature_results.get(&id) {
            Some(FeatureResult::Failed { error }) => Some(error.clone()),
            _ => None,
        },
        Err(ReplayError::FeatureFailed {
            label: failed,
            error,
        }) if failed == minted => Some(error),
        Err(other) => Some(other.to_string()),
    };
    if let Some(detail) = failure {
        // Transactional: the document is as it was before the call.
        let _ = tree.delete_feature(id, DeleteBehavior::Cascade);
        let _ = tree.replay(tolerance);
        return Err(ScriptError::FeatureFailed {
            label: minted,
            operation,
            detail,
        });
    }
    Ok(id)
}
