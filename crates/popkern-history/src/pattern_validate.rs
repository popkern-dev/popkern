//! Seed feature validation for pattern features.
//!
//! All four pattern variants (LinearPattern, CircularPattern, CurvePattern,
//! Mirror) share the same rules for which features may serve as seeds.  This
//! module centralises those checks so each pattern variant can call a single
//! function rather than duplicating logic.

use std::collections::HashSet;

use slotmap::SlotMap;

use crate::auto_merge::feature_mode;
use crate::error::HistoryError;
use crate::feature::{Feature, OperationMode};
use crate::node::{FeatureId, FeatureNode};

#[allow(dead_code)] // Called by pattern insertion helpers; callers land in a future task.
/// Validate that a proposed set of seed features is well-formed for use in a
/// pattern.
///
/// # Rules checked
///
/// 1. **Non-empty** — `seed_features` must contain at least one entry.
/// 2. **All exist** — every [`FeatureId`] must be present in `features`.
/// 3. **Same body** — all seed features must belong to the same [`crate::node::BodyId`].
/// 4. **Contiguous timeline** — the seeds must occupy consecutive positions in
///    `timeline` with no gaps between them.
/// 5. **Starts with Sketch** — the earliest seed (by timeline order) must be a
///    [`Feature::Sketch`].
/// 6. **Boss or Cut mode** — at least one seed must carry
///    [`OperationMode::Boss`] or [`OperationMode::Cut`].
/// 7. **No shared seeds** — no seed feature may already appear in
///    `existing_seeds` (i.e. be owned by a different pattern).
///
/// # Errors
///
/// Returns the first validation rule that is violated, as a [`HistoryError`]
/// variant.
pub(crate) fn validate_seed_features(
    seed_features: &[FeatureId],
    features: &SlotMap<FeatureId, FeatureNode>,
    timeline: &[FeatureId],
    existing_seeds: &HashSet<FeatureId>,
) -> Result<(), HistoryError> {
    // Rule 1: non-empty.
    if seed_features.is_empty() {
        return Err(HistoryError::EmptySeedFeatures);
    }

    // Rule 2: all exist; collect nodes at the same time.
    for &id in seed_features {
        if !features.contains_key(id) {
            return Err(HistoryError::FeatureNotFound);
        }
    }

    // Rule 3: all on the same body.
    let first_body = features[seed_features[0]].body;
    for &id in seed_features.iter().skip(1) {
        if features[id].body != first_body {
            return Err(HistoryError::SeedBodyMismatch);
        }
    }

    // Rule 4: contiguous in the timeline.
    // Find the timeline positions of all seed features.
    let mut positions: Vec<usize> = seed_features
        .iter()
        .filter_map(|id| timeline.iter().position(|t| t == id))
        .collect();

    // Every seed must actually appear in the timeline.
    if positions.len() != seed_features.len() {
        return Err(HistoryError::FeatureNotFound);
    }

    positions.sort_unstable();

    // Check there are no gaps.
    for window in positions.windows(2) {
        if window[1] != window[0] + 1 {
            return Err(HistoryError::SeedNotContiguous);
        }
    }

    // Rule 5: the first seed (by timeline position) must be a Sketch.
    let first_timeline_id = timeline[positions[0]];
    match &features[first_timeline_id].feature {
        Feature::Sketch { .. } => {}
        _ => return Err(HistoryError::SeedMissingSketch),
    }

    // Rule 6: at least one seed must have Boss or Cut mode.
    let has_mode = seed_features.iter().any(|&id| {
        matches!(
            feature_mode(&features[id].feature),
            OperationMode::Boss | OperationMode::Cut
        )
    });
    if !has_mode {
        return Err(HistoryError::SeedMissingMode);
    }

    // Rule 7: no seed may already be owned by another pattern.
    for &id in seed_features {
        if existing_seeds.contains(&id) {
            return Err(HistoryError::SeedAlreadyOwned);
        }
    }

    Ok(())
}

#[cfg(test)]
#[path = "pattern_validate_tests.rs"]
mod tests;
