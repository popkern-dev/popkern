// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! Tests for [`super::validate_seed_features`].
#[cfg(test)]
use super::validate_seed_features;

use std::collections::HashSet;

use popkern_core::naming::OperationId;
use slotmap::SlotMap;

use crate::error::HistoryError;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc};
use crate::node::{BodyId, FeatureId, FeatureNode};
use crate::parameter::ParameterId;

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

/// Insert a [`FeatureNode`] into `features` and append its key to `timeline`.
/// Returns the assigned [`FeatureId`].
fn insert(
    features: &mut SlotMap<FeatureId, FeatureNode>,
    timeline: &mut Vec<FeatureId>,
    feature: Feature,
    body: BodyId,
) -> FeatureId {
    let node = FeatureNode {
        feature,
        operation_id: OperationId(0),
        label: String::new(),
        depends_on: Vec::new(),
        suppressed: false,
        body,
        suppress_if: None,
        instance_op_ids: std::collections::HashMap::new(),
        last_result: None,
        selections: Vec::new(),
    };
    let id = features.insert(node);
    timeline.push(id);
    id
}

/// A [`Feature::Sketch`] placed on the XY plane (no dependencies).
fn sketch() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: crate::SketchPayload::Constraint(String::new()),
    }
}

/// A [`Feature::Extrude`] with `mode: Boss`.
///
/// `ParameterId` is a slotmap key so we fabricate a dummy one via a temporary
/// SlotMap.
fn extrude_boss() -> Feature {
    // Allocate a real ParameterId key.
    let mut p: slotmap::SlotMap<ParameterId, ()> = slotmap::SlotMap::with_key();
    let pid = p.insert(());
    Feature::Extrude {
        face: popkern_core::naming::EntityName {
            operation: OperationId(1),
            role: popkern_core::naming::EntityRole {
                category: popkern_core::naming::RoleCategory::Generated,
                index: 0,
            },
            ancestors: vec![popkern_core::naming::NameId(0)],
        },
        direction: None,
        end_condition: ExtrudeEnd::Blind { distance: pid },
        mode: OperationMode::Boss,
        draft_angle: None,
        second_direction: None,
    }
}

/// A [`Feature::Extrude`] with `mode: NewBody` (no Boss/Cut).
fn extrude_new_body() -> Feature {
    let mut p: slotmap::SlotMap<ParameterId, ()> = slotmap::SlotMap::with_key();
    let pid = p.insert(());
    Feature::Extrude {
        face: popkern_core::naming::EntityName {
            operation: OperationId(1),
            role: popkern_core::naming::EntityRole {
                category: popkern_core::naming::RoleCategory::Generated,
                index: 0,
            },
            ancestors: vec![popkern_core::naming::NameId(0)],
        },
        direction: None,
        end_condition: ExtrudeEnd::Blind { distance: pid },
        mode: OperationMode::NewBody,
        draft_angle: None,
        second_direction: None,
    }
}

// ---------------------------------------------------------------------------
// Rule 1: EmptySeedFeatures
// ---------------------------------------------------------------------------

#[test]
fn empty_seed_list_is_rejected() {
    let features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();

    let result = validate_seed_features(&[], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::EmptySeedFeatures)),
        "expected EmptySeedFeatures, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 2: FeatureNotFound (all must exist)
// ---------------------------------------------------------------------------

#[test]
fn missing_feature_id_is_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();

    let body = BodyId(0);
    let _real = insert(&mut features, &mut timeline, sketch(), body);

    // Remove a FeatureId from features so it no longer exists.
    let ghost_id = features.insert(FeatureNode {
        feature: sketch(),
        operation_id: OperationId(0),
        label: String::new(),
        depends_on: Vec::new(),
        suppressed: false,
        body,
        suppress_if: None,
        instance_op_ids: std::collections::HashMap::new(),
        last_result: None,
        selections: Vec::new(),
    });
    features.remove(ghost_id);

    let result = validate_seed_features(&[ghost_id], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::FeatureNotFound)),
        "expected FeatureNotFound, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 3: SeedBodyMismatch
// ---------------------------------------------------------------------------

#[test]
fn seeds_on_different_bodies_are_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();

    let id0 = insert(&mut features, &mut timeline, sketch(), BodyId(0));
    let id1 = insert(&mut features, &mut timeline, extrude_boss(), BodyId(1));

    let result = validate_seed_features(&[id0, id1], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::SeedBodyMismatch)),
        "expected SeedBodyMismatch, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 4: SeedNotContiguous
// ---------------------------------------------------------------------------

#[test]
fn non_contiguous_seeds_are_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    // Timeline: sketch, filler, extrude — sketch and extrude are not adjacent.
    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);
    let _filler = insert(&mut features, &mut timeline, extrude_new_body(), body);
    let extrude_id = insert(&mut features, &mut timeline, extrude_boss(), body);

    let result = validate_seed_features(&[sketch_id, extrude_id], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::SeedNotContiguous)),
        "expected SeedNotContiguous, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 5: SeedMissingSketch
// ---------------------------------------------------------------------------

#[test]
fn seed_chain_not_starting_with_sketch_is_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    // Only an extrude — no leading Sketch.
    let extrude_id = insert(&mut features, &mut timeline, extrude_boss(), body);

    let result = validate_seed_features(&[extrude_id], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::SeedMissingSketch)),
        "expected SeedMissingSketch, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 6: SeedMissingMode
// ---------------------------------------------------------------------------

#[test]
fn seed_chain_with_no_boss_or_cut_is_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);
    // extrude_new_body has mode=NewBody, not Boss/Cut.
    let extrude_id = insert(&mut features, &mut timeline, extrude_new_body(), body);

    let result = validate_seed_features(&[sketch_id, extrude_id], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::SeedMissingMode)),
        "expected SeedMissingMode, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Rule 7: SeedAlreadyOwned
// ---------------------------------------------------------------------------

#[test]
fn seed_already_owned_by_another_pattern_is_rejected() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let body = BodyId(0);

    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);
    let extrude_id = insert(&mut features, &mut timeline, extrude_boss(), body);

    // Mark extrude_id as already owned.
    let mut existing: HashSet<FeatureId> = HashSet::new();
    existing.insert(extrude_id);

    let result = validate_seed_features(&[sketch_id, extrude_id], &features, &timeline, &existing);
    assert!(
        matches!(result, Err(HistoryError::SeedAlreadyOwned)),
        "expected SeedAlreadyOwned, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Happy path
// ---------------------------------------------------------------------------

#[test]
fn valid_seed_chain_passes_all_rules() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);
    let extrude_id = insert(&mut features, &mut timeline, extrude_boss(), body);

    let result = validate_seed_features(&[sketch_id, extrude_id], &features, &timeline, &existing);
    assert!(result.is_ok(), "expected Ok, got {result:?}");
}

#[test]
fn single_seed_with_boss_extrude_starting_with_sketch_passes() {
    // The spec says "Starts with Sketch", meaning position[0] in timeline order
    // is a Sketch. If the seed list is just [sketch_id] it does NOT have a
    // Boss/Cut mode, so this tests that the sketch alone fails rule 6.
    // A single-feature seed must still include a Boss/Cut feature.
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);

    let result = validate_seed_features(&[sketch_id], &features, &timeline, &existing);
    // A lone Sketch has no Boss/Cut mode.
    assert!(
        matches!(result, Err(HistoryError::SeedMissingMode)),
        "expected SeedMissingMode for sketch-only seed, got {result:?}"
    );
}

#[test]
fn cut_mode_satisfies_boss_or_cut_rule() {
    let mut features: SlotMap<FeatureId, FeatureNode> = SlotMap::with_key();
    let mut timeline: Vec<FeatureId> = Vec::new();
    let existing: HashSet<FeatureId> = HashSet::new();
    let body = BodyId(0);

    let mut p: slotmap::SlotMap<ParameterId, ()> = slotmap::SlotMap::with_key();
    let pid = p.insert(());
    let extrude_cut = Feature::Extrude {
        face: popkern_core::naming::EntityName {
            operation: OperationId(1),
            role: popkern_core::naming::EntityRole {
                category: popkern_core::naming::RoleCategory::Generated,
                index: 0,
            },
            ancestors: vec![popkern_core::naming::NameId(0)],
        },
        direction: None,
        end_condition: ExtrudeEnd::Blind { distance: pid },
        mode: OperationMode::Cut,
        draft_angle: None,
        second_direction: None,
    };

    let sketch_id = insert(&mut features, &mut timeline, sketch(), body);
    let extrude_id = insert(&mut features, &mut timeline, extrude_cut, body);

    let result = validate_seed_features(&[sketch_id, extrude_id], &features, &timeline, &existing);
    assert!(
        result.is_ok(),
        "Cut mode should satisfy rule 6, got {result:?}"
    );
}
