// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use crate::FeatureTree;
use crate::feature::{Feature, SketchPlaneDesc};

/// Convenience: create a simple Sketch feature on the XY plane.
fn sketch_xy() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: String::new(),
    }
}

/// Convenience: create a simple Sketch feature on the XZ plane.
fn sketch_xz() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XZ,
        sketch_data: String::new(),
    }
}

// ---------------------------------------------------------------------------
// 1. Add → undo → feature removed from timeline → redo → feature back
// ---------------------------------------------------------------------------

#[test]
fn undo_add_feature_removes_from_timeline() {
    let mut tree = FeatureTree::new();
    assert!(!tree.can_undo());

    let id = tree.add_feature(sketch_xy(), "Sketch 1");
    assert_eq!(tree.timeline(), &[id]);
    assert!(tree.can_undo());
    assert!(!tree.can_redo());

    let undone = tree.undo();
    assert!(undone);
    assert!(
        tree.timeline().is_empty(),
        "timeline should be empty after undo"
    );
    assert!(!tree.can_undo());
    assert!(tree.can_redo());
}

#[test]
fn redo_add_feature_restores_timeline() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");

    tree.undo();
    assert!(tree.timeline().is_empty());

    let redone = tree.redo();
    assert!(redone);
    assert_eq!(
        tree.timeline(),
        &[id],
        "timeline should contain feature after redo"
    );
    assert!(!tree.can_redo());
    assert!(tree.can_undo());
}

// ---------------------------------------------------------------------------
// 2. Modify → undo → old params → redo → new params
// ---------------------------------------------------------------------------

#[test]
fn undo_modify_restores_old_params() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");

    // Clear the undo log so only the modify is tracked.
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.modify_feature(
        id,
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: "modified".to_string(),
        },
    );

    // After modify, the feature should have the new plane.
    assert!(matches!(
        tree.feature(id).unwrap().feature,
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            ..
        }
    ));

    tree.undo();

    assert!(
        matches!(
            tree.feature(id).unwrap().feature,
            Feature::Sketch {
                plane: SketchPlaneDesc::XY,
                ..
            }
        ),
        "undo should restore original XY plane"
    );
}

#[test]
fn redo_modify_restores_new_params() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.modify_feature(
        id,
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: "modified".to_string(),
        },
    );
    tree.undo();
    tree.redo();

    assert!(
        matches!(
            tree.feature(id).unwrap().feature,
            Feature::Sketch {
                plane: SketchPlaneDesc::XZ,
                ..
            }
        ),
        "redo should restore new XZ plane"
    );
}

// ---------------------------------------------------------------------------
// 3. Multiple undo/redo cycles
// ---------------------------------------------------------------------------

#[test]
fn multiple_undo_redo_cycles() {
    let mut tree = FeatureTree::new();

    let id0 = tree.add_feature(sketch_xy(), "A");
    let id1 = tree.add_feature(sketch_xz(), "B");

    // Undo both.
    assert!(tree.undo()); // undo B
    assert_eq!(tree.timeline(), &[id0]);
    assert!(tree.undo()); // undo A
    assert!(tree.timeline().is_empty());
    assert!(!tree.undo()); // nothing left

    // Redo both.
    assert!(tree.redo()); // redo A
    assert_eq!(tree.timeline(), &[id0]);
    assert!(tree.redo()); // redo B
    assert_eq!(tree.timeline(), &[id0, id1]);
    assert!(!tree.redo()); // nothing left
}

// ---------------------------------------------------------------------------
// 4. Suppress → undo → unsuppressed
// ---------------------------------------------------------------------------

#[test]
fn undo_suppress_restores_unsuppressed() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");
    // Clear log so we only test suppress.
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    assert!(!tree.feature(id).unwrap().suppressed);
    tree.suppress(id);
    assert!(tree.feature(id).unwrap().suppressed);

    tree.undo();
    assert!(
        !tree.feature(id).unwrap().suppressed,
        "undo should unsuppress"
    );
}

#[test]
fn redo_suppress_re_suppresses() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.suppress(id);
    tree.undo();
    assert!(!tree.feature(id).unwrap().suppressed);

    tree.redo();
    assert!(
        tree.feature(id).unwrap().suppressed,
        "redo should re-suppress"
    );
}

// ---------------------------------------------------------------------------
// 5. Rename → undo → old name
// ---------------------------------------------------------------------------

#[test]
fn undo_rename_restores_old_label() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Original");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.rename_feature(id, "New Name");
    assert_eq!(tree.feature(id).unwrap().label, "New Name");

    tree.undo();
    assert_eq!(
        tree.feature(id).unwrap().label,
        "Original",
        "undo should restore original label"
    );
}

#[test]
fn redo_rename_restores_new_label() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Original");
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.rename_feature(id, "New Name");
    tree.undo();
    tree.redo();
    assert_eq!(
        tree.feature(id).unwrap().label,
        "New Name",
        "redo should restore new label"
    );
}

// ---------------------------------------------------------------------------
// 6. New mutation truncates redo history
// ---------------------------------------------------------------------------

#[test]
fn new_mutation_truncates_redo_history() {
    let mut tree = FeatureTree::new();
    let id0 = tree.add_feature(sketch_xy(), "A");
    let _id1 = tree.add_feature(sketch_xz(), "B");

    // Undo B, then add a new feature — redo history for B should be gone.
    tree.undo();
    assert_eq!(tree.timeline(), &[id0]);

    let id2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::YZ,
            sketch_data: String::new(),
        },
        "C",
    );
    assert_eq!(tree.timeline(), &[id0, id2]);
    assert!(!tree.can_redo(), "redo history should have been truncated");
}

// ---------------------------------------------------------------------------
// 7. can_undo / can_redo correctness
// ---------------------------------------------------------------------------

#[test]
fn can_undo_can_redo_boundaries() {
    let mut tree = FeatureTree::new();
    assert!(!tree.can_undo());
    assert!(!tree.can_redo());

    tree.add_feature(sketch_xy(), "A");
    assert!(tree.can_undo());
    assert!(!tree.can_redo());

    tree.undo();
    assert!(!tree.can_undo());
    assert!(tree.can_redo());

    tree.redo();
    assert!(tree.can_undo());
    assert!(!tree.can_redo());
}

// ---------------------------------------------------------------------------
// 8. Unsuppress undo/redo
// ---------------------------------------------------------------------------

#[test]
fn undo_unsuppress_re_suppresses() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "Sketch 1");
    tree.suppress(id);
    // Now unsuppress and undo.
    tree.undo_log.clear();
    tree.undo_cursor = 0;

    tree.unsuppress(id);
    assert!(!tree.feature(id).unwrap().suppressed);

    tree.undo();
    assert!(
        tree.feature(id).unwrap().suppressed,
        "undo unsuppress should re-suppress"
    );

    tree.redo();
    assert!(
        !tree.feature(id).unwrap().suppressed,
        "redo unsuppress should unsuppress again"
    );
}

// ---------------------------------------------------------------------------
// 9. Reorder undo/redo
// ---------------------------------------------------------------------------

#[test]
fn undo_reorder_restores_original_position() {
    let mut tree = FeatureTree::new();
    let id0 = tree.add_feature(sketch_xy(), "A");
    let id1 = tree.add_feature(sketch_xz(), "B");
    let id2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::YZ,
            sketch_data: String::new(),
        },
        "C",
    );
    // Timeline: [id0, id1, id2]. Move id0 to position 2.
    tree.reorder(id0, 2).expect("reorder should succeed");
    assert_eq!(tree.timeline(), &[id1, id2, id0]);

    tree.undo(); // undo reorder
    tree.undo(); // undo add C
    tree.undo(); // undo add B
    // We undid back past the reorder, so just check we're consistent.
    // The reorder undo should have moved id0 back.
    // Let's do a more targeted test: undo only the reorder.
    // Reset.
    let mut tree2 = FeatureTree::new();
    let a = tree2.add_feature(sketch_xy(), "A");
    let b = tree2.add_feature(sketch_xz(), "B");
    // Clear add-A and add-B from log, only record the reorder.
    tree2.undo_log.clear();
    tree2.undo_cursor = 0;

    tree2.reorder(a, 1).expect("reorder should succeed");
    assert_eq!(tree2.timeline(), &[b, a]);

    tree2.undo();
    assert_eq!(
        tree2.timeline(),
        &[a, b],
        "undo reorder should restore original order"
    );

    tree2.redo();
    assert_eq!(
        tree2.timeline(),
        &[b, a],
        "redo reorder should re-apply move"
    );
}
