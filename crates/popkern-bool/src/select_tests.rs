//! Tests for Stage 4 fragment selection.
// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use slotmap::SlotMap;

use super::*;
use crate::error::{BoolOperation, Classification};

/// Helper: allocate a fresh FaceHandle backed by a throw-away SlotMap.
fn make_face() -> (SlotMap<FaceHandle, ()>, FaceHandle) {
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    (sm, h)
}

// ------------------------------------------------------------------
// Union -- 8 rows
// ------------------------------------------------------------------

#[test]
fn union_a_outside_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::A, Classification::Outside)],
    );
    assert!(result.is_kept(face, Operand::A));
    assert!(!result.is_flipped(face, Operand::A));
}

#[test]
fn union_a_inside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::A, Classification::Inside)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn union_a_cooriented_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::A, Classification::OnBoundaryCooriented)],
    );
    assert!(result.is_kept(face, Operand::A));
    assert!(!result.is_flipped(face, Operand::A));
}

#[test]
fn union_a_crossing_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::A, Classification::OnBoundaryCrossing)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn union_b_outside_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::B, Classification::Outside)],
    );
    assert!(result.is_kept(face, Operand::B));
    assert!(!result.is_flipped(face, Operand::B));
}

#[test]
fn union_b_inside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::B, Classification::Inside)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

#[test]
fn union_b_cooriented_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::B, Classification::OnBoundaryCooriented)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

#[test]
fn union_b_crossing_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Union,
        &[(face, Operand::B, Classification::OnBoundaryCrossing)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

// ------------------------------------------------------------------
// Intersection -- 8 rows
// ------------------------------------------------------------------

#[test]
fn intersection_a_outside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::A, Classification::Outside)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn intersection_a_inside_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::A, Classification::Inside)],
    );
    assert!(result.is_kept(face, Operand::A));
    assert!(!result.is_flipped(face, Operand::A));
}

#[test]
fn intersection_a_cooriented_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::A, Classification::OnBoundaryCooriented)],
    );
    assert!(result.is_kept(face, Operand::A));
    assert!(!result.is_flipped(face, Operand::A));
}

#[test]
fn intersection_a_crossing_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::A, Classification::OnBoundaryCrossing)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn intersection_b_outside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::B, Classification::Outside)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

#[test]
fn intersection_b_inside_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::B, Classification::Inside)],
    );
    assert!(result.is_kept(face, Operand::B));
    assert!(!result.is_flipped(face, Operand::B));
}

#[test]
fn intersection_b_cooriented_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::B, Classification::OnBoundaryCooriented)],
    );
    assert!(result.is_kept(face, Operand::B));
    assert!(!result.is_flipped(face, Operand::B));
}

#[test]
fn intersection_b_crossing_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Intersection,
        &[(face, Operand::B, Classification::OnBoundaryCrossing)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

// ------------------------------------------------------------------
// Subtraction -- 8 rows
// ------------------------------------------------------------------

#[test]
fn subtraction_a_outside_keep_no_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::A, Classification::Outside)],
    );
    assert!(result.is_kept(face, Operand::A));
    assert!(!result.is_flipped(face, Operand::A));
}

#[test]
fn subtraction_a_inside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::A, Classification::Inside)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn subtraction_a_cooriented_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::A, Classification::OnBoundaryCooriented)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn subtraction_a_crossing_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::A, Classification::OnBoundaryCrossing)],
    );
    assert!(!result.is_kept(face, Operand::A));
}

#[test]
fn subtraction_b_outside_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::B, Classification::Outside)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

#[test]
fn subtraction_b_inside_keep_and_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::B, Classification::Inside)],
    );
    assert!(result.is_kept(face, Operand::B));
    assert!(result.is_flipped(face, Operand::B));
}

#[test]
fn subtraction_b_cooriented_discard() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::B, Classification::OnBoundaryCooriented)],
    );
    assert!(!result.is_kept(face, Operand::B));
}

#[test]
fn subtraction_b_crossing_keep_and_flip() {
    let (_, face) = make_face();
    let result = select_faces(
        BoolOperation::Subtraction,
        &[(face, Operand::B, Classification::OnBoundaryCrossing)],
    );
    assert!(result.is_kept(face, Operand::B));
    assert!(result.is_flipped(face, Operand::B));
}

// ------------------------------------------------------------------
// Edge cases
// ------------------------------------------------------------------

#[test]
fn empty_input_gives_empty_output() {
    let result = select_faces(BoolOperation::Union, &[]);
    assert!(result.keep_a.is_empty());
    assert!(result.keep_b.is_empty());
    assert!(result.flip_a.is_empty());
    assert!(result.flip_b.is_empty());
}

#[test]
fn mixed_fragments_keep_discard_flip() {
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let f_keep = sm.insert(());
    let f_discard = sm.insert(());
    let f_flip = sm.insert(());

    let classifications = [
        // kept, no flip
        (f_keep, Operand::A, Classification::Outside),
        // discarded
        (f_discard, Operand::A, Classification::Inside),
        // kept + flipped
        (f_flip, Operand::B, Classification::Inside),
    ];

    let result = select_faces(BoolOperation::Subtraction, &classifications);

    assert!(result.is_kept(f_keep, Operand::A));
    assert!(!result.is_flipped(f_keep, Operand::A));

    assert!(!result.is_kept(f_discard, Operand::A));
    assert!(!result.is_flipped(f_discard, Operand::A));

    assert!(result.is_kept(f_flip, Operand::B));
    assert!(result.is_flipped(f_flip, Operand::B));
}

#[test]
fn flip_is_subset_of_keep() {
    // Generate all 24 combinations across all three operations and verify
    // that no face appears in flip without also appearing in keep.
    let ops = [
        BoolOperation::Union,
        BoolOperation::Intersection,
        BoolOperation::Subtraction,
    ];
    let operands = [Operand::A, Operand::B];
    let classes = [
        Classification::Outside,
        Classification::Inside,
        Classification::OnBoundaryCooriented,
        Classification::OnBoundaryCrossing,
    ];

    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();

    for op in ops {
        for operand in operands {
            for class in classes {
                let face = sm.insert(());
                let result = select_faces(op, &[(face, operand, class)]);
                if result.is_flipped(face, operand) {
                    assert!(
                        result.is_kept(face, operand),
                        "face in flip but not keep for {op:?} {operand:?} {class:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn handle_collision_between_operands() {
    // Verify that identical FaceHandle values from different operands
    // do not collide in the selection result.
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let shared_handle = sm.insert(());

    // A:Outside -> keep, B:Inside -> discard (for Union)
    let classifications = [
        (shared_handle, Operand::A, Classification::Outside),
        (shared_handle, Operand::B, Classification::Inside),
    ];

    let result = select_faces(BoolOperation::Union, &classifications);

    assert!(
        result.is_kept(shared_handle, Operand::A),
        "A:Outside should be kept in Union"
    );
    assert!(
        !result.is_kept(shared_handle, Operand::B),
        "B:Inside should NOT be kept in Union"
    );
}
