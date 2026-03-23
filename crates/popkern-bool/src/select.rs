//! Stage 4: Select face fragments by operation type.
//!
//! Given an operation type and a list of `(face, operand, classification)` tuples,
//! decides which fragments to keep in the result solid and which need their normals flipped.
//! This stage is pure table-driven logic — no geometry, no external dependencies beyond types.

// Items in this module are consumed by the pipeline orchestrator (not yet written).
#![allow(dead_code)]

use std::collections::HashSet;

use popkern_core::topology::FaceHandle;

use crate::error::{BoolOperation, Classification};

/// Which operand a face fragment belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Operand {
    /// Fragment from solid A.
    A,
    /// Fragment from solid B.
    B,
}

/// Result of fragment selection.
#[derive(Debug)]
pub(crate) struct SelectionResult {
    /// Face handles to keep in the result.
    pub keep: HashSet<FaceHandle>,
    /// Face handles that need their normals flipped.
    pub flip: HashSet<FaceHandle>,
}

/// Determine which face fragments to keep based on the operation type.
///
/// Each entry in `classifications` is a `(face, operand, class)` triple describing
/// a face fragment, which solid it came from, and how it was classified relative
/// to the opposite solid.  The returned [`SelectionResult`] contains:
///
/// - `keep`: all faces that belong in the result solid.
/// - `flip`: the subset of kept faces whose surface normals must be reversed before stitching.
pub(crate) fn select_faces(
    op: BoolOperation,
    classifications: &[(FaceHandle, Operand, Classification)],
) -> SelectionResult {
    let mut keep = HashSet::new();
    let mut flip = HashSet::new();

    for &(face, operand, class) in classifications {
        let (do_keep, do_flip) = should_keep(op, operand, class);
        if do_keep {
            keep.insert(face);
            if do_flip {
                flip.insert(face);
            }
        }
    }

    SelectionResult { keep, flip }
}

/// Consult the selection table from the design spec.
///
/// Returns `(keep, flip)` for a single fragment.
///
/// # Selection table
///
/// | Op           | Operand | Classification       | keep  | flip  |
/// |--------------|---------|----------------------|-------|-------|
/// | Union        | A       | Outside              | true  | false |
/// | Union        | A       | Inside               | false | —     |
/// | Union        | A       | OnBoundaryCooriented | true  | false |
/// | Union        | A       | OnBoundaryCrossing   | false | —     |
/// | Union        | B       | Outside              | true  | false |
/// | Union        | B       | Inside               | false | —     |
/// | Union        | B       | OnBoundaryCooriented | false | —     |
/// | Union        | B       | OnBoundaryCrossing   | false | —     |
/// | Intersection | A       | Outside              | false | —     |
/// | Intersection | A       | Inside               | true  | false |
/// | Intersection | A       | OnBoundaryCooriented | true  | false |
/// | Intersection | A       | OnBoundaryCrossing   | false | —     |
/// | Intersection | B       | Outside              | false | —     |
/// | Intersection | B       | Inside               | true  | false |
/// | Intersection | B       | OnBoundaryCooriented | true  | false |
/// | Intersection | B       | OnBoundaryCrossing   | false | —     |
/// | Subtraction  | A       | Outside              | true  | false |
/// | Subtraction  | A       | Inside               | false | —     |
/// | Subtraction  | A       | OnBoundaryCooriented | false | —     |
/// | Subtraction  | A       | OnBoundaryCrossing   | false | —     |
/// | Subtraction  | B       | Outside              | false | —     |
/// | Subtraction  | B       | Inside               | true  | true  |
/// | Subtraction  | B       | OnBoundaryCooriented | false | —     |
/// | Subtraction  | B       | OnBoundaryCrossing   | true  | true  |
fn should_keep(op: BoolOperation, operand: Operand, class: Classification) -> (bool, bool) {
    use BoolOperation::*;
    use Classification::*;
    use Operand::*;

    match (op, operand, class) {
        // Union ---------------------------------------------------------------
        (Union, A, Outside) => (true, false),
        (Union, A, Inside) => (false, false),
        (Union, A, OnBoundaryCooriented) => (true, false),
        (Union, A, OnBoundaryCrossing) => (false, false),
        (Union, B, Outside) => (true, false),
        (Union, B, Inside) => (false, false),
        (Union, B, OnBoundaryCooriented) => (false, false),
        (Union, B, OnBoundaryCrossing) => (false, false),

        // Intersection --------------------------------------------------------
        (Intersection, A, Outside) => (false, false),
        (Intersection, A, Inside) => (true, false),
        (Intersection, A, OnBoundaryCooriented) => (true, false),
        (Intersection, A, OnBoundaryCrossing) => (false, false),
        (Intersection, B, Outside) => (false, false),
        (Intersection, B, Inside) => (true, false),
        (Intersection, B, OnBoundaryCooriented) => (true, false),
        (Intersection, B, OnBoundaryCrossing) => (false, false),

        // Subtraction (A − B) ------------------------------------------------
        (Subtraction, A, Outside) => (true, false),
        (Subtraction, A, Inside) => (false, false),
        (Subtraction, A, OnBoundaryCooriented) => (false, false),
        (Subtraction, A, OnBoundaryCrossing) => (false, false),
        (Subtraction, B, Outside) => (false, false),
        (Subtraction, B, Inside) => (true, true),
        (Subtraction, B, OnBoundaryCooriented) => (false, false),
        (Subtraction, B, OnBoundaryCrossing) => (true, true),
    }
}

#[cfg(test)]
mod tests {
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
    // Union — 8 rows
    // ------------------------------------------------------------------

    #[test]
    fn union_a_outside_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::A, Classification::Outside)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn union_a_inside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::A, Classification::Inside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn union_a_cooriented_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::A, Classification::OnBoundaryCooriented)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn union_a_crossing_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::A, Classification::OnBoundaryCrossing)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn union_b_outside_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::B, Classification::Outside)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn union_b_inside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::B, Classification::Inside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn union_b_cooriented_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::B, Classification::OnBoundaryCooriented)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn union_b_crossing_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Union,
            &[(face, Operand::B, Classification::OnBoundaryCrossing)],
        );
        assert!(!result.keep.contains(&face));
    }

    // ------------------------------------------------------------------
    // Intersection — 8 rows
    // ------------------------------------------------------------------

    #[test]
    fn intersection_a_outside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::A, Classification::Outside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn intersection_a_inside_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::A, Classification::Inside)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn intersection_a_cooriented_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::A, Classification::OnBoundaryCooriented)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn intersection_a_crossing_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::A, Classification::OnBoundaryCrossing)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn intersection_b_outside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::B, Classification::Outside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn intersection_b_inside_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::B, Classification::Inside)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn intersection_b_cooriented_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::B, Classification::OnBoundaryCooriented)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn intersection_b_crossing_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Intersection,
            &[(face, Operand::B, Classification::OnBoundaryCrossing)],
        );
        assert!(!result.keep.contains(&face));
    }

    // ------------------------------------------------------------------
    // Subtraction — 8 rows
    // ------------------------------------------------------------------

    #[test]
    fn subtraction_a_outside_keep_no_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::A, Classification::Outside)],
        );
        assert!(result.keep.contains(&face));
        assert!(!result.flip.contains(&face));
    }

    #[test]
    fn subtraction_a_inside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::A, Classification::Inside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn subtraction_a_cooriented_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::A, Classification::OnBoundaryCooriented)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn subtraction_a_crossing_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::A, Classification::OnBoundaryCrossing)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn subtraction_b_outside_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::B, Classification::Outside)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn subtraction_b_inside_keep_and_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::B, Classification::Inside)],
        );
        assert!(result.keep.contains(&face));
        assert!(result.flip.contains(&face));
    }

    #[test]
    fn subtraction_b_cooriented_discard() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::B, Classification::OnBoundaryCooriented)],
        );
        assert!(!result.keep.contains(&face));
    }

    #[test]
    fn subtraction_b_crossing_keep_and_flip() {
        let (_, face) = make_face();
        let result = select_faces(
            BoolOperation::Subtraction,
            &[(face, Operand::B, Classification::OnBoundaryCrossing)],
        );
        assert!(result.keep.contains(&face));
        assert!(result.flip.contains(&face));
    }

    // ------------------------------------------------------------------
    // Edge cases
    // ------------------------------------------------------------------

    #[test]
    fn empty_input_gives_empty_output() {
        let result = select_faces(BoolOperation::Union, &[]);
        assert!(result.keep.is_empty());
        assert!(result.flip.is_empty());
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

        assert!(result.keep.contains(&f_keep));
        assert!(!result.flip.contains(&f_keep));

        assert!(!result.keep.contains(&f_discard));
        assert!(!result.flip.contains(&f_discard));

        assert!(result.keep.contains(&f_flip));
        assert!(result.flip.contains(&f_flip));
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
                    if result.flip.contains(&face) {
                        assert!(
                            result.keep.contains(&face),
                            "face in flip but not keep for {op:?} {operand:?} {class:?}"
                        );
                    }
                }
            }
        }
    }
}
