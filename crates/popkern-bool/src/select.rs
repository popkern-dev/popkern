//! Stage 4: Select face fragments by operation type.
//!
//! Given an operation type and a list of `(face, operand, classification)` tuples,
//! decides which fragments to keep in the result solid and which need their normals flipped.
//! This stage is pure table-driven logic — no geometry, no external dependencies beyond types.

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
///
/// Uses separate sets for each operand to avoid FaceHandle collisions
/// between different BReps (where handles from A and B may share the
/// same internal value).
#[derive(Debug)]
pub(crate) struct SelectionResult {
    /// Face handles from operand A to keep in the result.
    pub keep_a: HashSet<FaceHandle>,
    /// Face handles from operand B to keep in the result.
    pub keep_b: HashSet<FaceHandle>,
    /// Face handles from operand A that need their normals flipped.
    pub flip_a: HashSet<FaceHandle>,
    /// Face handles from operand B that need their normals flipped.
    pub flip_b: HashSet<FaceHandle>,
}

impl SelectionResult {
    /// Check if a face from operand A or B is kept, based on the operand.
    #[cfg(test)]
    pub(crate) fn is_kept(&self, face: FaceHandle, operand: Operand) -> bool {
        match operand {
            Operand::A => self.keep_a.contains(&face),
            Operand::B => self.keep_b.contains(&face),
        }
    }

    /// Check if a face from operand A or B should be flipped, based on the operand.
    #[cfg(test)]
    pub(crate) fn is_flipped(&self, face: FaceHandle, operand: Operand) -> bool {
        match operand {
            Operand::A => self.flip_a.contains(&face),
            Operand::B => self.flip_b.contains(&face),
        }
    }
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
    let mut keep_a = HashSet::new();
    let mut keep_b = HashSet::new();
    let mut flip_a = HashSet::new();
    let mut flip_b = HashSet::new();

    for &(face, operand, class) in classifications {
        let (do_keep, do_flip) = should_keep(op, operand, class);
        if do_keep {
            match operand {
                Operand::A => {
                    keep_a.insert(face);
                    if do_flip {
                        flip_a.insert(face);
                    }
                }
                Operand::B => {
                    keep_b.insert(face);
                    if do_flip {
                        flip_b.insert(face);
                    }
                }
            }
        }
    }

    SelectionResult {
        keep_a,
        keep_b,
        flip_a,
        flip_b,
    }
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
#[path = "select_tests.rs"]
mod tests;
