//! ON-face classification rules for coincident geometry.
//!
//! When faces from bodies A and B are coincident (ON), the boolean result
//! depends on the normal orientation relationship and the operation type.
//!
//! | Operation     | Same-oriented ON   | Opposite-oriented ON |
//! |---------------|--------------------|-----------------------|
//! | Union         | Keep one copy      | Remove both (cancel)  |
//! | Intersection  | Keep one copy      | Remove both (cancel)  |
//! | Subtraction   | Remove both        | Keep (from A)         |

use popkern_core::topology::FaceHandle;

use crate::coincidence::types::{CoincidenceMap, CoincidenceOrientation};
use crate::error::{BoolOperation, Classification};

/// Determine the classification of a face that lies ON the boundary of
/// the opposite body, given the operation and orientation.
///
/// Returns `None` if the face is not in the coincidence map (not ON).
pub fn classify_on_face(
    face: FaceHandle,
    operand_is_a: bool,
    operation: BoolOperation,
    coincidence: &CoincidenceMap,
) -> Option<Classification> {
    // Find this face in the coincidence map.
    // Only apply coincidence rules for full overlaps. For partial overlaps,
    // the face has been split into coincident and non-coincident fragments.
    // Ray-casting classification handles partial overlaps correctly at the
    // fragment level.
    let overlap = coincidence.face_overlaps.iter().find(|fo| {
        let matches = if operand_is_a {
            fo.face_a == face
        } else {
            fo.face_b == face
        };
        matches && fo.full_overlap
    })?;

    let orientation = overlap.orientation;

    Some(match (operation, orientation) {
        // Union: same-oriented → keep one (from A); opposite → cancel both
        (BoolOperation::Union, CoincidenceOrientation::Same) => {
            if operand_is_a {
                Classification::OnBoundaryCooriented
            } else {
                // B's copy is redundant — classify as Inside to exclude
                Classification::Inside
            }
        }
        (BoolOperation::Union, CoincidenceOrientation::Opposite) => {
            // Both cancel — internal sheet, remove
            Classification::Inside
        }

        // Intersection: same-oriented → keep one (from A); opposite → cancel
        (BoolOperation::Intersection, CoincidenceOrientation::Same) => {
            if operand_is_a {
                Classification::OnBoundaryCooriented
            } else {
                Classification::Outside
            }
        }
        (BoolOperation::Intersection, CoincidenceOrientation::Opposite) => {
            // No common material at overlap
            Classification::Outside
        }

        // Subtraction: same-oriented → cancel; opposite → keep from A
        (BoolOperation::Subtraction, CoincidenceOrientation::Same) => Classification::Inside,
        (BoolOperation::Subtraction, CoincidenceOrientation::Opposite) => {
            if operand_is_a {
                Classification::OnBoundaryCrossing
            } else {
                Classification::Inside
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coincidence::types::{CoincidenceOrientation, FaceOverlap};
    use popkern_core::geometry::Point3;

    fn make_map(
        face_a: FaceHandle,
        face_b: FaceHandle,
        orientation: CoincidenceOrientation,
    ) -> CoincidenceMap {
        let mut map = CoincidenceMap::default();
        map.face_overlaps.push(FaceOverlap {
            face_a,
            face_b,
            orientation,
            full_overlap: true,
            overlap_center: Point3::new(0.0, 0.0, 0.0),
        });
        map
    }

    #[test]
    fn union_same_oriented_keeps_a() {
        let fa = FaceHandle::default();
        let fb = FaceHandle::default();
        let map = make_map(fa, fb, CoincidenceOrientation::Same);

        let class_a = classify_on_face(fa, true, BoolOperation::Union, &map);
        assert_eq!(class_a, Some(Classification::OnBoundaryCooriented));

        // B's copy should be excluded
        let class_b = classify_on_face(fb, false, BoolOperation::Union, &map);
        assert_eq!(class_b, Some(Classification::Inside));
    }

    #[test]
    fn union_opposite_oriented_cancels() {
        let fa = FaceHandle::default();
        let fb = FaceHandle::default();
        let map = make_map(fa, fb, CoincidenceOrientation::Opposite);

        let class = classify_on_face(fa, true, BoolOperation::Union, &map);
        assert_eq!(class, Some(Classification::Inside));
    }

    #[test]
    fn subtraction_same_cancels() {
        let fa = FaceHandle::default();
        let fb = FaceHandle::default();
        let map = make_map(fa, fb, CoincidenceOrientation::Same);

        let class = classify_on_face(fa, true, BoolOperation::Subtraction, &map);
        assert_eq!(class, Some(Classification::Inside));
    }

    #[test]
    fn subtraction_opposite_keeps_a() {
        let fa = FaceHandle::default();
        let fb = FaceHandle::default();
        let map = make_map(fa, fb, CoincidenceOrientation::Opposite);

        let class = classify_on_face(fa, true, BoolOperation::Subtraction, &map);
        assert_eq!(class, Some(Classification::OnBoundaryCrossing));
    }

    #[test]
    fn non_coincident_face_returns_none() {
        let map = CoincidenceMap::default();
        let result = classify_on_face(FaceHandle::default(), true, BoolOperation::Union, &map);
        assert!(result.is_none());
    }
}
