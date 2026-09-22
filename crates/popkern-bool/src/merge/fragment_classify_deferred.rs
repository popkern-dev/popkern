//! §6.4 deferred ON-resolution overrides — split from
//! `fragment_classify.rs` under the 500-line budget.

use popkern_core::topology::BRep;

use super::super::graph::Operand;
use super::super::walk::Fragment;
use super::{FragmentSelection, apply_selection_table, classification_to_point};
use crate::coincidence::CoincidenceMap;
use crate::error::{BoolOperation, Classification};

/// Walk fragments whose face is in an ambiguous-on pair and apply
/// `disambiguate_with_neighbors` to refine the classification.
pub(super) fn apply_deferred_on_overrides(
    fragments: &[Fragment],
    selections: &mut [FragmentSelection],
    brep_a: &BRep,
    brep_b: &BRep,
    op: BoolOperation,
    coincidence: &CoincidenceMap,
) {
    use crate::classify_on_deferred::find_ambiguous_pairs;

    let pairs = find_ambiguous_pairs(coincidence);
    if pairs.is_empty() {
        return;
    }

    // Build the map (face_handle, classification) for both operands by
    // reading back the existing selections. The classification is
    // synthesized from the FragmentSelection: keep=true → Inside,
    // keep=false → Outside (for the disambiguator's vote logic).
    let mut classes_a: Vec<(popkern_core::topology::FaceHandle, Classification)> = Vec::new();
    let mut classes_b: Vec<(popkern_core::topology::FaceHandle, Classification)> = Vec::new();
    for (frag, sel) in fragments.iter().zip(selections.iter()) {
        let class = if sel.keep {
            Classification::Inside
        } else {
            Classification::Outside
        };
        match frag.face.operand {
            Operand::A => classes_a.push((frag.face.face, class)),
            Operand::B => classes_b.push((frag.face.face, class)),
        }
    }

    // For each ambiguous pair, query the disambiguator. The pair
    // contains face_a (operand A) and face_b (operand B).
    for pair in &pairs {
        // Operand A side.
        if let Some((_, original)) = classes_a.iter().find(|(f, _)| *f == pair.face_a) {
            let suggested = crate::classify_on_deferred::disambiguate_with_neighbors(
                pair.face_a,
                *original,
                brep_a,
                &classes_a,
            );
            if suggested != *original {
                let point_class = classification_to_point(suggested);
                for (i, frag) in fragments.iter().enumerate() {
                    if frag.face.operand == Operand::A && frag.face.face == pair.face_a {
                        selections[i] = apply_selection_table(op, Operand::A, point_class);
                    }
                }
            }
        }

        // Operand B side.
        if let Some((_, original)) = classes_b.iter().find(|(f, _)| *f == pair.face_b) {
            let suggested = crate::classify_on_deferred::disambiguate_with_neighbors(
                pair.face_b,
                *original,
                brep_b,
                &classes_b,
            );
            if suggested != *original {
                let point_class = classification_to_point(suggested);
                for (i, frag) in fragments.iter().enumerate() {
                    if frag.face.operand == Operand::B && frag.face.face == pair.face_b {
                        selections[i] = apply_selection_table(op, Operand::B, point_class);
                    }
                }
            }
        }
    }
}
