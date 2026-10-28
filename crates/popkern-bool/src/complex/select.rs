//! Region selection and the fragment selections it implies (PK-975).
//!
//! A boolean is a predicate on region membership (∪: in A or B; ∩: in
//! both; A − B: in A, not B). A complex face is on the result's
//! boundary iff exactly one of its sides bounds a selected region; its
//! output normal points away from the selected side. In the
//! assembler's terms: the face's primary fragment is kept, flipped iff
//! its folded normal would otherwise disagree with operand A's
//! folded-normal convention for the result (so a mixed-sense pair
//! needs no winding vote). The other fragments of a same-domain face
//! are never kept — one skin per coincident patch.
//!
//! Every coincident-face rule the per-fragment classifier spells out
//! case by case (keep A's copy of a cooriented skin, drop an
//! anti-oriented contact in a union, open a flush subtraction's mouth)
//! follows from this one rule applied to the shared face's two sides.

use super::ids::{EdgeId, VertexId};
use super::lift::Lift;
use super::model::SideIndex;
use super::operand_set::OperandSet;
use crate::error::BoolOperation;
use crate::merge::fragment_classify::FragmentSelection;

/// Whether a region with `membership` is in the result of `op` over
/// operands A (0) and B (1).
pub(crate) fn selects(op: BoolOperation, membership: &OperandSet) -> bool {
    let (a, b) = (membership.contains(0), membership.contains(1));
    match op {
        BoolOperation::Union => a || b,
        BoolOperation::Intersection => a && b,
        BoolOperation::Subtraction => a && !b,
    }
}

/// The complex's verdict for one boolean.
#[derive(Debug, Clone)]
pub(crate) struct ComplexSelection {
    /// Keep/flip per fragment, in fragment order.
    pub selections: Vec<FragmentSelection>,
    /// Edges more than two kept fins use: the result is non-manifold
    /// there.
    pub nonmanifold_edges: Vec<EdgeId>,
    /// Vertices where the kept faces form more than one fan.
    pub nonmanifold_vertices: Vec<VertexId>,
}

/// Select the fragments of `lift` for `op`. Shells and membership must
/// be assigned.
pub(crate) fn select_fragments(
    lift: &Lift,
    op: BoolOperation,
    fragment_count: usize,
) -> ComplexSelection {
    let c = &lift.complex;
    let selected_shell: Vec<bool> = c
        .shells
        .iter()
        .map(|s| s.membership.as_ref().is_some_and(|m| selects(op, m)))
        .collect();
    let mut selections = vec![
        FragmentSelection {
            keep: false,
            flip: false,
        };
        fragment_count
    ];
    let mut kept_face = vec![false; c.faces.len()];
    for (fi, face) in c.faces.iter().enumerate() {
        let side_selected = |side: SideIndex| {
            face.sides[side as usize]
                .shell
                .is_some_and(|s| selected_shell[s.index()])
        };
        let ahead = side_selected(SideIndex::Ahead);
        let behind = side_selected(SideIndex::Behind);
        if ahead == behind {
            continue;
        }
        kept_face[fi] = true;
        // The assembler keeps a fragment's folded-normal sense unless
        // flipped, so every kept face must agree on ONE convention:
        // operand A's (its folded normals point out of material iff
        // `output_folded_outward`). The face's result-outward normal
        // is its owner's outward normal when the selected material is
        // behind it, the reverse when ahead.
        let side_sign = if behind { 1.0 } else { -1.0 };
        let face_sense = if face.folded_is_outward { 1.0 } else { -1.0 };
        let output_sense = if lift.output_folded_outward {
            1.0
        } else {
            -1.0
        };
        if let Some(&primary) = face.source_fragments.first()
            && let Some(sel) = selections.get_mut(primary)
        {
            *sel = FragmentSelection {
                keep: true,
                flip: output_sense * side_sign * face_sense < 0.0,
            };
        }
    }
    let nonmanifold_edges = c
        .edges
        .iter()
        .enumerate()
        .filter(|(_, e)| {
            e.radial
                .iter()
                .filter(|f| kept_face[c.fins[f.index()].face.index()])
                .count()
                > 2
        })
        .map(|(i, _)| EdgeId::from_index(i))
        .collect();
    let nonmanifold_vertices = super::fans::split_vertices(c, &kept_face);
    ComplexSelection {
        selections,
        nonmanifold_edges,
        nonmanifold_vertices,
    }
}
