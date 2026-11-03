//! Region selection and the fragment selections it implies (PK-975).
//!
//! A boolean is a predicate on region membership (∪: in A or B; ∩: in
//! both; A − B: in A, not B). A complex face is on the result's
//! boundary iff exactly one of its sides bounds a selected region; its
//! output normal points away from the selected side. In the
//! assembler's terms: the face's primary fragment is kept, flipped iff
//! its folded normal would otherwise disagree with its output
//! component's folded-normal convention (so a mixed-sense pair needs
//! no winding vote). A component holding faces of A takes the sense of
//! its first A face's operand component: operand A's faces keep their
//! flags, as the cellular assembler keeps them, even when A's
//! components fold different ways (a cavity, two solids built apart).
//! A component of B faces alone keeps B's: a disjoint tool kept by a
//! union, or a piece meeting A only at a non-manifold edge. The other fragments of a same-domain face
//! are never kept — one skin per coincident patch.
//!
//! Every coincident-face rule the per-fragment classifier spells out
//! case by case (keep A's copy of a cooriented skin, drop an
//! anti-oriented contact in a union, open a flush subtraction's mouth)
//! follows from this one rule applied to the shared face's two sides.

use std::collections::BTreeMap;

use popkern_core::topology::FaceHandle;

use super::ids::{EdgeId, VertexId};
use super::lift::Lift;
use super::model::SideIndex;
use super::operand_set::OperandSet;
use super::regions::Forest;
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
    /// Per complex face: on the result's boundary.
    pub kept_face: Vec<bool>,
    /// Per shell: its region is in the result.
    pub selected_shell: Vec<bool>,
}

/// Per complex face: the folded-normal convention of its output
/// component (kept faces joined through manifold edges): the sense of
/// its A face of least handle, else of its B face of least handle.
fn component_conventions(lift: &Lift, kept: &[bool]) -> Vec<bool> {
    let c = &lift.complex;
    let mut forest = Forest::new(c.faces.len());
    for e in &c.edges {
        let faces: Vec<usize> = e
            .radial
            .iter()
            .map(|f| c.fins[f.index()].face.index())
            .filter(|&fi| kept[fi])
            .collect();
        // Through manifold edges only: the result splits at a
        // non-manifold edge, and the shells either side keep their own
        // operands' conventions.
        if let [a, b] = faces[..] {
            forest.union(a, b);
        }
    }
    // The least (operand, face handle) primary owner of each component
    // decides: an A face where it has one.
    let mut chosen: BTreeMap<usize, ((usize, FaceHandle), bool)> = BTreeMap::new();
    for (fi, face) in c.faces.iter().enumerate() {
        let Some(owner) = face.owners.first().filter(|_| kept[fi]) else {
            continue;
        };
        let key = (owner.operand, owner.face);
        let entry = chosen
            .entry(forest.find(fi))
            .or_insert((key, face.folded_is_outward));
        if key < entry.0 {
            *entry = (key, face.folded_is_outward);
        }
    }
    (0..c.faces.len())
        .map(|fi| {
            chosen
                .get(&forest.find(fi))
                .map_or(lift.output_folded_outward, |&(_, sense)| sense)
        })
        .collect()
}

/// Select the fragments of `lift` for `op`. Shells and membership must
/// be assigned.
pub(crate) fn select_fragments(
    lift: &Lift,
    op: BoolOperation,
    fragment_count: usize,
) -> ComplexSelection {
    select_fragments_by(lift, |m| selects(op, m), fragment_count)
}

/// Select the fragments of `lift` bounding the regions whose membership
/// `keep` accepts (an N-operand partition's set expression; the
/// two-operand booleans are three such predicates).
pub(crate) fn select_fragments_by(
    lift: &Lift,
    keep: impl Fn(&OperandSet) -> bool,
    fragment_count: usize,
) -> ComplexSelection {
    let selected_shell: Vec<bool> = lift
        .complex
        .shells
        .iter()
        .map(|s| s.membership.as_ref().is_some_and(&keep))
        .collect();
    select_shells(lift, selected_shell, fragment_count)
}

/// Select the fragments of `lift` bounding exactly the regions
/// `selected_shell` marks (one region alone: a partition's region
/// solid, which may share its membership with a disjoint region).
pub(crate) fn select_shells(
    lift: &Lift,
    selected_shell: Vec<bool>,
    fragment_count: usize,
) -> ComplexSelection {
    let c = &lift.complex;
    let mut selections = vec![
        FragmentSelection {
            keep: false,
            flip: false,
        };
        fragment_count
    ];
    // Which faces are kept, and which side of each holds the result.
    let mut kept_behind: Vec<Option<bool>> = vec![None; c.faces.len()];
    for (fi, face) in c.faces.iter().enumerate() {
        let side_selected = |side: SideIndex| {
            face.sides[side as usize]
                .shell
                .is_some_and(|s| selected_shell[s.index()])
        };
        let ahead = side_selected(SideIndex::Ahead);
        let behind = side_selected(SideIndex::Behind);
        if ahead != behind {
            kept_behind[fi] = Some(behind);
        }
    }
    let kept_face: Vec<bool> = kept_behind.iter().map(Option::is_some).collect();
    let convention = component_conventions(lift, &kept_face);
    for (fi, face) in c.faces.iter().enumerate() {
        let Some(behind) = kept_behind[fi] else {
            continue;
        };
        // The assembler keeps a fragment's folded-normal sense unless
        // flipped, so every kept face of a component must agree on its
        // convention. The face's result-outward normal is its owner's
        // outward normal when the selected material is behind it, the
        // reverse when ahead.
        let side_sign = if behind { 1.0 } else { -1.0 };
        let face_sense = if face.folded_is_outward { 1.0 } else { -1.0 };
        let output_sense = if convention[fi] { 1.0 } else { -1.0 };
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
        kept_face,
        selected_shell,
    }
}
