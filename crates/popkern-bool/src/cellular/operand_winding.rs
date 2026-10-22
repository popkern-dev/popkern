//! Operand winding agreement before cellular assembly (PK-913).
//!
//! Each operand's fragments wind as its own faces do. Two operands of
//! the same material sense meet along every shared record in opposite
//! directions (each face's loop keeps its own material on the left),
//! and the assembly twins those uses into one manifold edge. A
//! MIXED-sense pair (one shell wound inward, as `BRep::prism` and
//! every blend rebuild are, one outward) meets along every shared
//! record in the SAME direction: the assembly built a non-manifold
//! edge and declined (measured: a disc ∩ inward boss, record
//! half-edges 43, 49, 55 and 61), and the pair served only through the
//! PK-478 rescue, which reruns the whole pipeline with the inward
//! operand reversed.
//!
//! The shared records are the witness: when every shared edge that a
//! kept A fragment and a kept B fragment both use is used in the same
//! direction, B's kept fragments are wound opposite to A's, and
//! toggling their flips makes the pair wire as a same-sense pair. The
//! result keeps A's sense. Any shared edge used in opposite directions
//! (the matched case, or a genuine local conflict) leaves the
//! selections alone, so a same-sense pair is untouched. When B's far
//! faces are carried (a scoped union, PK-816), the assembler reverses
//! them in the same pass (PK-916), so B's sense never splits between
//! its fragments and its far faces.

use std::collections::BTreeMap;

use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// How the kept fragments of A and B use their shared record edges.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct SharedWinding {
    /// Shared edges both operands use in the same direction.
    pub(crate) same: usize,
    /// Shared edges the operands use in opposite directions.
    pub(crate) opposite: usize,
}

impl SharedWinding {
    /// Every jointly used shared edge runs the same way: B is wound
    /// opposite to A.
    pub(crate) fn operands_disagree(&self) -> bool {
        self.same > 0 && self.opposite == 0
    }
}

/// Tally the shared-edge directions of the kept fragments. An edge
/// one operand uses in both directions (two of its kept fragments on
/// either side) is not counted.
pub(crate) fn shared_winding(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    selections: &[FragmentSelection],
) -> SharedWinding {
    // Per shared edge (keyed by its lower half-edge): the direction
    // each operand uses, `Some(true)` = along the key half-edge,
    // `None` once that operand has used both.
    let mut uses: BTreeMap<usize, [Option<Option<bool>>; 2]> = BTreeMap::new();
    for (frag, sel) in fragments.iter().zip(selections) {
        if !sel.keep {
            continue;
        }
        let side = usize::from(frag.face.operand == Operand::B);
        for &raw in frag
            .outer_loop
            .iter()
            .chain(frag.inner_loops.iter().flatten())
        {
            let he = if sel.flip {
                graph.half_edges[raw].twin
            } else {
                raw
            };
            if !matches!(graph.half_edges[he].source, EdgeSource::Shared { .. }) {
                continue;
            }
            let key = he.min(graph.half_edges[he].twin);
            let along = he == key;
            let slot = &mut uses.entry(key).or_default()[side];
            *slot = match *slot {
                None => Some(Some(along)),
                Some(Some(d)) if d == along => Some(Some(d)),
                _ => Some(None),
            };
        }
    }
    let mut tally = SharedWinding::default();
    for [a, b] in uses.into_values() {
        if let (Some(Some(a)), Some(Some(b))) = (a, b) {
            if a == b {
                tally.same += 1;
            } else {
                tally.opposite += 1;
            }
        }
    }
    tally
}

/// Toggle every kept B fragment's flip when the operands' windings
/// disagree on every jointly used shared edge. Returns the tally and
/// whether B was flipped (the assembler then reverses B's carried far
/// faces too).
pub(crate) fn agree_operand_winding(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    selections: &mut [FragmentSelection],
) -> (SharedWinding, bool) {
    let tally = shared_winding(graph, fragments, selections);
    let flip_b = tally.operands_disagree();
    if flip_b {
        for (frag, sel) in fragments.iter().zip(selections.iter_mut()) {
            if sel.keep && frag.face.operand == Operand::B {
                sel.flip = !sel.flip;
            }
        }
    }
    (tally, flip_b)
}

#[cfg(test)]
#[path = "operand_winding_tests.rs"]
mod tests;
