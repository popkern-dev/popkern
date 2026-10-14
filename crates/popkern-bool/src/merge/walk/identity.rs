//! Identity fragments: a face no intersection touches walks to its own
//! loops (PK-802).
//!
//! The chart walk over an UNTOUCHED face — every graph half-edge on it
//! is a copy of one of its own original half-edges, none split, none a
//! record piece — can only reproduce the face's stored loops: there is
//! no junction at which it could route differently. Building the chart
//! subdivision, its embeddings and orbits to rediscover them cost
//! ~0.05 ms per face, so a 64-slot disc's late cut spent 12 of its
//! 24 ms merge re-walking the 250 faces the slot never meets.
//!
//! The shortcut is decided by topology alone (no geometry, no
//! tolerance): the face's graph half-edge set must equal, exactly, the
//! images of its original loop half-edges under the operand's
//! half-edge map, every one still `Original*`-sourced and registered to
//! this face. Anything else — a split boundary, a record piece, a seam
//! substitute, a coincident-overlap face, a doubly-periodic surface
//! (whose fragments carry chart views the torus sampler reads) — walks
//! as before.
//!
//! Each loop starts at its lowest graph half-edge index, as the walk
//! seeds it. Inner loops keep the face's STORED order; the walk lists
//! them in discovery order with the outer swap-removed, which rotated
//! an untouched face's holes on every boolean (the same loops — the
//! compare mode checks them as a set). `POPKERN_BOOL_SCOPE=compare` walks every
//! identity face too and reports any difference (see
//! [`compare_identity`]).

use std::collections::HashMap;

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use super::super::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};
use super::Fragment;

/// Original-half-edge → graph-half-edge maps of both operands, as
/// `insert_operand_topology` returned them.
pub(crate) struct OperandHalfEdgeMaps<'a> {
    /// Operand A's map.
    pub a: &'a HashMap<HalfEdgeHandle, usize>,
    /// Operand B's map.
    pub b: &'a HashMap<HalfEdgeHandle, usize>,
}

/// The graph half-edges of one stored loop, in loop order; `None` when
/// any half-edge is unmapped or the loop does not close.
pub(super) fn mapped_loop(
    brep: &BRep,
    loop_h: LoopHandle,
    map: &HashMap<HalfEdgeHandle, usize>,
) -> Option<Vec<usize>> {
    let start = brep.loops.get(loop_h)?.half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(*map.get(&he)?);
        he = brep.half_edges.get(he)?.next;
        if he == start {
            // The walk seeds each loop at its lowest-index half-edge.
            let k = (0..out.len()).min_by_key(|&i| out[i])?;
            out.rotate_left(k);
            return Some(out);
        }
    }
    None
}

/// The identity fragment of `face_ref`, when the face is untouched (see
/// the module docs); `None` sends the face through the walk.
pub(crate) fn identity_fragment(
    graph: &IntersectionGraph,
    face_ref: &FaceRef,
    he_indices: &[usize],
    brep_a: &BRep,
    brep_b: &BRep,
    maps: &OperandHalfEdgeMaps<'_>,
    overlap_pairs: &[(FaceHandle, FaceHandle)],
) -> Option<Fragment> {
    let (brep, map) = match face_ref.operand {
        Operand::A => (brep_a, maps.a),
        Operand::B => (brep_b, maps.b),
    };
    let in_overlap = overlap_pairs
        .iter()
        .any(|&(fa, fb)| match face_ref.operand {
            Operand::A => fa == face_ref.face,
            Operand::B => fb == face_ref.face,
        });
    if in_overlap || graph.seam_substitutes.iter().any(|(f, _)| f == face_ref) {
        return None;
    }
    let face = brep.faces.get(face_ref.face)?;
    if face.surface.as_ref().is_some_and(|s| {
        let (u, v) = s.geometry.is_periodic();
        u && v
    }) {
        return None;
    }
    let outer_loop = mapped_loop(brep, face.outer_loop, map)?;
    let inner_loops = face
        .inner_loops
        .iter()
        .map(|&lh| mapped_loop(brep, lh, map))
        .collect::<Option<Vec<_>>>()?;
    if !std::iter::once(&outer_loop)
        .chain(&inner_loops)
        .all(|lp| closes_in_graph(graph, lp))
    {
        return None;
    }
    let mapped = outer_loop.len() + inner_loops.iter().map(Vec::len).sum::<usize>();
    if mapped != he_indices.len() {
        return None;
    }
    let mut sorted: Vec<usize> = he_indices.to_vec();
    sorted.sort_unstable();
    let original = |h: usize| {
        let he = &graph.half_edges[h];
        he.face.as_ref() == Some(face_ref)
            && matches!(
                he.source,
                EdgeSource::OriginalA(_) | EdgeSource::OriginalB(_)
            )
            && sorted.binary_search(&h).is_ok()
    };
    let all = outer_loop.iter().chain(inner_loops.iter().flatten());
    let mut seen = std::collections::HashSet::with_capacity(mapped);
    for &h in all {
        if !original(h) || !seen.insert(h) {
            return None;
        }
    }
    Some(Fragment {
        face: *face_ref,
        outer_loop,
        inner_loops,
        chart: None,
    })
}

/// Whether `POPKERN_BOOL_SCOPE=compare` asks the walk to check every
/// identity fragment against the full walk.
pub(crate) fn compare_mode() -> bool {
    std::env::var_os("POPKERN_BOOL_SCOPE").is_some_and(|v| v == "compare")
}

/// Compare an identity fragment against what the full walk produced
/// for the same face; a difference is a defect in the identity
/// precondition and panics with both loop sets (compare mode is a
/// test-harness switch, never a serving path).
pub(crate) fn compare_identity(identity: &Fragment, walked: &[&Fragment]) {
    let sorted = |inner: &[Vec<usize>]| {
        let mut inner = inner.to_vec();
        inner.sort();
        inner
    };
    let same = walked.len() == 1
        && walked[0].outer_loop == identity.outer_loop
        && sorted(&walked[0].inner_loops) == sorted(&identity.inner_loops);
    assert!(
        same,
        "identity fragment differs from the walk on {:?}: identity outer {:?} inner {:?}; \
         walk {:?}",
        identity.face,
        identity.outer_loop,
        identity.inner_loops,
        walked
            .iter()
            .map(|f| (&f.outer_loop, &f.inner_loops))
            .collect::<Vec<_>>()
    );
}

/// Consecutive half-edges of `lp` meet head to tail in the graph (no
/// split node interrupts the stored loop).
pub(super) fn closes_in_graph(graph: &IntersectionGraph, lp: &[usize]) -> bool {
    let end = |h: usize| graph.half_edges[graph.half_edges[h].twin].origin;
    lp.iter()
        .zip(lp.iter().cycle().skip(1))
        .all(|(&h, &n)| end(h) == graph.half_edges[n].origin)
}
