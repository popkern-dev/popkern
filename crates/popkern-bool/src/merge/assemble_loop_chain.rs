//! Loop half-edge re-sequencing for result assembly — split from
//! `assemble.rs` (500-line cap ratchet). The manifold-assembly repair
//! for multi-edge shared loops lives here.

use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Reorder resolved BRep half-edges into a single closed directed chain
/// (each half-edge's end vertex is the next's start vertex).
///
/// Returns `Some(reordered)` only when the whole handle set links into
/// exactly one closed cycle; `None` otherwise (a caller keeps its
/// original order and declines downstream). This is the manifold-assembly
/// repair for multi-edge shared loops: when a fragment's primary
/// half-edge on a shared edge was already claimed by the adjacent
/// (oppositely-oriented) fragment, `resolve_loop_hes` falls back to the
/// twin half-edge — the geometrically correct opposite direction, but
/// emitted in the walker's original sequence, which no longer chains
/// head-to-tail. The resolved twins still form one closed cycle (they are
/// the same edges, reversed); re-sequencing by vertex adjacency recovers
/// the loop without fabricating or re-choosing any edge. A plain circular
/// hole never needs this (its rim is a single self-loop, order-immune);
/// only order-sensitive multi-edge shared loops (a keyed bore rim over a
/// union-derived cap) do.
pub(super) fn reorder_into_closed_chain(
    brep: &BRep,
    handles: &[HalfEdgeHandle],
) -> Option<Vec<HalfEdgeHandle>> {
    if handles.len() < 2 {
        return None;
    }
    let start_of = |h: HalfEdgeHandle| brep.half_edges[h].origin;
    let end_of = |h: HalfEdgeHandle| brep.half_edges[brep.half_edges[h].twin].origin;
    // Leave an already-connected chain untouched (byte-identical wiring
    // for the currently-passing corpus): repair activates only on loops
    // whose resolved handles do not link head-to-tail.
    let n = handles.len();
    if (0..n).all(|i| end_of(handles[i]) == start_of(handles[(i + 1) % n])) {
        return None;
    }
    // A self-loop half-edge (start == end) is its own mini-cycle and
    // cannot be sequenced against neighbours by vertex adjacency; leave
    // any loop containing one to its original order.
    if handles.iter().any(|&h| start_of(h) == end_of(h)) {
        return None;
    }
    let mut remaining: Vec<HalfEdgeHandle> = handles.to_vec();
    let first = remaining.remove(0);
    let chain_start = start_of(first);
    let mut chain = vec![first];
    let mut tip = end_of(first);
    while !remaining.is_empty() {
        // Exactly one remaining half-edge should start at the tip; an
        // ambiguous junction (two candidates) means the handle set is
        // not a simple cycle — decline the repair.
        let matches: Vec<usize> = remaining
            .iter()
            .enumerate()
            .filter(|(_, h)| start_of(**h) == tip)
            .map(|(i, _)| i)
            .collect();
        if matches.len() != 1 {
            trace_reorder_decline(brep, handles, "ambiguous-or-missing-successor", tip);
            return None;
        }
        let next = remaining.remove(matches[0]);
        tip = end_of(next);
        chain.push(next);
    }
    if tip == chain_start {
        Some(chain)
    } else {
        trace_reorder_decline(brep, handles, "open-chain", tip);
        None
    }
}

/// Why a loop's resolved handles could not be re-sequenced — the
/// invalid-loop class is only diagnosable here, where the set's
/// orientation mix is still visible.
fn trace_reorder_decline(
    brep: &BRep,
    handles: &[HalfEdgeHandle],
    reason: &str,
    tip: popkern_core::topology::VertexHandle,
) {
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [loop-chain] decline ({reason}) at tip {tip:?}; {} handles:",
            handles.len()
        );
        for &h in handles {
            let s = brep.half_edges[h].origin;
            let e = brep.half_edges[brep.half_edges[h].twin].origin;
            let sp = brep.vertices.get(s).map(|v| v.position);
            let ep = brep.vertices.get(e).map(|v| v.position);
            if let (Some(sp), Some(ep)) = (sp, ep) {
                eprintln!(
                    "    he {h:?} {s:?}({:.3},{:.3},{:.3}) -> {e:?}({:.3},{:.3},{:.3})",
                    sp.x, sp.y, sp.z, ep.x, ep.y, ep.z
                );
            }
        }
    }
}
