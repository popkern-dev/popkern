//! Fragments that fill one of their own face's holes.
//!
//! A coplanar partner can bring a region into a holed face's walk that
//! the face's own hole rim closes off: the walk then emits the hole's
//! interior as a fragment of the holed face (measured: a hub standing
//! in a bearing race's bore, their bottoms on one plane; the race's
//! bottom walked a disc over its own bore, which tied radially with
//! the hub's coplanar bottom on the bore rim). Such a fragment bounds
//! none of its face's material. It is recognisable by topology alone:
//! every edge of every loop is a span of the face's own hole rims, and
//! the face's material always reaches its outer loop.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, EdgeHandle};

use super::lift::operand_brep;
use super::spans::SpanRep;
use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// The original operand edge a graph half-edge copies or splits, on
/// `operand`'s side.
fn original_edge(source: &EdgeSource, operand: Operand) -> Option<EdgeHandle> {
    match (source, operand) {
        (EdgeSource::OriginalA(e), Operand::A)
        | (EdgeSource::SplitA { original: e }, Operand::A)
        | (EdgeSource::OriginalB(e), Operand::B)
        | (EdgeSource::SplitB { original: e }, Operand::B) => Some(*e),
        _ => None,
    }
}

/// The edges of `brep`'s face's inner loops.
fn hole_rim_edges(brep: &BRep, face: popkern_core::topology::FaceHandle) -> BTreeSet<EdgeHandle> {
    let mut rims = BTreeSet::new();
    let Some(f) = brep.faces.get(face) else {
        return rims;
    };
    for &lp in &f.inner_loops {
        let Some(start) = brep.loops.get(lp).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        while let Some(h) = brep.half_edges.get(he) {
            rims.insert(h.edge);
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    rims
}

/// Per fragment, whether it fills a hole of its own face.
pub(crate) fn hole_filling_fragments(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    operands: [&BRep; 2],
    reps: &BTreeMap<usize, SpanRep>,
) -> Vec<bool> {
    // Graph keys per span representative.
    let mut members: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (&key, rep) in reps {
        members.entry(rep.key).or_default().push(key);
    }
    fragments
        .iter()
        .map(|frag| {
            let brep = operand_brep(operands, frag.face);
            if brep
                .faces
                .get(frag.face.face)
                .is_none_or(|f| f.inner_loops.is_empty())
            {
                return false;
            }
            let rims = hole_rim_edges(brep, frag.face.face);
            let on_rim = |h: usize| {
                let raw = h.min(graph.half_edges[h].twin);
                let Some(rep) = reps.get(&raw) else {
                    // Degenerate keys bound nothing.
                    return true;
                };
                members.get(&rep.key).is_some_and(|keys| {
                    keys.iter().any(|&k| {
                        original_edge(&graph.half_edges[k].source, frag.face.operand)
                            .is_some_and(|e| rims.contains(&e))
                    })
                })
            };
            !frag.outer_loop.is_empty()
                && frag
                    .outer_loop
                    .iter()
                    .chain(frag.inner_loops.iter().flatten())
                    .all(|&h| on_rim(h))
        })
        .collect()
}

#[cfg(test)]
#[path = "hole_phantoms_tests.rs"]
mod tests;
