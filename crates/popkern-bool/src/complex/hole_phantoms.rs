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
//!
//! A rim-only fragment can also be the hole itself, filed apart by the
//! walk: it then winds against the material, and rejoins the fragment
//! around it ([`rejoin_own_holes`], applied after the walk for every
//! classifier).

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, EdgeHandle};

use super::lift::operand_brep;
use super::spans::SpanRep;
use crate::merge::graph::{EdgeSource, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// The original operand edge a graph half-edge copies or splits, on
/// `operand`'s side.
fn original_edge(source: &EdgeSource, operand: Operand) -> Option<EdgeHandle> {
    match *source {
        EdgeSource::Original(o, e)
        | EdgeSource::Split {
            operand: o,
            original: e,
        } if o == operand => Some(e),
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
    operands: &[&BRep],
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

/// Per rim-only fragment, the fragment whose hole it is, if any.
///
/// The walk can file a face's own hole ring as a separate piece: where
/// a coincident interface splits the face, a ring lying outside the
/// elected outer loop is never offered to hole nesting (measured: a
/// print-in-place hinge box's back face, split by the swung lid, lost
/// its three knuckle footprints, and the knuckle walls' rims were left
/// free). Such a ring walks its rim AGAINST the face's material, where
/// a filler walks it with the material; its host is the smallest
/// fragment of the same face winding the other way that strictly
/// contains it. Containment is the walk's own positive proof, read in
/// the plane facing the host's Newell normal: a ring it cannot place
/// (a ring around a full turn of a curved face) keeps no host, and is
/// left out as before (the complex then declines, recorded, instead of
/// guessing).
pub(crate) fn own_hole_hosts(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    rim_only: &[bool],
) -> Vec<Option<usize>> {
    use crate::merge::walk::hole_nesting::loop_strictly_inside_outer;
    use crate::merge::walk::walk_geometry::loop_newell;
    let newell: Vec<_> = fragments
        .iter()
        .map(|f| loop_newell(graph, &f.outer_loop))
        .collect();
    fragments
        .iter()
        .enumerate()
        .map(|(k, ring)| {
            if !rim_only[k] {
                return None;
            }
            let mut best: Option<(f64, usize)> = None;
            for (j, host) in fragments.iter().enumerate() {
                if j == k || rim_only[j] || host.face != ring.face {
                    continue;
                }
                if newell[j].dot(&newell[k]) >= 0.0
                    || !loop_strictly_inside_outer(
                        graph,
                        &ring.outer_loop,
                        &host.outer_loop,
                        &newell[j],
                    )
                {
                    continue;
                }
                let area = newell[j].length();
                if best.is_none_or(|(a, _)| area < a) {
                    best = Some((area, j));
                }
            }
            best.map(|(_, j)| j)
        })
        .collect()
}

/// The fragments with every own hole filed apart rejoined to its host
/// and the ring itself removed. `None` when no ring was filed apart.
pub(crate) fn rejoin_own_holes(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    operands: &[&BRep],
    linear: f64,
) -> Option<Vec<Fragment>> {
    let reps = super::spans::unify_spans(graph, fragments, linear);
    let rim_only = hole_filling_fragments(graph, fragments, operands, &reps);
    let hosts = own_hole_hosts(graph, fragments, &rim_only);
    if hosts.iter().all(Option::is_none) {
        return None;
    }
    let mut out = fragments.to_vec();
    for (k, host) in hosts.iter().enumerate() {
        let Some(j) = *host else { continue };
        out[j].inner_loops.push(fragments[k].outer_loop.clone());
        // The host's chart view gains the ring's run, or is dropped
        // when the ring has none (a view must cover every loop).
        let ring_run = fragments[k]
            .chart
            .as_ref()
            .filter(|c| out[j].chart.as_ref().is_some_and(|h| h.sigma == c.sigma))
            .and_then(|c| c.loops.first().cloned());
        match (out[j].chart.as_mut(), ring_run) {
            (Some(view), Some(run)) => view.loops.push(run),
            _ => out[j].chart = None,
        }
    }
    Some(
        out.into_iter()
            .zip(&hosts)
            .filter(|(_, h)| h.is_none())
            .map(|(f, _)| f)
            .collect(),
    )
}

#[cfg(test)]
#[path = "hole_phantoms_tests.rs"]
mod tests;
