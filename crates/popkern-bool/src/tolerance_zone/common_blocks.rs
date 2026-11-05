//! Common blocks (plan §4.3): a record lying along an operand edge for
//! its whole length is that edge, not a second curve beside it.
//!
//! Records of a face with a coplanar neighbour of the other operand
//! trace the neighbour's boundary: a gear flank meets the mating gear's
//! bottom plane exactly along its own bottom outline. The graph inserts
//! both the record and the edge, and each operand's walk reconciles them
//! on its own (the duplicate-cover dedup on one face, nothing on the
//! other), so the two operands partition their common overlap
//! differently and the complex cannot merge it as one face.
//!
//! A record is a common block with operand edge `e` when both its ends
//! resolve onto `e` (an interior incidence, or an operand vertex `e`
//! ends at) and its interior samples lie on `e`'s curve within the sum
//! of the record's and the edge's tolerances.

use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
use std::collections::BTreeMap;

use popkern_core::topology::{BRep, EdgeHandle};
use slotmap::Key;

use super::candidates::{CandidateVertex, Incidence};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// Interior samples of a record checked against the edge.
const SAMPLES: usize = 5;

/// Each operand's edges by the vertices they end at.
type EdgesAt = [BTreeMap<u64, Vec<EdgeHandle>>; 2];

/// `brep`'s edges by end vertex, built once: a scan of every edge per
/// record end costs a boolean its records times its operands' edges.
fn edges_at(brep: &BRep) -> BTreeMap<u64, Vec<EdgeHandle>> {
    let mut out: BTreeMap<u64, Vec<EdgeHandle>> = BTreeMap::new();
    for (eh, e) in &brep.edges {
        let Some(h) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let twin_origin = h
            .mate()
            .and_then(|m| brep.half_edges.get(m))
            .map(|t| t.origin);
        for v in std::iter::once(h.origin).chain(twin_origin) {
            let list = out.entry(v.data().as_ffi()).or_default();
            if !list.contains(&eh) {
                list.push(eh);
            }
        }
    }
    out
}

/// The operand edges the group of candidate `i` lies on.
fn group_edges(
    cands: &[CandidateVertex],
    res: &Resolution,
    edges_at: &EdgesAt,
    i: usize,
) -> Vec<(usize, EdgeHandle)> {
    let mut out = Vec::new();
    for &m in &res.groups[res.group_of[i]].members {
        for inc in &cands[m].incidence {
            match *inc {
                Incidence::Edge { operand, edge, .. } => out.push((operand, edge)),
                Incidence::Vertex { operand, vertex } => {
                    if let Some(edges) = edges_at[operand].get(&vertex.data().as_ffi()) {
                        out.extend(edges.iter().map(|&eh| (operand, eh)));
                    }
                }
                _ => {}
            }
        }
    }
    out.sort_by_key(|(o, e)| (*o, e.data().as_ffi()));
    out.dedup();
    out
}

/// Every record of `records` that is a common block, with each edge it
/// lies along.
pub(crate) fn common_block_records(
    cands: &[CandidateVertex],
    res: &Resolution,
    operands: [&BRep; 2],
    records: &[IntersectionRecord],
    linear: f64,
) -> Vec<(usize, usize, EdgeHandle)> {
    let ends = super::consequences::record_end_candidates(cands, records.len());
    let edges_at: EdgesAt = [edges_at(operands[0]), edges_at(operands[1])];
    let mut out = Vec::new();
    for (ri, [s, e]) in ends.into_iter().enumerate() {
        let (Some(s), Some(e)) = (s, e) else {
            continue;
        };
        let at_end = group_edges(cands, res, &edges_at, e);
        let shared: Vec<_> = group_edges(cands, res, &edges_at, s)
            .into_iter()
            .filter(|x| at_end.contains(x))
            .collect();
        let curve = &records[ri].curve.curve_3d;
        let (t0, t1) = curve.parameter_range();
        for (operand, eh) in shared {
            let Some(edge) = operands[operand].edges.get(eh) else {
                continue;
            };
            let Some(c) = edge.curve.as_ref() else {
                continue;
            };
            let reach = 2.0 * linear + edge.tolerance.unwrap_or(linear);
            let along = (1..=SAMPLES).all(|k| {
                let p = curve.point_at(t0 + (t1 - t0) * k as f64 / (SAMPLES + 1) as f64);
                let t = closest_point_in_window(&c.geometry, &p, c.t_start, c.t_end);
                c.geometry.point_at(t).distance_to(&p) < reach
            });
            // Every edge it lies along: where the two operands touch
            // along an edge, the record is a block with an edge of each.
            if along {
                out.push((ri, operand, eh));
            }
        }
    }
    out
}
