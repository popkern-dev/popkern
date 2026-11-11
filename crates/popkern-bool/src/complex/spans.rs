//! Coincident graph spans are one complex edge.
//!
//! Two graph edge pairs connecting the same two points along the same
//! curve — a boundary piece and its overlap record on a flush
//! interface, or a record running along an edge an operand already has
//! — are one edge of the result. The assembler unifies them with its
//! span registry; the complex must too, or each copy carries only some
//! of the edge's fins and the complex has free edges where the result
//! has none (measured: the shadow harness's first suite run, 210 of 253
//! complex declines, every flush family).

use popkern_core::topology::BRep;
use std::collections::{BTreeMap, BTreeSet};

use crate::cellular::assemble_edge_geometry::span_sample_points;
use crate::cellular::assemble_span_unify::{SpanRegistry, runs_same_way};
use crate::merge::graph::IntersectionGraph;
use crate::merge::walk::Fragment;

use super::lift::is_degenerate_key;

/// The edge a graph edge key lifts to: the representative key of its
/// coincident span, and whether the key runs the representative's way.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SpanRep {
    /// Representative graph edge key.
    pub key: usize,
    /// Whether the key's forward half-edge runs as the representative's.
    pub same_dir: bool,
}

/// Every non-degenerate key the fragments use, mapped to its span's
/// representative. Keys are visited in canonical (ascending) order, as
/// the assembler visits them, so both pick the same representative.
pub(crate) fn unify_spans(
    operands: &[&BRep],
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    linear: f64,
) -> BTreeMap<usize, SpanRep> {
    let keys: BTreeSet<usize> = fragments
        .iter()
        .flat_map(|f| f.outer_loop.iter().chain(f.inner_loops.iter().flatten()))
        .map(|&h| h.min(graph.half_edges[h].twin))
        .filter(|&k| !is_degenerate_key(graph, k))
        .collect();
    let reach = popkern_core::tolerance::Tolerance::default()
        .with_linear(linear)
        .linear_reach();
    let mut spans = SpanRegistry::new(linear * 10.0, reach);
    let mut out = BTreeMap::new();
    for key in keys {
        let sample = span_sample_points(graph, key);
        let rep = match spans.prior_for(graph, operands, key, sample) {
            // A registered span is its own representative.
            Some(prior) => SpanRep {
                key: prior,
                same_dir: runs_same_way(graph, key, prior, linear * 10.0),
            },
            None => {
                let (nodes, a, b, mid) = sample;
                spans.register(nodes, a, b, mid, key);
                SpanRep {
                    key,
                    same_dir: true,
                }
            }
        };
        out.insert(key, rep);
    }
    follow_duplicate_removals(graph, &mut out, linear);
    out
}

/// An operand piece the duplicate-cover pass removed from one face (in
/// favour of the record it duplicated there) keeps its twin on its
/// other face, and the graph judged the two pairs one curve
/// (`IntersectionGraph::duplicate_of`): the surviving half lifts to the
/// record's edge. Measured: a boss cap 1e-6 past a ring's extreme; the
/// boss wall's rim arc was removed for the ring-wall record 1e-6 off it
/// (inside the duplicate band, outside the span registry's linear
/// band), the cap kept its rim half, and the complex held the record
/// with three fins and the rim with one. A closed pair (coincident
/// ends) gives no direction by its ends and keeps its own edge.
fn follow_duplicate_removals(
    graph: &IntersectionGraph,
    out: &mut BTreeMap<usize, SpanRep>,
    linear: f64,
) {
    let start = |k: usize| graph.nodes[graph.half_edges[k].origin].position;
    let end = |k: usize| graph.nodes[graph.half_edges[graph.half_edges[k].twin].origin].position;
    for (&removed, &record) in &graph.duplicate_of {
        let key = removed.min(graph.half_edges[removed].twin);
        let record_key = record.min(graph.half_edges[record].twin);
        let (Some(&from), Some(&to)) = (out.get(&key), out.get(&record_key)) else {
            continue;
        };
        if from.key == to.key || start(to.key).distance_to(&end(to.key)) <= linear {
            continue;
        }
        for rep in out.values_mut().filter(|r| r.key == from.key) {
            let here = start(from.key);
            let along = here.distance_to(&start(to.key)) <= here.distance_to(&end(to.key));
            *rep = SpanRep {
                key: to.key,
                same_dir: rep.same_dir == along,
            };
        }
    }
}
