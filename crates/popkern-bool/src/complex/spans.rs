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

use std::collections::{BTreeMap, BTreeSet};

use crate::cellular::assemble_edge_geometry::span_sample_points;
use crate::cellular::assemble_span_unify::SpanRegistry;
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
    let mut spans = SpanRegistry::new(linear * 10.0, linear);
    let mut out = BTreeMap::new();
    for key in keys {
        let sample = span_sample_points(graph, key);
        let rep = match spans.prior_for(graph, key, sample) {
            // A registered span is its own representative.
            Some(prior) => {
                let start = |k: usize| graph.nodes[graph.half_edges[k].origin].position;
                let end = |k: usize| {
                    graph.nodes[graph.half_edges[graph.half_edges[k].twin].origin].position
                };
                let here = start(key);
                SpanRep {
                    key: prior,
                    same_dir: here.distance_to(&start(prior)) <= here.distance_to(&end(prior)),
                }
            }
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
    out
}
