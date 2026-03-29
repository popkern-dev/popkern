// Items are pub(crate) for use by future modules (self-intersection, containment).
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! Loop tracing from an edge graph using smallest-angle planar traversal.
//!
//! Converts an [`EdgeGraph`] into a set of closed [`TracedLoop`]s, each
//! with a signed area (positive = CCW, negative = CW). The outermost
//! boundary loop (largest-magnitude CW area) is discarded.

use crate::edge_graph::EdgeGraph;
use popkern_constraint::ElementHandle;
use std::collections::HashSet;

/// A sketch edge with its traversal direction.
#[derive(Debug, Clone)]
pub(crate) struct DirectedEdge {
    /// Handle to the source sketch element.
    pub element: ElementHandle,
    /// 2D position of the start endpoint (in traversal direction).
    pub start: (f64, f64),
    /// 2D position of the end endpoint (in traversal direction).
    pub end: (f64, f64),
    /// `true` if this edge is traversed opposite to its graph definition.
    pub reversed: bool,
}

/// A closed loop of directed edges with its signed area.
#[derive(Debug, Clone)]
pub(crate) struct TracedLoop {
    /// Ordered directed edges forming the closed loop.
    pub edges: Vec<DirectedEdge>,
    /// Signed area via shoelace formula: positive = CCW, negative = CW.
    pub signed_area: f64,
}

/// A directed half-edge identifier: (edge_index, forward).
///
/// `forward = true` means traversal from `start_node` to `end_node`,
/// `forward = false` means traversal from `end_node` to `start_node`.
type HalfEdgeId = (usize, bool);

/// An entry in the per-node adjacency list.
struct AdjEntry {
    /// Index into `EdgeGraph::edges`.
    edge_index: usize,
    /// Whether this is the forward direction of the edge.
    forward: bool,
    /// The node at the other end of this directed half-edge.
    other_node: usize,
    /// Outgoing angle from this node toward `other_node`.
    outgoing_angle: f64,
}

/// Trace all closed loops from the edge graph.
///
/// Uses smallest-angle planar traversal:
/// 1. Build adjacency list with outgoing angles
/// 2. For each unvisited directed half-edge, trace a loop by always
///    picking the smallest CW turn (rightmost turn)
/// 3. Discard the outermost boundary (largest |negative area|)
pub(crate) fn trace_loops(graph: &EdgeGraph) -> Vec<TracedLoop> {
    if graph.edges.is_empty() {
        return Vec::new();
    }

    let adjacency = build_adjacency(graph);
    let mut visited: HashSet<HalfEdgeId> = HashSet::new();
    let mut loops = Vec::new();

    // Each undirected edge gives two directed half-edges.
    for edge_idx in 0..graph.edges.len() {
        for forward in [true, false] {
            if visited.contains(&(edge_idx, forward)) {
                continue;
            }
            if let Some(traced) = trace_one_loop(graph, &adjacency, &mut visited, edge_idx, forward)
            {
                loops.push(traced);
            }
        }
    }

    discard_outer_boundary(loops)
}

/// Build the adjacency list: for each node, list of outgoing directed
/// half-edges sorted by outgoing angle.
fn build_adjacency(graph: &EdgeGraph) -> Vec<Vec<AdjEntry>> {
    let n = graph.nodes.len();
    let mut adj: Vec<Vec<AdjEntry>> = (0..n).map(|_| Vec::new()).collect();

    for (idx, edge) in graph.edges.iter().enumerate() {
        let s = edge.start_node;
        let e = edge.end_node;
        let s_pos = graph.nodes[s];
        let e_pos = graph.nodes[e];

        // Forward half-edge: leaves from start_node toward end_node.
        let angle_fwd = (e_pos.1 - s_pos.1).atan2(e_pos.0 - s_pos.0);
        adj[s].push(AdjEntry {
            edge_index: idx,
            forward: true,
            other_node: e,
            outgoing_angle: angle_fwd,
        });

        // Reverse half-edge: leaves from end_node toward start_node.
        let angle_rev = (s_pos.1 - e_pos.1).atan2(s_pos.0 - e_pos.0);
        adj[e].push(AdjEntry {
            edge_index: idx,
            forward: false,
            other_node: s,
            outgoing_angle: angle_rev,
        });
    }

    // Sort each node's adjacency by outgoing angle.
    for entries in &mut adj {
        entries.sort_by(|a, b| {
            a.outgoing_angle
                .partial_cmp(&b.outgoing_angle)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    adj
}

/// Trace a single loop starting from a specific directed half-edge.
///
/// Returns `None` if the loop does not close or is degenerate.
fn trace_one_loop(
    graph: &EdgeGraph,
    adjacency: &[Vec<AdjEntry>],
    visited: &mut HashSet<HalfEdgeId>,
    start_edge_idx: usize,
    start_forward: bool,
) -> Option<TracedLoop> {
    let start_edge = &graph.edges[start_edge_idx];
    let (start_from, start_to) = if start_forward {
        (start_edge.start_node, start_edge.end_node)
    } else {
        (start_edge.end_node, start_edge.start_node)
    };

    let loop_start_node = start_from;
    let mut directed_edges: Vec<DirectedEdge> = Vec::new();

    let mut cur_edge_idx = start_edge_idx;
    let mut cur_forward = start_forward;
    let mut cur_from = start_from;
    let mut cur_to = start_to;

    // Safety limit to prevent infinite loops on malformed graphs.
    let max_steps = graph.edges.len() * 2 + 1;

    for _ in 0..max_steps {
        // Mark this half-edge as visited.
        visited.insert((cur_edge_idx, cur_forward));

        // Record this directed edge.
        let ge = &graph.edges[cur_edge_idx];
        let (start_pos, end_pos) = if cur_forward {
            (ge.start_pos, ge.end_pos)
        } else {
            (ge.end_pos, ge.start_pos)
        };

        directed_edges.push(DirectedEdge {
            element: ge.element,
            start: start_pos,
            end: end_pos,
            reversed: !cur_forward,
        });

        // Have we closed the loop?
        if cur_to == loop_start_node {
            let signed_area = compute_signed_area(&directed_edges);
            return Some(TracedLoop {
                edges: directed_edges,
                signed_area,
            });
        }

        // Pick the next half-edge at cur_to.
        // The "incoming angle" at cur_to is the angle from cur_from to cur_to.
        let from_pos = graph.nodes[cur_from];
        let to_pos = graph.nodes[cur_to];
        let incoming_angle = (to_pos.1 - from_pos.1).atan2(to_pos.0 - from_pos.0);

        // The twin of our incoming half-edge arrives at cur_to with
        // outgoing_angle = incoming_angle + PI (pointing back to cur_from).
        // We want the *next* edge CW from this twin's outgoing angle.
        let twin_outgoing = normalize_angle(incoming_angle + std::f64::consts::PI);

        // Among all outgoing half-edges at cur_to (excluding the twin of
        // the edge we just traversed), find the one whose outgoing angle
        // is the largest CW rotation from twin_outgoing. Equivalently,
        // it's the first edge CCW *before* the twin direction, i.e., the
        // one with the smallest positive (twin_outgoing - outgoing_angle)
        // mod 2*pi, but > 0.
        //
        // In practice: for each candidate, compute
        //   delta = (twin_outgoing - candidate_angle) mod 2*pi
        // and pick the smallest delta > 0. This gives the rightmost turn.

        let next = pick_next_edge(adjacency, cur_to, cur_edge_idx, cur_forward, twin_outgoing);

        match next {
            Some((next_idx, next_forward, next_other)) => {
                cur_from = cur_to;
                cur_to = next_other;
                cur_edge_idx = next_idx;
                cur_forward = next_forward;
            }
            None => {
                // Dead end -- not a closed loop.
                return None;
            }
        }
    }

    // Did not close within max steps.
    None
}

/// Pick the next half-edge at `node`, excluding the twin of the incoming
/// half-edge. Returns the candidate with the smallest positive CW delta
/// from `twin_outgoing`.
///
/// Returns `(edge_index, forward, other_node)`.
fn pick_next_edge(
    adjacency: &[Vec<AdjEntry>],
    node: usize,
    incoming_edge_idx: usize,
    incoming_forward: bool,
    twin_outgoing: f64,
) -> Option<(usize, bool, usize)> {
    let entries = &adjacency[node];
    if entries.is_empty() {
        return None;
    }

    // The twin of the incoming half-edge: same edge, opposite direction.
    let twin_forward = !incoming_forward;

    let mut best: Option<(f64, usize, &AdjEntry)> = None;
    let two_pi = 2.0 * std::f64::consts::PI;

    for (i, entry) in entries.iter().enumerate() {
        // Skip the twin of the edge we just traversed.
        if entry.edge_index == incoming_edge_idx && entry.forward == twin_forward {
            continue;
        }

        // CW delta from twin_outgoing to this candidate's outgoing angle.
        // We want: delta = (twin_outgoing - candidate_angle) mod 2*pi,
        // in range (0, 2*pi].
        let mut delta = twin_outgoing - entry.outgoing_angle;
        delta = delta.rem_euclid(two_pi);
        if delta < 1e-12 {
            delta = two_pi; // Treat ~0 as full rotation (same direction as twin)
        }

        match &best {
            Some((best_delta, best_idx, _)) => {
                if delta < *best_delta - 1e-9 {
                    best = Some((delta, i, entry));
                } else if (delta - *best_delta).abs() < 1e-9 && i < *best_idx {
                    // Tie-break: lower adjacency index (deterministic).
                    best = Some((delta, i, entry));
                }
            }
            None => {
                best = Some((delta, i, entry));
            }
        }
    }

    best.map(|(_, _, entry)| (entry.edge_index, entry.forward, entry.other_node))
}

/// Normalize an angle to `(-PI, PI]`.
fn normalize_angle(a: f64) -> f64 {
    let mut a = a.rem_euclid(2.0 * std::f64::consts::PI);
    if a > std::f64::consts::PI {
        a -= 2.0 * std::f64::consts::PI;
    }
    a
}

/// Compute signed area of a loop via the shoelace formula on edge
/// endpoint positions.
fn compute_signed_area(edges: &[DirectedEdge]) -> f64 {
    let mut sum = 0.0;
    for edge in edges {
        let (x_i, y_i) = edge.start;
        let (x_next, y_next) = edge.end;
        sum += x_i * y_next - x_next * y_i;
    }
    0.5 * sum
}

/// Discard the outermost boundary loop of each connected component.
///
/// Each connected component in a planar graph produces exactly one
/// unbounded face whose loop has the largest-magnitude negative (CW)
/// signed area within that component. Those loops are discarded.
fn discard_outer_boundary(loops: Vec<TracedLoop>) -> Vec<TracedLoop> {
    if loops.is_empty() {
        return loops;
    }

    // Group loops by connected component. Two loops are in the same
    // component if they share an edge (by ElementHandle identity).
    // Build a union-find over loop indices.
    let n = loops.len();
    let mut parent: Vec<usize> = (0..n).collect();

    // Collect all element handles per loop for intersection testing.
    let loop_elements: Vec<HashSet<usize>> = loops
        .iter()
        .map(|l| l.edges.iter().map(|e| element_hash(&e.element)).collect())
        .collect();

    for i in 0..n {
        for j in (i + 1)..n {
            if loop_elements[i]
                .iter()
                .any(|e| loop_elements[j].contains(e))
            {
                union(&mut parent, i, j);
            }
        }
    }

    // For each component, find the outer boundary (largest |negative area|).
    let mut outer_indices: HashSet<usize> = HashSet::new();
    let mut component_outer: std::collections::HashMap<usize, (usize, f64)> =
        std::collections::HashMap::new();

    for (i, l) in loops.iter().enumerate() {
        if l.signed_area < 0.0 {
            let root = find(&mut parent, i);
            let mag = l.signed_area.abs();
            let entry = component_outer.entry(root).or_insert((i, mag));
            if mag > entry.1 {
                *entry = (i, mag);
            }
        }
    }

    for (idx, _) in component_outer.values() {
        outer_indices.insert(*idx);
    }

    loops
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !outer_indices.contains(i))
        .map(|(_, l)| l)
        .collect()
}

/// Hash an element handle to a usize for set-based comparison.
fn element_hash(element: &ElementHandle) -> usize {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    element.hash(&mut hasher);
    hasher.finish() as usize
}

/// Union-find: find root with path compression.
fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Union-find: merge two sets.
fn union(parent: &mut [usize], a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[rb] = ra;
    }
}

#[cfg(test)]
#[path = "loop_tracer_tests.rs"]
mod tests;
