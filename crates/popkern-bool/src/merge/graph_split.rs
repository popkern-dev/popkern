//! Endpoint classification and graph edge splitting.
//!
//! Provides tools to determine where a 3D point falls relative to
//! existing graph topology (on a vertex, on an edge, or in a face
//! interior), and to split a graph half-edge pair at a given point.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::geometry::Point3;

use super::graph::{EdgeSource, GraphHalfEdge, IntersectionGraph, NodeSource};

/// Classification of an intersection curve endpoint relative to the graph.
#[derive(Debug)]
#[allow(clippy::enum_variant_names)]
pub(crate) enum EndpointClass {
    /// The point coincides with an existing graph node.
    OnVertex(usize),
    /// The point lies on a graph edge. Carries the half-edge index and the
    /// parametric position along the segment (0.0 = origin, 1.0 = end).
    OnEdge(usize, f64),
    /// The point is interior to a face (not on any existing edge or vertex).
    OnFaceInterior,
}

/// Classify where a point falls in the intersection graph.
///
/// Checks nodes first (within `tolerance`), then edges (point-to-segment
/// distance), and falls back to [`EndpointClass::OnFaceInterior`].
pub(crate) fn classify_endpoint(
    graph: &IntersectionGraph,
    point: &Point3,
    tolerance: f64,
) -> EndpointClass {
    // 1. Check all nodes.
    if let Some(idx) = graph.find_node(point, tolerance) {
        return EndpointClass::OnVertex(idx);
    }

    // 2. Check all half-edges (only check forward half-edges, i.e. even indices).
    //    For each edge pair (i, i+1), check if point is on the line segment.
    let mut i = 0;
    while i < graph.half_edges.len() {
        let he = &graph.half_edges[i];
        let twin = &graph.half_edges[he.twin];

        let p0 = graph.nodes[he.origin].position;
        let p1 = graph.nodes[twin.origin].position;

        if let Some(t) = point_on_segment(point, &p0, &p1, tolerance) {
            return EndpointClass::OnEdge(i, t);
        }

        // Skip the twin — we only need to check each edge once.
        i += 2;
    }

    // 3. Not on any vertex or edge.
    EndpointClass::OnFaceInterior
}

/// Check whether `point` lies on segment `p0..p1` within `tolerance`.
///
/// Returns `Some(t)` where t in (0,1) is the parameter, or `None`.
/// Endpoints (t near 0 or 1) are excluded — they would be caught by
/// the vertex check.
fn point_on_segment(point: &Point3, p0: &Point3, p1: &Point3, tolerance: f64) -> Option<f64> {
    let seg = *p1 - *p0;
    let seg_len_sq = seg.length_squared();
    if seg_len_sq < tolerance * tolerance {
        return None; // degenerate edge
    }

    let to_point = *point - *p0;
    let t = to_point.dot(&seg) / seg_len_sq;

    // Exclude endpoints (those should match as vertices).
    if t <= tolerance || t >= 1.0 - tolerance {
        return None;
    }

    // Check perpendicular distance.
    let projected = *p0 + seg * t;
    let dist = point.distance_to(&projected);
    if dist < tolerance { Some(t) } else { None }
}

/// Split a graph edge at a point, creating a new node.
///
/// The original half-edge pair (`he_idx` and its twin) is shortened to end
/// at the new node. A second pair is created from the new node to the
/// original endpoint. Face assignments are preserved on all affected
/// half-edges.
///
/// Returns the index of the new node.
pub(crate) fn split_graph_edge(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    split_point: Point3,
    _t_param: f64,
) -> usize {
    let twin_idx = graph.half_edges[he_idx].twin;

    // Original endpoints.
    let _start_node = graph.half_edges[he_idx].origin;
    let end_node = graph.half_edges[twin_idx].origin;

    // Preserve face assignments and source before mutation.
    let fwd_face = graph.half_edges[he_idx].face;
    let twin_face = graph.half_edges[twin_idx].face;
    let source = graph.half_edges[he_idx].source.clone();

    // Derive a split source from the original.
    let split_source = match &source {
        EdgeSource::OriginalA(eh) => EdgeSource::SplitA { original: *eh },
        EdgeSource::OriginalB(eh) => EdgeSource::SplitB { original: *eh },
        other => other.clone(),
    };

    // Create the new node at the split point.
    let new_node = graph.add_node(split_point, NodeSource::Intersection);

    // Create the second segment: new_node → end_node (and twin end_node → new_node).
    let new_fwd_idx = graph.half_edges.len();
    let new_twin_idx = new_fwd_idx + 1;

    graph.half_edges.push(GraphHalfEdge {
        origin: new_node,
        twin: new_twin_idx,
        face: fwd_face,
        curve: None,  // Curve splitting handled separately.
        pcurve: None, // Pcurve splitting handled separately.
        source: split_source.clone(),
    });
    graph.half_edges.push(GraphHalfEdge {
        origin: end_node,
        twin: new_fwd_idx,
        face: twin_face,
        curve: None,
        pcurve: None,
        source: split_source,
    });

    // Update outgoing lists.
    graph.nodes[new_node].outgoing.push(new_fwd_idx);
    graph.nodes[end_node].outgoing.push(new_twin_idx);

    // Remove old outgoing entry for end_node → twin_idx, since twin now
    // originates at end_node but goes to new_node instead of start_node.
    // Actually twin_idx still originates at end_node, but it now should go
    // to new_node. We update its target by re-pointing.

    // Update the original pair to end at new_node instead of end_node.
    // The forward half-edge (he_idx) keeps its origin (start_node).
    // Its twin (twin_idx) now originates at new_node instead of end_node.
    graph.half_edges[twin_idx].origin = new_node;

    // Fix outgoing: remove twin_idx from end_node, add to new_node.
    graph.nodes[end_node]
        .outgoing
        .retain(|&idx| idx != twin_idx);
    graph.nodes[new_node].outgoing.push(twin_idx);

    // Wire twins correctly:
    // he_idx ↔ twin_idx  (start_node → new_node, new_node → start_node)
    // new_fwd_idx ↔ new_twin_idx (new_node → end_node, end_node → new_node)
    // These are already set correctly by the push and original values.

    // Update face_edges: add new half-edges to the same faces.
    if let Some(face) = fwd_face {
        graph.face_edges.entry(face).or_default().push(new_fwd_idx);
    }
    if let Some(face) = twin_face {
        graph.face_edges.entry(face).or_default().push(new_twin_idx);
    }

    new_node
}

/// Expand self-loop edges — currently a no-op stub.
///
/// Self-loop expansion is deferred: shared self-loops work correctly
/// as inner loops, and original self-loops are handled as boundary
/// edges by the walk's self-loop preference.
pub(crate) fn expand_self_loops(
    _graph: &mut IntersectionGraph,
    _num_segments: usize,
    _snap_tol: f64,
) {
}

/// Split non-shared edges at shared intersection curve endpoint nodes.
///
/// After `insert_intersection_curves` (Phase 3), shared edges exist in the
/// graph with their endpoint nodes. Some of these nodes may lie in the
/// interior of existing non-shared edges — typically when an operand edge
/// lies ON a coincident face and the intersection curve endpoint is at
/// the boundary crossing. The normal edge-face hit detection misses these
/// because the edge lies on the face, not through it.
///
/// This function finds non-shared edges that pass through shared edge
/// endpoint nodes and splits them at those nodes (reusing the existing
/// node, no duplicates). Must be called after Phase 3.
pub(crate) fn split_edges_at_shared_nodes(graph: &mut IntersectionGraph, tolerance: f64) {
    // Collect shared edge endpoint node indices.
    let shared_nodes: Vec<usize> = {
        let set: std::collections::HashSet<usize> = graph
            .half_edges
            .iter()
            .filter(|he| matches!(he.source, EdgeSource::Shared { .. }))
            .flat_map(|he| [he.origin, graph.half_edges[he.twin].origin])
            .collect();
        set.into_iter().collect()
    };

    // For each shared node, check if any non-shared edge passes through it.
    // Restart after each split since indices change.
    let mut did_split = true;
    while did_split {
        did_split = false;
        for &node_idx in &shared_nodes {
            let pos = graph.nodes[node_idx].position;
            let mut i = 0;
            while i < graph.half_edges.len() {
                let he = &graph.half_edges[i];
                if i >= he.twin || matches!(he.source, EdgeSource::Shared { .. }) {
                    i += 2;
                    continue;
                }
                let twin_origin = graph.half_edges[he.twin].origin;
                if he.origin == node_idx || twin_origin == node_idx {
                    i += 2;
                    continue;
                }
                let p0 = graph.nodes[he.origin].position;
                let p1 = graph.nodes[twin_origin].position;
                if point_on_segment(&pos, &p0, &p1, tolerance).is_some() {
                    split_at_existing_node(graph, i, node_idx);
                    did_split = true;
                    break;
                }
                i += 2;
            }
            if did_split {
                break;
            }
        }
    }
}

/// Remove non-shared edges that duplicate a shared edge on the same face.
///
/// After `split_edges_at_shared_nodes`, split segments may connect the
/// same pair of nodes as a shared intersection curve. The walk can't
/// distinguish them, so remove the non-shared duplicate.
pub(crate) fn remove_collinear_duplicates(graph: &mut IntersectionGraph) {
    // For each face, find non-shared half-edges that connect the same
    // pair of nodes as a shared half-edge on the same face. Compare by
    // 3D position (immune to node index issues from split pointer chains).
    let tol = 1e-10;
    let mut to_remove: Vec<usize> = Vec::new();

    for he_indices in graph.face_edges.values() {
        // Collect shared edge endpoint positions on this face.
        let shared_segs: Vec<(Point3, Point3)> = he_indices
            .iter()
            .filter(|&&idx| matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }))
            .map(|&idx| {
                let he = &graph.half_edges[idx];
                let p0 = graph.nodes[he.origin].position;
                let p1 = graph.nodes[graph.half_edges[he.twin].origin].position;
                (p0, p1)
            })
            .collect();

        if shared_segs.is_empty() {
            continue;
        }

        // Check non-shared half-edges on this face.
        for &idx in he_indices {
            if matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }) {
                continue;
            }
            let he = &graph.half_edges[idx];
            let p0 = graph.nodes[he.origin].position;
            let p1 = graph.nodes[graph.half_edges[he.twin].origin].position;

            // Use the adaptive squared-distance predicate so the
            // duplicate-detection decision is order-independent at
            // the tolerance boundary.
            let tol_sq = tol * tol;
            let close =
                |a: &popkern_core::geometry::Point3, b: &popkern_core::geometry::Point3| -> bool {
                    popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                        (a.x, a.y, a.z),
                        (b.x, b.y, b.z),
                        tol_sq,
                    )
                };
            let is_dup = shared_segs.iter().any(|(sp0, sp1)| {
                (close(sp0, &p0) && close(sp1, &p1)) || (close(sp0, &p1) && close(sp1, &p0))
            });

            if is_dup {
                to_remove.push(idx);
            }
        }
    }

    for &idx in &to_remove {
        if let Some(f) = graph.half_edges[idx].face
            && let Some(edges) = graph.face_edges.get_mut(&f)
        {
            edges.retain(|&e| e != idx);
        }
        let origin = graph.half_edges[idx].origin;
        graph.nodes[origin].outgoing.retain(|&e| e != idx);
        graph.half_edges[idx].face = None;
    }
}

/// Split a graph edge at an existing node (no new node creation).
/// Split SHARED (intersection-curve) edges at the given operand
/// nodes: an edge-face hit puts a node on an operand edge, and any
/// intersection curve passing THROUGH that node (a boss seam
/// crossing a boss/plate intersection circle) must be split there or
/// the loop walk cannot route through the junction. This is the
/// inverse of [`split_edges_at_shared_nodes`].
pub(crate) fn split_shared_edges_at_nodes(
    graph: &mut IntersectionGraph,
    nodes: &[usize],
    tolerance: f64,
) {
    use popkern_core::geometry::curve::CurveEval as _;
    let mut did_split = true;
    while did_split {
        did_split = false;
        for &node_idx in nodes {
            let pos = graph.nodes[node_idx].position;
            let mut i = 0;
            while i < graph.half_edges.len() {
                let he = &graph.half_edges[i];
                if i >= he.twin || !matches!(he.source, EdgeSource::Shared { .. }) {
                    i += 2;
                    continue;
                }
                let twin_origin = graph.half_edges[he.twin].origin;
                if he.origin == node_idx || twin_origin == node_idx {
                    i += 2;
                    continue;
                }
                // Distance from the node to the edge's CURVE (shared
                // edges carry true intersection geometry — circles
                // included; the chord test would miss them).
                let on_curve = if let Some(curve) = &he.curve {
                    let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
                    let (t0, t1) = (curve.t_start, curve.t_end);
                    let mut best = f64::INFINITY;
                    let mut best_t = t0;
                    for k in 0..=64 {
                        let t = t0 + (t1 - t0) * k as f64 / 64.0;
                        let d = eval.point_at(t).distance_to(&pos);
                        if d < best {
                            best = d;
                            best_t = t;
                        }
                    }
                    let (mut lo, mut hi) = (
                        (best_t - (t1 - t0) / 64.0).max(t0.min(t1)),
                        (best_t + (t1 - t0) / 64.0).min(t0.max(t1)),
                    );
                    for _ in 0..40 {
                        let m1 = lo + (hi - lo) / 3.0;
                        let m2 = hi - (hi - lo) / 3.0;
                        if eval.point_at(m1).distance_to(&pos) < eval.point_at(m2).distance_to(&pos)
                        {
                            hi = m2;
                        } else {
                            lo = m1;
                        }
                    }
                    let t = 0.5 * (lo + hi);
                    // Interior only — endpoint touches are already nodes.
                    let interior = (t - t0).abs() > (t1 - t0).abs() * 1e-3
                        && (t - t1).abs() > (t1 - t0).abs() * 1e-3;
                    interior && eval.point_at(t).distance_to(&pos) < tolerance
                } else {
                    let p0 = graph.nodes[he.origin].position;
                    let p1 = graph.nodes[twin_origin].position;
                    point_on_segment(&pos, &p0, &p1, tolerance).is_some()
                };
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [shared-split] node {node_idx} vs shared he {i}: on_curve={on_curve}"
                    );
                }
                if on_curve {
                    split_at_existing_node(graph, i, node_idx);
                    did_split = true;
                    break;
                }
                i += 2;
            }
            if did_split {
                break;
            }
        }
    }
}

fn split_at_existing_node(graph: &mut IntersectionGraph, he_idx: usize, node: usize) {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let fwd_face = graph.half_edges[he_idx].face;
    let twin_face = graph.half_edges[twin_idx].face;
    let source = graph.half_edges[he_idx].source.clone();
    let split_source = match &source {
        EdgeSource::OriginalA(eh) => EdgeSource::SplitA { original: *eh },
        EdgeSource::OriginalB(eh) => EdgeSource::SplitB { original: *eh },
        other => other.clone(),
    };

    // Split the curve at the node's parameter: the halves must keep
    // TRUE geometry (same curve, sub-ranges). Dropping the curve on
    // a circular edge left the second half with a straight chord
    // ACROSS the circle — radial sorting then mis-routes the loop
    // walk. The node keeps the parent range's tail; the new half
    // takes [t_split, t_end].
    let (parent_curve_trimmed, new_half_curve) = match &graph.half_edges[he_idx].curve {
        Some(curve) => {
            use popkern_core::geometry::curve::CurveEval as _;
            let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
            let pos = graph.nodes[node].position;
            let (t0, t1) = (curve.t_start, curve.t_end);
            let mut best = (f64::INFINITY, t0);
            for k in 0..=128 {
                let t = t0 + (t1 - t0) * k as f64 / 128.0;
                let d = eval.point_at(t).distance_to(&pos);
                if d < best.0 {
                    best = (d, t);
                }
            }
            let (mut lo, mut hi) = (
                (best.1 - (t1 - t0) / 128.0).max(t0.min(t1)),
                (best.1 + (t1 - t0) / 128.0).min(t0.max(t1)),
            );
            for _ in 0..40 {
                let m1 = lo + (hi - lo) / 3.0;
                let m2 = hi - (hi - lo) / 3.0;
                if eval.point_at(m1).distance_to(&pos) < eval.point_at(m2).distance_to(&pos) {
                    hi = m2;
                } else {
                    lo = m1;
                }
            }
            let t_split = 0.5 * (lo + hi);
            (
                Some(popkern_core::geometry::CurveOnEdge {
                    geometry: curve.geometry.clone(),
                    t_start: t0,
                    t_end: t_split,
                }),
                Some(popkern_core::geometry::CurveOnEdge {
                    geometry: curve.geometry.clone(),
                    t_start: t_split,
                    t_end: t1,
                }),
            )
        }
        None => (None, None),
    };

    let new_fwd = graph.half_edges.len();
    let new_twin = new_fwd + 1;
    let new_half_curve_twin = new_half_curve.clone();

    graph.half_edges.push(GraphHalfEdge {
        origin: node,
        twin: new_twin,
        face: fwd_face,
        curve: new_half_curve,
        pcurve: None,
        source: split_source.clone(),
    });
    graph.half_edges.push(GraphHalfEdge {
        origin: end_node,
        twin: new_fwd,
        face: twin_face,
        curve: new_half_curve_twin,
        pcurve: None,
        source: split_source,
    });
    if let Some(trimmed) = parent_curve_trimmed {
        graph.half_edges[he_idx].curve = Some(trimmed.clone());
        // The parent's twin walks the same sub-range.
        graph.half_edges[twin_idx].curve = Some(trimmed);
    }

    graph.nodes[node].outgoing.push(new_fwd);
    graph.nodes[end_node].outgoing.push(new_twin);

    // Shorten original: twin now starts at the split node.
    graph.half_edges[twin_idx].origin = node;
    graph.nodes[end_node]
        .outgoing
        .retain(|&idx| idx != twin_idx);
    graph.nodes[node].outgoing.push(twin_idx);

    if let Some(face) = fwd_face {
        graph.face_edges.entry(face).or_default().push(new_fwd);
    }
    if let Some(face) = twin_face {
        graph.face_edges.entry(face).or_default().push(new_twin);
    }
}

#[cfg(test)]
#[path = "graph_split_tests.rs"]
mod tests;
