//! PK-295: split SHARED intersection-record edges at pre-existing
//! pierce nodes interior to their curves. Split out of
//! `graph_split.rs` (500-line ratchet).

use super::graph::{EdgeSource, IntersectionGraph};
use super::graph_split::{
    point_on_curved_graph_edge, point_on_graph_edge_curve, point_on_segment, split_at_existing_node,
};

/// Split SHARED intersection-record edges at existing nodes lying on
/// their curve INTERIOR. Phase 2 splits operand edges where they
/// pierce opposite faces (creating pierce nodes), and phase 3 then
/// inserts records whose curves can pass THROUGH such a node
/// mid-span — a cutter wall's spiric crossing a periodic face's seam
/// runs through the seam's pierce node unsplit, radial routing at
/// that node cannot turn onto the record, and the two fragments
/// beside the seam can never close their loops (PK-295, measured:
/// the θ=0 torus bulge of a rim-filleted cylinder dropped as open
/// chains, stranding its wall spiric and both flanking rim arcs as
/// self-twinned half-edges). Mirror of
/// [`split_edges_at_shared_nodes`], which splits NON-shared edges at
/// shared-record endpoints.
pub(crate) fn split_shared_edges_at_interior_nodes(graph: &mut IntersectionGraph, tolerance: f64) {
    let mut did_split = true;
    while did_split {
        did_split = false;
        for node_idx in 0..graph.nodes.len() {
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
                let p0 = graph.nodes[he.origin].position;
                let p1 = graph.nodes[twin_origin].position;
                let on_edge = match point_on_segment(&pos, &p0, &p1, tolerance) {
                    Some(t) => point_on_graph_edge_curve(graph, i, &pos, t, tolerance),
                    None => point_on_curved_graph_edge(graph, i, &pos, tolerance).is_some(),
                };
                if !on_edge
                    && std::env::var_os("PK_SPLIT_TRACE").is_some()
                    && point_on_curved_graph_edge(graph, i, &pos, tolerance * 1.0e4).is_some()
                {
                    eprintln!(
                        "[nsplit-nearmiss] node {node_idx} ({:.4},{:.4},{:.4}) near shared he {i} within {:.1e}",
                        pos.x,
                        pos.y,
                        pos.z,
                        tolerance * 1.0e4
                    );
                }
                if on_edge {
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

    // COMPLETE-THE-CUT pass (PK-305): a record chain legitimately
    // dead-ends ON a face boundary edge, but its endpoint carries the
    // marching/window slop class (~1e-3) — far beyond the snap
    // tolerance the passes above use — so the boundary edge never
    // split, the fragment walk carried the record as a slit inside
    // the un-cut face, and the served loop's UV chain kept a
    // half-period gap (measured: rotated torus × torus faces 18v1 /
    // 29v1, STEP refused by OCC). Scope: the node must END a Shared
    // record he, and the candidate edge must belong to one of that
    // record's OWN faces — the record's registration is the warrant
    // for cutting that face's boundary. Verification stays
    // curve-true at the widened gate.
    let loose = tolerance * 1.0e3;
    let mut did_split = true;
    while did_split {
        did_split = false;
        // (node, faces-of-record) pairs for shared endpoints.
        let mut shared_ends: Vec<(
            usize,
            popkern_core::topology::FaceHandle,
            popkern_core::topology::FaceHandle,
        )> = Vec::new();
        for he in graph.half_edges.iter() {
            if let EdgeSource::Shared { face_a, face_b } = &he.source {
                shared_ends.push((he.origin, *face_a, *face_b));
            }
        }
        'outer: for &(node_idx, rec_fa, rec_fb) in &shared_ends {
            let pos = graph.nodes[node_idx].position;
            let mut i = 0;
            while i < graph.half_edges.len() {
                let he = &graph.half_edges[i];
                if i >= he.twin || matches!(he.source, EdgeSource::Shared { .. }) {
                    i += 2;
                    continue;
                }
                let edge_face_matches = he
                    .face
                    .as_ref()
                    .is_some_and(|f| f.face == rec_fa || f.face == rec_fb)
                    || graph.half_edges[he.twin]
                        .face
                        .as_ref()
                        .is_some_and(|f| f.face == rec_fa || f.face == rec_fb);
                if !edge_face_matches {
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
                // Endpoint guard at the widened gate: a node within
                // the slop class of an edge VERTEX is a vertex touch,
                // not a mid-edge crossing.
                if pos.distance_to(&p0) < loose || pos.distance_to(&p1) < loose {
                    i += 2;
                    continue;
                }
                let on_edge = match point_on_segment(&pos, &p0, &p1, loose) {
                    Some(t) => point_on_graph_edge_curve(graph, i, &pos, t, loose),
                    None => point_on_curved_graph_edge(graph, i, &pos, loose).is_some(),
                };
                if on_edge {
                    split_at_existing_node(graph, i, node_idx);
                    did_split = true;
                    break 'outer;
                }
                i += 2;
            }
        }
    }
}

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
                if let Some(t) = point_on_segment(&pos, &p0, &p1, tolerance) {
                    // Chord test only — verify against the real
                    // curve for curved edges (same defect as
                    // classify_endpoint: a 180° arc was split at
                    // points on its chord, off-surface).
                    if point_on_graph_edge_curve(graph, i, &pos, t, tolerance) {
                        split_at_existing_node(graph, i, node_idx);
                        did_split = true;
                        break;
                    }
                } else if point_on_curved_graph_edge(graph, i, &pos, tolerance).is_some() {
                    // Mirror case: the node sits ON the curve but
                    // off its chord.
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

    // COMPLETE-THE-CUT pass (PK-305): a record chain legitimately
    // dead-ends ON a face boundary edge, but its endpoint carries the
    // marching/window slop class (~1e-3) — far beyond the snap
    // tolerance the passes above use — so the boundary edge never
    // split, the fragment walk carried the record as a slit inside
    // the un-cut face, and the served loop's UV chain kept a
    // half-period gap (measured: rotated torus × torus faces 18v1 /
    // 29v1, STEP refused by OCC). Scope: the node must END a Shared
    // record he, and the candidate edge must belong to one of that
    // record's OWN faces — the record's registration is the warrant
    // for cutting that face's boundary. Verification stays
    // curve-true at the widened gate.
    let loose = tolerance * 1.0e3;
    let mut did_split = true;
    while did_split {
        did_split = false;
        // (node, faces-of-record) pairs for shared endpoints.
        let mut shared_ends: Vec<(
            usize,
            popkern_core::topology::FaceHandle,
            popkern_core::topology::FaceHandle,
        )> = Vec::new();
        for he in graph.half_edges.iter() {
            if let EdgeSource::Shared { face_a, face_b } = &he.source {
                shared_ends.push((he.origin, *face_a, *face_b));
            }
        }
        'outer: for &(node_idx, rec_fa, rec_fb) in &shared_ends {
            let pos = graph.nodes[node_idx].position;
            let mut i = 0;
            while i < graph.half_edges.len() {
                let he = &graph.half_edges[i];
                if i >= he.twin || matches!(he.source, EdgeSource::Shared { .. }) {
                    i += 2;
                    continue;
                }
                let edge_face_matches = he
                    .face
                    .as_ref()
                    .is_some_and(|f| f.face == rec_fa || f.face == rec_fb)
                    || graph.half_edges[he.twin]
                        .face
                        .as_ref()
                        .is_some_and(|f| f.face == rec_fa || f.face == rec_fb);
                if !edge_face_matches {
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
                // Endpoint guard at the widened gate: a node within
                // the slop class of an edge VERTEX is a vertex touch,
                // not a mid-edge crossing.
                if pos.distance_to(&p0) < loose || pos.distance_to(&p1) < loose {
                    i += 2;
                    continue;
                }
                let on_edge = match point_on_segment(&pos, &p0, &p1, loose) {
                    Some(t) => point_on_graph_edge_curve(graph, i, &pos, t, loose),
                    None => point_on_curved_graph_edge(graph, i, &pos, loose).is_some(),
                };
                if on_edge {
                    split_at_existing_node(graph, i, node_idx);
                    did_split = true;
                    break 'outer;
                }
                i += 2;
            }
        }
    }
}
