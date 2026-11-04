//! Splitting operand edges at shared-record endpoints, and completing
//! dead-end cuts. Split out of `graph_split.rs` (500-line ratchet).

use super::graph::{EdgeSource, IntersectionGraph};
use super::graph_split::{point_on_curved_graph_edge, point_on_graph_edge_curve, point_on_segment};
use super::graph_split_node::split_at_existing_node_with_pcurves;

pub(crate) fn split_edges_at_shared_nodes(
    graph: &mut IntersectionGraph,
    tolerance: f64,
    breps: Option<(&[&popkern_core::topology::BRep], f64)>,
) {
    // Collect shared edge endpoint node indices, in ASCENDING order:
    // the split order below decides every new graph index, and a
    // HashSet's per-process seed made the same boolean mint its split
    // pieces in different orders run to run (measured: a torus-groove
    // cut's output swapped edge handles 6/7 between runs, and the
    // follow-on drill served or declined by run — the PK-512 class).
    let shared_nodes: Vec<usize> = {
        let set: std::collections::BTreeSet<usize> = graph
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
                if node_on_edge(graph, i, pos, tolerance) {
                    did_split = split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
                    if did_split {
                        break;
                    }
                }
                i += 2;
            }
            if did_split {
                break;
            }
        }
    }

    complete_dead_end_cuts(graph, tolerance, breps);
}

/// Whether `pos` lies on graph edge `i`, by its chord then its true
/// curve, or on its curve off the chord: where
/// [`split_edges_at_shared_nodes`] splits an operand edge.
pub(crate) fn node_on_edge(
    graph: &IntersectionGraph,
    i: usize,
    pos: popkern_core::geometry::Point3,
    tolerance: f64,
) -> bool {
    let he = &graph.half_edges[i];
    let p0 = graph.nodes[he.origin].position;
    let p1 = graph.nodes[graph.half_edges[he.twin].origin].position;
    match point_on_segment(&pos, &p0, &p1, tolerance) {
        // Chord test only — verify against the real curve for curved
        // edges (same defect as classify_endpoint: a 180° arc was
        // split at points on its chord, off-surface).
        Some(t) => point_on_graph_edge_curve(graph, i, &pos, t, tolerance),
        // Mirror case: the node sits ON the curve but off its chord.
        None => point_on_curved_graph_edge(graph, i, &pos, tolerance).is_some(),
    }
}

/// Every (record-end node, operand edge) junction
/// [`split_edges_at_shared_nodes`] would still split in `graph`: zero
/// once whatever split the operand edges found every junction the pass
/// finds.
pub(crate) fn unsplit_operand_junctions(graph: &IntersectionGraph, tolerance: f64) -> usize {
    let shared_nodes: std::collections::BTreeSet<usize> = graph
        .half_edges
        .iter()
        .filter(|he| matches!(he.source, EdgeSource::Shared { .. }))
        .flat_map(|he| [he.origin, graph.half_edges[he.twin].origin])
        .collect();
    let mut count = 0;
    for i in (0..graph.half_edges.len()).filter(|&i| {
        let he = &graph.half_edges[i];
        i < he.twin && !matches!(he.source, EdgeSource::Shared { .. })
    }) {
        let ends = [
            graph.half_edges[i].origin,
            graph.half_edges[graph.half_edges[i].twin].origin,
        ];
        let bound = super::graph_split_at_nodes::shared_edge_bounds(graph, i);
        count += shared_nodes
            .iter()
            .filter(|&&n| {
                let p = graph.nodes[n].position;
                !ends.contains(&n)
                    && super::graph_split_at_nodes::within_bounds(bound, &p, tolerance)
                    && node_on_edge(graph, i, p, tolerance)
            })
            .count();
    }
    count
}

/// COMPLETE-THE-CUT (PK-305), shared by both split passes: a record
/// chain that dead-ends within the marching slop of a boundary edge of
/// one of its own faces splits that edge at the dead end.
pub(crate) fn complete_dead_end_cuts(
    graph: &mut IntersectionGraph,
    tolerance: f64,
    breps: Option<(&[&popkern_core::topology::BRep], f64)>,
) {
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
    let loose_base = tolerance * 1.0e3;
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
            if let EdgeSource::Shared { face_a, face_b, .. } = &he.source {
                shared_ends.push((he.origin, *face_a, *face_b));
            }
        }
        // (record face pair, boundary original) combinations already
        // JOINED at some node: the chain has its junction with that
        // boundary edge, and a SECOND loose split of the same pair
        // fuses geometrically distinct points (measured, PK-491:
        // the M5 12/5 groove-wall chain meets the tool's
        // crest-flank corner at their true junction, while the
        // chain's cone-seam crossing node 2.4e-3 away sat 3.3e-4
        // off the corner curve — the loose split there wired the
        // crest's boundary through an off-crest vertex and the
        // transactional validate declined every assembly). A chain
        // nowhere joined to the candidate edge keeps the PK-305
        // completion.
        let original_of = |src: &EdgeSource| match *src {
            EdgeSource::Original(o, e)
            | EdgeSource::Split {
                operand: o,
                original: e,
            } => Some((o.index() as u8, e)),
            _ => None,
        };
        let joined: std::collections::HashSet<(
            popkern_core::topology::FaceHandle,
            popkern_core::topology::FaceHandle,
            u8,
            popkern_core::topology::EdgeHandle,
        )> = {
            let mut set = std::collections::HashSet::new();
            for node in &graph.nodes {
                let mut pairs = Vec::new();
                let mut origs = Vec::new();
                for &h in &node.outgoing {
                    match &graph.half_edges[h].source {
                        EdgeSource::Shared { face_a, face_b, .. } => pairs.push((*face_a, *face_b)),
                        s => {
                            if let Some(o) = original_of(s) {
                                origs.push(o);
                            }
                        }
                    }
                }
                for &(fa, fb) in &pairs {
                    for &(op, e) in &origs {
                        set.insert((fa, fb, op, e));
                    }
                }
            }
            set
        };
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
                // Only a DEAD END completes: where two or more records
                // on the candidate edge's face meet, the chain turns a
                // corner inside the face and goes on — pulling that
                // corner onto a boundary edge within the slop gate
                // fuses two distinct points (PK-790: a stadium box's
                // corner sat 1e-3 inside the rod rim, inside the 2e-3
                // gate; the rim was split at the corner and its arc
                // piece served 1e-3 off the wall).
                let twin = &graph.half_edges[he.twin];
                if [he.face.as_ref(), twin.face.as_ref()]
                    .into_iter()
                    .flatten()
                    .any(|f| {
                        (f.face == rec_fa || f.face == rec_fb)
                            && chain_degree(graph, node_idx, f) >= 2
                    })
                {
                    i += 2;
                    continue;
                }
                // Already joined elsewhere: skip (see `joined`).
                let Some((cand_op, cand_orig)) = original_of(&he.source) else {
                    i += 2;
                    continue;
                };
                if joined.contains(&(rec_fa, rec_fb, cand_op, cand_orig)) {
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
                // The widened gate is capped GEOMETRY-RELATIVE (5% of
                // the candidate edge's chord): tolerance × 1e3 is an
                // absolute distance, and at linear 1e-4 it reached
                // 0.2 — exactly a bore circle's clearance from the
                // cube rim — so the rim was split at a node 0.2 OFF
                // it and every coarse-tolerance boolean in that
                // resonance served garbage (PK-342; at 1e-3 the same
                // runaway gate tripped its own endpoint guard, which
                // masked the defect). The torus half-period seam slop
                // this gate exists for is ~1e-3 on unit-scale edges —
                // comfortably inside the cap.
                let loose = loose_base.min(p0.distance_to(&p1) * 0.05);
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
                    did_split = split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
                    if !did_split {
                        i += 2;
                        continue;
                    }
                    break 'outer;
                }
                i += 2;
            }
        }
    }
}

/// Intersection records at `node` with `face` on one side (a record's
/// `face_a` is operand A's, `face_b` operand B's).
fn chain_degree(graph: &IntersectionGraph, node: usize, face: &super::graph::FaceRef) -> usize {
    graph.nodes[node]
        .outgoing
        .iter()
        .filter(|&&h| match &graph.half_edges[h].source {
            EdgeSource::Shared {
                face_a,
                face_b,
                operands,
            } => {
                (operands[0] == face.operand && *face_a == face.face)
                    || (operands[1] == face.operand && *face_b == face.face)
            }
            _ => false,
        })
        .count()
}
