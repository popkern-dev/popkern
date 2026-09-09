//! PK-295: split SHARED intersection-record edges at pre-existing
//! pierce nodes interior to their curves. Split out of
//! `graph_split.rs` (500-line ratchet).

use super::graph::{EdgeSource, IntersectionGraph};
use super::graph_split::{point_on_curved_graph_edge, point_on_graph_edge_curve, point_on_segment};
use super::graph_split_node::split_at_existing_node_with_pcurves;

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
pub(crate) fn split_shared_edges_at_interior_nodes(
    graph: &mut IntersectionGraph,
    tolerance: f64,
    breps: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) {
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
                    && popkern_core::trace::trace_enabled("face-split")
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
                    split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
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
                    split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
                    did_split = true;
                    break 'outer;
                }
                i += 2;
            }
        }
    }
}

pub(crate) fn split_edges_at_shared_nodes(
    graph: &mut IntersectionGraph,
    tolerance: f64,
    breps: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) {
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
                        split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
                        did_split = true;
                        break;
                    }
                } else if point_on_curved_graph_edge(graph, i, &pos, tolerance).is_some() {
                    // Mirror case: the node sits ON the curve but
                    // off its chord.
                    split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
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
                    split_at_existing_node_with_pcurves(graph, i, node_idx, breps);
                    did_split = true;
                    break 'outer;
                }
                i += 2;
            }
        }
    }
}

/// Weld DEAD-END shared-record endpoint nodes onto the boundary node
/// they missed.
///
/// A clipped record's endpoint is evaluated on the FITTED record
/// curve and carries its lateral slop (~1e-4..1e-3 measured), while
/// the boundary edge it crosses was split at the edge-face hit's
/// exact Newton solve — two nodes for one junction, up to the fit
/// slop apart, far beyond the graph snap (2× linear). The record
/// then dead-ends at its own node, the chart walk's graze filter
/// deletes the pair, the face never partitions, and the merge
/// strands (measured: a chamfered rod's bottom-face groove chain
/// ended at n43 = (−1.239468, 2.731961) while the rim split sat at
/// n18 = (−1.239313, 2.732050), 1.8e-4 apart — every headed
/// fastener declined, PK-460).
///
/// Scope is the graze filter's own trigger: only endpoint nodes whose
/// ENTIRE incidence is shared-record half-edges are candidates, and
/// only a node carrying non-record incidence (a boundary split or
/// original vertex) within the marching-slop band may receive them.
/// A genuine interior graze has no such neighbour and keeps its
/// honest drop; fully-wired junction clusters (the 5.6e-4 exact
/// connectors of PK-441) are out of scope by the dead-end condition.
pub(crate) fn weld_dead_end_shared_endpoints(graph: &mut IntersectionGraph, band: f64) {
    let n_nodes = graph.nodes.len();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (i, he) in graph.half_edges.iter().enumerate() {
        incident[he.origin].push(i);
    }
    let is_shared = |h: usize| {
        matches!(
            graph.half_edges[h].source,
            crate::merge::graph::EdgeSource::Shared { .. }
        )
    };
    // Beyond the tuned unconditional band, a VERIFIED annulus (up to
    // 10× further): weld only when the candidate junction lies on the
    // dead-end record curve's extension — perpendicular distance from
    // the node to curve samples just past the endpoint. A shallow
    // record × boundary crossing AMPLIFIES pcurve slop into the
    // along-curve endpoint miss (measured, PK-464: the mouth-chamfer
    // circle crossed the nut's crest-band helix boundary at slope
    // 0.16, a sub-1e-3 boundary-pcurve wobble moved the clip window
    // end 6.2e-4 from the true junction node, the ring segment
    // dead-ended and the graze cascade unravelled the whole chamfer
    // ring). The perpendicular gate keeps wrong-junction welds out:
    // an off-curve node at the same spatial distance stays unwelded.
    let on_extension = |n: usize, m_pos: &popkern_core::geometry::Point3| -> bool {
        use popkern_core::geometry::curve::CurveEval as _;
        let perp_band = band * 0.2; // 100× linear at the tuned 500× band
        for &h in &incident[n] {
            let he = &graph.half_edges[h];
            let Some(c) = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref())
            else {
                continue;
            };
            let span = c.t_end - c.t_start;
            if !span.is_finite() || span.abs() < 1e-30 {
                continue;
            }
            // Which curve end sits at `n`; outward tangent by an
            // inside finite difference (a clamped NURBS does not
            // extrapolate past its knots, so sampling beyond the end
            // is not an option).
            let p_n = graph.nodes[n].position;
            let d_start = c.geometry.point_at(c.t_start).distance_to(&p_n);
            let d_end = c.geometry.point_at(c.t_end).distance_to(&p_n);
            let (t_at, t_in) = if d_start <= d_end {
                (c.t_start, c.t_start + span * 1e-3)
            } else {
                (c.t_end, c.t_end - span * 1e-3)
            };
            let p_end = c.geometry.point_at(t_at);
            let inward = c.geometry.point_at(t_in) - p_end;
            let Some(outward) = (inward * -1.0).normalized() else {
                continue;
            };
            let to_m = *m_pos - p_end;
            let along = to_m.dot(&outward);
            if along < -perp_band {
                continue; // candidate sits back along the curve, not past the end
            }
            let perp = (to_m - outward * along).length();
            if perp < perp_band {
                return true;
            }
        }
        false
    };
    let mut welds: Vec<(usize, usize)> = Vec::new();
    for n in 0..n_nodes {
        if incident[n].is_empty() {
            continue;
        }
        let all_record = incident[n].iter().all(|&h| is_shared(h));
        if !all_record {
            continue;
        }
        let p = graph.nodes[n].position;
        let mut best: Option<(f64, usize)> = None;
        for (m, inc) in incident.iter().enumerate() {
            if m == n || inc.is_empty() {
                continue;
            }
            let has_boundary = inc.iter().any(|&h| !is_shared(h));
            if !has_boundary {
                continue;
            }
            let d = graph.nodes[m].position.distance_to(&p);
            let accept = d < band || (d < band * 10.0 && on_extension(n, &graph.nodes[m].position));
            if accept && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, m));
            }
        }
        if let Some((d, m)) = best {
            if popkern_core::trace::trace_enabled("boolean") {
                let q = graph.nodes[m].position;
                eprintln!(
                    "  [weld-dead-end] node {n} ({:.6},{:.6},{:.6}) -> {m} ({:.6},{:.6},{:.6}) d={d:.2e}",
                    p.x, p.y, p.z, q.x, q.y, q.z
                );
            }
            welds.push((n, m));
        }
    }
    for (n, m) in welds {
        // Snap the welded record curves' end control points onto the
        // target node (the closed-record rotation weld's convention):
        // the graph re-origin alone leaves the CURVE ending at the
        // old position, and the assembled body then carries the weld
        // distance as a real shell gap between neighbours (measured:
        // 6.2e-4 at the nut's chamfer-ring junction — the serve
        // failed validation with exactly the welded distance).
        use popkern_core::geometry::curve::CurveEval as _;
        let old_pos = graph.nodes[n].position;
        let new_pos = graph.nodes[m].position;
        let d_weld = old_pos.distance_to(&new_pos);
        let moved: Vec<usize> = graph.nodes[n].outgoing.clone();
        for &h in &moved {
            for hh in [h, graph.half_edges[h].twin] {
                if let Some(c) = graph.half_edges[hh].curve.as_mut() {
                    let snap_band = d_weld * 2.0 + band * 0.2;
                    let d_start = c.geometry.point_at(c.t_start).distance_to(&old_pos);
                    let d_end = c.geometry.point_at(c.t_end).distance_to(&old_pos);
                    // A clamped NURBS's FIRST control point sits at
                    // the parametric MINIMUM — map the near half-edge
                    // end through the traversal direction.
                    let near_is_t_start = d_start <= d_end;
                    let near_is_param_min = near_is_t_start == (c.t_start <= c.t_end);
                    let near_d = if near_is_t_start { d_start } else { d_end };
                    if near_d < snap_band
                        && let popkern_core::geometry::curve::CurveGeometry::Nurbs(nc) =
                            &mut c.geometry
                    {
                        let cp = if near_is_param_min {
                            nc.control_points.first_mut()
                        } else {
                            nc.control_points.last_mut()
                        };
                        if let Some(cp) = cp {
                            *cp = new_pos;
                        }
                    }
                }
            }
        }
        for h in moved {
            graph.half_edges[h].origin = m;
            graph.nodes[m].outgoing.push(h);
        }
        graph.nodes[n].outgoing.clear();
    }
}
