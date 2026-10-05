//! Splitting a graph edge at an EXISTING node (extracted from
//! `graph_split.rs` under the 500-line file budget). See
//! [`split_at_existing_node_with_pcurves`].

use super::graph::{EdgeSource, GraphHalfEdge, IntersectionGraph};

/// Returns `true` when the edge was actually split — a refused
/// split (window-boundary minimizer) MUST NOT count as progress in
/// a caller's restart loop (measured, PK-491/PK-342: at linear 1e-2
/// a node re-matched the same piece every restart while the guard
/// refused the split; the tolerance-envelope sweep spun for two
/// hours in `split_shared_edges_at_nodes`).
pub(super) fn split_at_existing_node(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    node: usize,
) -> bool {
    split_at_existing_node_with_pcurves(graph, he_idx, node, None)
}

/// Closest curve-window parameter to `pos`, or `None` when the
/// minimizer sits ON the window boundary (relative slop 1e-9) — a
/// split there would produce a zero-span half.
///
/// Winding-aware (PK-670): the historical 128-station scan refined
/// around its single nearest station aliases onto an adjacent turn of
/// a many-turn edge (past ~25 turns at radius/pitch 2.2 the stations
/// fall below about two per turn), splitting the edge one pitch away
/// from the node. The search refines every scan basin at a density
/// scaled to the curve's knot spans and keeps the nearest
/// ([`nearest_param_all_basins`]); the historical answer stands
/// whenever it is already on the curve or no basin is nearer.
pub(super) fn curve_window_interior_split(
    curve: &popkern_core::geometry::CurveOnEdge,
    pos: &popkern_core::geometry::Point3,
) -> Option<f64> {
    use super::param_basin_search::{curve_scan_count, nearest_param_all_basins};
    const LEGACY_STATIONS: usize = 128;
    let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
    let (t0, t1) = (curve.t_start, curve.t_end);
    let (w_lo, w_hi) = (t0.min(t1), t0.max(t1));
    let dist = |t: f64| eval.point_at(t).distance_to(pos);
    let refine = |center: f64, half: f64| {
        let (mut lo, mut hi) = ((center - half).max(w_lo), (center + half).min(w_hi));
        for _ in 0..40 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            if dist(m1) < dist(m2) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        0.5 * (lo + hi)
    };
    let dense = curve_scan_count(&curve.geometry, t0, t1, LEGACY_STATIONS);
    // "On the curve" for the early exit: far below any turn spacing,
    // above the evaluation noise of the window.
    let exact = 1.0e-9 * (1.0 + pos.x.abs() + pos.y.abs() + pos.z.abs());
    let (t_split, _) =
        nearest_param_all_basins(&dist, (t0, t1), (LEGACY_STATIONS, dense), exact, &refine);
    let span = (t1 - t0).abs().max(f64::MIN_POSITIVE);
    if (t_split - t0).abs() < span * 1.0e-9 || (t_split - t1).abs() < span * 1.0e-9 {
        return None;
    }
    Some(t_split)
}

/// [`split_at_existing_node`] that also CUTS both sides' pcurves at
/// the node's image when the operand breps are supplied (PK-441 —
/// see `graph_split_pcurve`); without them the pcurves drop, the
/// historical behavior.
pub(super) fn split_at_existing_node_with_pcurves(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    node: usize,
    breps: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) -> bool {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    if popkern_core::trace::trace_enabled("face-split") {
        let sp = graph.nodes[graph.half_edges[he_idx].origin].position;
        let ep = graph.nodes[end_node].position;
        let np = graph.nodes[node].position;
        let w = graph.half_edges[he_idx]
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[twin_idx].curve.as_ref())
            .map(|c| format!("[{:.3},{:.3}]", c.t_start, c.t_end))
            .unwrap_or_else(|| "none".into());
        eprintln!(
            "[nsplit] he {he_idx} ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1}) w{w} at node ({:.1},{:.1},{:.1})",
            sp.x, sp.y, sp.z, ep.x, ep.y, ep.z, np.x, np.y, np.z
        );
    }
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
    // The curve may live on EITHER half of the pair — reading only
    // he_idx left a twin-carried curve with its stale FULL window,
    // which the curved-edge matcher then re-matched at the same
    // node, splitting the edge again and again (duplicate n↔n
    // pairs on a gusset cap's arc).
    let pair_curve = graph.half_edges[he_idx]
        .curve
        .clone()
        .or_else(|| graph.half_edges[twin_idx].curve.clone());
    let (parent_curve_trimmed, new_half_curve) = match &pair_curve {
        Some(curve) => {
            // A closest point ON the window boundary would mint a
            // zero-span half whose endpoints are nonetheless
            // distinct nodes — geometry that lies (measured,
            // PK-491: the groove tool's crest-flank edge split at
            // the cone-seam crossing node 2.4e-3 from its own
            // junction endpoint; the [0.90,0.90] stub read as a
            // second curve through the node and derailed the flank
            // walk's wedge fan). Refuse such a split.
            let Some(t_split) = curve_window_interior_split(curve, &graph.nodes[node].position)
            else {
                if popkern_core::trace::trace_enabled("face-split") {
                    let np = graph.nodes[node].position;
                    eprintln!(
                        "[nsplit] REFUSED window-boundary split of he {he_idx} at node ({:.4},{:.4},{:.4})",
                        np.x, np.y, np.z
                    );
                }
                return false;
            };
            let (t0, t1) = (curve.t_start, curve.t_end);
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
        // Stale full-span pcurves poison embeddings and the exported
        // result (see split_graph_edge) — drop them; consumers
        // re-derive from the trimmed curve. With operand breps in
        // hand, cut them honestly instead (PK-441).
        let fwd_pc = graph.half_edges[he_idx].pcurve.take();
        let twin_pc = graph.half_edges[twin_idx].pcurve.take();
        if let Some((brep_a, brep_b, tol)) = breps {
            let split_pos = graph.nodes[node].position;
            let origin_pos = graph.nodes[graph.half_edges[he_idx].origin].position;
            let end_pos = graph.nodes[end_node].position;
            super::graph_split_pcurve::distribute_split_pcurves(
                graph,
                he_idx,
                twin_idx,
                new_fwd,
                new_twin,
                fwd_pc,
                twin_pc,
                &split_pos,
                &origin_pos,
                &end_pos,
                brep_a,
                brep_b,
                tol,
            );
        }
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
    // SHARED records are cross-registered: BOTH directions live in
    // BOTH faces' walk sets (each face's two fragments consume one
    // direction each). The per-face registration above covers only
    // one direction per face, so a walk arriving at the split node
    // could not turn onto the new piece and sailed across the
    // junction (measured: a joint-shoulder record split at a rail
    // flank's T-point still walked as if unsplit, merging the lap
    // ring with the region beyond the shoulder).
    if let EdgeSource::Shared { face_a, face_b } = graph.half_edges[new_fwd].source {
        for face in [
            crate::merge::graph::FaceRef {
                operand: crate::merge::graph::Operand::A,
                face: face_a,
            },
            crate::merge::graph::FaceRef {
                operand: crate::merge::graph::Operand::B,
                face: face_b,
            },
        ] {
            let edges = graph.face_edges.entry(face).or_default();
            for idx in [new_fwd, new_twin] {
                if !edges.contains(&idx) {
                    edges.push(idx);
                }
            }
        }
    }
    true
}

/// Split original edges at edge-face pierce hits, preserving pcurves,
/// and return the created/adopted node indices. Shared by the
/// cellular and graph-merge pipelines.
pub(crate) fn split_edges_at_face_hits(
    graph: &mut IntersectionGraph,
    hits: &[crate::intersect_edge_face::EdgeFaceHit],
    snap_tol: f64,
    breps: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) -> Vec<usize> {
    let mut nodes = Vec::new();
    for hit in hits {
        let cls = super::graph_split::classify_endpoint(graph, &hit.point, snap_tol);
        if let super::graph_split::EndpointClass::OnEdge(he_idx, _t) = cls {
            nodes.push(super::graph_split::split_graph_edge_with_pcurves(
                graph, he_idx, hit.point, breps,
            ));
        }
    }
    nodes
}
