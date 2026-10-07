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
            // The segment test uses the node-to-node CHORD. For a
            // CURVED edge a point can sit on the chord while being
            // nowhere near the curve (a star tip at (12,0) lay on a
            // 180° arc's chord and the hit machinery split the arc
            // there, off-surface). Verify against the real curve.
            if !point_on_graph_edge_curve(graph, i, point, t, tolerance) {
                i += 2;
                continue;
            }
            return EndpointClass::OnEdge(i, t);
        }
        // The MIRROR failure: a point ON a curved edge can be far
        // from its chord (an arc's step-face crossing sat 5.5 mm off
        // the chord), so the segment test never fires. Test curved
        // edges against the real curve directly.
        if let Some(t) = point_on_curved_graph_edge(graph, i, point, tolerance) {
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
pub(super) fn point_on_segment(
    point: &Point3,
    p0: &Point3,
    p1: &Point3,
    tolerance: f64,
) -> Option<f64> {
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
    split_graph_edge_with_pcurves(graph, he_idx, split_point, None)
}

/// [`split_graph_edge`] that also CUTS both sides' pcurves at the
/// split point's image when the operand breps are supplied (PK-441 —
/// see `graph_split_pcurve`); without them the pcurves drop, the
/// historical behavior.
pub(crate) fn split_graph_edge_with_pcurves(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    split_point: Point3,
    breps: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) -> usize {
    let twin_idx = graph.half_edges[he_idx].twin;

    // Original endpoints.
    let _start_node = graph.half_edges[he_idx].origin;
    let end_node = graph.half_edges[twin_idx].origin;

    // Split the CURVE window at the split point. The original pair
    // keeps the window on the start side, the new fragment gets the
    // tail. Leaving the original with its FULL window (and the
    // fragment with none) meant a hit-split arc kept a curve
    // overrunning its nodes — assembly's endpoint-matched adoption
    // failed and the arc degraded to chords (a gusset cap's arc
    // boundary, split where a step face crosses it, walked as
    // garbage bigons). Inversion of the split point on the curve
    // gives the split parameter; on inversion failure the old
    // behavior is preserved.
    let start_pos = graph.nodes[graph.half_edges[he_idx].origin].position;
    if popkern_core::trace::trace_enabled("face-split") {
        let ep = graph.nodes[end_node].position;
        let w = graph.half_edges[he_idx]
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[twin_idx].curve.as_ref())
            .map(|c| format!("[{:.3},{:.3}]", c.t_start, c.t_end))
            .unwrap_or_else(|| "none".into());
        eprintln!(
            "[gsplit] he {he_idx} ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1}) w{w} at ({:.1},{:.1},{:.1})",
            start_pos.x,
            start_pos.y,
            start_pos.z,
            ep.x,
            ep.y,
            ep.z,
            split_point.x,
            split_point.y,
            split_point.z
        );
    }
    let curve_windows = graph.half_edges[he_idx]
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[twin_idx].curve.as_ref())
        .and_then(|c| {
            use popkern_core::geometry::curve::CurveEval;
            let (t_split, dist) = nearest_param_on_curve(c, &split_point);
            // Defensive: a split point far off the curve means the
            // caller's bookkeeping is wrong — keep the old (no
            // window split) behavior rather than corrupt windows.
            // This decline is otherwise silent (the edge keeps its
            // full window and downstream adoption degrades to
            // chords), so it traces under the face-split key.
            if dist > 1.0e-3 {
                if popkern_core::trace::trace_enabled("face-split") {
                    eprintln!(
                        "[gsplit-decline] he {he_idx} window split declined: nearest t {t_split:.6} dist {dist:.6} on w[{:.3},{:.3}]",
                        c.t_start, c.t_end
                    );
                }
                return None;
            }
            let p_ts = c.geometry.point_at(c.t_start);
            let d_start = p_ts.distance_to(&start_pos);
            let d_end = p_ts.distance_to(&graph.nodes[end_node].position);
            // A CLOSED edge (both curve ends at one vertex — a
            // boolean-promoted rim circle) ties this test; the pcurve
            // then carries the fwd side's true traversal direction
            // (PK-441 — see `curve_co_directed_with_fwd`).
            let start_first = if (d_start - d_end).abs() > 1.0e-9 {
                d_start <= d_end
            } else {
                breps
                    .and_then(|(a, b, _)| {
                        super::graph_split_pcurve::curve_co_directed_with_fwd(graph, he_idx, a, b)
                    })
                    .unwrap_or(true)
            };
            if start_first {
                Some(((c.t_start, t_split), (t_split, c.t_end)))
            } else {
                Some(((t_split, c.t_end), (c.t_start, t_split)))
            }
        });

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

    let tail_curve = curve_windows.and_then(|(_, (ta, tb))| {
        let src = graph.half_edges[he_idx]
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[twin_idx].curve.as_ref())?;
        let mut c = src.clone();
        c.t_start = ta;
        c.t_end = tb;
        Some(c)
    });
    graph.half_edges.push(GraphHalfEdge {
        origin: new_node,
        twin: new_twin_idx,
        face: fwd_face,
        curve: tail_curve,
        pcurve: None, // Pcurves are regenerated downstream.
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
    // Trim the ORIGINAL pair's curve to the head window, and DROP its
    // pcurves: a split parent's stored pcurve still spans the WHOLE
    // original edge, and every consumer that trusted it inherited the
    // lie — the chart embedding ran the full original extent (one
    // wrong-branch run flipped a face's boundary shoelace, σ
    // mis-calibrated, and the role assembly turned region-separating
    // records into holes), and the assembled result carried the
    // full-span pcurve so the STEP export's branch normalization saw
    // whole-period UV gaps. With no pcurve the embedding falls back
    // to projection of the trimmed 3D curve and the result pcurve
    // recompute re-derives the correct sub-span.
    let mut split_side_pcurves: Option<(
        Option<popkern_core::geometry::pcurve::Pcurve>,
        Option<popkern_core::geometry::pcurve::Pcurve>,
    )> = None;
    if let Some(((ha, hb), _)) = curve_windows {
        split_side_pcurves = Some((
            graph.half_edges[he_idx].pcurve.take(),
            graph.half_edges[twin_idx].pcurve.take(),
        ));
        for idx in [he_idx, twin_idx] {
            if let Some(c) = graph.half_edges[idx].curve.as_mut() {
                c.t_start = ha;
                c.t_end = hb;
            }
            graph.half_edges[idx].pcurve = None;
        }
    }
    if popkern_core::trace::trace_enabled("face-split")
        && let Some(((ha, hb), (na2, nb2))) = curve_windows
    {
        eprintln!("[gsplit]   orig w[{ha:.3},{hb:.3}] frag w[{na2:.3},{nb2:.3}]");
    }

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

    if let (Some((fwd_pc, twin_pc)), Some((brep_a, brep_b, tol))) = (split_side_pcurves, breps) {
        let end_pos = graph.nodes[end_node].position;
        super::graph_split_pcurve::distribute_split_pcurves(
            graph,
            he_idx,
            twin_idx,
            new_fwd_idx,
            new_twin_idx,
            fwd_pc,
            twin_pc,
            &split_point,
            &start_pos,
            &end_pos,
            brep_a,
            brep_b,
            tol,
        );
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
/// True when `point` (already on the node-to-node CHORD of graph
/// edge `he_idx` at chord parameter `t_chord`) also lies on the
/// edge's REAL curve within `tolerance`, NEAR the chord position.
/// Straight/curveless edges are trivially true; curved edges get a
/// local ternary-search refinement around the chord-parameter
/// estimate. Locality is load-bearing: a GLOBAL on-curve test
/// matches every node lying anywhere on a closed ring against every
/// chord segment of that ring, and the spurious OnEdge splits
/// stranded the filleted-torus assembly (measured during PK-449 —
/// do not widen this to the knot-seeded global inversion).
pub(super) fn point_on_graph_edge_curve(
    graph: &IntersectionGraph,
    he_idx: usize,
    point: &Point3,
    t_chord: f64,
    tolerance: f64,
) -> bool {
    let he = &graph.half_edges[he_idx];
    let curve = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref());
    let curved_curve = curve.filter(|c| match &c.geometry {
        popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
        popkern_core::geometry::curve::CurveGeometry::Nurbs(n) => n.degree > 1,
        _ => true,
    });
    use popkern_core::geometry::curve::CurveEval;
    let Some(c) = curved_curve else {
        return true;
    };
    let tc = c.t_start + (c.t_end - c.t_start) * t_chord;
    let mut best = c.geometry.point_at(tc).distance_to(point);
    let span = (c.t_end - c.t_start).abs().max(1e-12);
    let mut lo = (tc - 0.1 * span).max(c.t_start.min(c.t_end));
    let mut hi = (tc + 0.1 * span).min(c.t_start.max(c.t_end));
    for _ in 0..24 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = c.geometry.point_at(m1).distance_to(point);
        let d2 = c.geometry.point_at(m2).distance_to(point);
        if d1 < d2 {
            hi = m2;
            best = best.min(d1);
        } else {
            lo = m1;
            best = best.min(d2);
        }
    }
    best <= tolerance
}

/// Nearest parameter on a curve window to `point`: coarse scan +
/// local ternary refinement. Returns `(t, distance)`. This is the
/// split-parameter workhorse — Newton-based inversion failed
/// outright on mixed-path rational arcs (measured: both hit-splits
/// of a gusset cap arc lost their windows, the stale full window
/// kept re-matching, and the edge was split five ways with
/// duplicate pairs). Its first-minimum single-bracket semantics are
/// LOAD-BEARING for the split-window and curve-adoption sites:
/// delegating to the adaptive multi-sheet inversion cracked a
/// cap-cut sphere's rim windows during PK-449 — multi-sheet
/// robustness lives at the classification site
/// (`point_on_curved_graph_edge`), which is where the density
/// defect was first measured. Only the coarse DENSITY adapts here
/// (PK-451): 32 stays the floor — every window with ≤ 16 NURBS
/// knot spans scans bit-identically to the historical 32-sample
/// body — but a many-turn curve (the 6.5-turn coil seam put 4.9
/// samples per turn at n = 32) gets 2 samples per span so the
/// coarse minimum cannot land on the adjacent turn's basin.
pub(crate) fn nearest_param_on_curve(
    curve: &popkern_core::geometry::CurveOnEdge,
    point: &Point3,
) -> (f64, f64) {
    use popkern_core::geometry::curve::CurveEval;
    let (ta, tb) = (curve.t_start, curve.t_end);
    let n: usize = match &curve.geometry {
        popkern_core::geometry::curve::CurveGeometry::Nurbs(nurbs) => {
            let (lo, hi) = (ta.min(tb), ta.max(tb));
            let mut spans = 0usize;
            let mut last = f64::NEG_INFINITY;
            for &k in &nurbs.knots {
                if k >= lo && k <= hi && k > last {
                    spans += 1;
                    last = k;
                }
            }
            (2 * spans).clamp(32, 4096)
        }
        _ => 32,
    };
    let mut best = (0usize, f64::MAX);
    for i in 0..=n {
        let t = ta + (tb - ta) * i as f64 / n as f64;
        let d = curve.geometry.point_at(t).distance_to(point);
        if d < best.1 {
            best = (i, d);
        }
    }
    let (mut lo, mut hi) = (
        ta + (tb - ta) * (best.0.saturating_sub(1)) as f64 / n as f64,
        ta + (tb - ta) * (best.0 + 1).min(n) as f64 / n as f64,
    );
    let mut best_t = ta + (tb - ta) * best.0 as f64 / n as f64;
    let mut best_d = best.1;
    for _ in 0..48 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = curve.geometry.point_at(m1).distance_to(point);
        let d2 = curve.geometry.point_at(m2).distance_to(point);
        if d1 < d2 {
            hi = m2;
            if d1 < best_d {
                best_d = d1;
                best_t = m1;
            }
        } else {
            lo = m1;
            if d2 < best_d {
                best_d = d2;
                best_t = m2;
            }
        }
    }
    (best_t, best_d)
}

/// Whether any knot span of `n` overlapping `[ta, tb]` has `p`
/// inside its LOCAL control hull's AABB inflated by `pad`. Each
/// span's curve piece is a convex combination of its `degree + 1`
/// control points, so a point beyond every span box provably cannot
/// lie within `pad` of the queried curve segment — an EXACT reject
/// for on-curve decisions, at three float compares per span instead
/// of a knot-seeded global inversion. Multi-wrap tool boundary coils
/// carry hundreds of spans, and the splitters' node × edge scans
/// paid the full inversion per pair — O(nodes × spans) at 80+
/// thread turns, the cellular merge stage's residual quadratic term
/// after PK-496 verdict propagation (measured by stack sampling:
/// `nearest_param_on_curve` under `split_edges_at_shared_nodes`).
fn span_hull_within(
    n: &popkern_core::geometry::NurbsCurve,
    ta: f64,
    tb: f64,
    p: &Point3,
    pad: f64,
) -> bool {
    let d = n.degree;
    let (lo_t, hi_t) = if ta <= tb { (ta, tb) } else { (tb, ta) };
    if n.knots.len() < 2 * d + 2 || n.control_points.len() + d + 1 != n.knots.len() {
        return true; // malformed shape: no opinion, run the exact test
    }
    for j in d..(n.knots.len() - 1 - d) {
        let (k0, k1) = (n.knots[j], n.knots[j + 1]);
        if k1 <= k0 || k1 < lo_t || k0 > hi_t {
            continue;
        }
        let mut inside = true;
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for cp in &n.control_points[j - d..=j] {
            for (k, c) in [cp.x, cp.y, cp.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
        for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
            if c < lo[k] - pad || c > hi[k] + pad {
                inside = false;
                break;
            }
        }
        if inside {
            return true;
        }
    }
    false
}

/// If `point` lies on the REAL curve of graph edge `he_idx` (which
/// must be CURVED — straight edges are fully served by the chord
/// test), return its chord-equivalent fraction. Coarse scan +
/// local ternary refinement; interior hits only (endpoint
/// matches are the caller's node-snap job).
pub(super) fn point_on_curved_graph_edge(
    graph: &IntersectionGraph,
    he_idx: usize,
    point: &Point3,
    tolerance: f64,
) -> Option<f64> {
    let he = &graph.half_edges[he_idx];
    let curve = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())?;
    let curved = match &curve.geometry {
        popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
        popkern_core::geometry::curve::CurveGeometry::Nurbs(n) => n.degree > 1,
        _ => true,
    };
    if !curved {
        return None;
    }
    // Knot-seeded global inversion: the private 16-sample scan here
    // ranked the ADJACENT COIL TURN's basin above the true one on a
    // 6.5-turn helical seam (0.209 vs 0.219), so the seam × plate
    // pierce point classified OnFaceInterior, the seam never split at
    // the crossing, and the whole wall dropped to a knife-edge sample
    // (PK-449 — the station-aligned deep grind's real mechanism).
    let (ta, tb) = (curve.t_start, curve.t_end);
    // Span-hull pre-gate (PK-496): a point beyond every overlapping
    // span's local control hull cannot be an interior on-curve hit —
    // reject before paying the global inversion. A hit outside the
    // queried [ta, tb] sub-range would clamp to an endpoint fraction
    // and be rejected below anyway, so scoping the gate to
    // overlapping spans preserves the answer exactly.
    if let popkern_core::geometry::curve::CurveGeometry::Nurbs(n) = &curve.geometry
        && !span_hull_within(n, ta, tb, point, tolerance)
    {
        return None;
    }
    let (best_t, best_d) =
        popkern_core::geometry::curve_inversion::nearest_param_on_curve(curve, point);
    if best_d > tolerance {
        return None;
    }
    let f = ((best_t - ta) / (tb - ta)).clamp(0.0, 1.0);
    // Interior only: endpoints belong to node snapping.
    if !(1e-6..=(1.0 - 1e-6)).contains(&f) {
        return None;
    }
    Some(f)
}

pub(crate) use super::graph_split_at_nodes::split_shared_edges_at_nodes;

#[cfg(test)]
#[path = "graph_split_tests.rs"]
mod tests;
