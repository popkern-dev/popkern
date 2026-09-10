//! Dead-end record-endpoint welding for the merge graph (split from
//! `graph_split_shared.rs` under the 500-line file cap). See
//! [`weld_dead_end_shared_endpoints`].

use super::graph::{EdgeSource, IntersectionGraph};

#[path = "graph_weld_certify.rs"]
mod graph_weld_certify;
use graph_weld_certify::{on_carrier, on_extension};

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
pub(crate) fn weld_dead_end_shared_endpoints(
    graph: &mut IntersectionGraph,
    band: f64,
    carriers: Option<(
        &popkern_core::topology::BRep,
        &popkern_core::topology::BRep,
        f64,
    )>,
) {
    let n_nodes = graph.nodes.len();
    // The band is tolerance-scaled by the caller (linear × 500) —
    // meaningful at fine tolerances, MODEL-SIZED at coarse ones (a
    // 1e-2-tolerance unit model got a 5.0 band; the dead-end welds
    // swallowed its record network and the subtraction SERVED an
    // empty body). Cap by the node cloud's diagonal.
    let band = {
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for node in &graph.nodes {
            let p = node.position;
            for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[i] = lo[i].min(c);
                hi[i] = hi[i].max(c);
            }
        }
        let diag =
            ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
        if diag.is_finite() && diag > 0.0 {
            band.min(diag * 1.0e-3)
        } else {
            band
        }
    };
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
    let mut welds: Vec<(usize, usize)> = Vec::new();
    for n in 0..n_nodes {
        if incident[n].is_empty() {
            continue;
        }
        let all_record = incident[n].iter().all(|&h| is_shared(h));
        if !all_record {
            continue;
        }
        // The incident-degree (lone-endpoint) tiering restriction was
        // measured NOT load-bearing for the tap's chamfer (6 strands
        // either way) while it blocked the far-tier haul the m6
        // nut's mouth chamfer needs — the far tier keeps its own
        // carrier+cone certificates as the acceptance bar.
        let lone_endpoint = true;

        let p = graph.nodes[n].position;
        // Carrier pair of this dead end's record, for the close tier.
        let node_carriers = carriers.and_then(|(brep_a, brep_b, _)| {
            incident[n].iter().find_map(|&h| {
                if let EdgeSource::Shared { face_a, face_b } = graph.half_edges[h].source {
                    let sa = &brep_a.faces.get(face_a)?.surface.as_ref()?.geometry;
                    let sb = &brep_b.faces.get(face_b)?.surface.as_ref()?.geometry;
                    Some((sa, sb))
                } else {
                    None
                }
            })
        });
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
            // Never weld across an existing record connection: two
            // junctions JOINED by a record piece are distinct by
            // construction (the tuned exact-connector class — the
            // m6 nut's mouth-chamfer chain carries real 1.5e-3
            // connectors between junctions 1.0e-3 apart; welding
            // their ends collapsed the connector into a self-loop,
            // the stub cleanup retired it, and the chain gapped into
            // a graze cascade). A weld target must be an UNCONNECTED
            // duplicate of this node.
            let connected = incident[n]
                .iter()
                .any(|&h| graph.half_edges[graph.half_edges[h].twin].origin == m);
            if connected {
                continue;
            }
            // CLOSE TIER (band..3×band): a record endpoint stopping a
            // hair past the unconditional band next to an EXISTING
            // junction is the same missed-junction class, but at that
            // range the tangent-ray direction is noise (measured,
            // PK-465: two mouth junctions 6.1e-4 and 7.3e-4 from
            // their record ends on a 5e-4 band — unwelded, each
            // rooted a graze cascade that unwalked a whole helix
            // ring). Direction cannot certify there; lying on BOTH
            // carrier surfaces of the record's face pair can.
            let close_certified = || -> bool {
                let Some((sa, sb)) = node_carriers else {
                    return false;
                };
                let Some((_, _, lin)) = carriers else {
                    return false;
                };
                let q = graph.nodes[m].position;
                let gate = (d * 0.02).max(lin * 100.0);
                on_carrier(sa, &q, gate) && on_carrier(sb, &q, gate)
            };
            let accept = d < band
                || (d < band * 3.0 && close_certified())
                || (lone_endpoint
                    && d < band * 100.0
                    && on_extension(graph, &incident, band, n, &graph.nodes[m].position));
            if accept && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, m));
            }
        }
        // FAR TIER: a record crossing the face boundary NEAR-TANGENT
        // amplifies clip slop by 1/sin(crossing angle) into an
        // along-boundary endpoint miss far past the tangent-ray
        // annulus (measured, PK-465: the angled tap's thread exit
        // crosses the mouth rim at helix pitch angle — every exit
        // record stopped 0.04–0.15 short of an EXISTING rim/seam/
        // pierce junction node, the graze cascade unwalked the exit
        // ring, and the merge stranded 20 sides). At this scale the
        // tangent-ray perp test is meaningless (curvature sag alone
        // exceeds the band), so the certificate changes shape: the
        // junction must lie on BOTH carrier surfaces of the record's
        // face pair (within a slop-scaled gate), sit inside a 25°
        // forward cone of the record's endpoint tangent, be within
        // 20% of the record's own chord, and be UNAMBIGUOUS (second
        // candidate at least twice as far).
        if best.is_none()
            && lone_endpoint
            && let Some((brep_a, brep_b, lin)) = carriers
        {
            use popkern_core::geometry::curve::CurveEval as _;
            let rec_scale = incident[n]
                .iter()
                .map(|&h| {
                    let q = graph.half_edges[graph.half_edges[h].twin].origin;
                    graph.nodes[q].position.distance_to(&p)
                })
                .fold(0.0f64, f64::max);
            let far_cap = (0.2 * rec_scale).max(band * 100.0);
            let outward_dirs: Vec<popkern_core::geometry::Vec3> = incident[n]
                .iter()
                .filter_map(|&h| {
                    let he = &graph.half_edges[h];
                    let c = he
                        .curve
                        .as_ref()
                        .or_else(|| graph.half_edges[he.twin].curve.as_ref())?;
                    let span = c.t_end - c.t_start;
                    if !span.is_finite() || span.abs() < 1e-30 {
                        return None;
                    }
                    let d_start = c.geometry.point_at(c.t_start).distance_to(&p);
                    let d_end = c.geometry.point_at(c.t_end).distance_to(&p);
                    let (t_at, t_in) = if d_start <= d_end {
                        (c.t_start, c.t_start + span * 1e-3)
                    } else {
                        (c.t_end, c.t_end - span * 1e-3)
                    };
                    (c.geometry.point_at(t_at) - c.geometry.point_at(t_in)).normalized()
                })
                .collect();
            let pair_surfaces = incident[n].iter().find_map(|&h| {
                if let EdgeSource::Shared { face_a, face_b } = graph.half_edges[h].source {
                    let sa = &brep_a.faces.get(face_a)?.surface.as_ref()?.geometry;
                    let sb = &brep_b.faces.get(face_b)?.surface.as_ref()?.geometry;
                    Some((sa, sb))
                } else {
                    None
                }
            });
            let far_trace = popkern_core::trace::trace_enabled("weld-far");
            if far_trace {
                eprintln!(
                    "  [weld-far?] node {n} ({:.4},{:.4},{:.4}) rec_scale {rec_scale:.3} far_cap {far_cap:.3} dirs {} pair {}",
                    p.x,
                    p.y,
                    p.z,
                    outward_dirs.len(),
                    pair_surfaces.is_some()
                );
            }
            if let Some((surf_a, surf_b)) = pair_surfaces {
                let mut cands: Vec<(f64, usize)> = Vec::new();
                for (m, inc) in incident.iter().enumerate() {
                    if m == n || inc.is_empty() || !inc.iter().any(|&h| !is_shared(h)) {
                        continue;
                    }
                    let q = graph.nodes[m].position;
                    let d = q.distance_to(&p);
                    // Anything the unconditional tier owns is done;
                    // the tangent-ray annulus overlaps this range but
                    // ran first — the far tier only sees its refusals.
                    if d < band || d > far_cap {
                        continue;
                    }
                    let Some(dir_to_m) = (q - p).normalized() else {
                        continue;
                    };
                    let cone = outward_dirs
                        .iter()
                        .map(|o| o.dot(&dir_to_m))
                        .fold(f64::NEG_INFINITY, f64::max);
                    let gate = (d * 0.02).max(lin * 100.0);
                    let ca = on_carrier(surf_a, &q, gate);
                    let cb = on_carrier(surf_b, &q, gate);
                    if far_trace {
                        eprintln!(
                            "  [weld-far?]   cand {m} ({:.4},{:.4},{:.4}) d {d:.3e} cone {cone:.3} on_a {ca} on_b {cb}",
                            q.x, q.y, q.z
                        );
                    }
                    if cone <= 0.9 || !ca || !cb {
                        continue;
                    }
                    cands.push((d, m));
                }
                cands.sort_by(|x, y| x.0.total_cmp(&y.0));
                if let [first, rest @ ..] = cands.as_slice()
                    && rest.first().is_none_or(|s| s.0 > first.0 * 2.0)
                {
                    best = Some((first.0, first.1));
                    if popkern_core::trace::trace_enabled("boolean") {
                        eprintln!("  [weld-far] node {n} certified at d={:.2e}", first.0);
                    }
                }
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
    // SPLIT-JUNCTION COALESCING: a record endpoint landing fit-slop
    // off its true boundary junction can also SPLIT the boundary edge
    // at its own sloppy position — the graph then carries TWO
    // junction nodes a sliver apart on one boundary edge, and on the
    // OTHER face's chart (where the boundary edge does not exist)
    // the ring gaps and the graze cascade unwalks it (measured,
    // PK-464: the chamfered nut's crest arc minted a junction 7.1e-3
    // from the flank's exact one). Weld the record-lighter node onto
    // its DIRECTLY-ADJACENT junction sibling when the record's
    // tangent ray points at it; the boundary sliver between them
    // collapses to zero length for the degenerate cleanup.
    for (n, inc_n) in incident.iter().enumerate().take(n_nodes) {
        if inc_n.is_empty() {
            continue;
        }
        let records_n: Vec<usize> = inc_n.iter().copied().filter(|&h| is_shared(h)).collect();
        let boundary_n: Vec<usize> = inc_n.iter().copied().filter(|&h| !is_shared(h)).collect();
        if records_n.len() != 1 || boundary_n.is_empty() {
            continue;
        }
        let p = graph.nodes[n].position;
        // Directly-adjacent nodes along boundary half-edges.
        for &bh in &boundary_n {
            let m = graph.half_edges[graph.half_edges[bh].twin].origin;
            if m == n {
                continue;
            }
            let d = graph.nodes[m].position.distance_to(&p);
            // Never weld across an existing record connection: two
            // junctions JOINED by a record piece are distinct by
            // construction (the tuned exact-connector class — the
            // m6 nut's mouth-chamfer chain carries real 1.5e-3
            // connectors between junctions 1.0e-3 apart; welding
            // their ends collapsed the connector into a self-loop,
            // the stub cleanup retired it, and the chain gapped into
            // a graze cascade). A weld target must be an UNCONNECTED
            // duplicate of this node.
            let connected = incident[n]
                .iter()
                .any(|&h| graph.half_edges[graph.half_edges[h].twin].origin == m);
            if connected {
                continue;
            }
            let m_has_record = graph
                .half_edges
                .iter()
                .any(|hd| hd.origin == m && matches!(hd.source, EdgeSource::Shared { .. }));
            // Verified by the RECORD's own tangent ray — the boundary
            // hes trivially point along the boundary at any sibling.
            let record_points_at_m = {
                use popkern_core::geometry::curve::CurveEval as _;
                let rh = records_n[0];
                let he = &graph.half_edges[rh];
                let q = graph.nodes[graph.half_edges[he.twin].origin].position;
                let tangent_from = |outward: popkern_core::geometry::Vec3| {
                    let to_m = graph.nodes[m].position - p;
                    let along = to_m.dot(&outward);
                    let perp = (to_m - outward * along).length();
                    along >= -band && perp < band
                };
                match he
                    .curve
                    .as_ref()
                    .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                {
                    Some(c) => {
                        let d_start = c.geometry.point_at(c.t_start).distance_to(&p);
                        let d_end = c.geometry.point_at(c.t_end).distance_to(&p);
                        let span = c.t_end - c.t_start;
                        let (t_at, t_in) = if d_start <= d_end {
                            (c.t_start, c.t_start + span * 1e-3)
                        } else {
                            (c.t_end, c.t_end - span * 1e-3)
                        };
                        let dvec = c.geometry.point_at(t_at) - c.geometry.point_at(t_in);
                        dvec.normalized().map(tangent_from).unwrap_or(false)
                    }
                    None => (p - q).normalized().map(tangent_from).unwrap_or(false),
                }
            };
            if d > band * 100.0 || !m_has_record || !record_points_at_m {
                continue;
            }
            if popkern_core::trace::trace_enabled("boolean") {
                let q = graph.nodes[m].position;
                eprintln!(
                    "  [weld-junction] node {n} ({:.6},{:.6},{:.6}) -> {m} ({:.6},{:.6},{:.6}) d={d:.2e}",
                    p.x, p.y, p.z, q.x, q.y, q.z
                );
            }
            welds.push((n, m));
            break;
        }
    }
    let mut max_weld_d = 0.0_f64;
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
        max_weld_d = max_weld_d.max(d_weld);
        // Move EVERY half-edge that originates at the welded node —
        // scanned from the half-edge table, not the node's outgoing
        // list: pieces minted by splits/bridges can carry an origin
        // the outgoing bookkeeping never recorded, and a weld that
        // trusts the list leaves them pointing at the emptied node,
        // which then reads as a degree-1 dead end and roots a graze
        // cascade that unwalks the whole chain (measured, PK-465:
        // every welded mouth junction re-appeared as a cascade root
        // and the coil flank lost its helix column).
        let moved: Vec<usize> = graph
            .half_edges
            .iter()
            .enumerate()
            .filter(|(_, he)| he.origin == n)
            .map(|(i, _)| i)
            .collect();
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
                    if near_d >= snap_band {
                        continue;
                    }
                    // LARGE welds re-author the curve by interpolation
                    // through on-curve samples with the corrected end
                    // sample, in the original parameter domain — a
                    // bare end-CP move shifts the whole final span by
                    // the endpoint delta, whose direction rotates off
                    // the local tangent along the span (measured
                    // 4.0e-5 off-surface on a 6.2e-4 move; a 3.7e-2
                    // weld would gouge visibly). Small welds keep the
                    // cheap CP snap.
                    if d_weld > band {
                        let n_s = 17usize;
                        let (t0, t1) = (c.t_start, c.t_end);
                        let mut pts: Vec<popkern_core::geometry::Point3> = (0..n_s)
                            .map(|k| {
                                c.geometry
                                    .point_at(t0 + (t1 - t0) * k as f64 / (n_s - 1) as f64)
                            })
                            .collect();
                        let idx = if near_is_t_start { 0 } else { n_s - 1 };
                        pts[idx] = new_pos;
                        let params: Vec<f64> = (0..n_s)
                            .map(|k| t0 + (t1 - t0) * k as f64 / (n_s - 1) as f64)
                            .collect();
                        if let Ok(fitted) =
                            popkern_core::geometry::fitting::fit_curve(&pts, &params, 3, n_s)
                        {
                            c.geometry =
                                popkern_core::geometry::curve::CurveGeometry::Nurbs(fitted);
                            continue;
                        }
                    }
                    if let popkern_core::geometry::curve::CurveGeometry::Nurbs(nc) = &mut c.geometry
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
    graph_weld_certify::retire_self_loop_stubs(graph, max_weld_d, band);
}
