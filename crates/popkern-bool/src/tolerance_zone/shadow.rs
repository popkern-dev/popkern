//! Shadow comparison of candidate resolution against the intersection
//! graph's node identity (plan Stage 2 gate).
//!
//! Each candidate is mapped to the graph node the graph build gave it
//! (the node within the snap tolerance of its position; a record end,
//! the node within the wide adoption radius as the build adopts it).
//! Two partitions of the same candidates then disagree two ways:
//! - **split**: a resolved vertex whose members sit on several graph
//!   nodes (the graph keeps apart what resolution identifies);
//! - **merged**: a graph node holding members of several resolved
//!   vertices (the graph identifies what a compatibility rule, or the
//!   growth cap, keeps apart).
//!
//! The counts go to `POPKERN_COMPLEX_SHADOW_LOG` beside the complex
//! shadow events, one line per boolean.

use std::collections::{BTreeMap, BTreeSet};

use super::candidates::{CandidateVertex, SourceKind};
use super::resolve::{Resolution, ResolutionDecline};
use crate::merge::graph::IntersectionGraph;

/// One boolean's comparison.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct ResolutionShadow {
    /// Candidates collected.
    pub candidates: usize,
    /// Proximity components with more than one candidate.
    pub components: usize,
    /// Resolved vertices.
    pub vertices: usize,
    /// Pairs kept apart by a compatibility rule.
    pub separations: usize,
    /// Resolved vertices whose members sit on several graph nodes.
    pub split: usize,
    /// Graph nodes holding members of several resolved vertices.
    pub merged: usize,
    /// The widest such merge: the largest distance between the anchors
    /// of two resolved vertices one graph node holds.
    pub merged_worst: f64,
    /// Candidates with no graph node, by source: operand vertices
    /// (outside a scoped build's active faces), edge-face hits, record
    /// ends (a record the graph build filtered or re-cut).
    pub unmapped: [usize; 3],
    /// Resolved vertices with a T-junction member (on another edge's
    /// or record's interior): paves the table commits.
    pub t_junctions: usize,
    /// Of those, the ones some graph edge still passes through
    /// unsplit: a junction the graph's split passes missed.
    pub unsplit: usize,
    /// Record ends whose vertex lies on no boundary: no operand vertex
    /// or edge and no other record among its group's incidences, and not
    /// the record's own other end. Each is an `OpenIntersectionLoop`
    /// (plan §4.3a endpoint closure): today's gap closure stitches it
    /// after the walk.
    pub open_ends: usize,
    /// Vertices joining ends of two records of one face pair away from
    /// every boundary: two branches of one intersection meeting in the
    /// faces' interiors (a branch point, or two branches passing within
    /// tolerance, which the branch rule must keep apart).
    pub branch_meets: usize,
    /// Records lying along an operand edge for their whole length
    /// (both ends and the middle on one edge): common blocks, one edge
    /// that today enters the graph twice (the record and the edge).
    pub common_blocks: usize,
    /// Records whose ends resolved into one vertex and that lie within
    /// it: micro edges, collapsed into the vertex by construction.
    pub micro_edges: usize,
    /// Records whose ends resolved into one vertex and that run out
    /// and back along one segment: spurs, which bound nothing.
    pub spurs: usize,
    /// Records whose ends resolved into one vertex but that stray from
    /// it: loops through the vertex, with the largest extent among them.
    pub vertex_loops: (usize, f64),
    /// Intersection nodes of the finished graph (after its split
    /// passes) that carry a record edge and that no candidate maps to,
    /// with the largest distance from one to its nearest resolved
    /// vertex: the junctions the passes find that the table does not
    /// predict, which must be none before the table's incidences
    /// replace the passes. Measured (suite 22, table on): 21 nodes in
    /// 12 booleans, every one in a scene with common blocks (same-domain
    /// overlaps, plan §4.4); none in the scenes without.
    pub unpredicted: (usize, f64),
    /// (node, record edge) junctions of the finished graph the graph's
    /// own splitting pass would still split: with the table deciding,
    /// the junctions its paves missed.
    pub junctions_unpaved: usize,
    /// (record-end node, operand edge) junctions of the finished graph
    /// the graph's operand-edge pass would still split: with the table
    /// deciding, the operand-edge junctions its paves missed.
    pub operand_junctions_unpaved: usize,
    /// The resolution's decline, if it declined.
    pub declined: Option<String>,
}

/// Whether record end `i` resolved onto a boundary.
pub(crate) fn end_closes(cands: &[CandidateVertex], res: &Resolution, i: usize) -> bool {
    use super::candidates::Incidence;
    let Some(&Incidence::Record { record: own, .. }) = cands[i].incidence.first() else {
        return true;
    };
    let group = &res.groups[res.group_of[i]];
    let own_ends = group
        .members
        .iter()
        .filter(|&&m| {
            matches!(cands[m].incidence.first(), Some(Incidence::Record { record, end: Some(_), .. }) if *record == own)
        })
        .count();
    own_ends > 1
        || group.members.iter().any(|&m| {
            cands[m].incidence.iter().any(|inc| match inc {
                Incidence::Vertex { .. } | Incidence::Edge { .. } => true,
                Incidence::Record { record, .. } => *record != own,
                Incidence::Face { .. } => false,
            })
        })
}

/// For each `(p, reach)` of `queries`, whether a graph edge's curve
/// passes through `p` (within `reach`) at its interior, clear of both
/// its end nodes. One sweep: the queries sorted by x, each graph edge
/// tests only those in its conservative box (a per-query scan of
/// every edge, rejecting by distance from the edge's midpoint, cost a
/// 100-knuckle hinge's 893 junctions against 2,310 edges 0.55 s: the
/// knuckle edges run the hinge's length).
fn passes_unsplit(
    graph: &IntersectionGraph,
    queries: &[(popkern_core::geometry::Point3, f64)],
) -> Vec<bool> {
    use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};
    let mut out = vec![false; queries.len()];
    let reach_max = queries.iter().map(|q| q.1).fold(0.0_f64, f64::max);
    let mut by_x: Vec<(f64, usize)> = queries
        .iter()
        .enumerate()
        .map(|(i, q)| (q.0.x, i))
        .collect();
    by_x.sort_by(|p, q| p.0.total_cmp(&q.0));
    for (i, he) in graph.half_edges.iter().enumerate() {
        if he.twin < i {
            continue;
        }
        let Some(c) = he.curve.as_ref() else {
            continue;
        };
        let a = graph.nodes[he.origin].position;
        let b = graph.nodes[graph.half_edges[he.twin].origin].position;
        // The curve's conservative box: a query outside it grown by
        // the reach cannot lie on the curve.
        let mut hull = vec![c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)];
        crate::face_box::conservative_curve_extent(c, &mut hull);
        let Some(bounds) = super::region::Bounds::of(&hull, reach_max) else {
            continue;
        };
        let lo = by_x.partition_point(|(x, _)| *x < bounds.lo[0]);
        let hi = by_x.partition_point(|(x, _)| *x <= bounds.hi[0]);
        for &(_, qi) in &by_x[lo..hi] {
            let (p, reach) = queries[qi];
            if out[qi]
                || !bounds.holds(&p)
                || a.distance_to(&p) < reach
                || b.distance_to(&p) < reach
            {
                continue;
            }
            let t = closest_point_in_window(&c.geometry, &p, c.t_start, c.t_end);
            if c.geometry.point_at(t).distance_to(&p) < reach {
                out[qi] = true;
            }
        }
    }
    out
}

/// Compare `resolution` of `cands` against `graph`'s nodes.
pub(crate) fn compare(
    graph: &IntersectionGraph,
    cands: &[CandidateVertex],
    resolution: &Result<Resolution, ResolutionDecline>,
    snap_tol: f64,
    adoption_tol: f64,
) -> ResolutionShadow {
    let mut out = ResolutionShadow {
        candidates: cands.len(),
        ..ResolutionShadow::default()
    };
    let res = match resolution {
        Ok(r) => r,
        Err(d) => {
            out.declined = Some(format!("{d:?}"));
            return out;
        }
    };
    out.components = res.components;
    out.vertices = res.groups.len();
    out.separations = res.separations.len();
    let node_of: Vec<Option<usize>> = cands
        .iter()
        .map(|c| {
            graph.find_node(&c.position, snap_tol).or_else(|| {
                (c.source == SourceKind::RecordEndpoint)
                    .then(|| graph.find_node(&c.position, adoption_tol))
                    .flatten()
            })
        })
        .collect();
    for (c, n) in cands.iter().zip(&node_of) {
        if n.is_none() {
            out.unmapped[c.source as usize] += 1;
        }
    }
    for g in &res.groups {
        let nodes: BTreeSet<usize> = g.members.iter().filter_map(|&m| node_of[m]).collect();
        if nodes.len() > 1 {
            out.split += 1;
        }
    }
    let mut groups_at: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    for (c, node) in node_of.iter().enumerate() {
        if let Some(n) = node {
            groups_at.entry(*n).or_default().insert(res.group_of[c]);
        }
    }
    out.merged = groups_at.values().filter(|gs| gs.len() > 1).count();
    for (n, node) in graph.nodes.iter().enumerate() {
        // A node splitting one operand edge with no record at it (a hit
        // outside the hit face's trim) is a no-op split that
        // post-processing dissolves; only junctions count.
        let junction = graph.half_edges.iter().any(|h| {
            h.origin == n && matches!(h.source, crate::merge::graph::EdgeSource::Shared { .. })
        });
        if groups_at.contains_key(&n)
            || !junction
            || !matches!(node.source, crate::merge::graph::NodeSource::Intersection)
        {
            continue;
        }
        let nearest = res
            .groups
            .iter()
            .map(|g| cands[g.anchor].position.distance_to(&node.position))
            .fold(f64::INFINITY, f64::min);
        out.unpredicted.0 += 1;
        out.unpredicted.1 = out.unpredicted.1.max(nearest);
    }
    for gs in groups_at.values().filter(|gs| gs.len() > 1) {
        let anchors: Vec<_> = gs
            .iter()
            .map(|&g| cands[res.groups[g].anchor].position)
            .collect();
        for (i, p) in anchors.iter().enumerate() {
            for q in &anchors[i + 1..] {
                out.merged_worst = out.merged_worst.max(p.distance_to(q));
            }
        }
    }
    out.branch_meets = res
        .groups
        .iter()
        .filter(|g| {
            use super::candidates::Incidence;
            let on_boundary = g.members.iter().any(|&m| {
                cands[m]
                    .incidence
                    .iter()
                    .any(|i| matches!(i, Incidence::Vertex { .. } | Incidence::Edge { .. }))
            });
            let mut pairs = std::collections::BTreeMap::<(u64, u64), BTreeSet<usize>>::new();
            for &m in &g.members {
                let inc = &cands[m].incidence;
                let (
                    Some(Incidence::Record { record, .. }),
                    Some(Incidence::Face { face: fa, .. }),
                    Some(Incidence::Face { face: fb, .. }),
                ) = (inc.first(), inc.get(1), inc.get(2))
                else {
                    continue;
                };
                use slotmap::Key;
                pairs
                    .entry((fa.data().as_ffi(), fb.data().as_ffi()))
                    .or_default()
                    .insert(*record);
            }
            !on_boundary && pairs.values().any(|rs| rs.len() > 1)
        })
        .count();
    out.open_ends = (0..cands.len())
        .filter(|&i| cands[i].source == SourceKind::RecordEndpoint && !end_closes(cands, res, i))
        .count();
    let junctions: Vec<_> = res
        .groups
        .iter()
        .filter(|g| {
            g.members
                .iter()
                .any(|&m| super::incidence::is_t_junction(&cands[m]))
        })
        .map(|g| (cands[g.anchor].position, g.tolerance))
        .collect();
    out.t_junctions = junctions.len();
    out.unsplit = passes_unsplit(graph, &junctions)
        .into_iter()
        .filter(|&u| u)
        .count();
    out
}

/// Build the boolean's vertex table and record its comparison against
/// `graph` (built with the cellular snap tolerance, 2 × linear, and,
/// unless the table decided them, record ends adopted within 20 ×
/// linear).
pub(crate) fn run(
    graph: &IntersectionGraph,
    operands: (&popkern_core::topology::BRep, &popkern_core::topology::BRep),
    ix: &crate::intersect::IntersectionResult,
    records: &[crate::intersect::IntersectionRecord],
    tolerance: &popkern_core::tolerance::Tolerance,
) {
    let linear = tolerance.linear();
    let table = super::table::build(operands, ix, records, linear);
    let mut shadow = compare(
        graph,
        &table.cands,
        &table.resolution,
        2.0 * linear,
        20.0 * linear,
    );
    shadow.junctions_unpaved =
        crate::merge::graph_split_at_nodes::unsplit_junctions(graph, 2.0 * linear);
    shadow.operand_junctions_unpaved =
        crate::merge::graph_split_shared::unsplit_operand_junctions(graph, 2.0 * linear);
    if let Ok(res) = &table.resolution {
        shadow.common_blocks = super::common_blocks::common_block_records(
            &table.cands,
            res,
            [operands.0, operands.1],
            records,
            linear,
        )
        .len();
        for c in super::consequences::collapsed_ends(&table.cands, res, records) {
            match c.kind() {
                super::consequences::Collapse::MicroEdge => shadow.micro_edges += 1,
                super::consequences::Collapse::Spur => shadow.spurs += 1,
                super::consequences::Collapse::Loop => {
                    shadow.vertex_loops.0 += 1;
                    shadow.vertex_loops.1 = shadow.vertex_loops.1.max(c.extent);
                }
            }
        }
    }
    record(&shadow);
}

/// Append `shadow` to `POPKERN_COMPLEX_SHADOW_LOG`, when set.
pub(crate) fn record(shadow: &ResolutionShadow) {
    if popkern_core::trace::trace_enabled("complex") {
        eprintln!("[complex] resolution shadow {shadow:?}");
    }
    crate::complex_shadow::append_log_line(&format!("{shadow:?}"));
}
