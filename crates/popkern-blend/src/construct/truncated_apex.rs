//! Truncated apex corners for planar chamfers (PK-997).
//!
//! **The class.** A vertex `V` where at least one edge is chamfered
//! with a planar bevel and at least one edge stays sharp, every face at
//! `V` planar, and the bevels cut the WHOLE cone at `V`: a bevel is
//! parallel to its own edge, so it crosses every other edge of a convex
//! vertex, not only its two neighbours (a square pyramid's apex with
//! one slant edge chamfered loses the apex; the bevel meets the far
//! slant edge too). The per-edge pipeline builds the exact faces
//! around `V` and every edge leaving the corner to its exact far end,
//! but leaves a stale cluster at the corner: `V` itself, off the bevel
//! planes (measured on 4- to 6-gonal pyramids: off 1 to `m` of its `m`
//! incident planes, by 0.04–0.06 at d = 0.1), plus, where a sharp edge
//! lies between two bevels of unequal reach, that edge's end at the
//! shallower crossing, off the deeper bevel (irregular cones; an
//! earlier settle may already have merged `V` into it).
//!
//! **The exact corner** is the half-space intersection of the `m`
//! planes around the collapsed cluster: for `m = 3` their common point
//! (forced: three planar faces meet in one vertex); for `m ≥ 4` a
//! ridge tree whose leaves are the far ends of the corner's edges,
//! each on two consecutive fan planes, and whose interior nodes are
//! feasible triple points (`plane_ridge_tree`). A 5-gonal pyramid
//! with one slant edge chamfered: the bevel meets the two far slant
//! edges at two triple points joined by a ridge on the far face.
//!
//! **Construction, by Euler operators.** One `kve` per edge inside the
//! stale cluster collapses it into `V` (those edges are pipeline
//! pieces of the corner, not input edges); `V` becomes the tree's
//! root; every interior ridge on planes `(i, j)` separates the fan
//! edges `i..j` (those between the two planes) from the rest, so it is
//! one `split_vertex` moving that contiguous sector to the new node.
//! Then every tree vertex takes its node position and its edges are
//! re-authored as exact segments. Faces keep their handles and the
//! edges leaving the corner are re-homed, not re-made, so names are
//! untouched.
//!
//! **Invariants and postconditions:** every tree vertex lies on all its
//! incident planes and every edge keeps its direction (no node beyond
//! the far end of an edge it carries). **Fallback:** a cluster that does
//! not read as above (a curved face, a face twice in the fan, a
//! non-convex fan of four or more planes, a cluster `kve` cannot
//! collapse, a tree with other leaves, a node that overshoots) is
//! restored and left for the serve gate's structured check-9 decline;
//! an operator failure restores the corner's snapshot and declines the
//! operation naming the vertex. Bevels that consume a whole face
//! decline earlier, in the face trims, and are served by growth from a
//! smaller size (`bevel_growth`, PK-1001). **Locality:** the
//! cluster is grown only through off-plane vertices, seeded at the
//! corner (or the off-plane vertex nearest it within half its shortest
//! edge). **Tolerance:** the caller's blend tolerance; positions gate
//! at `linear × 10`.

#[path = "truncated_apex_fan.rs"]
mod fan;

use super::EdgeGeometry;
use super::mitre_apex_corner::{destination, outgoing, plane_distance};
use super::mitre_apex_scaffold::three_plane_point;
use super::plane_ridge_tree::{RidgeNode, ridge_tree, split_fan};
use crate::error::BlendFailure;
use fan::{Corner, Fan, corner_vertex, off_plane, plane_of, read_fan, stale_cluster};
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};
use std::collections::BTreeSet;

/// Outcome of the stage.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct TruncatedApexReport {
    /// Corners rebuilt exactly.
    pub corners_settled: usize,
}

/// Settle every truncated apex corner of `input` in `out`.
pub(crate) fn settle_truncated_apexes(
    input: &BRep,
    out: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<TruncatedApexReport, (EdgeHandle, BlendFailure)> {
    let trace = popkern_core::trace::trace_enabled("truncated-apex");
    let gate = (tolerance.linear() * 10.0).max(1e-9);
    let mut report = TruncatedApexReport::default();
    for corner in candidate_corners(input, edge_geometries) {
        let Some(vertex) = corner_vertex(out, &corner, gate) else {
            continue;
        };
        let Some(cluster) = stale_cluster(out, vertex, gate) else {
            continue;
        };
        if cluster.len() == 1 && !off_plane(out, vertex, gate).unwrap_or(false) {
            continue; // the planes concur at V: already exact
        }
        let snapshot = out.clone();
        match settle_one(out, vertex, &cluster, gate) {
            Ok(Ok((planes, nodes))) => {
                report.corners_settled += 1;
                if trace {
                    eprintln!(
                        "[truncated-apex] corner {vertex:?}: cluster of {}, {planes} planes, \
                         {nodes} tree nodes",
                        cluster.len()
                    );
                }
            }
            Ok(Err(why)) => {
                *out = snapshot;
                if trace {
                    eprintln!("[truncated-apex] corner {vertex:?}: not recognised: {why}");
                }
            }
            Err(reason) => {
                *out = snapshot;
                return Err((
                    corner.edge,
                    BlendFailure::InvalidTopology {
                        reason: format!("truncated apex corner at {vertex:?}: {reason}"),
                        location: Some(corner.position),
                    },
                ));
            }
        }
    }
    Ok(report)
}

/// Vertices of `input` with at least one edge chamfered by a planar
/// bevel and at least one edge not chamfered (all-chamfered corners
/// are the mitre apex stage's).
fn candidate_corners(input: &BRep, edge_geometries: &[EdgeGeometry]) -> Vec<Corner> {
    let mut out = Vec::new();
    for (vh, vertex) in input.vertices.iter() {
        let edges: BTreeSet<EdgeHandle> = outgoing(input, vh)
            .into_iter()
            .filter_map(|he| input.half_edges.get(he).map(|h| h.edge))
            .collect();
        let bevelled: Vec<EdgeHandle> = edges
            .iter()
            .copied()
            .filter(|e| {
                edge_geometries
                    .iter()
                    .any(|g| g.analysis.edge == *e && plane_of(&g.surface).is_some())
            })
            .collect();
        let chamfered = edges
            .iter()
            .filter(|e| edge_geometries.iter().any(|g| g.analysis.edge == **e))
            .count();
        let shortest_edge = outgoing(input, vh)
            .into_iter()
            .filter_map(|he| destination(input, he))
            .filter_map(|d| input.vertices.get(d))
            .map(|d| d.position.distance_to(&vertex.position))
            .fold(f64::INFINITY, f64::min);
        if let Some(&edge) = bevelled.first()
            && chamfered < edges.len()
        {
            out.push(Corner {
                position: vertex.position,
                edge,
                shortest_edge,
            });
        }
    }
    out
}

/// Rebuild one corner. `Ok(Err)` = the fan does not read as a
/// truncated apex (left alone); `Err` = an operator failed.
fn settle_one(
    out: &mut BRep,
    vertex: VertexHandle,
    cluster: &[VertexHandle],
    gate: f64,
) -> Result<Result<(usize, usize), String>, String> {
    let corner = out
        .vertices
        .get(vertex)
        .map(|v| v.position)
        .ok_or("stale corner vertex")?;
    // Collapse the stale cluster into V: one `kve` per internal edge.
    while let Some(inner) = outgoing(out, vertex)
        .into_iter()
        .find(|&he| destination(out, he).is_some_and(|d| d != vertex && cluster.contains(&d)))
    {
        if let Err(e) = out.kve(inner) {
            return Ok(Err(format!("the stale cluster does not collapse: {e}")));
        }
    }
    let Some(fan) = read_fan(out, vertex) else {
        return Ok(Err(
            "the collapsed fan has a curved face, a repeated face or a curved edge".into(),
        ));
    };
    let fan = &fan;
    let m = fan.planes.len();
    // Edge k lies on planes k and k + 1: its far end must too.
    for (k, leaf) in fan.leaves.iter().enumerate() {
        let (a, b) = (&fan.planes[k], &fan.planes[(k + 1) % m]);
        if plane_distance(a, leaf) > gate || plane_distance(b, leaf) > gate {
            return Ok(Err(format!("edge {k}'s far end is off its face planes")));
        }
    }
    let (nodes, ridges) = if m == 3 {
        let Some(apex) = three_plane_point(&fan.planes[0], &fan.planes[1], &fan.planes[2]) else {
            return Ok(Err("the three face planes share a line".into()));
        };
        let mut nodes: Vec<RidgeNode> = (0..3)
            .map(|k| RidgeNode {
                position: fan.leaves[k],
                planes: [k, (k + 1) % 3].into_iter().collect(),
                leaf: Some(k),
            })
            .collect();
        nodes.push(RidgeNode {
            position: apex,
            planes: (0..3).collect(),
            leaf: None,
        });
        (nodes, Vec::new())
    } else {
        let Some(removed_side) = removed_sides(fan, gate) else {
            return Ok(Err(format!("the {m}-plane fan is not convex")));
        };
        let leaves: Vec<(Point3, (usize, usize))> = fan
            .leaves
            .iter()
            .enumerate()
            .map(|(k, &p)| (p, (k, (k + 1) % m)))
            .collect();
        match ridge_tree(corner, &fan.planes, &removed_side, &leaves, gate) {
            Ok(tree) => tree,
            Err(why) => return Ok(Err(why)),
        }
    };
    // Each leaf's interior neighbour carries that edge: it must lie
    // before the far end (an overshoot needs more topology than a fan).
    let carrier = |k: usize| -> Option<usize> {
        if m == 3 {
            return Some(3);
        }
        ridges.iter().find_map(|&(a, b, _)| {
            if a == k {
                Some(b)
            } else if b == k {
                Some(a)
            } else {
                None
            }
        })
    };
    for k in 0..m {
        let Some(c) = carrier(k) else {
            return Ok(Err(format!("leaf {k} has no tree neighbour")));
        };
        let along = (fan.leaves[k] - corner).dot(&(fan.leaves[k] - nodes[c].position));
        if along <= gate * (fan.leaves[k] - corner).length() {
            return Ok(Err(format!(
                "the tree node carrying edge {k} lies beyond its far end"
            )));
        }
    }
    let interior: Vec<usize> = (m..nodes.len()).collect();
    let positions: Vec<Point3> = nodes.iter().map(|n| n.position).collect();
    let built = match split_fan(out, &fan.half_edges, fan.vertex, &positions, &ridges)? {
        Ok(built) => built,
        Err(why) => return Ok(Err(why)),
    };
    for &k in &interior {
        let v = built[k].ok_or("an interior node was never built")?;
        out.vertices.get_mut(v).ok_or("stale tree vertex")?.position = nodes[k].position;
    }
    for &k in &interior {
        let v = built[k].ok_or("an interior node was never built")?;
        for he in outgoing(out, v) {
            author_segment(out, he, gate / 10.0)?;
        }
        // Postcondition: the vertex lies on every face it now touches
        // (a ridge between the wrong face pair would not).
        if off_plane(out, v, gate) != Some(false) {
            return Ok(Err(format!(
                "tree node {k} is off a face it touches after the splits"
            )));
        }
    }
    Ok(Ok((m, interior.len())))
}

/// The removed side's sign for every fan plane, read from the far ends
/// of the edges not on it (all on the material side of a convex fan),
/// or `None` when they disagree.
fn removed_sides(fan: &Fan, gate: f64) -> Option<Vec<f64>> {
    let m = fan.planes.len();
    let mut sides = Vec::with_capacity(m);
    for (k, (o, n)) in fan.planes.iter().enumerate() {
        // Edges k − 1 and k lie on plane k.
        let mut sign = 0.0f64;
        for (j, leaf) in fan.leaves.iter().enumerate() {
            if j == k || (j + 1) % m == k {
                continue;
            }
            let s = (*leaf - *o).dot(n);
            if s.abs() <= gate {
                return None;
            }
            if sign != 0.0 && sign != s.signum() {
                return None;
            }
            sign = s.signum();
        }
        sides.push(-sign);
    }
    Some(sides)
}

/// Author the edge of `he` as the straight segment between its
/// endpoints (both faces planar, so the edge is their line).
fn author_segment(out: &mut BRep, he: HalfEdgeHandle, tol: f64) -> Result<(), String> {
    let edge = out.half_edges.get(he).ok_or("stale half-edge")?.edge;
    let primary = out.edges.get(edge).ok_or("stale edge")?.half_edge;
    let o = out.half_edges.get(primary).ok_or("stale primary")?.origin;
    let d = destination(out, primary).ok_or("edge without destination")?;
    let (po, pd) = (
        out.vertices.get(o).ok_or("stale origin")?.position,
        out.vertices.get(d).ok_or("stale destination")?.position,
    );
    let segment = CurveOnEdge::line_segment(po, pd, tol).ok_or("degenerate segment")?;
    if let Some(e) = out.edges.get_mut(edge) {
        e.curve = Some(segment);
    }
    Ok(())
}
