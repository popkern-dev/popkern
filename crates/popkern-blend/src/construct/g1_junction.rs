//! Tangent-chain junctions (PK-733): two filleted edges meeting at a
//! vertex with one tangent line there. Their strips meet on the
//! junction's cross-section — the ball at the vertex is the same ball
//! for both — so the junction is a continuation, never a corner.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// Two unit tangents closer than this (sine) are one tangent line.
/// Fitted intersection edges carry ~1e-6 of tangent noise at their
/// split points (measured: an r40 × r20 tee saddle split at a seam
/// crossing, sine 1.2e-6); a genuine corner is orders above 1e-4, and
/// a near-tangent one is still caught by the chain construction's
/// section-agreement check.
const TANGENT_SINE: f64 = 1e-4;

/// The unit tangent of `edge` at its endpoint `vertex` (along the
/// stored curve, or the chord for an edge without one).
pub(crate) fn edge_tangent(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> Option<Vec3> {
    let data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let (from, to) = (he.origin, twin.origin);
    if from == to {
        return None; // a closed edge has no junction end
    }
    let tangent = match data.curve.as_ref() {
        Some(curve) => {
            let t = if vertex == from {
                curve.t_start
            } else {
                curve.t_end
            };
            curve.geometry.tangent_at(t)
        }
        None => brep.vertices.get(to)?.position - brep.vertices.get(from)?.position,
    };
    tangent.normalized()
}

/// Whether filleted edges `a` and `b` share one tangent line at
/// `vertex` (either sense: one runs into the vertex, one out).
pub(crate) fn tangent_at_vertex(
    brep: &BRep,
    a: EdgeHandle,
    b: EdgeHandle,
    vertex: VertexHandle,
) -> bool {
    match (edge_tangent(brep, a, vertex), edge_tangent(brep, b, vertex)) {
        (Some(ta), Some(tb)) => ta.cross(&tb).length() < TANGENT_SINE,
        _ => false,
    }
}

/// The tangent-chain junctions among the fillet edges `ends`
/// (`(edge, start, end)` per edge): vertices where exactly two distinct
/// fillet edges meet on one tangent line. They are continuations, so
/// the corner machinery must not see them as corners.
pub(crate) fn junction_vertices(
    brep: &BRep,
    ends: impl Iterator<Item = (EdgeHandle, VertexHandle, VertexHandle)>,
) -> std::collections::BTreeSet<VertexHandle> {
    let mut by_vertex: std::collections::BTreeMap<VertexHandle, Vec<EdgeHandle>> =
        std::collections::BTreeMap::new();
    for (edge, start, end) in ends {
        by_vertex.entry(start).or_default().push(edge);
        if end != start {
            by_vertex.entry(end).or_default().push(edge);
        }
    }
    by_vertex
        .into_iter()
        .filter(|(v, edges)| {
            edges.len() == 2
                && edges[0] != edges[1]
                && tangent_at_vertex(brep, edges[0], edges[1], *v)
        })
        .map(|(v, _)| v)
        .collect()
}

/// The CLOSED tangent chains among `edge_geometries`: groups of fillet
/// edges whose every endpoint is a tangent-chain junction, joined into
/// cycles (each junction shared by exactly two members). Indices in
/// `edge_geometries`, one group per chain, ascending.
pub(crate) fn closed_tangent_chains(
    brep: &BRep,
    edge_geometries: &[super::EdgeGeometry],
) -> Vec<Vec<usize>> {
    let junctions = junction_vertices(
        brep,
        edge_geometries.iter().map(|eg| {
            (
                eg.analysis.edge,
                eg.analysis.vertex_start,
                eg.analysis.vertex_end,
            )
        }),
    );
    let open: Vec<bool> = edge_geometries
        .iter()
        .map(|eg| {
            eg.analysis.vertex_start != eg.analysis.vertex_end
                && junctions.contains(&eg.analysis.vertex_start)
                && junctions.contains(&eg.analysis.vertex_end)
        })
        .collect();
    let mut seen = vec![false; edge_geometries.len()];
    let mut chains = Vec::new();
    for start in 0..edge_geometries.len() {
        if seen[start] || !open[start] {
            continue;
        }
        let mut group = vec![start];
        seen[start] = true;
        let mut grew = true;
        while grew {
            grew = false;
            for j in 0..edge_geometries.len() {
                if seen[j] || !open[j] {
                    continue;
                }
                let (a, b) = (
                    edge_geometries[j].analysis.vertex_start,
                    edge_geometries[j].analysis.vertex_end,
                );
                let touches = group.iter().any(|&g| {
                    let (c, d) = (
                        edge_geometries[g].analysis.vertex_start,
                        edge_geometries[g].analysis.vertex_end,
                    );
                    a == c || a == d || b == c || b == d
                });
                if touches {
                    seen[j] = true;
                    group.push(j);
                    grew = true;
                }
            }
        }
        group.sort_unstable();
        // Closed: every member endpoint is shared by exactly two members.
        let closed = group.iter().all(|&g| {
            [
                edge_geometries[g].analysis.vertex_start,
                edge_geometries[g].analysis.vertex_end,
            ]
            .iter()
            .all(|v| {
                group
                    .iter()
                    .filter(|&&o| {
                        edge_geometries[o].analysis.vertex_start == *v
                            || edge_geometries[o].analysis.vertex_end == *v
                    })
                    .count()
                    == 2
            })
        });
        if closed && group.len() >= 2 {
            chains.push(group);
        }
    }
    chains
}

/// At a tangent-chain junction `vertex` between fillet edge `fillet`
/// and its filleted neighbour `other`, the anchor is the section
/// contact — the trim's point nearest the junction, interior to the
/// face — not a corner walk `r` along the tangent neighbour, which
/// lands off the section entirely (measured: (27, −15) for the
/// section (27, −12) on a rounded-rectangle rim). `None` when the
/// neighbour does not continue `fillet` tangentially.
pub(crate) fn junction_section_contact(
    brep: &BRep,
    fillet: EdgeHandle,
    other: EdgeHandle,
    vertex: VertexHandle,
    trim: &popkern_core::geometry::NurbsCurve,
) -> Option<popkern_core::geometry::Point3> {
    if !tangent_at_vertex(brep, fillet, other, vertex) {
        return None;
    }
    let p = brep.vertices.get(vertex)?.position;
    let s = trim.closest_point(&p)?;
    Some(trim.point_at(s))
}
