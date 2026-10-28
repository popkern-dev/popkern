//! Reading a truncated apex corner's served cluster (PK-997): the
//! vertex standing for the input corner, the stale vertices around it
//! and the fan they collapse to.

use super::super::mitre_apex_corner::{Plane, destination, outgoing, plane_distance};
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};
use std::collections::BTreeSet;

/// A candidate corner of the input: its position, one chamfered edge
/// (the failure's subject) and its shortest edge's length (the
/// corner's own scale).
pub(super) struct Corner {
    pub(super) position: Point3,
    pub(super) edge: EdgeHandle,
    pub(super) shortest_edge: f64,
}

/// The served vertex of `corner`: the original vertex when it is still
/// there, else the off-plane vertex nearest it (an earlier settle can
/// merge the original into a stale neighbour), within half the
/// corner's shortest edge so a cluster never reaches another corner.
pub(super) fn corner_vertex(out: &BRep, corner: &Corner, gate: f64) -> Option<VertexHandle> {
    if let Some((h, _)) = out
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&corner.position) <= gate)
    {
        return Some(h);
    }
    out.vertices
        .iter()
        .map(|(h, v)| (h, v.position.distance_to(&corner.position)))
        .filter(|&(h, d)| {
            d < 0.5 * corner.shortest_edge && off_plane(out, h, gate).unwrap_or(false)
        })
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(h, _)| h)
}

/// The served fan at a corner: outgoing half-edges in fan order
/// (`twin.next`) and the plane of the face left of each. Edge `k`
/// lies between planes `k` and `k + 1`.
pub(super) struct Fan {
    pub(super) vertex: VertexHandle,
    pub(super) half_edges: Vec<HalfEdgeHandle>,
    pub(super) planes: Vec<Plane>,
    /// The far end of each edge.
    pub(super) leaves: Vec<Point3>,
}

pub(super) fn plane_of(surface: &SurfaceGeometry) -> Option<Plane> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

/// The fan of `vertex`, when every face around it is planar and
/// appears once, and every edge is a line.
pub(super) fn read_fan(out: &BRep, vertex: VertexHandle) -> Option<Fan> {
    let first = outgoing(out, vertex).into_iter().next()?;
    let mut half_edges = vec![first];
    loop {
        let last = *half_edges.last()?;
        let twin = out.half_edges.get(last)?.twin;
        let next = out.half_edges.get(twin)?.next;
        if next == first {
            break;
        }
        if half_edges.len() > out.half_edges.len() {
            return None;
        }
        half_edges.push(next);
    }
    if half_edges.len() != outgoing(out, vertex).len() || half_edges.len() < 3 {
        return None;
    }
    let mut faces = BTreeSet::new();
    let mut planes = Vec::new();
    let mut leaves = Vec::new();
    for &he in &half_edges {
        let h = out.half_edges.get(he)?;
        let face = out.loops.get(h.loop_)?.face;
        if !faces.insert(face) {
            return None;
        }
        planes.push(plane_of(&out.faces.get(face)?.surface.as_ref()?.geometry)?);
        let curve = out.edges.get(h.edge)?.curve.as_ref()?;
        if !matches!(curve.geometry, CurveGeometry::Line { .. }) {
            return None;
        }
        leaves.push(out.vertices.get(destination(out, he)?)?.position);
    }
    Some(Fan {
        vertex,
        half_edges,
        planes,
        leaves,
    })
}

/// Whether `vertex` lies off one of its face planes (`None`: a face
/// around it is not planar).
pub(super) fn off_plane(out: &BRep, vertex: VertexHandle, gate: f64) -> Option<bool> {
    let position = out.vertices.get(vertex)?.position;
    let mut off = false;
    for he in outgoing(out, vertex) {
        let face = out.loops.get(out.half_edges.get(he)?.loop_)?.face;
        let plane = plane_of(&out.faces.get(face)?.surface.as_ref()?.geometry)?;
        off |= plane_distance(&plane, &position) > gate;
    }
    Some(off)
}

/// `vertex` and the vertices off their face planes connected to it
/// through each other: the per-edge pipeline's stale corner (a sharp
/// edge between two bevels gets its end at the shallower crossing,
/// off the deeper bevel). `None` when a face around it is curved.
pub(super) fn stale_cluster(
    out: &BRep,
    vertex: VertexHandle,
    gate: f64,
) -> Option<Vec<VertexHandle>> {
    off_plane(out, vertex, gate)?;
    let mut cluster = vec![vertex];
    let mut frontier = vec![vertex];
    while let Some(at) = frontier.pop() {
        for he in outgoing(out, at) {
            let Some(next) = destination(out, he) else {
                continue;
            };
            if cluster.contains(&next) || !off_plane(out, next, gate).unwrap_or(false) {
                continue;
            }
            cluster.push(next);
            frontier.push(next);
        }
    }
    Some(cluster)
}
