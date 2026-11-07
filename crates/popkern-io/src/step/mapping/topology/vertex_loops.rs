//! A face's `VERTEX_LOOP` bounds (PK-1116): a vertex on the face's
//! interior becomes a vertex loop of the face, holding the file's vertex
//! (shared with the edges using it, else made at its point); a vertex
//! on the face's trim (a cone's apex bound, which the face's other loop
//! already closes) carries no topology of its own and is skipped.

use popkern_core::topology::brep::BRep;
use popkern_core::topology::{FaceHandle, Loop, Vertex, VertexHandle};
use popkern_core::validation::FaceUvContainment;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::topology::VertexPoint;

/// A face's vertex loop bound awaiting its built face: the face, the
/// bound's `VERTEX_POINT`, and its vertex when the face's edges made it.
pub(super) struct PendingVertexLoop {
    pub face: FaceHandle,
    pub vertex_point: EntityId,
    pub vertex: Option<VertexHandle>,
}

/// Attach each pending vertex loop standing on its face's interior.
pub(super) fn attach_vertex_loops(
    brep: &mut BRep,
    pending: Vec<PendingVertexLoop>,
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
) {
    for p in pending {
        let point = match p.vertex {
            Some(v) => brep.vertices.get(v).map(|v| v.position),
            None => graph
                .get(p.vertex_point)
                .and_then(|e| VertexPoint::from_params(&e.params).ok())
                .and_then(|vp| geo.resolve_point(vp.point, log)),
        };
        let Some(point) = point else {
            log.warn(
                Phase::MapTopology,
                Some(p.vertex_point),
                "VERTEX_LOOP vertex unread",
            );
            continue;
        };
        let tolerance = brep.tolerance.linear();
        let inside = brep
            .face_point_trim_containment(p.face, &point)
            .is_some_and(|c| {
                c.distance <= tolerance && c.containment == FaceUvContainment::Interior
            });
        if !inside {
            log.warn(
                Phase::MapTopology,
                Some(p.vertex_point),
                "VERTEX_LOOP bound on the face's trim (degenerate apex) skipped",
            );
            continue;
        }
        let vertex = p.vertex.unwrap_or_else(|| {
            brep.vertices.insert(Vertex {
                position: point,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            })
        });
        brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: p.face,
            vertex: Some(vertex),
        });
    }
}
