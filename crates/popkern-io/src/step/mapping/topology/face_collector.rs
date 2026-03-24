//! Face data collection: walks STEP topology to gather unique vertices, edges,
//! and per-face bound/loop structure before BRep construction.

use std::collections::HashMap;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::topology::{
    AdvancedFace, EdgeCurve, EdgeLoop, FaceBound, FaceOuterBound, OrientedEdge,
};

/// Pre-collected face topology data for a shell.
pub(super) struct FaceData {
    /// All unique vertex entity IDs in discovery order.
    pub vertex_ids: Vec<EntityId>,
    /// All unique edge curve entity IDs in discovery order.
    pub edge_curve_ids: Vec<EntityId>,
    /// Per-face topology information.
    pub faces: Vec<FaceInfo>,
}

/// Topology of a single face.
pub(super) struct FaceInfo {
    /// Entity ID of the ADVANCED_FACE.
    pub face_id: EntityId,
    /// Entity ID of the face surface geometry.
    pub surface_id: EntityId,
    /// Whether the face normal agrees with the surface normal.
    pub same_sense: bool,
    /// Face bounds (loops), with outer bound first if present.
    pub bounds: Vec<BoundInfo>,
}

/// Topology of a single face bound (loop).
pub(super) struct BoundInfo {
    /// Whether this is the outer bound of the face.
    pub is_outer: bool,
    /// Oriented edges composing this loop.
    pub oriented_edges: Vec<OrientedEdgeInfo>,
}

/// Resolved oriented edge reference.
pub(super) struct OrientedEdgeInfo {
    /// Entity ID of the underlying EDGE_CURVE.
    pub edge_curve_id: EntityId,
    /// Effective orientation (composed from ORIENTED_EDGE and FACE_BOUND).
    pub orientation: bool,
}

/// Walk all faces in a shell and collect unique vertex/edge IDs and face topology.
pub(super) fn collect_face_data(
    graph: &EntityGraph,
    face_ids: &[EntityId],
    log: &mut DiagnosticLog,
) -> Option<FaceData> {
    let mut vertex_ids = Vec::new();
    let mut edge_curve_ids = Vec::new();
    let mut faces = Vec::new();
    let mut seen_vertices: HashMap<EntityId, ()> = HashMap::new();
    let mut seen_edges: HashMap<EntityId, ()> = HashMap::new();

    for &face_id in face_ids {
        let face_entity = graph.get(face_id)?;
        let af = match AdvancedFace::from_params(&face_entity.params) {
            Ok(f) => f,
            Err(e) => {
                log.error(Phase::MapTopology, Some(face_id), format!("{e}"));
                return None;
            }
        };

        let mut bounds = Vec::new();
        for &bound_id in &af.bounds {
            let bound_entity = graph.get(bound_id)?;
            let (loop_id, is_outer, bound_orientation) = match bound_entity.type_name.as_str() {
                "FACE_OUTER_BOUND" => {
                    let fob = match FaceOuterBound::from_params(&bound_entity.params) {
                        Ok(b) => b,
                        Err(e) => {
                            log.error(Phase::MapTopology, Some(bound_id), format!("{e}"));
                            return None;
                        }
                    };
                    (fob.bound, true, fob.orientation)
                }
                "FACE_BOUND" => {
                    let fb = match FaceBound::from_params(&bound_entity.params) {
                        Ok(b) => b,
                        Err(e) => {
                            log.error(Phase::MapTopology, Some(bound_id), format!("{e}"));
                            return None;
                        }
                    };
                    (fb.bound, false, fb.orientation)
                }
                other => {
                    log.warn(
                        Phase::MapTopology,
                        Some(bound_id),
                        format!("unknown face bound type: {other}"),
                    );
                    continue;
                }
            };

            let loop_entity = graph.get(loop_id)?;
            let el = match EdgeLoop::from_params(&loop_entity.params) {
                Ok(l) => l,
                Err(e) => {
                    log.error(Phase::MapTopology, Some(loop_id), format!("{e}"));
                    return None;
                }
            };

            let mut oriented_edges = Vec::new();
            for &oe_id in &el.edge_list {
                let oe_entity = graph.get(oe_id)?;
                let oe = match OrientedEdge::from_params(&oe_entity.params) {
                    Ok(o) => o,
                    Err(e) => {
                        log.error(Phase::MapTopology, Some(oe_id), format!("{e}"));
                        return None;
                    }
                };

                let ec_entity = graph.get(oe.edge_element)?;
                let ec = match EdgeCurve::from_params(&ec_entity.params) {
                    Ok(e) => e,
                    Err(e) => {
                        log.error(Phase::MapTopology, Some(oe.edge_element), format!("{e}"));
                        return None;
                    }
                };

                // Track unique vertices
                for &vp_id in &[ec.start, ec.end] {
                    if seen_vertices.insert(vp_id, ()).is_none() {
                        vertex_ids.push(vp_id);
                    }
                }

                // Track unique edges
                if seen_edges.insert(oe.edge_element, ()).is_none() {
                    edge_curve_ids.push(oe.edge_element);
                }

                // Compose ORIENTED_EDGE and FACE_BOUND orientations
                let effective_orientation = oe.orientation == bound_orientation;
                oriented_edges.push(OrientedEdgeInfo {
                    edge_curve_id: oe.edge_element,
                    orientation: effective_orientation,
                });
            }

            bounds.push(BoundInfo {
                is_outer,
                oriented_edges,
            });
        }

        // Sort bounds so outer bound comes first
        bounds.sort_by_key(|b| if b.is_outer { 0 } else { 1 });

        faces.push(FaceInfo {
            face_id,
            surface_id: af.face_geometry,
            same_sense: af.same_sense,
            bounds,
        });
    }

    Some(FaceData {
        vertex_ids,
        edge_curve_ids,
        faces,
    })
}
