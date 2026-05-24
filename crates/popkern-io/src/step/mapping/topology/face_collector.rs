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
    /// Entity ID of the ORIENTED_EDGE.
    pub oriented_edge_id: EntityId,
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
                    oriented_edge_id: oe_id,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::{EntityInstance, Parameter};

    fn ent(id: EntityId, ty: &str, p: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id,
            type_name: ty.into(),
            params: p,
            complex_parts: vec![],
        }
    }
    fn s(v: &str) -> Parameter {
        Parameter::String(v.into())
    }
    #[allow(dead_code)]
    fn r(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn er(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn l(items: Vec<Parameter>) -> Parameter {
        Parameter::List(items)
    }
    fn b(val: bool) -> Parameter {
        Parameter::Enum(if val { "T" } else { "F" }.into())
    }

    fn graph(entities: Vec<EntityInstance>) -> EntityGraph {
        EntityGraph::build(entities, &mut DiagnosticLog::new())
    }

    /// Build a minimal face topology graph:
    /// #1,#2 VERTEX_POINT → #10 EDGE_CURVE → #20 ORIENTED_EDGE
    /// → #30 EDGE_LOOP → #40 FACE_OUTER_BOUND → #50 ADVANCED_FACE
    fn minimal_face_entities() -> Vec<EntityInstance> {
        vec![
            // Vertex points
            ent(1, "VERTEX_POINT", vec![s(""), er(100)]),
            ent(2, "VERTEX_POINT", vec![s(""), er(101)]),
            // Edge curve: start=#1, end=#2, curve=#200, same_sense=T
            ent(
                10,
                "EDGE_CURVE",
                vec![s(""), er(1), er(2), er(200), b(true)],
            ),
            // Oriented edge: *, *, #10, orientation=T
            ent(
                20,
                "ORIENTED_EDGE",
                vec![
                    s(""),
                    Parameter::Derived,
                    Parameter::Derived,
                    er(10),
                    b(true),
                ],
            ),
            // Edge loop with one oriented edge
            ent(30, "EDGE_LOOP", vec![s(""), l(vec![er(20)])]),
            // Face outer bound: loop=#30, orientation=T
            ent(40, "FACE_OUTER_BOUND", vec![s(""), er(30), b(true)]),
            // Advanced face: bounds=[#40], surface=#300, same_sense=T
            ent(
                50,
                "ADVANCED_FACE",
                vec![s(""), l(vec![er(40)]), er(300), b(true)],
            ),
        ]
    }

    #[test]
    fn collect_minimal_face() {
        let ents = minimal_face_entities();
        let g = graph(ents);
        let mut log = DiagnosticLog::new();
        let data = collect_face_data(&g, &[50], &mut log).unwrap();
        assert_eq!(data.vertex_ids, vec![1, 2]);
        assert_eq!(data.edge_curve_ids, vec![10]);
        assert_eq!(data.faces.len(), 1);
        assert_eq!(data.faces[0].surface_id, 300);
        assert!(data.faces[0].same_sense);
        assert_eq!(data.faces[0].bounds.len(), 1);
        assert!(data.faces[0].bounds[0].is_outer);
    }

    #[test]
    fn collect_face_data_missing_face_returns_none() {
        let g = graph(vec![]);
        let mut log = DiagnosticLog::new();
        assert!(collect_face_data(&g, &[999], &mut log).is_none());
    }

    #[test]
    fn collect_face_data_bad_advanced_face_returns_none() {
        // Entity exists but params are invalid for AdvancedFace
        let ents = vec![ent(1, "ADVANCED_FACE", vec![s("bad")])];
        let g = graph(ents);
        let mut log = DiagnosticLog::new();
        assert!(collect_face_data(&g, &[1], &mut log).is_none());
    }

    #[test]
    fn collect_face_data_unknown_bound_type_skipped() {
        let ents = vec![
            // A bound entity that's neither FACE_OUTER_BOUND nor FACE_BOUND
            ent(40, "SOME_OTHER_BOUND", vec![s(""), er(30), b(true)]),
            ent(
                50,
                "ADVANCED_FACE",
                vec![s(""), l(vec![er(40)]), er(300), b(true)],
            ),
        ];
        let g = graph(ents);
        let mut log = DiagnosticLog::new();
        let data = collect_face_data(&g, &[50], &mut log);
        // The unknown bound is skipped, face should still be present but with empty bounds
        assert!(data.is_some());
        assert_eq!(data.unwrap().faces[0].bounds.len(), 0);
    }

    #[test]
    fn collect_face_with_inner_bound() {
        let mut ents = minimal_face_entities();
        // Add second edge and oriented edge for inner loop
        ents.push(ent(
            11,
            "EDGE_CURVE",
            vec![s(""), er(1), er(2), er(201), b(true)],
        ));
        ents.push(ent(
            21,
            "ORIENTED_EDGE",
            vec![
                s(""),
                Parameter::Derived,
                Parameter::Derived,
                er(11),
                b(false),
            ],
        ));
        ents.push(ent(31, "EDGE_LOOP", vec![s(""), l(vec![er(21)])]));
        ents.push(ent(41, "FACE_BOUND", vec![s(""), er(31), b(true)]));
        // Update ADVANCED_FACE to include inner bound
        // Remove old face entity and re-add with both bounds
        ents.retain(|e| e.id != 50);
        ents.push(ent(
            50,
            "ADVANCED_FACE",
            vec![s(""), l(vec![er(40), er(41)]), er(300), b(true)],
        ));
        let g = graph(ents);
        let mut log = DiagnosticLog::new();
        let data = collect_face_data(&g, &[50], &mut log).unwrap();
        assert_eq!(data.faces[0].bounds.len(), 2);
        // Outer bound should be sorted first
        assert!(data.faces[0].bounds[0].is_outer);
        assert!(!data.faces[0].bounds[1].is_outer);
    }
}
