//! Face data collection: walks STEP topology to gather unique vertices, edges,
//! and per-face bound/loop structure before BRep construction.

use std::collections::HashMap;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, Parameter};
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
#[cfg(test)]
pub(super) fn collect_face_data(
    graph: &EntityGraph,
    face_ids: &[EntityId],
    log: &mut DiagnosticLog,
) -> Option<FaceData> {
    collect_face_data_reversing(graph, face_ids, &std::collections::HashSet::new(), log)
}

/// Face data collection with the shell-listed faces in `reversed_uses`
/// used reversed on top of their own `ORIENTED_FACE` sense: the faces
/// of a void shell whose `ORIENTED_CLOSED_SHELL` is `.F.` (PK-880).
pub(super) fn collect_face_data_reversing(
    graph: &EntityGraph,
    face_ids: &[EntityId],
    reversed_uses: &std::collections::HashSet<EntityId>,
    log: &mut DiagnosticLog,
) -> Option<FaceData> {
    let mut vertex_ids = Vec::new();
    let mut edge_curve_ids = Vec::new();
    let mut faces = Vec::new();
    let mut seen_vertices: HashMap<EntityId, ()> = HashMap::new();
    let mut seen_edges: HashMap<EntityId, ()> = HashMap::new();

    for &shell_face_id in face_ids {
        // A shell may list an ORIENTED_FACE('', *, face, orientation):
        // the face it wraps, used reversed when the orientation is .F.
        // — its sense and every bound's orientation flip (PK-635,
        // measured: piston.stp's closed shell).
        let (face_id, reversed) = match graph.get(shell_face_id) {
            Some(e) if e.type_name == "ORIENTED_FACE" => match (e.params.get(2), e.params.get(3)) {
                (Some(Parameter::EntityRef(inner)), orientation) => (
                    *inner,
                    matches!(orientation, Some(Parameter::Enum(v)) if v == "F"),
                ),
                _ => {
                    log.error(
                        Phase::MapTopology,
                        Some(shell_face_id),
                        "ORIENTED_FACE without a face element",
                    );
                    return None;
                }
            },
            _ => (shell_face_id, false),
        };
        let reversed = reversed != reversed_uses.contains(&shell_face_id);
        let face_entity = graph.get(face_id)?;
        let mut af = match AdvancedFace::from_params(&face_entity.params) {
            Ok(f) => f,
            Err(e) => {
                log.error(Phase::MapTopology, Some(face_id), format!("{e}"));
                return None;
            }
        };
        if reversed {
            af.same_sense = !af.same_sense;
        }

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

            // A reversed face reverses each of its bounds.
            let bound_orientation = bound_orientation != reversed;
            let loop_entity = graph.get(loop_id)?;
            if loop_entity.type_name == "VERTEX_LOOP" {
                // A degenerate bound — one vertex, a cone's apex: the
                // face's other loop already closes the region, so the
                // point loop carries no boundary to wire (PK-635,
                // measured: team.stp).
                log.warn(
                    Phase::MapTopology,
                    Some(loop_id),
                    "VERTEX_LOOP bound (degenerate apex) skipped",
                );
                continue;
            }
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

            // A reversed bound (`orientation .F.`) is the loop traversed
            // backwards: every edge flipped (composed above) AND the
            // order reversed. Flipping alone leaves a 3+-edge loop open
            // (measured: every face of an OpenCASCADE-written box, whose
            // writer marks its outer loops `FACE_BOUND(.F.)`, refused
            // as "loop is not closed").
            if !bound_orientation {
                oriented_edges.reverse();
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

/// The bound whose sampled boundary has the largest extent (bounding
/// box diagonal): the outer loop of a face whose file marks none. A
/// single bound is its own outer loop.
pub(super) fn widest_bound(
    bounds: &[BoundInfo],
    edge_map: &HashMap<
        EntityId,
        (
            popkern_core::topology::EdgeHandle,
            popkern_core::topology::HalfEdgeHandle,
            popkern_core::topology::HalfEdgeHandle,
        ),
    >,
    brep: &popkern_core::topology::BRep,
) -> Option<usize> {
    use popkern_core::geometry::curve::CurveEval;
    if bounds.len() <= 1 {
        return (!bounds.is_empty()).then_some(0);
    }
    let extent = |bound: &BoundInfo| {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for oe in &bound.oriented_edges {
            let Some(curve) = edge_map
                .get(&oe.edge_curve_id)
                .and_then(|(eh, _, _)| brep.edges.get(*eh))
                .and_then(|e| e.curve.as_ref())
            else {
                continue;
            };
            for k in 0..=8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * f64::from(k) / 8.0;
                let p = curve.geometry.point_at(t);
                for (i, x) in [p.x, p.y, p.z].into_iter().enumerate() {
                    lo[i] = lo[i].min(x);
                    hi[i] = hi[i].max(x);
                }
            }
        }
        (0..3)
            .map(|i| (hi[i] - lo[i]).max(0.0).powi(2))
            .sum::<f64>()
    };
    (0..bounds.len()).max_by(|&a, &b| extent(&bounds[a]).total_cmp(&extent(&bounds[b])))
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
