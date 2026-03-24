//! Topology mapping: STEP topology entities to popkern-core [`BRep`] via [`BRepBuilder`].
//!
//! Walks `MANIFOLD_SOLID_BREP` → `CLOSED_SHELL` → `ADVANCED_FACE` → edge loops
//! and constructs a fully-wired BRep for each solid.
//!
//! Face data collection is in the [`face_collector`] submodule.

mod face_collector;

use std::collections::HashMap;

use popkern_core::constructors::BRepBuilder;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::brep::BRep;
use popkern_core::topology::{EdgeHandle, HalfEdgeHandle, VertexHandle};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::topology::{ClosedShell, EdgeCurve, ManifoldSolidBrep, VertexPoint};

use face_collector::collect_face_data;

/// Find all `MANIFOLD_SOLID_BREP` entities and build a [`BRep`] for each.
///
/// Returns a list of `(entity_id, BRep)` pairs.
pub fn map_brep_solids(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Vec<(EntityId, BRep)> {
    let mut results = Vec::new();
    let msb_ids: Vec<EntityId> = graph
        .iter()
        .filter(|(_, e)| e.type_name == "MANIFOLD_SOLID_BREP")
        .map(|(&id, _)| id)
        .collect();

    for msb_id in msb_ids {
        match build_solid(graph, geo, msb_id, log, tolerance) {
            Some(brep) => results.push((msb_id, brep)),
            None => {
                log.error(
                    Phase::MapTopology,
                    Some(msb_id),
                    "failed to build solid from MANIFOLD_SOLID_BREP",
                );
            }
        }
    }
    results
}

/// Build a single BRep solid from a `MANIFOLD_SOLID_BREP` entity.
fn build_solid(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    msb_id: EntityId,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Option<BRep> {
    // 1. Resolve the ManifoldSolidBrep
    let entity = graph.get(msb_id)?;
    let msb = match ManifoldSolidBrep::from_params(&entity.params) {
        Ok(m) => m,
        Err(e) => {
            log.error(Phase::MapTopology, Some(msb_id), format!("{e}"));
            return None;
        }
    };

    // 2. Resolve the ClosedShell
    let shell_entity = graph.get(msb.outer)?;
    let shell = match ClosedShell::from_params(&shell_entity.params) {
        Ok(s) => s,
        Err(e) => {
            log.error(Phase::MapTopology, Some(msb.outer), format!("{e}"));
            return None;
        }
    };

    // 3. Gather all unique vertex and edge entities across all faces
    let mut vertex_map: HashMap<EntityId, VertexHandle> = HashMap::new();
    let mut edge_map: HashMap<EntityId, (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle)> =
        HashMap::new();

    let face_data = collect_face_data(graph, &shell.faces, log)?;

    let mut builder = BRepBuilder::new(tolerance);

    // Create vertices
    for &vp_id in &face_data.vertex_ids {
        if vertex_map.contains_key(&vp_id) {
            continue;
        }
        let vp_entity = graph.get(vp_id)?;
        let vp = match VertexPoint::from_params(&vp_entity.params) {
            Ok(v) => v,
            Err(e) => {
                log.error(Phase::MapTopology, Some(vp_id), format!("{e}"));
                return None;
            }
        };
        let pt = geo.resolve_point(vp.point, log)?;
        let vh = builder.add_vertex(pt);
        vertex_map.insert(vp_id, vh);
    }

    // Create edges
    for &ec_id in &face_data.edge_curve_ids {
        if edge_map.contains_key(&ec_id) {
            continue;
        }
        let ec_entity = graph.get(ec_id)?;
        let ec = match EdgeCurve::from_params(&ec_entity.params) {
            Ok(e) => e,
            Err(e) => {
                log.error(Phase::MapTopology, Some(ec_id), format!("{e}"));
                return None;
            }
        };
        let v1 = *vertex_map.get(&ec.start)?;
        let v2 = *vertex_map.get(&ec.end)?;
        let (eh, he1, he2) = builder.add_edge(v1, v2);
        edge_map.insert(ec_id, (eh, he1, he2));

        // Attach edge curve geometry
        if let Some(curve) = geo.resolve_curve(ec.edge_geometry, log) {
            builder.set_edge_curve(eh, curve);
        }
    }

    // Create solid and shell
    let solid = builder.add_solid();
    let shell_h = builder.add_shell(solid);

    // 4. Create faces with loops
    for face_info in &face_data.faces {
        let face_h = builder.add_face(shell_h);

        // Attach surface geometry
        if let Some(sof) = geo.resolve_surface(face_info.surface_id, face_info.same_sense, log) {
            builder.set_face_surface(face_h, sof);
        }

        // Wire each bound as a loop
        for (bound_idx, bound) in face_info.bounds.iter().enumerate() {
            let mut half_edges = Vec::with_capacity(bound.oriented_edges.len());
            for oe_info in &bound.oriented_edges {
                let (_eh, he_fwd, he_rev) = edge_map.get(&oe_info.edge_curve_id)?;
                let he = if oe_info.orientation {
                    *he_fwd
                } else {
                    *he_rev
                };
                half_edges.push(he);
            }

            if half_edges.is_empty() {
                log.warn(
                    Phase::MapTopology,
                    Some(face_info.face_id),
                    "empty edge loop, skipping",
                );
                continue;
            }

            let loop_result = if bound_idx == 0 && bound.is_outer {
                builder.add_loop(&half_edges, face_h)
            } else {
                builder.add_inner_loop(&half_edges, face_h)
            };

            if let Err(e) = loop_result {
                log.error(
                    Phase::MapTopology,
                    Some(face_info.face_id),
                    format!("failed to create loop: {e}"),
                );
                return None;
            }
        }
    }

    Some(builder.build_unchecked())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::{EntityInstance, Parameter};

    fn make_entity(id: EntityId, type_name: &str, params: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id,
            type_name: type_name.into(),
            params,
            complex_parts: vec![],
        }
    }

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn p_list(items: Vec<Parameter>) -> Parameter {
        Parameter::List(items)
    }

    /// Build a minimal triangle solid from STEP entities and verify topology.
    #[test]
    fn build_triangle_solid() {
        let entities = vec![
            make_entity(
                1,
                "CARTESIAN_POINT",
                vec![
                    p_str(""),
                    p_list(vec![p_real(0.0), p_real(0.0), p_real(0.0)]),
                ],
            ),
            make_entity(
                2,
                "CARTESIAN_POINT",
                vec![
                    p_str(""),
                    p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
                ],
            ),
            make_entity(
                3,
                "CARTESIAN_POINT",
                vec![
                    p_str(""),
                    p_list(vec![p_real(0.0), p_real(1.0), p_real(0.0)]),
                ],
            ),
            make_entity(
                4,
                "DIRECTION",
                vec![
                    p_str(""),
                    p_list(vec![p_real(0.0), p_real(0.0), p_real(1.0)]),
                ],
            ),
            make_entity(
                5,
                "DIRECTION",
                vec![
                    p_str(""),
                    p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
                ],
            ),
            make_entity(
                6,
                "AXIS2_PLACEMENT_3D",
                vec![p_str(""), p_ref(1), p_ref(4), p_ref(5)],
            ),
            make_entity(7, "PLANE", vec![p_str(""), p_ref(6)]),
            make_entity(
                8,
                "DIRECTION",
                vec![
                    p_str(""),
                    p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
                ],
            ),
            make_entity(9, "VECTOR", vec![p_str(""), p_ref(8), p_real(1.0)]),
            make_entity(
                10,
                "DIRECTION",
                vec![
                    p_str(""),
                    p_list(vec![p_real(-1.0), p_real(1.0), p_real(0.0)]),
                ],
            ),
            make_entity(11, "VECTOR", vec![p_str(""), p_ref(10), p_real(1.0)]),
            make_entity(
                12,
                "DIRECTION",
                vec![
                    p_str(""),
                    p_list(vec![p_real(0.0), p_real(-1.0), p_real(0.0)]),
                ],
            ),
            make_entity(13, "VECTOR", vec![p_str(""), p_ref(12), p_real(1.0)]),
            make_entity(14, "LINE", vec![p_str(""), p_ref(1), p_ref(9)]),
            make_entity(15, "LINE", vec![p_str(""), p_ref(2), p_ref(11)]),
            make_entity(16, "LINE", vec![p_str(""), p_ref(3), p_ref(13)]),
            make_entity(20, "VERTEX_POINT", vec![p_str(""), p_ref(1)]),
            make_entity(21, "VERTEX_POINT", vec![p_str(""), p_ref(2)]),
            make_entity(22, "VERTEX_POINT", vec![p_str(""), p_ref(3)]),
            make_entity(
                30,
                "EDGE_CURVE",
                vec![p_str(""), p_ref(20), p_ref(21), p_ref(14), p_enum("T")],
            ),
            make_entity(
                31,
                "EDGE_CURVE",
                vec![p_str(""), p_ref(21), p_ref(22), p_ref(15), p_enum("T")],
            ),
            make_entity(
                32,
                "EDGE_CURVE",
                vec![p_str(""), p_ref(22), p_ref(20), p_ref(16), p_enum("T")],
            ),
            make_entity(
                40,
                "ORIENTED_EDGE",
                vec![
                    p_str(""),
                    Parameter::Derived,
                    Parameter::Derived,
                    p_ref(30),
                    p_enum("T"),
                ],
            ),
            make_entity(
                41,
                "ORIENTED_EDGE",
                vec![
                    p_str(""),
                    Parameter::Derived,
                    Parameter::Derived,
                    p_ref(31),
                    p_enum("T"),
                ],
            ),
            make_entity(
                42,
                "ORIENTED_EDGE",
                vec![
                    p_str(""),
                    Parameter::Derived,
                    Parameter::Derived,
                    p_ref(32),
                    p_enum("T"),
                ],
            ),
            make_entity(
                50,
                "EDGE_LOOP",
                vec![p_str(""), p_list(vec![p_ref(40), p_ref(41), p_ref(42)])],
            ),
            make_entity(
                60,
                "FACE_OUTER_BOUND",
                vec![p_str(""), p_ref(50), p_enum("T")],
            ),
            make_entity(
                70,
                "ADVANCED_FACE",
                vec![p_str(""), p_list(vec![p_ref(60)]), p_ref(7), p_enum("T")],
            ),
            make_entity(80, "CLOSED_SHELL", vec![p_str(""), p_list(vec![p_ref(70)])]),
            make_entity(90, "MANIFOLD_SOLID_BREP", vec![p_str("body"), p_ref(80)]),
        ];

        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(entities, &mut log);
        assert!(!log.has_errors(), "graph build errors: {:?}", log.entries());

        let mut geo_cache = GeometryCache::new(&graph);
        let tolerance = Tolerance::default();

        let solids = map_brep_solids(&graph, &mut geo_cache, &mut log, tolerance);
        assert_eq!(solids.len(), 1, "expected 1 solid, got {}", solids.len());

        let (entity_id, brep) = &solids[0];
        assert_eq!(*entity_id, 90);
        assert_eq!(brep.vertices.len(), 3);
        assert_eq!(brep.edges.len(), 3);
        assert_eq!(brep.faces.len(), 1);
        assert_eq!(brep.solids.len(), 1);
        assert_eq!(brep.shells.len(), 1);
    }
}
