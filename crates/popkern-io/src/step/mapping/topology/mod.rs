//! Topology mapping: STEP topology entities to popkern-core [`BRep`] via [`BRepBuilder`].
//!
//! Walks `MANIFOLD_SOLID_BREP` / `SHELL_BASED_SURFACE_MODEL` → shell →
//! `ADVANCED_FACE` → edge loops and constructs a fully-wired BRep for each body.
//!
//! Face data collection is in the [`face_collector`] submodule.

mod face_collector;
mod open_shell;
mod pcurve_orientation;

use std::collections::HashMap;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::brep::BRep;
use popkern_core::topology::{EdgeHandle, HalfEdgeHandle, VertexHandle};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::topology::{
    ClosedShell, EdgeCurve, ManifoldSolidBrep, OpenShell, ShellBasedSurfaceModel, VertexPoint,
};

use face_collector::collect_face_data;
use open_shell::normalize_open_shell_boundaries;
use pcurve_orientation::select_orient_pcurve_to_half_edge;

/// Find BRep body representation items and build a [`BRep`] for each.
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

    let surface_model_ids: Vec<EntityId> = graph
        .iter()
        .filter(|(_, e)| e.type_name == "SHELL_BASED_SURFACE_MODEL")
        .map(|(&id, _)| id)
        .collect();

    for item_id in surface_model_ids {
        let Some(entity) = graph.get(item_id) else {
            continue;
        };
        let model = match ShellBasedSurfaceModel::from_params(&entity.params) {
            Ok(m) => m,
            Err(e) => {
                log.error(Phase::MapTopology, Some(item_id), format!("{e}"));
                continue;
            }
        };

        for shell_id in model.shells {
            match build_shell_body(graph, geo, shell_id, true, log, tolerance) {
                Some(brep) => results.push((item_id, brep)),
                None => log.error(
                    Phase::MapTopology,
                    Some(shell_id),
                    "failed to build sheet body from SHELL_BASED_SURFACE_MODEL",
                ),
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

    build_shell_body(graph, geo, msb.outer, false, log, tolerance)
}

fn build_shell_body(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    shell_id: EntityId,
    is_open_shell: bool,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Option<BRep> {
    let shell_entity = graph.get(shell_id)?;
    let faces = match shell_entity.type_name.as_str() {
        "CLOSED_SHELL" => match ClosedShell::from_params(&shell_entity.params) {
            Ok(shell) => shell.faces,
            Err(e) => {
                log.error(Phase::MapTopology, Some(shell_id), format!("{e}"));
                return None;
            }
        },
        "OPEN_SHELL" => match OpenShell::from_params(&shell_entity.params) {
            Ok(shell) => shell.faces,
            Err(e) => {
                log.error(Phase::MapTopology, Some(shell_id), format!("{e}"));
                return None;
            }
        },
        other => {
            log.error(
                Phase::MapTopology,
                Some(shell_id),
                format!("expected shell entity, got {other}"),
            );
            return None;
        }
    };

    build_brep_from_faces(graph, geo, &faces, is_open_shell, log, tolerance)
}

fn build_brep_from_faces(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    faces: &[EntityId],
    is_open_shell: bool,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Option<BRep> {
    let mut vertex_map: HashMap<EntityId, VertexHandle> = HashMap::new();
    let mut edge_map: HashMap<EntityId, (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle)> =
        HashMap::new();

    let face_data = collect_face_data(graph, faces, log)?;

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

        // STEP EDGE_CURVE geometry is bounded by the topological start/end
        // vertices. Resolve the native curve, then recover the edge parameter
        // range from those vertices so arcs and NURBS trims do not import as
        // their full underlying curve domains.
        if let Some(curve) = geo.resolve_curve(ec.edge_geometry, log) {
            let start = builder.brep().vertices[v1].position;
            let end = builder.brep().vertices[v2].position;
            let curve =
                recover_edge_curve_trim(ec_id, curve, start, end, ec.same_sense, tolerance, log);
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
                if builder.brep().half_edges[he].pcurve.is_none()
                    && let Some((edge_handle, _he_fwd, _he_rev)) =
                        edge_map.get(&oe_info.edge_curve_id)
                    && let Some(edge) = builder.brep().edges.get(*edge_handle)
                    && let Some(curve) = edge.curve.as_ref()
                    && let Some(surface) = builder.brep().faces[face_h].surface.as_ref()
                    && let Some(edge_curve) = graph
                        .get(oe_info.edge_curve_id)
                        .and_then(|entity| EdgeCurve::from_params(&entity.params).ok())
                {
                    let pcurves = geo.resolve_surface_curve_pcurves(
                        edge_curve.edge_geometry,
                        face_info.surface_id,
                        log,
                    );
                    if let Some(oriented) = select_orient_pcurve_to_half_edge(
                        pcurves,
                        curve,
                        surface,
                        oe_info.orientation,
                        face_info.face_id,
                        tolerance,
                        log,
                    ) {
                        builder.brep_mut().half_edges[he].pcurve = Some(oriented);
                    }
                }
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

    if is_open_shell {
        normalize_open_shell_boundaries(&mut builder);
    }

    Some(builder.build_unchecked())
}

fn recover_edge_curve_trim(
    edge_curve_id: EntityId,
    mut curve: CurveOnEdge,
    start: Point3,
    end: Point3,
    same_sense: bool,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) -> CurveOnEdge {
    let Some(mut start_t) = curve.geometry.closest_point(&start) else {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            "could not recover EDGE_CURVE start parameter from topological vertex",
        );
        return curve;
    };
    let Some(mut end_t) = curve.geometry.closest_point(&end) else {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            "could not recover EDGE_CURVE end parameter from topological vertex",
        );
        return curve;
    };

    let start_gap = curve.geometry.point_at(start_t).distance_to(&start);
    let end_gap = curve.geometry.point_at(end_t).distance_to(&end);
    let limit = tolerance.linear().max(1e-9);
    if start_gap > limit || end_gap > limit {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            format!(
                "declined EDGE_CURVE trim recovery: endpoint gaps {:.3e}, {:.3e} exceed {:.3e}",
                start_gap, end_gap, limit
            ),
        );
        return curve;
    }

    normalize_periodic_parameter_order(&curve.geometry, same_sense, &mut start_t, &mut end_t);

    curve.t_start = start_t;
    curve.t_end = end_t;

    curve
}

fn normalize_periodic_parameter_order(
    geometry: &CurveGeometry,
    same_sense: bool,
    start_t: &mut f64,
    end_t: &mut f64,
) {
    let is_periodic = matches!(
        geometry,
        CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. } | CurveGeometry::Ellipse { .. }
    );
    if !is_periodic {
        return;
    }

    if same_sense {
        while *end_t < *start_t {
            *end_t += std::f64::consts::TAU;
        }
    } else {
        while *start_t < *end_t {
            *start_t += std::f64::consts::TAU;
        }
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
