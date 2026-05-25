//! Face-local STEP edge export for incoherent shared face uses.

use std::collections::HashMap;

use popkern_core::geometry::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::error::IoError;
use crate::step::export::{ExportContext, bool_, ref_, str_param};
use crate::step::export_pcurve_synthesis::{
    export_pcurve_trace_matches_edge, step_export_pcurve_association_tolerance,
    synthesize_export_face_use_curve_and_pcurve,
};
use crate::step::export_topology::ExportedFaceSurface;
use crate::step::export_topology_face_uses::{EdgeFaceUse, edge_face_uses};
use crate::step::export_topology_helpers::pcurve_endpoints_match_edge_segment;
use crate::step::parser::types::EntityId;

pub(super) fn pre_export_face_use_edges(
    export: &mut ExportContext,
    brep: &BRep,
    vertex_map: &HashMap<VertexHandle, EntityId>,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
    edge_map: &HashMap<EdgeHandle, EntityId>,
    local_edge_map: &mut HashMap<EdgeFaceUse, EntityId>,
) -> Result<(), IoError> {
    for edge in brep.edges.keys() {
        for face_use in edge_face_uses(brep, edge) {
            if pcurve_exportable_for_face_use(brep, edge, face_use, surface_map).is_some() {
                continue;
            }
            let Some(local_edge_id) =
                export_local_face_use_edge(export, brep, edge, face_use, vertex_map, surface_map)?
            else {
                continue;
            };
            if edge_map
                .get(&edge)
                .is_some_and(|global| *global == local_edge_id)
            {
                continue;
            }
            local_edge_map.insert(face_use, local_edge_id);
        }
    }
    Ok(())
}

fn pcurve_exportable_for_face_use(
    brep: &BRep,
    edge_handle: EdgeHandle,
    face_use: EdgeFaceUse,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
) -> Option<()> {
    let he = brep.half_edges.get(face_use.half_edge)?;
    let surface = surface_map.get(&face_use.face)?;
    let pcurve = he.pcurve.as_ref()?;
    (pcurve_endpoints_match_edge_segment(
        brep,
        edge_handle,
        pcurve,
        &surface.surface.geometry,
        step_export_pcurve_association_tolerance(brep),
    ) || export_pcurve_trace_matches_edge(
        brep,
        edge_handle,
        face_use.half_edge,
        pcurve,
        &surface.surface.geometry,
        step_export_pcurve_association_tolerance(brep),
    ))
    .then_some(())
}

fn export_local_face_use_edge(
    export: &mut ExportContext,
    brep: &BRep,
    edge_handle: EdgeHandle,
    face_use: EdgeFaceUse,
    _vertex_map: &HashMap<VertexHandle, EntityId>,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
) -> Result<Option<EntityId>, IoError> {
    if !brep.half_edges.contains_key(face_use.half_edge) {
        return Ok(None);
    }
    let Some(surface) = surface_map.get(&face_use.face) else {
        return Ok(None);
    };
    let Some((curve, pcurve)) = synthesize_export_face_use_curve_and_pcurve(
        brep,
        edge_handle,
        face_use.half_edge,
        &surface.surface.geometry,
    ) else {
        return Ok(None);
    };
    let start_point = curve.geometry.point_at(curve.t_start);
    let end_point = curve.geometry.point_at(curve.t_end);
    let start_point_id = export.export_point(&start_point);
    let end_point_id = export.export_point(&end_point);
    let start_vp = export.emit("VERTEX_POINT", vec![str_param(""), ref_(start_point_id)]);
    let end_vp = export.emit("VERTEX_POINT", vec![str_param(""), ref_(end_point_id)]);
    let curve_3d_id = export.export_curve(&curve);
    let pcurve_id = export.export_pcurve(&pcurve.to_nurbs(), surface.id);
    let surface_curve_id = export.export_surface_curve(curve_3d_id, &[pcurve_id]);
    Ok(Some(export.emit(
        "EDGE_CURVE",
        vec![
            str_param(""),
            ref_(start_vp),
            ref_(end_vp),
            ref_(surface_curve_id),
            bool_(true),
        ],
    )))
}
