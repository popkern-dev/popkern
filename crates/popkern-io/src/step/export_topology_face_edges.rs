//! Face-local STEP edge export for incoherent shared face uses.

use std::collections::HashMap;

use popkern_core::geometry::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::error::IoError;
use crate::step::export::{ExportContext, bool_, ref_, str_param};
use crate::step::export_pcurve_synthesis::{
    export_pcurve_trace_matches_edge, step_export_pcurve_association_tolerance,
    synthesize_export_face_use_curve_and_pcurve, synthesize_export_half_edge_pcurve_on_surface,
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
            if popkern_core::trace::trace_enabled("step-io") {
                explain_face_use_fork(brep, edge, face_use, surface_map);
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
    // Mirror `export_pcurve_for_use` exactly: a stored pcurve that
    // fails validation must NOT shadow a synthesized one that passes.
    // (A cylinder cap's stored rim pcurve is a phase-shifted ring —
    // valid for tessellation, but not co-parameterized with the edge
    // curve as STEP requires; synthesis projects the edge curve and
    // is co-parameterized by construction.)
    let tolerance = step_export_pcurve_association_tolerance(brep);
    let matches = |pcurve: &popkern_core::geometry::pcurve::Pcurve| {
        pcurve_endpoints_match_edge_segment(
            brep,
            edge_handle,
            pcurve,
            &surface.surface.geometry,
            tolerance,
        ) || export_pcurve_trace_matches_edge(
            brep,
            edge_handle,
            face_use.half_edge,
            pcurve,
            &surface.surface.geometry,
            tolerance,
        )
    };
    brep.face_use_pcurve(face_use.face, face_use.half_edge)
        .or(he.pcurve.as_ref())
        .cloned()
        .filter(|pcurve| matches(pcurve))
        .or_else(|| {
            synthesize_export_half_edge_pcurve_on_surface(
                brep,
                edge_handle,
                face_use.half_edge,
                &surface.surface.geometry,
            )
            .filter(|pcurve| matches(pcurve))
        })
        // Mirror export_pcurve_for_use's final fallback: a faithful
        // ungated projection keeps the edge SHARED (sewable) — the
        // face-local fork guarantees a once-used edge (free bound).
        .or_else(|| {
            crate::step::export_pcurve_synthesis::synthesize_half_edge_pcurve_ungated(
                brep,
                edge_handle,
                face_use.half_edge,
                &surface.surface.geometry,
            )
        })
        .map(|_| ())
}

fn explain_face_use_fork(
    brep: &BRep,
    edge_handle: EdgeHandle,
    face_use: EdgeFaceUse,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
) {
    let Some(he) = brep.half_edges.get(face_use.half_edge) else {
        eprintln!(
            "  [step-fork] {edge_handle:?} on {:?}: he missing",
            face_use.face
        );
        return;
    };
    let Some(surface) = surface_map.get(&face_use.face) else {
        eprintln!(
            "  [step-fork] {edge_handle:?} on {:?}: no surface",
            face_use.face
        );
        return;
    };
    let stored = brep
        .face_use_pcurve(face_use.face, face_use.half_edge)
        .or(he.pcurve.as_ref())
        .cloned();
    let tolerance = step_export_pcurve_association_tolerance(brep);
    let describe = |label: &str, pc: &popkern_core::geometry::pcurve::Pcurve| {
        let ep = pcurve_endpoints_match_edge_segment(
            brep,
            edge_handle,
            pc,
            &surface.surface.geometry,
            tolerance,
        );
        let tr = export_pcurve_trace_matches_edge(
            brep,
            edge_handle,
            face_use.half_edge,
            pc,
            &surface.surface.geometry,
            tolerance,
        );
        let (t0, t1) = pc.parameter_range();
        let a = pc.point_at(t0);
        let b = pc.point_at(t1);
        eprintln!(
            "  [step-fork] {edge_handle:?} on {:?}: {label} endpoints_ok={ep} trace_ok={tr} pc ({:.2},{:.2})->({:.2},{:.2}) t[{t0:.3},{t1:.3}]",
            face_use.face, a.x, a.y, b.x, b.y
        );
    };
    match &stored {
        Some(pc) => describe("stored", pc),
        None => eprintln!(
            "  [step-fork] {edge_handle:?} on {:?}: NO stored pcurve",
            face_use.face
        ),
    }
    if let Some(pc) = synthesize_export_half_edge_pcurve_on_surface(
        brep,
        edge_handle,
        face_use.half_edge,
        &surface.surface.geometry,
    ) {
        describe("synth ", &pc);
    } else {
        eprintln!(
            "  [step-fork] {edge_handle:?} on {:?}: synthesis returned None",
            face_use.face
        );
    }
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
    let step_nurbs = crate::step::export_geometry::pcurve_to_step_chart(
        &surface.surface.geometry,
        crate::step::export_pcurve::export_pcurve_nurbs(
            &pcurve,
            &surface.surface.geometry,
            brep.tolerance.linear(),
        ),
    );
    let pcurve_id = export.export_pcurve(&step_nurbs, surface.id);
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
