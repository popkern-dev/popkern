//! STEP topology import fidelity diagnostics with entity provenance.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, EdgeHandle, FaceHandle, HalfEdgeHandle,
    analyze_brep_edge_surface_compatibility,
};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

const STEP_IMPORT_FIDELITY_SAMPLES: usize = 9;

#[derive(Debug, Default)]
pub(super) struct StepTopologyImportProvenance {
    pub edge_curves: HashMap<EdgeHandle, EntityId>,
    pub half_edges: HashMap<HalfEdgeHandle, StepHalfEdgeImportProvenance>,
}

#[derive(Debug, Clone)]
pub(super) struct StepHalfEdgeImportProvenance {
    pub edge_curve_id: EntityId,
    pub oriented_edge_id: EntityId,
    pub face_id: EntityId,
    pub surface_id: EntityId,
    pub pcurve_id: Option<EntityId>,
    pub reference_curve_id: Option<EntityId>,
    pub selected_pcurve_gap: Option<f64>,
}

pub(super) fn report_step_import_fidelity_diagnostics(
    brep: &BRep,
    provenance: &StepTopologyImportProvenance,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) {
    report_pcurve_fidelity_diagnostics(brep, provenance, tolerance, log);
    report_edge_surface_fidelity_diagnostics(brep, provenance, tolerance, log);
}

fn report_pcurve_fidelity_diagnostics(
    brep: &BRep,
    provenance: &StepTopologyImportProvenance,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) {
    let limit = tolerance.linear();
    for edge in brep.edges.keys() {
        let Some(report) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        if report.max_deviation <= limit {
            continue;
        }
        for half_edge in report
            .half_edges
            .iter()
            .filter(|half_edge| half_edge.max_deviation > limit)
        {
            let imported = provenance.half_edges.get(&half_edge.half_edge);
            let edge_curve_id = imported
                .map(|entry| entry.edge_curve_id)
                .or_else(|| provenance.edge_curves.get(&edge).copied());
            let face_id = imported.map(|entry| entry.face_id);
            let surface_id = imported.map(|entry| entry.surface_id);
            let pcurve_id = imported.and_then(|entry| entry.pcurve_id);
            let reference_curve_id = imported.and_then(|entry| entry.reference_curve_id);
            let selected_gap = imported.and_then(|entry| entry.selected_pcurve_gap);
            log.warn(
                Phase::MapTopology,
                edge_curve_id,
                format!(
                    "STEP import pcurve fidelity blocker: edge={edge:?} half_edge={:?} face={:?} surface_kind={} edge_curve={} oriented_edge={} advanced_face={} surface={} pcurve={} reference_curve={} selected_gap={} max_deviation={:.3e} same_fraction={:.3e} trace={:.3e} worst_fraction={:.3e}",
                    half_edge.half_edge,
                    half_edge.face,
                    face_surface_kind(brep, half_edge.face),
                    entity_label(edge_curve_id),
                    entity_label(imported.map(|entry| entry.oriented_edge_id)),
                    entity_label(face_id),
                    entity_label(surface_id),
                    entity_label(pcurve_id),
                    entity_label(reference_curve_id),
                    selected_gap
                        .map(|gap| format!("{gap:.3e}"))
                        .unwrap_or_else(|| "none".to_string()),
                    half_edge.max_deviation,
                    half_edge.max_same_fraction_deviation,
                    half_edge.max_trace_deviation,
                    half_edge.worst_fraction,
                ),
            );
        }
    }
}

fn report_edge_surface_fidelity_diagnostics(
    brep: &BRep,
    provenance: &StepTopologyImportProvenance,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) {
    let limit = tolerance.linear();
    let report = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: limit,
            max_anchor_motion: 0.0,
            samples: STEP_IMPORT_FIDELITY_SAMPLES,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    for finding in report
        .findings
        .iter()
        .filter(|finding| finding.max_deviation > limit)
    {
        let imported = provenance.half_edges.get(&finding.half_edge);
        let edge_curve_id = imported
            .map(|entry| entry.edge_curve_id)
            .or_else(|| provenance.edge_curves.get(&finding.edge).copied());
        log.warn(
            Phase::MapTopology,
            edge_curve_id,
            format!(
                "STEP import edge/surface fidelity blocker: edge={:?} half_edge={:?} face={:?} surface_kind={} edge_curve={} oriented_edge={} advanced_face={} surface={} max_deviation={:.3e} curve_parameter={:.3e}",
                finding.edge,
                finding.half_edge,
                finding.face,
                face_surface_kind(brep, finding.face),
                entity_label(edge_curve_id),
                entity_label(imported.map(|entry| entry.oriented_edge_id)),
                entity_label(imported.map(|entry| entry.face_id)),
                entity_label(imported.map(|entry| entry.surface_id)),
                finding.max_deviation,
                finding.max_curve_parameter,
            ),
        );
    }
}

fn face_surface_kind(brep: &BRep, face: FaceHandle) -> &'static str {
    match brep
        .faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "cone",
        Some(SurfaceGeometry::Sphere { .. }) => "sphere",
        Some(SurfaceGeometry::Torus { .. }) => "torus",
        Some(SurfaceGeometry::Nurbs(_)) => "nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "gregory",
        None => "missing",
    }
}

fn entity_label(entity: Option<EntityId>) -> String {
    entity
        .map(|id| format!("#{id}"))
        .unwrap_or_else(|| "none".to_string())
}
