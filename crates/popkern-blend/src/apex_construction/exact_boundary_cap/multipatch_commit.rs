//! Transactional commit boundary for validated multi-patch candidates.
#![allow(clippy::result_large_err)]

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::ExactBoundaryMultiPatchCandidate;
use super::report::ExactBoundaryMultiPatchDiagnostic;
use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

#[path = "multipatch_commit_topology.rs"]
mod commit_topology;
use commit_topology::materialize_patch_faces;

/// Candidate-only topology delta that a later transaction may materialize.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryMultiPatchCommitPreflightReport {
    /// One logical cap feature represented by multiple generated faces.
    pub logical_feature_name: String,
    /// Generated cap face names in deterministic commit order.
    pub generated_face_names: Vec<String>,
    /// Internal interface edge names in deterministic commit order.
    pub internal_interface_names: Vec<String>,
    /// Source outer edges reused by the candidate.
    pub reused_outer_edges: Vec<EdgeHandle>,
    /// Number of pcurves required by the candidate commit.
    pub pcurve_count: usize,
    /// Number of generated pole spokes required by the candidate commit.
    pub generated_spoke_count: usize,
}

/// Live topology committed for one validated multi-patch sector split.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryMultiPatchCommitReport {
    /// Candidate-only preflight evidence used by this transaction.
    pub preflight: ExactBoundaryMultiPatchCommitPreflightReport,
    /// Source cap face reused as the outer strip patch face.
    pub strip_face: FaceHandle,
    /// Generated cap face for the inner pole patch.
    pub inner_face: FaceHandle,
    /// Generated strip/inner interface edge.
    pub internal_interface_edges: Vec<EdgeHandle>,
    /// Generated connector edges closing the strip and inner loops.
    pub generated_connector_edges: Vec<EdgeHandle>,
    /// Number of pcurves authored on live half-edges by this commit.
    pub committed_pcurve_count: usize,
}

/// Verify commit readiness without mutating live topology.
pub(crate) fn preflight_multipatch_commit(
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
) -> Result<ExactBoundaryMultiPatchCommitPreflightReport, ExactBoundaryMultiPatchDiagnostic> {
    let strip = &candidate.strip_inventory;
    if !surface_layout_is_valid(&candidate.strip_candidate.surface)
        || !surface_layout_is_valid(&candidate.inner_patch.surface)
    {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::commit_preflight_missing_surface(
                sector,
                strip,
                first_invalid_surface_name(candidate),
            ),
        );
    }
    if !pcurve_is_valid(&candidate.strip_candidate.outer_pcurve)
        || !pcurve_is_valid(&candidate.interface.strip_pcurve)
        || !pcurve_is_valid(&candidate.inner_patch.interface_pcurve)
    {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::commit_preflight_missing_pcurve(
                sector,
                strip,
                candidate.strip_inventory.patch_name.clone(),
                Some(candidate.interface.interface_name.clone()),
            ),
        );
    }
    if !curve_layout_is_valid(&candidate.interface.curve)
        || candidate
            .support_authority_boundary
            .as_ref()
            .is_some_and(|curve| !curve_layout_is_valid(curve))
        || candidate.interface.interface_name != candidate.inner_patch.interface_name
        || candidate.interface.interface_name != candidate.inventory.interfaces[0].interface_name
    {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::commit_preflight_missing_interface(
                sector,
                strip,
                candidate.interface.interface_name.clone(),
            ),
        );
    }
    if generated_spoke_count(&candidate.inner_patch.surface) < 2 {
        return Err(
            ExactBoundaryMultiPatchDiagnostic::commit_preflight_missing_generated_spoke(
                sector,
                strip,
                candidate.inner_patch.patch_name.clone(),
            ),
        );
    }
    Ok(ExactBoundaryMultiPatchCommitPreflightReport {
        logical_feature_name: format!("{}.multipatch_cap", sector.sector_name),
        generated_face_names: vec![
            candidate.strip_inventory.patch_name.clone(),
            candidate.inner_patch.patch_name.clone(),
        ],
        internal_interface_names: vec![candidate.interface.interface_name.clone()],
        reused_outer_edges: vec![candidate.strip_inventory.edge],
        pcurve_count: 3,
        generated_spoke_count: generated_spoke_count(&candidate.inner_patch.surface),
    })
}

/// Materialize a validated multi-patch candidate into live BRep topology.
///
/// The transaction reuses the source sector cap face as the outer strip face,
/// creates one generated inner face, wires one generated interface edge between
/// them, and restores the full BRep snapshot on any materialization or
/// validation failure.
pub(crate) fn commit_multipatch_candidate_transactionally(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
) -> Result<ExactBoundaryMultiPatchCommitReport, ExactBoundaryMultiPatchDiagnostic> {
    let preflight = preflight_multipatch_commit(sector, candidate)?;
    let snapshot = brep.clone();
    let result = commit_preflighted_candidate(brep, sector, candidate, preflight);
    match result {
        Ok(report) => Ok(report),
        Err(diagnostic) => {
            *brep = snapshot;
            Err(diagnostic)
        }
    }
}

fn commit_preflighted_candidate(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
    preflight: ExactBoundaryMultiPatchCommitPreflightReport,
) -> Result<ExactBoundaryMultiPatchCommitReport, ExactBoundaryMultiPatchDiagnostic> {
    let strip = &candidate.strip_inventory;
    let split = materialize_patch_faces(brep, sector, candidate)?;
    brep.faces[split.strip_face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(candidate.strip_candidate.surface.clone()),
        orientation: true,
    });
    brep.faces[split.inner_face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(candidate.inner_patch.surface.clone()),
        orientation: true,
    });
    brep.half_edges[strip.cap_half_edge].pcurve =
        Some(candidate.strip_candidate.outer_pcurve.clone());
    if let Some(curve) = &candidate.support_authority_boundary {
        brep.edges[strip.edge].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve.clone()),
            t_start: curve.parameter_range().0,
            t_end: curve.parameter_range().1,
        });
    }
    brep.half_edges[split.interface_strip_half_edge].pcurve =
        Some(reversed_pcurve(&candidate.interface.strip_pcurve));
    brep.half_edges[split.interface_inner_half_edge].pcurve =
        Some(candidate.inner_patch.interface_pcurve.clone());

    if let Err(err) = brep.validate() {
        return Err(ExactBoundaryMultiPatchDiagnostic::commit_validation_failed(
            sector,
            strip,
            Some(err.to_string().len() as f64),
        ));
    }

    Ok(ExactBoundaryMultiPatchCommitReport {
        preflight,
        strip_face: split.strip_face,
        inner_face: split.inner_face,
        internal_interface_edges: vec![split.interface_edge],
        generated_connector_edges: vec![split.right_connector_edge, split.left_connector_edge],
        committed_pcurve_count: 3,
    })
}

fn first_invalid_surface_name(candidate: &ExactBoundaryMultiPatchCandidate) -> String {
    if !surface_layout_is_valid(&candidate.strip_candidate.surface) {
        candidate.strip_inventory.patch_name.clone()
    } else {
        candidate.inner_patch.patch_name.clone()
    }
}

fn reversed_pcurve(pcurve: &Pcurve) -> Pcurve {
    let start = pcurve.point_at(pcurve.t_end);
    let end = pcurve.point_at(pcurve.t_start);
    Pcurve::linear(
        (start.x, start.y),
        (end.x, end.y),
        pcurve.t_start,
        pcurve.t_end,
    )
}

fn surface_layout_is_valid(surface: &NurbsSurface) -> bool {
    let expected = surface.num_u.saturating_mul(surface.num_v);
    expected > 0
        && surface.control_points.len() == expected
        && surface.weights.len() == expected
        && surface.knots_u.len() == surface.num_u + surface.degree_u + 1
        && surface.knots_v.len() == surface.num_v + surface.degree_v + 1
        && surface
            .weights
            .iter()
            .all(|weight| weight.is_finite() && *weight > 0.0)
        && surface
            .control_points
            .iter()
            .all(|point| point.x.is_finite() && point.y.is_finite() && point.z.is_finite())
}

fn curve_layout_is_valid(curve: &NurbsCurve) -> bool {
    !curve.control_points.is_empty()
        && curve.weights.len() == curve.control_points.len()
        && curve.knots.len() == curve.control_points.len() + curve.degree + 1
        && curve
            .weights
            .iter()
            .all(|weight| weight.is_finite() && *weight > 0.0)
        && curve
            .control_points
            .iter()
            .all(|point| point.x.is_finite() && point.y.is_finite() && point.z.is_finite())
}

fn pcurve_is_valid(pcurve: &Pcurve) -> bool {
    pcurve.t_start.is_finite()
        && pcurve.t_end.is_finite()
        && pcurve.t_end > pcurve.t_start
        && pcurve.point_at(pcurve.t_start).x.is_finite()
        && pcurve.point_at(pcurve.t_end).x.is_finite()
        && pcurve.point_at(pcurve.t_start).y.is_finite()
        && pcurve.point_at(pcurve.t_end).y.is_finite()
}

fn generated_spoke_count(surface: &NurbsSurface) -> usize {
    if surface.num_v < 2 || surface.control_points.len() < surface.num_u.saturating_mul(2) {
        return 0;
    }
    let pole_row_start = surface.num_v - 1;
    let Some(reference) = surface.control_points.get(pole_row_start).copied() else {
        return 0;
    };
    let mut count = 0;
    for index in 0..surface.num_u {
        let point_index = index * surface.num_v + pole_row_start;
        let Some(point) = surface.control_points.get(point_index) else {
            return 0;
        };
        if point.distance_to(&reference) <= 1.0e-7 {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
#[path = "multipatch_commit_tests.rs"]
mod tests;
