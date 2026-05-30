//! Live commit bridge for exact-boundary multipatch route attempts.

use popkern_core::error::ValidationError;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryRouteDecisionStep};

use super::super::commit::{
    commit_exact_boundary_sector_surfaces, post_commit_validation_diagnostic,
};
use super::super::multipatch::ExactBoundaryMultiPatchCandidate;
use super::super::multipatch::commit::commit_multipatch_candidate_transactionally;
use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::sector::ExactBoundaryCapSectorInput;
use super::super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};
use super::candidate_continuity::candidate_selected_strip_continuity_diagnostic;
use super::multipatch_route::ExactBoundaryMultipatchRouteAttempt;
use crate::construct::cap_outer_acceptance::cap_support_continuity_diagnostic;

pub(super) enum ExactBoundaryMultipatchCommitOutcome {
    Accepted {
        faces: Vec<FaceHandle>,
    },
    Declined {
        step: ExactBoundaryRouteDecisionStep,
        diagnostic: ApexConstructionDiagnostic,
    },
}

pub(super) fn commit_multipatch_routes_transactionally(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    baseline_validation_errors: &[ValidationError],
    attempts: &[(
        &ExactBoundaryMultipatchRouteAttempt,
        &ExactBoundaryMultiPatchCandidate,
    )],
) -> ExactBoundaryMultipatchCommitOutcome {
    let mut inner_faces = Vec::with_capacity(attempts.len());
    let mut strip_faces = Vec::with_capacity(attempts.len());
    for (attempt, candidate) in attempts {
        if let Some(diagnostic) =
            candidate_selected_strip_continuity_diagnostic(&attempt.sector, candidate, pcurves)
        {
            return declined(apex_diagnostic_code(&diagnostic), diagnostic);
        }
        let multipatch_commit =
            match commit_multipatch_candidate_transactionally(brep, &attempt.sector, candidate) {
                Ok(report) => report,
                Err(diagnostic) => {
                    return declined(diagnostic.code.clone(), diagnostic.into_apex_diagnostic());
                }
            };
        bind_inner_face_retained_boundaries(
            brep,
            &attempt.sector,
            candidate,
            solved_surfaces,
            pcurves,
            multipatch_commit.inner_face,
        );
        bind_strip_face_support_boundary(brep, &attempt.sector, candidate, pcurves);
        strip_faces.push(multipatch_commit.strip_face);
        inner_faces.push(multipatch_commit.inner_face);
    }
    if let Err(diagnostic) = commit_unaffected_one_patch_sectors(
        brep,
        cap_face,
        sectors,
        solved_surfaces,
        pcurves,
        attempts
            .iter()
            .map(|(attempt, _)| attempt.sector.sector_name.as_str()),
    ) {
        return declined(
            diagnostic.code.to_string(),
            diagnostic.into_apex_diagnostic(),
        );
    }
    if let Some(diagnostic) = post_commit_validation_diagnostic(
        baseline_validation_errors,
        &brep.validate_all(),
        brep.tolerance.linear(),
    ) {
        return declined(
            diagnostic.code.to_string(),
            diagnostic.into_apex_diagnostic(),
        );
    }
    let mut faces = sectors
        .iter()
        .map(|sector| sector.cap_face)
        .collect::<Vec<_>>();
    faces.extend(inner_faces.iter().copied());
    faces.extend(strip_faces.iter().copied());
    if let Some(diagnostic) = cap_support_continuity_diagnostic(brep, &faces, 5.0e-3, 5.0e-2) {
        return declined(apex_diagnostic_code(&diagnostic), diagnostic);
    }
    ExactBoundaryMultipatchCommitOutcome::Accepted { faces }
}

fn bind_strip_face_support_boundary(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
    pcurves: &PreparedExactBoundaryPcurves,
) {
    let Some(pair) = pcurves.pairs.iter().find(|pair| {
        pair.sector_name == sector.sector_name && pair.side == candidate.strip_inventory.side
    }) else {
        return;
    };
    let side_index = match candidate.strip_inventory.side {
        super::super::sector::ExactBoundaryOuterSide::FilletEnd => 0,
        super::super::sector::ExactBoundaryOuterSide::FaceEdge => 1,
    };
    let cap_half_edge = sector.cap_half_edges[side_index];
    let cap_surface = SurfaceGeometry::Nurbs(candidate.strip_candidate.surface.clone());
    brep.half_edges[cap_half_edge].pcurve = Some(orient_pcurve_to_half_edge_vertices(
        brep,
        cap_half_edge,
        candidate.strip_candidate.outer_pcurve.clone(),
        &cap_surface,
    ));
    if brep.half_edges.contains_key(pair.support_half_edge) {
        brep.half_edges[pair.support_half_edge].pcurve = Some(orient_pcurve_to_half_edge_vertices(
            brep,
            pair.support_half_edge,
            pair.support_pcurve.clone(),
            &pair.support_surface,
        ));
    }
    if brep.faces.contains_key(pair.support_face) {
        brep.faces[pair.support_face].surface = Some(SurfaceOnFace {
            geometry: pair.support_surface.clone(),
            orientation: true,
        });
    }
}

fn orient_pcurve_to_half_edge_vertices(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: Pcurve,
    surface: &SurfaceGeometry,
) -> Pcurve {
    let Some(direct_gap) = pcurve_endpoint_gap(brep, half_edge, &pcurve, surface) else {
        return pcurve;
    };
    let reversed = reversed_pcurve(&pcurve);
    let Some(reversed_gap) = pcurve_endpoint_gap(brep, half_edge, &reversed, surface) else {
        return pcurve;
    };
    if reversed_gap < direct_gap {
        reversed
    } else {
        pcurve
    }
}

fn pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let destination = brep
        .vertices
        .get(half_edge_destination_vertex(brep, half_edge)?)?
        .position;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    Some(
        start
            .distance_to(&origin)
            .max(end.distance_to(&destination)),
    )
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin| twin.origin)
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

fn bind_inner_face_retained_boundaries(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    inner_face: FaceHandle,
) {
    let Some(surface) = solved_surfaces
        .surfaces
        .iter()
        .find(|surface| surface.sector_name == sector.sector_name)
    else {
        return;
    };
    brep.faces[inner_face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(surface.surface.clone()),
        orientation: true,
    });
    for pair in pcurves
        .pairs
        .iter()
        .filter(|pair| pair.sector_name == sector.sector_name)
        .filter(|pair| pair.side != candidate.strip_inventory.side)
    {
        let side_index = match pair.side {
            super::super::sector::ExactBoundaryOuterSide::FilletEnd => 0,
            super::super::sector::ExactBoundaryOuterSide::FaceEdge => 1,
        };
        brep.half_edges[sector.cap_half_edges[side_index]].pcurve = Some(pair.cap_pcurve.clone());
        if brep.half_edges.contains_key(pair.support_half_edge) {
            brep.half_edges[pair.support_half_edge].pcurve = Some(pair.support_pcurve.clone());
        }
    }
}

fn commit_unaffected_one_patch_sectors<'a>(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved_surfaces: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    multipatch_sector_names: impl Iterator<Item = &'a str>,
) -> Result<(), super::super::commit::ExactBoundaryCommitDiagnostic> {
    let multipatch_sector_names = multipatch_sector_names.collect::<Vec<_>>();
    let mut remaining_sectors = sectors
        .iter()
        .filter(|sector| {
            !multipatch_sector_names
                .iter()
                .any(|name| *name == sector.sector_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    for (index, sector) in remaining_sectors.iter_mut().enumerate() {
        sector.sector_index = index;
    }
    if remaining_sectors.is_empty() {
        return Ok(());
    }
    let remaining_solved = ExactBoundarySolvedSurfaces {
        surfaces: filter_surfaces_for_sectors(&solved_surfaces.surfaces, &remaining_sectors),
    };
    let remaining_pcurves = PreparedExactBoundaryPcurves {
        pairs: filter_pcurves_for_sectors(&pcurves.pairs, &remaining_sectors),
    };
    let commit_cap_face = validation_cap_face(cap_face, sectors, &remaining_sectors);
    commit_exact_boundary_sector_surfaces(
        brep,
        commit_cap_face,
        &remaining_sectors,
        &remaining_solved,
        &remaining_pcurves,
    )
    .map(|_| ())
}

fn validation_cap_face(
    cap_face: FaceHandle,
    all_sectors: &[ExactBoundaryCapSectorInput],
    committed_sectors: &[ExactBoundaryCapSectorInput],
) -> FaceHandle {
    if committed_sectors.len() == 1 && committed_sectors[0].cap_face == cap_face {
        all_sectors
            .iter()
            .find(|sector| sector.cap_face != cap_face)
            .map(|sector| sector.cap_face)
            .unwrap_or(cap_face)
    } else {
        cap_face
    }
}

fn filter_surfaces_for_sectors(
    surfaces: &[ExactBoundarySectorSurface],
    sectors: &[ExactBoundaryCapSectorInput],
) -> Vec<ExactBoundarySectorSurface> {
    let mut filtered = surfaces
        .iter()
        .filter(|surface| {
            sectors
                .iter()
                .any(|sector| sector.sector_name == surface.sector_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    for surface in &mut filtered {
        if let Some(index) = sectors
            .iter()
            .position(|sector| sector.sector_name == surface.sector_name)
        {
            surface.sector_index = index;
        }
    }
    filtered
}

fn filter_pcurves_for_sectors(
    pairs: &[PreparedExactBoundaryPcurvePair],
    sectors: &[ExactBoundaryCapSectorInput],
) -> Vec<PreparedExactBoundaryPcurvePair> {
    let mut filtered = pairs
        .iter()
        .filter(|pair| {
            sectors
                .iter()
                .any(|sector| sector.sector_name == pair.sector_name)
        })
        .cloned()
        .collect::<Vec<_>>();
    for pair in &mut filtered {
        if let Some(index) = sectors
            .iter()
            .position(|sector| sector.sector_name == pair.sector_name)
        {
            pair.sector_index = index;
        }
    }
    filtered
}

fn declined(
    code: String,
    diagnostic: ApexConstructionDiagnostic,
) -> ExactBoundaryMultipatchCommitOutcome {
    ExactBoundaryMultipatchCommitOutcome::Declined {
        step: ExactBoundaryRouteDecisionStep::declined("exact_boundary.multipatch", Some(code)),
        diagnostic,
    }
}

fn apex_diagnostic_code(diagnostic: &ApexConstructionDiagnostic) -> String {
    match diagnostic {
        ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) => diagnostic.code.clone(),
        ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => diagnostic
            .continuity
            .decline_code
            .clone()
            .unwrap_or_else(|| "cap_outer_g1_exceeds_gate".to_string()),
        _ => "apex_construction_declined".to_string(),
    }
}
