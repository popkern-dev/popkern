//! Transactional BRep commit boundary for exact-boundary apex cap sectors.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{
    BRep, FaceHandle, SharedTrimEdgeInput, TopologyRole, commit_shared_trim_edge_transactionally,
};

#[path = "commit_helpers.rs"]
mod commit_helpers;
#[path = "commit_post_validation.rs"]
mod commit_post_validation;
use super::cap_spoke_commit::commit_cap_spoke_shared_trims;
pub(crate) use super::commit_diagnostics::{
    ExactBoundaryCommitDeclineReason, ExactBoundaryCommitDiagnostic,
};
use super::pcurves::PreparedExactBoundaryPcurves;
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::solve::ExactBoundarySolvedSurfaces;
use super::support_subdivision::materialize_support_subdivision;
use commit_helpers::{
    SharedTrimCurveSelection, canonical_shared_trim_curve, curve_from_pcurve_in_edge_direction,
    normalize_semantic_line_edge_curve, project_curve_to_half_edge_pcurve,
    resolved_support_side_for_commit, reversed_pcurve, solved_surface,
};
pub(crate) use commit_post_validation::post_commit_validation_diagnostic;

/// Result of committing exact-boundary sector candidates into live topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ExactBoundaryCommitReport {
    /// Number of sector faces assigned a solved NURBS surface.
    pub committed_sector_count: usize,
    /// Number of cap/support pcurve pairs written to live half-edges.
    pub committed_pcurve_pair_count: usize,
}

/// Commit solved exact-boundary sector surfaces and paired pcurves transactionally.
pub(crate) fn commit_exact_boundary_sector_surfaces(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
) -> Result<ExactBoundaryCommitReport, ExactBoundaryCommitDiagnostic> {
    validate_commit_inputs(brep, cap_face, sectors, solved, pcurves)?;
    let snapshot = brep.clone();
    let result = commit_validated_inputs(brep, sectors, solved, pcurves);
    match result {
        Ok(report) => Ok(report),
        Err(diagnostic) => {
            *brep = snapshot;
            Err(diagnostic)
        }
    }
}

fn validate_commit_inputs(
    brep: &BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
) -> Result<(), ExactBoundaryCommitDiagnostic> {
    if !brep.faces.contains_key(cap_face) {
        return Err(ExactBoundaryCommitDiagnostic::new(
            ExactBoundaryCommitDeclineReason::CapFaceMissing,
            None,
            Some(cap_face),
        ));
    }
    if solved.validated_surface_count() != sectors.len() {
        return Err(ExactBoundaryCommitDiagnostic::new(
            ExactBoundaryCommitDeclineReason::SurfaceInventoryMismatch,
            None,
            Some(cap_face),
        ));
    }
    if pcurves.validated_pair_count() != sectors.len() * 2 {
        return Err(ExactBoundaryCommitDiagnostic::new(
            ExactBoundaryCommitDeclineReason::PcurveInventoryMismatch,
            None,
            Some(cap_face),
        ));
    }
    for sector in sectors {
        if !has_unique_surface_for_sector(solved, sector) {
            return Err(ExactBoundaryCommitDiagnostic::new(
                ExactBoundaryCommitDeclineReason::SurfaceInventoryMismatch,
                Some(sector.sector_name.clone()),
                Some(sector.cap_face),
            ));
        }
        if !has_required_pcurve_pairs_for_sector(pcurves, sector) {
            return Err(ExactBoundaryCommitDiagnostic::new(
                ExactBoundaryCommitDeclineReason::PcurveInventoryMismatch,
                Some(sector.sector_name.clone()),
                Some(sector.cap_face),
            ));
        }
        if (sectors.len() == 1 && sector.cap_face == cap_face)
            || !sector_faces_are_unique(sectors, sector)
        {
            return Err(ExactBoundaryCommitDiagnostic::new(
                ExactBoundaryCommitDeclineReason::SectorTopologyMissing,
                Some(sector.sector_name.clone()),
                Some(sector.cap_face),
            ));
        }
        if !brep.faces.contains_key(sector.cap_face) {
            return Err(ExactBoundaryCommitDiagnostic::new(
                ExactBoundaryCommitDeclineReason::SectorFaceMissing,
                Some(sector.sector_name.clone()),
                Some(sector.cap_face),
            ));
        }
        for half_edge in sector.cap_half_edges {
            if !brep.half_edges.contains_key(half_edge) {
                return Err(ExactBoundaryCommitDiagnostic::new(
                    ExactBoundaryCommitDeclineReason::CapHalfEdgeMissing,
                    Some(sector.sector_name.clone()),
                    Some(sector.cap_face),
                ));
            }
        }
        for half_edge in sector.support_half_edges {
            if !brep.half_edges.contains_key(half_edge) {
                return Err(ExactBoundaryCommitDiagnostic::new(
                    ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing,
                    Some(sector.sector_name.clone()),
                    Some(sector.cap_face),
                ));
            }
        }
    }
    Ok(())
}

fn sector_faces_are_unique(
    sectors: &[ExactBoundaryCapSectorInput],
    sector: &ExactBoundaryCapSectorInput,
) -> bool {
    sectors
        .iter()
        .filter(|candidate| candidate.cap_face == sector.cap_face)
        .count()
        == 1
}

fn has_unique_surface_for_sector(
    solved: &ExactBoundarySolvedSurfaces,
    sector: &ExactBoundaryCapSectorInput,
) -> bool {
    solved
        .surfaces
        .iter()
        .filter(|surface| {
            surface.sector_index == sector.sector_index && surface.sector_name == sector.sector_name
        })
        .count()
        == 1
}

fn has_required_pcurve_pairs_for_sector(
    pcurves: &PreparedExactBoundaryPcurves,
    sector: &ExactBoundaryCapSectorInput,
) -> bool {
    [
        ExactBoundaryOuterSide::FilletEnd,
        ExactBoundaryOuterSide::FaceEdge,
    ]
    .into_iter()
    .all(|side| {
        pcurves
            .pairs
            .iter()
            .filter(|pair| {
                pair.sector_index == sector.sector_index
                    && pair.sector_name == sector.sector_name
                    && pair.side == side
            })
            .count()
            == 1
    })
}

fn commit_validated_inputs(
    brep: &mut BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
) -> Result<ExactBoundaryCommitReport, ExactBoundaryCommitDiagnostic> {
    for sector in sectors {
        let surface = solved_surface(solved, sector).ok_or_else(|| {
            ExactBoundaryCommitDiagnostic::new(
                ExactBoundaryCommitDeclineReason::SurfaceInventoryMismatch,
                Some(sector.sector_name.clone()),
                Some(sector.cap_face),
            )
        })?;
        brep.faces[sector.cap_face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(surface.surface.clone()),
            orientation: true,
        });

        for pair in pcurves
            .pairs
            .iter()
            .filter(|pair| pair.sector_index == sector.sector_index)
        {
            let side_index = match pair.side {
                ExactBoundaryOuterSide::FilletEnd => 0,
                ExactBoundaryOuterSide::FaceEdge => 1,
            };
            let cap_half_edge = sector.cap_half_edges[side_index];
            let support_half_edge = pair.support_half_edge;
            normalize_semantic_line_edge_curve(brep, pair.edge);
            brep.half_edges[cap_half_edge].pcurve = Some(pair.cap_pcurve.clone());
            let cap_surface = brep.faces[sector.cap_face]
                .surface
                .as_ref()
                .map(|surface| surface.geometry.clone())
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SectorFaceMissing,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
                })?;
            let (support_half_edge, support_face) = resolved_support_side_for_commit(
                brep,
                pair.side,
                cap_half_edge,
                support_half_edge,
                &pair.support_surface,
            )
            .ok_or_else(|| {
                ExactBoundaryCommitDiagnostic::support_side_missing(
                    sector,
                    side_index,
                    cap_half_edge,
                    pair.support_half_edge,
                    None,
                )
            })?;
            let mut support_half_edge = support_half_edge;
            let mut support_face = support_face;
            let mut support_surface = brep
                .faces
                .get(support_face)
                .and_then(|face| face.surface.as_ref())
                .map(|surface| surface.geometry.clone())
                .or_else(|| {
                    assign_loop_valid_support_surface(
                        brep,
                        support_face,
                        pair.support_surface.clone(),
                    )
                });
            let support_residual =
                support_loop_surface_residual(brep, support_face, &pair.support_surface);
            if (support_surface.is_none()
                || support_residual
                    .is_some_and(|residual| residual > brep.tolerance.linear().max(1.0e-6)))
                && let Some(subdivision) = materialize_support_subdivision(
                    brep,
                    cap_half_edge,
                    support_half_edge,
                    support_face,
                    pair.support_surface.clone(),
                )
            {
                support_half_edge = subdivision.support_half_edge;
                support_face = subdivision.support_face;
                support_surface = Some(pair.support_surface.clone());
            }
            let support_surface = support_surface.ok_or_else(|| {
                if let Some(residual) = support_residual
                    && residual > brep.tolerance.linear().max(1.0e-6)
                {
                    ExactBoundaryCommitDiagnostic::support_face_subdivision_required(
                        sector,
                        side_index,
                        cap_half_edge,
                        support_half_edge,
                        support_face,
                        residual,
                    )
                } else {
                    ExactBoundaryCommitDiagnostic::support_surface_missing(
                        sector,
                        side_index,
                        cap_half_edge,
                        support_half_edge,
                        support_face,
                        support_residual,
                    )
                }
            })?;
            let cap_edge = brep.half_edges[cap_half_edge].edge;
            let support_edge = brep.half_edges[support_half_edge].edge;
            let support_is_cap_twin = brep.half_edges[cap_half_edge].twin == support_half_edge;
            let support_is_residual_patch = brep
                .half_edges
                .get(support_half_edge)
                .and_then(|half_edge| brep.loops.get(half_edge.loop_))
                .and_then(|loop_| brep.faces.get(loop_.face))
                .is_some_and(|face| face.role == TopologyRole::ConstructionResidual);
            let support_is_cap_face = support_face == sector.cap_face;
            let support_pcurve = if support_half_edge == pair.support_half_edge {
                pair.support_pcurve.clone()
            } else {
                let cap_curve = curve_from_pcurve_in_edge_direction(
                    brep,
                    cap_edge,
                    cap_half_edge,
                    &pair.cap_pcurve,
                    &cap_surface,
                )
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
                })?;
                project_curve_to_half_edge_pcurve(
                    brep,
                    support_half_edge,
                    &cap_curve,
                    &support_surface,
                )
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
                })?
            };
            brep.half_edges[cap_half_edge].pcurve = Some(pair.cap_pcurve.clone());
            brep.half_edges[support_half_edge].pcurve = Some(support_pcurve.clone());
            if cap_edge != support_edge
                || !support_is_cap_twin
                || support_is_residual_patch
                || support_is_cap_face
                || pair.side == ExactBoundaryOuterSide::FaceEdge
            {
                if let Some(cap_curve) = curve_from_pcurve_in_edge_direction(
                    brep,
                    cap_edge,
                    cap_half_edge,
                    &pair.cap_pcurve,
                    &cap_surface,
                ) {
                    brep.edges[cap_edge].curve = Some(cap_curve);
                }
                if support_edge != cap_edge
                    && let Some(support_curve) = curve_from_pcurve_in_edge_direction(
                        brep,
                        support_edge,
                        support_half_edge,
                        &support_pcurve,
                        &support_surface,
                    )
                {
                    brep.edges[support_edge].curve = Some(support_curve);
                }
                continue;
            }
            let curve = canonical_shared_trim_curve(
                brep,
                SharedTrimCurveSelection {
                    edge: pair.edge,
                    cap_half_edge,
                    cap_pcurve: &pair.cap_pcurve,
                    cap_surface: &cap_surface,
                    support_half_edge,
                    support_pcurve: &support_pcurve,
                    support_surface: &support_surface,
                },
            )
            .ok_or_else(|| {
                ExactBoundaryCommitDiagnostic::new(
                    ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed,
                    Some(sector.sector_name.clone()),
                    Some(sector.cap_face),
                )
            })?;
            if let Err(failure) = commit_shared_trim_edge_transactionally(
                brep,
                SharedTrimEdgeInput {
                    edge: pair.edge,
                    half_edge_a: cap_half_edge,
                    half_edge_b: support_half_edge,
                    curve: curve.clone(),
                    pcurve_a: pair.cap_pcurve.clone(),
                    pcurve_b: support_pcurve.clone(),
                    tolerance: 1.0e-6,
                    validate_brep: false,
                    provenance: "exact_boundary_cap_commit_shared_trim",
                },
            ) {
                let reversed_support_pcurve = reversed_pcurve(&support_pcurve);
                if commit_shared_trim_edge_transactionally(
                    brep,
                    SharedTrimEdgeInput {
                        edge: pair.edge,
                        half_edge_a: cap_half_edge,
                        half_edge_b: support_half_edge,
                        curve,
                        pcurve_a: pair.cap_pcurve.clone(),
                        pcurve_b: reversed_support_pcurve,
                        tolerance: 1.0e-6,
                        validate_brep: false,
                        provenance: "exact_boundary_cap_commit_shared_trim",
                    },
                )
                .is_ok()
                {
                    continue;
                }
                return Err(
                    ExactBoundaryCommitDiagnostic::shared_trim_validation_failed(
                        sector, side_index, failure,
                    ),
                );
            }
        }
    }
    commit_cap_spoke_shared_trims(brep, sectors)?;
    Ok(ExactBoundaryCommitReport {
        committed_sector_count: sectors.len(),
        committed_pcurve_pair_count: pcurves.pairs.len(),
    })
}

fn assign_loop_valid_support_surface(
    brep: &mut BRep,
    support_face: FaceHandle,
    support_surface: SurfaceGeometry,
) -> Option<SurfaceGeometry> {
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let residual = support_loop_surface_residual(brep, support_face, &support_surface)?;
    if residual > tolerance {
        return None;
    }
    let face = brep.faces.get_mut(support_face)?;
    face.surface = Some(SurfaceOnFace {
        geometry: support_surface.clone(),
        orientation: true,
    });
    Some(support_surface)
}

fn support_loop_surface_residual(
    brep: &BRep,
    support_face: FaceHandle,
    support_surface: &SurfaceGeometry,
) -> Option<f64> {
    let loop_handle = brep.faces.get(support_face)?.outer_loop;
    collect_loop_half_edges(brep, loop_handle)?
        .into_iter()
        .filter_map(|half_edge| {
            let edge = brep.edges.get(brep.half_edges.get(half_edge)?.edge)?;
            let curve = edge.curve.as_ref()?;
            Some(
                edge_surface_deviation(
                    &curve.geometry,
                    (curve.t_start, curve.t_end),
                    support_surface,
                    9,
                )
                .max_deviation,
            )
        })
        .reduce(f64::max)
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<popkern_core::topology::HalfEdgeHandle>> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}
