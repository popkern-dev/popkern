//! Input inventory validation for the exact-boundary sector commit.

use popkern_core::topology::{BRep, FaceHandle};

use super::super::pcurves::PreparedExactBoundaryPcurves;
use super::super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::super::solve::ExactBoundarySolvedSurfaces;
use super::{ExactBoundaryCommitDeclineReason, ExactBoundaryCommitDiagnostic};

#[allow(clippy::result_large_err)]
pub(super) fn validate_commit_inputs(
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
