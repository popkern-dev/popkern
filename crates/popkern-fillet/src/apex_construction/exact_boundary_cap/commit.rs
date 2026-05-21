//! Transactional BRep commit boundary for exact-boundary apex cap sectors.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle};

use super::pcurves::PreparedExactBoundaryPcurves;
use super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};

/// Result of committing exact-boundary sector candidates into live topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ExactBoundaryCommitReport {
    /// Number of sector faces assigned a solved NURBS surface.
    pub committed_sector_count: usize,
    /// Number of cap/support pcurve pairs written to live half-edges.
    pub committed_pcurve_pair_count: usize,
}

/// Stable decline reason for the BRep commit boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryCommitDeclineReason {
    /// The aggregate cap face handle is stale.
    CapFaceMissing,
    /// Candidate sector inventory does not match solved surfaces.
    SurfaceInventoryMismatch,
    /// Candidate sector inventory does not match prepared pcurve pairs.
    PcurveInventoryMismatch,
    /// Candidate topology still has one aggregate cap face instead of sector faces.
    SectorTopologyMissing,
    /// A sector face handle is stale.
    SectorFaceMissing,
    /// A cap-side half-edge to receive a pcurve is stale.
    CapHalfEdgeMissing,
    /// A support-side half-edge to receive a pcurve is stale.
    SupportHalfEdgeMissing,
}

impl ExactBoundaryCommitDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::CapFaceMissing => "exact_boundary_cap.commit.cap_face_missing",
            Self::SurfaceInventoryMismatch => {
                "exact_boundary_cap.commit.surface_inventory_mismatch"
            }
            Self::PcurveInventoryMismatch => "exact_boundary_cap.commit.pcurve_inventory_mismatch",
            Self::SectorTopologyMissing => "exact_boundary_cap.commit.sector_topology_missing",
            Self::SectorFaceMissing => "exact_boundary_cap.commit.sector_face_missing",
            Self::CapHalfEdgeMissing => "exact_boundary_cap.commit.cap_half_edge_missing",
            Self::SupportHalfEdgeMissing => "exact_boundary_cap.commit.support_half_edge_missing",
        }
    }
}

/// Structured exact-boundary commit decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryCommitDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryCommitDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Stable sector name when the failure is sector-scoped.
    pub sector_name: Option<String>,
    /// Cap or sector face being inspected.
    pub face: Option<FaceHandle>,
}

impl ExactBoundaryCommitDiagnostic {
    fn new(
        reason: ExactBoundaryCommitDeclineReason,
        sector_name: Option<String>,
        face: Option<FaceHandle>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_name,
            face,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.commit".to_string(),
                code: self.code.to_string(),
                sector_name: self
                    .sector_name
                    .unwrap_or_else(|| "exact_boundary_cap".to_string()),
                edge: None,
                support_half_edge: None,
                support_face: self.face,
                residual: None,
            },
        )
    }
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
    if result.is_err() {
        *brep = snapshot;
    }
    result
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
            brep.half_edges[sector.cap_half_edges[side_index]].pcurve =
                Some(pair.cap_pcurve.clone());
            brep.half_edges[sector.support_half_edges[side_index]].pcurve =
                Some(pair.support_pcurve.clone());
        }
    }
    Ok(ExactBoundaryCommitReport {
        committed_sector_count: sectors.len(),
        committed_pcurve_pair_count: pcurves.pairs.len(),
    })
}

fn solved_surface<'a>(
    solved: &'a ExactBoundarySolvedSurfaces,
    sector: &ExactBoundaryCapSectorInput,
) -> Option<&'a ExactBoundarySectorSurface> {
    solved
        .surfaces
        .iter()
        .find(|surface| surface.sector_index == sector.sector_index)
}
