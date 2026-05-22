//! Transactional BRep commit boundary for exact-boundary apex cap sectors.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, SharedTrimEdgeFailureReason, SharedTrimEdgeInput,
    commit_shared_trim_edge_transactionally,
};

use super::cap_spoke_commit::commit_cap_spoke_shared_trims;
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
    /// Atomic shared trim validation rejected the committed edge/pcurve triple.
    SharedTrimValidationFailed,
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
            Self::SharedTrimValidationFailed => {
                "exact_boundary_cap.commit.shared_trim_validation_failed"
            }
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
    /// Edge being inspected or committed.
    pub edge: Option<EdgeHandle>,
    /// Cap-side half-edge being committed.
    pub cap_half_edge: Option<HalfEdgeHandle>,
    /// Support-side half-edge being committed.
    pub support_half_edge: Option<HalfEdgeHandle>,
    /// Shared trim failure reason, when available.
    pub shared_trim_failure: Option<SharedTrimEdgeFailureReason>,
    /// Measured validation residual, when available.
    pub residual: Option<f64>,
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
            edge: None,
            cap_half_edge: None,
            support_half_edge: None,
            shared_trim_failure: None,
            residual: None,
        }
    }

    fn shared_trim_validation_failed(
        sector: &ExactBoundaryCapSectorInput,
        side_index: usize,
        failure: popkern_core::topology::SharedTrimEdgeFailure,
    ) -> Self {
        Self {
            reason: ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed,
            code: ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed.code(),
            sector_name: Some(sector.sector_name.clone()),
            face: Some(sector.cap_face),
            edge: Some(sector.outer_seams[side_index].edge),
            cap_half_edge: Some(sector.cap_half_edges[side_index]),
            support_half_edge: Some(sector.support_half_edges[side_index]),
            shared_trim_failure: Some(failure.reason),
            residual: Some(
                failure
                    .validation
                    .max_curve_surface_deviation
                    .max(failure.validation.max_side_agreement_deviation)
                    .max(failure.validation.max_endpoint_gap),
            ),
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
                edge: self.edge,
                support_half_edge: self.support_half_edge.or(self.cap_half_edge),
                support_face: self.face,
                residual: self.residual,
                shared_fact_ids: Vec::new(),
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
            normalize_semantic_line_edge_curve(brep, pair.edge);
            brep.half_edges[sector.cap_half_edges[side_index]].pcurve =
                Some(pair.cap_pcurve.clone());
            brep.half_edges[sector.support_half_edges[side_index]].pcurve =
                Some(pair.support_pcurve.clone());
            let cap_surface = brep.faces[sector.cap_face]
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SectorFaceMissing,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
                })?;
            let support_surface = brep
                .half_edges
                .get(sector.support_half_edges[side_index])
                .and_then(|half_edge| brep.loops.get(half_edge.loop_))
                .and_then(|loop_| brep.faces.get(loop_.face))
                .and_then(|face| face.surface.as_ref())
                .map(|surface| &surface.geometry)
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
                })?;
            let curve = curve_from_pcurve_in_edge_direction(
                brep,
                pair.edge,
                sector.support_half_edges[side_index],
                &pair.support_pcurve,
                support_surface,
            )
            .or_else(|| {
                curve_from_pcurve_in_edge_direction(
                    brep,
                    pair.edge,
                    sector.cap_half_edges[side_index],
                    &pair.cap_pcurve,
                    cap_surface,
                )
            })
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
                    half_edge_a: sector.cap_half_edges[side_index],
                    half_edge_b: sector.support_half_edges[side_index],
                    curve,
                    pcurve_a: pair.cap_pcurve.clone(),
                    pcurve_b: pair.support_pcurve.clone(),
                    tolerance: 1.0e-6,
                    validate_brep: false,
                    provenance: "exact_boundary_cap_commit_shared_trim",
                },
            ) {
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

fn curve_from_pcurve_in_edge_direction(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let edge_data = brep.edges.get(edge)?;
    let forward = half_edge == edge_data.half_edge;
    let (pc_start, pc_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let pcurve_fraction = if forward { fraction } else { 1.0 - fraction };
        let pcurve_t = pc_start + (pc_end - pc_start) * pcurve_fraction;
        let uv = pcurve.point_at(pcurve_t);
        points.push(surface.point_at(uv.x, uv.y));
        parameters.push(fraction);
    }
    let degree = 9.min(points.len().saturating_sub(1));
    let curve =
        popkern_core::geometry::fitting::fit_curve(&points, &parameters, degree, points.len())
            .ok()
            .map(CurveGeometry::Nurbs)
            .or_else(|| line_curve_from_endpoints(&points))?;
    Some(CurveOnEdge {
        geometry: curve,
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn line_curve_from_endpoints(points: &[Point3]) -> Option<CurveGeometry> {
    let start = *points.first()?;
    let end = *points.last()?;
    CurveOnEdge::line_segment(start, end, 1.0e-12).map(|curve| curve.geometry)
}

fn normalize_semantic_line_edge_curve(brep: &mut BRep, edge: popkern_core::topology::EdgeHandle) {
    let Some(curve) = brep
        .edges
        .get_mut(edge)
        .and_then(|edge| edge.curve.as_mut())
    else {
        return;
    };
    if let Some(line) = curve
        .geometry
        .canonical_line_in_range((curve.t_start, curve.t_end), 1.0e-8)
    {
        curve.geometry = line;
    }
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
