//! Transactional BRep commit boundary for exact-boundary apex cap sectors.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, SharedTrimEdgeFailureReason, SharedTrimEdgeInput,
    TopologyRole, commit_shared_trim_edge_transactionally,
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
            let (support_half_edge, support_face) =
                resolved_support_side_for_commit(brep, pair.side, cap_half_edge, support_half_edge)
                    .ok_or_else(|| {
                        ExactBoundaryCommitDiagnostic::new(
                            ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing,
                            Some(sector.sector_name.clone()),
                            Some(sector.cap_face),
                        )
                    })?;
            let support_surface = brep
                .faces
                .get(support_face)
                .and_then(|face| face.surface.as_ref())
                .map(|surface| surface.geometry.clone())
                .ok_or_else(|| {
                    ExactBoundaryCommitDiagnostic::new(
                        ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing,
                        Some(sector.sector_name.clone()),
                        Some(sector.cap_face),
                    )
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
                if support_edge != cap_edge {
                    if let Some(support_curve) = curve_from_pcurve_in_edge_direction(
                        brep,
                        support_edge,
                        support_half_edge,
                        &support_pcurve,
                        &support_surface,
                    ) {
                        brep.edges[support_edge].curve = Some(support_curve);
                    }
                }
                continue;
            }
            let curve = canonical_shared_trim_curve(
                brep,
                pair.edge,
                cap_half_edge,
                &pair.cap_pcurve,
                &cap_surface,
                support_half_edge,
                &support_pcurve,
                &support_surface,
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

fn resolved_support_side_for_commit(
    brep: &BRep,
    side: ExactBoundaryOuterSide,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    let support_face = face_for_half_edge(brep, support_half_edge)?;
    if side == ExactBoundaryOuterSide::FilletEnd
        && !half_edge_has_cylinder_surface(brep, support_half_edge)
        && let Some(cylinder_side) = same_oriented_cylinder_support_side(brep, cap_half_edge)
    {
        return Some(cylinder_side);
    }
    Some((support_half_edge, support_face))
}

fn half_edge_has_cylinder_surface(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    face_for_half_edge(brep, half_edge)
        .and_then(|face| brep.faces.get(face))
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
}

fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if same_oriented_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
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

fn same_oriented_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != support_start {
        return false;
    }
    half_edge_destination_vertex(brep, cap_side)
        .zip(half_edge_destination_vertex(brep, support_side))
        .is_some_and(|(cap_end, support_end)| cap_end == support_end)
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let loop_ = brep.half_edges.get(half_edge)?.loop_;
    brep.loops.get(loop_).map(|loop_| loop_.face)
}

fn half_edge_destination_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<popkern_core::topology::VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin| twin.origin)
}

fn project_curve_to_half_edge_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = 33usize;
    let mut uv_points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let point = curve.geometry.point_at(curve_t);
        let (u, v) = surface.closest_point(&point)?;
        uv_points.push(Point3::new(u, v, 0.0));
        parameters.push(fraction);
    }
    let degree = 5.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &parameters, degree)
        .ok()
        .map(Pcurve::from_nurbs)
}

fn half_edge_curve_forward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Option<bool> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&curve_start) <= origin.distance_to(&curve_end))
}

fn canonical_shared_trim_curve(
    brep: &BRep,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    cap_pcurve: &Pcurve,
    cap_surface: &SurfaceGeometry,
    support_half_edge: HalfEdgeHandle,
    support_pcurve: &Pcurve,
    support_surface: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let cap_curve =
        curve_from_pcurve_in_edge_direction(brep, edge, cap_half_edge, cap_pcurve, cap_surface)
            .and_then(|curve| orient_curve_to_edge_vertices(brep, edge, curve));
    let support_curve = curve_from_pcurve_in_edge_direction(
        brep,
        edge,
        support_half_edge,
        support_pcurve,
        support_surface,
    )
    .and_then(|curve| orient_curve_to_edge_vertices(brep, edge, curve));
    match (cap_curve, support_curve) {
        (Some(cap_curve), Some(support_curve)) => {
            let cap_endpoint_gap = edge_endpoint_gap(brep, edge, &cap_curve)?;
            let support_endpoint_gap = edge_endpoint_gap(brep, edge, &support_curve)?;
            if matches!(support_surface, SurfaceGeometry::Cylinder { .. })
                && support_endpoint_gap <= cap_endpoint_gap.max(1.0e-6)
            {
                return Some(support_curve);
            }
            if cap_endpoint_gap <= support_endpoint_gap {
                Some(cap_curve)
            } else {
                Some(support_curve)
            }
        }
        (Some(curve), None) | (None, Some(curve)) => Some(curve),
        (None, None) => None,
    }
}

fn reversed_pcurve(pcurve: &Pcurve) -> Pcurve {
    let (t_start, t_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut uv_points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = t_end + (t_start - t_end) * fraction;
        uv_points.push(pcurve.point_at(t));
        parameters.push(fraction);
    }
    let degree = 5.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &parameters, degree)
        .ok()
        .map(Pcurve::from_nurbs)
        .unwrap_or_else(|| pcurve.clone())
}

fn orient_curve_to_edge_vertices(
    brep: &BRep,
    edge: EdgeHandle,
    curve: CurveOnEdge,
) -> Option<CurveOnEdge> {
    let direct = edge_endpoint_gap(brep, edge, &curve)?;
    let reversed = reversed_curve_on_edge(&curve)?;
    let reverse = edge_endpoint_gap(brep, edge, &reversed)?;
    Some(if reverse < direct { reversed } else { curve })
}

fn reversed_curve_on_edge(curve: &CurveOnEdge) -> Option<CurveOnEdge> {
    let samples = 17usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = curve.t_end + (curve.t_start - curve.t_end) * fraction;
        points.push(curve.geometry.point_at(t));
        parameters.push(fraction);
    }
    let degree = 5.min(points.len().saturating_sub(1));
    let geometry = popkern_core::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)
        .or_else(|| line_curve_from_endpoints(&points))?;
    Some(CurveOnEdge {
        geometry,
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn edge_endpoint_gap(brep: &BRep, edge: EdgeHandle, curve: &CurveOnEdge) -> Option<f64> {
    let canonical = brep.edges.get(edge)?.half_edge;
    let start = brep
        .vertices
        .get(brep.half_edges.get(canonical)?.origin)?
        .position;
    let end = brep
        .vertices
        .get(
            brep.half_edges
                .get(brep.half_edges.get(canonical)?.twin)?
                .origin,
        )?
        .position;
    Some(
        curve
            .geometry
            .point_at(curve.t_start)
            .distance_to(&start)
            .max(curve.geometry.point_at(curve.t_end).distance_to(&end)),
    )
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
    let degree = 5.min(points.len().saturating_sub(1));
    let curve = popkern_core::geometry::fitting::interpolate_curve(&points, &parameters, degree)
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
