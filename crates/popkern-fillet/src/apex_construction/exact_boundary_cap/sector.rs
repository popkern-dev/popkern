//! Exact-boundary sector extraction for direct-apex Gregory caps.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::gregory_cap_boundary::half_edge_for_face_on_edge;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{TargetCapFace, TargetEdgeRole};

/// Sector-local side occupied by an immutable cap/support seam.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryOuterSide {
    /// Sector `u_min`: cap-vs-fillet support seam.
    FilletEnd,
    /// Sector `v_max`: cap-vs-bridge support seam.
    FaceEdge,
}

/// Live BRep ownership facts for one immutable outer seam.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryOuterSeam {
    /// Sector-local side.
    pub side: ExactBoundaryOuterSide,
    /// Shared topological edge.
    pub edge: EdgeHandle,
    /// Half-edge owned by the cap side.
    pub cap_half_edge: HalfEdgeHandle,
    /// Twin half-edge owned by the support side.
    pub support_half_edge: HalfEdgeHandle,
    /// Support face across the shared edge.
    pub support_face: FaceHandle,
    /// Whether the cap half-edge already carries a pcurve.
    pub cap_pcurve_present: bool,
    /// Whether the support half-edge already carries a pcurve.
    pub support_pcurve_present: bool,
}

/// Exact-boundary sector input extracted from already-authored topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryCapSectorInput {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable semantic name for diagnostics and future UI highlighting.
    pub sector_name: String,
    /// Fillet index owned by this sector.
    pub fillet_index: usize,
    /// Current cap face that owns the cap-side half-edges.
    pub cap_face: FaceHandle,
    /// Fillet-end edge for this sector.
    pub fillet_end_edge: EdgeHandle,
    /// Plane/support bridge edge for this sector.
    pub face_edge: EdgeHandle,
    /// Support face across [`Self::fillet_end_edge`].
    pub fillet_support_face: FaceHandle,
    /// Support face across [`Self::face_edge`].
    pub bridge_support_face: FaceHandle,
    /// Cap-side half-edges, `[FilletEnd, FaceEdge]`.
    pub cap_half_edges: [HalfEdgeHandle; 2],
    /// Support-side half-edges, `[FilletEnd, FaceEdge]`.
    pub support_half_edges: [HalfEdgeHandle; 2],
    /// Immutable outer seam ownership facts.
    pub outer_seams: [ExactBoundaryOuterSeam; 2],
}

/// Stable decline reason for sector extraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundarySectorDeclineReason {
    /// The plan is not a Gregory cap topology.
    NonGregoryCap,
    /// Planner order and sector inventory disagree.
    SectorInventoryMismatch,
    /// A sector edge handle is stale.
    StaleEdge,
    /// A sector edge has no exact 3D curve.
    MissingEdgeCurve,
    /// The edge is not bounded by the cap face.
    MissingCapHalfEdge,
    /// The cap half-edge has no live twin.
    MissingSupportHalfEdge,
    /// The support half-edge does not resolve to a non-cap face.
    MissingSupportFace,
    /// The two immutable outer seams collapsed to the same edge.
    DuplicateOuterSeam,
}

impl ExactBoundarySectorDeclineReason {
    /// Stable diagnostic code for API, GUI, and test assertions.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::NonGregoryCap => "exact_boundary_cap.sector.non_gregory_cap",
            Self::SectorInventoryMismatch => "exact_boundary_cap.sector.inventory_mismatch",
            Self::StaleEdge => "exact_boundary_cap.sector.stale_edge",
            Self::MissingEdgeCurve => "exact_boundary_cap.sector.boundary_missing_curve",
            Self::MissingCapHalfEdge => "exact_boundary_cap.sector.cap_half_edge_missing",
            Self::MissingSupportHalfEdge => "exact_boundary_cap.sector.support_half_edge_missing",
            Self::MissingSupportFace => "exact_boundary_cap.sector.support_face_missing",
            Self::DuplicateOuterSeam => "exact_boundary_cap.sector.duplicate_outer_seam",
        }
    }
}

/// Structured extraction decline with stable handles where available.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundarySectorDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundarySectorDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Planner-order sector index, if known.
    pub sector_index: Option<usize>,
    /// Stable sector name, if known.
    pub sector_name: Option<String>,
    /// Cap face being inspected.
    pub cap_face: Option<FaceHandle>,
    /// Edge being inspected.
    pub edge: Option<EdgeHandle>,
    /// Half-edge being inspected.
    pub half_edge: Option<HalfEdgeHandle>,
}

impl ExactBoundarySectorDiagnostic {
    fn new(
        reason: ExactBoundarySectorDeclineReason,
        sector_index: Option<usize>,
        sector_name: Option<String>,
        cap_face: Option<FaceHandle>,
        edge: Option<EdgeHandle>,
        half_edge: Option<HalfEdgeHandle>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_index,
            sector_name,
            cap_face,
            edge,
            half_edge,
        }
    }

    pub(crate) fn into_build_error(self) -> ApexBuildError {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "{} sector={:?} edge={:?} half_edge={:?}",
                self.code, self.sector_name, self.edge, self.half_edge
            ),
        }
    }
}

/// Extract planner-order exact-boundary sector facts from live cap topology.
pub(crate) fn extract_exact_boundary_cap_sectors(
    brep: &BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) -> Result<Vec<ExactBoundaryCapSectorInput>, ExactBoundarySectorDiagnostic> {
    let outer_indices = match &plan.target_topology.cap_face {
        TargetCapFace::Gregory {
            outer_loop_edge_indices,
        } => outer_loop_edge_indices,
        _ => {
            return Err(ExactBoundarySectorDiagnostic::new(
                ExactBoundarySectorDeclineReason::NonGregoryCap,
                None,
                None,
                Some(cap_face),
                None,
                None,
            ));
        }
    };
    let n = outer_indices.len() / 2;
    if outer_indices.len() < 6 || !outer_indices.len().is_multiple_of(2) || sectors.len() != n {
        return Err(ExactBoundarySectorDiagnostic::new(
            ExactBoundarySectorDeclineReason::SectorInventoryMismatch,
            None,
            None,
            Some(cap_face),
            None,
            None,
        ));
    }

    let mut out = Vec::with_capacity(n);
    for sector_index in 0..n {
        let sector_name = sector_name(sector_index);
        let sector = sectors
            .iter()
            .find(|candidate| candidate.sector_index == sector_index)
            .ok_or_else(|| {
                ExactBoundarySectorDiagnostic::new(
                    ExactBoundarySectorDeclineReason::SectorInventoryMismatch,
                    Some(sector_index),
                    Some(sector_name.clone()),
                    Some(cap_face),
                    None,
                    None,
                )
            })?;
        validate_planner_slots(plan, outer_indices, sector, &sector_name, cap_face)?;
        if sector.fillet_end_edge == sector.face_edge {
            return Err(ExactBoundarySectorDiagnostic::new(
                ExactBoundarySectorDeclineReason::DuplicateOuterSeam,
                Some(sector_index),
                Some(sector_name),
                Some(cap_face),
                Some(sector.fillet_end_edge),
                None,
            ));
        }

        let fillet = extract_outer_seam(
            brep,
            cap_face,
            sector,
            &sector_name,
            ExactBoundaryOuterSide::FilletEnd,
            sector.fillet_end_edge,
        )?;
        let bridge = extract_outer_seam(
            brep,
            cap_face,
            sector,
            &sector_name,
            ExactBoundaryOuterSide::FaceEdge,
            sector.face_edge,
        )?;

        out.push(ExactBoundaryCapSectorInput {
            sector_index,
            sector_name,
            fillet_index: sector.fillet_index,
            cap_face,
            fillet_end_edge: sector.fillet_end_edge,
            face_edge: sector.face_edge,
            fillet_support_face: fillet.support_face,
            bridge_support_face: bridge.support_face,
            cap_half_edges: [fillet.cap_half_edge, bridge.cap_half_edge],
            support_half_edges: [fillet.support_half_edge, bridge.support_half_edge],
            outer_seams: [fillet, bridge],
        });
    }
    Ok(out)
}

fn validate_planner_slots(
    plan: &ApexPlan,
    outer_indices: &[usize],
    sector: &GregoryCapSector,
    sector_name: &str,
    cap_face: FaceHandle,
) -> Result<(), ExactBoundarySectorDiagnostic> {
    let fillet_slot = 2 * sector.sector_index;
    let face_slot = fillet_slot + 1;
    let fillet_edge = plan
        .target_topology
        .apex_edges
        .get(outer_indices[fillet_slot]);
    let face_edge = plan
        .target_topology
        .apex_edges
        .get(outer_indices[face_slot]);
    let roles_match = matches!(
        fillet_edge.map(|edge| edge.role),
        Some(TargetEdgeRole::CapVsCyl { fillet_idx }) if fillet_idx == sector.fillet_index
    ) && matches!(
        face_edge.map(|edge| edge.role),
        Some(TargetEdgeRole::CapVsPlane { .. })
    );
    if roles_match {
        Ok(())
    } else {
        Err(ExactBoundarySectorDiagnostic::new(
            ExactBoundarySectorDeclineReason::SectorInventoryMismatch,
            Some(sector.sector_index),
            Some(sector_name.to_string()),
            Some(cap_face),
            None,
            None,
        ))
    }
}

fn extract_outer_seam(
    brep: &BRep,
    cap_face: FaceHandle,
    sector: &GregoryCapSector,
    sector_name: &str,
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
) -> Result<ExactBoundaryOuterSeam, ExactBoundarySectorDiagnostic> {
    let edge_data = brep.edges.get(edge).ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::StaleEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            None,
        )
    })?;
    if edge_data.curve.is_none() {
        return Err(diagnostic(
            ExactBoundarySectorDeclineReason::MissingEdgeCurve,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            None,
        ));
    }
    let cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face).ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingCapHalfEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            None,
        )
    })?;
    let cap_half_edge_data = brep.half_edges.get(cap_half_edge).ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingCapHalfEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            Some(cap_half_edge),
        )
    })?;
    let mut support_half_edge = cap_half_edge_data.twin;
    let mut support_half_edge_data = brep.half_edges.get(support_half_edge).ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingSupportHalfEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            Some(support_half_edge),
        )
    })?;
    let mut support_face = brep
        .loops
        .get(support_half_edge_data.loop_)
        .map(|loop_data| loop_data.face)
        .filter(|face| *face != cap_face)
        .ok_or_else(|| {
            diagnostic(
                ExactBoundarySectorDeclineReason::MissingSupportFace,
                sector,
                sector_name,
                cap_face,
                Some(edge),
                Some(support_half_edge),
            )
        })?;
    if side == ExactBoundaryOuterSide::FilletEnd
        && brep
            .faces
            .get(support_face)
            .and_then(|face| face.surface.as_ref())
            .is_none()
        && let Some((semantic_support_half_edge, semantic_support_face)) =
            same_oriented_cylinder_support_side(brep, cap_face, cap_half_edge)
    {
        support_half_edge = semantic_support_half_edge;
        support_half_edge_data = brep.half_edges.get(support_half_edge).ok_or_else(|| {
            diagnostic(
                ExactBoundarySectorDeclineReason::MissingSupportHalfEdge,
                sector,
                sector_name,
                cap_face,
                Some(edge),
                Some(support_half_edge),
            )
        })?;
        support_face = semantic_support_face;
    }

    Ok(ExactBoundaryOuterSeam {
        side,
        edge,
        cap_half_edge,
        support_half_edge,
        support_face,
        cap_pcurve_present: cap_half_edge_data.pcurve.is_some(),
        support_pcurve_present: support_half_edge_data.pcurve.is_some(),
    })
}

fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
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
    let Some(cap_end) = half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    cap_end == support_end
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(cycle);
        }
    }
    None
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    if let Some(twin) = brep.half_edges.get(data.twin) {
        return Some(twin.origin);
    }
    if !brep.boundary_edges.contains(&data.edge) {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    }
    let edge = brep.edges.get(data.edge)?;
    let Some(curve) = edge.curve.as_ref() else {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    };
    let start_position = brep.vertices.get(data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let destination = if start_position.distance_to(&p0) <= start_position.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&destination);
            (distance <= tolerance).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
}

fn diagnostic(
    reason: ExactBoundarySectorDeclineReason,
    sector: &GregoryCapSector,
    sector_name: &str,
    cap_face: FaceHandle,
    edge: Option<EdgeHandle>,
    half_edge: Option<HalfEdgeHandle>,
) -> ExactBoundarySectorDiagnostic {
    ExactBoundarySectorDiagnostic::new(
        reason,
        Some(sector.sector_index),
        Some(sector_name.to_string()),
        Some(cap_face),
        edge,
        half_edge,
    )
}

fn sector_name(sector_index: usize) -> String {
    format!("direct_apex.cap_sector.{sector_index}")
}
