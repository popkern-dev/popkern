//! Exact-boundary sector extraction for direct-apex Gregory caps.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

#[path = "sector_support_side.rs"]
mod sector_support_side;

use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::gregory_cap_boundary::half_edge_for_face_on_edge;
use crate::apex_construction::plan::ApexPlan;
use crate::apex_construction::target_topology::{TargetCapFace, TargetEdgeRole};
use sector_support_side::{
    coincident_support_side, same_oriented_cylinder_support_side, same_oriented_support_side,
};

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

/// Exact-boundary sector input from already-authored topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryCapSectorInput {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable semantic name for diagnostics and UI highlighting.
    pub sector_name: String,
    /// Fillet index owned by this sector.
    pub fillet_index: usize,
    /// Current cap face that owns the cap-side half-edges.
    pub cap_face: FaceHandle,
    /// Fillet-end edge for this sector.
    pub fillet_end_edge: EdgeHandle,
    /// Plane/support bridge edge; `None` when the trailing corner is
    /// a welded shared corner (PK-219) — the sector's outer boundary
    /// is the fillet seam alone.
    pub face_edge: Option<EdgeHandle>,
    /// Support face across [`Self::fillet_end_edge`].
    pub fillet_support_face: FaceHandle,
    /// Support face across [`Self::face_edge`] (bridged only).
    pub bridge_support_face: Option<FaceHandle>,
    /// Fillet-side outer seam ownership facts.
    pub fillet_seam: ExactBoundaryOuterSeam,
    /// Bridge-side outer seam ownership facts (bridged sectors only).
    pub bridge_seam: Option<ExactBoundaryOuterSeam>,
}

impl ExactBoundaryCapSectorInput {
    /// Outer seams, `[FilletEnd, FaceEdge?]` (1 welded / 2 bridged).
    pub fn outer_seams(&self) -> impl Iterator<Item = &ExactBoundaryOuterSeam> {
        std::iter::once(&self.fillet_seam).chain(self.bridge_seam.as_ref())
    }

    /// Number of outer seams (1 welded, 2 bridged).
    #[allow(dead_code)]
    pub fn outer_seam_count(&self) -> usize {
        1 + usize::from(self.bridge_seam.is_some())
    }

    /// True when the trailing corner is a welded shared corner.
    #[allow(dead_code)]
    pub fn is_welded(&self) -> bool {
        self.bridge_seam.is_none()
    }

    /// Seam by side index (0 = FilletEnd, 1 = FaceEdge / welded None).
    pub fn seam_at(&self, side_index: usize) -> Option<&ExactBoundaryOuterSeam> {
        match side_index {
            0 => Some(&self.fillet_seam),
            1 => self.bridge_seam.as_ref(),
            _ => None,
        }
    }

    /// Seam for a side; FaceEdge side of a welded sector is `None`.
    pub fn seam_for_side(&self, side: ExactBoundaryOuterSide) -> Option<&ExactBoundaryOuterSeam> {
        match side {
            ExactBoundaryOuterSide::FilletEnd => Some(&self.fillet_seam),
            ExactBoundaryOuterSide::FaceEdge => self.bridge_seam.as_ref(),
        }
    }
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
    // With the shared-corner contract the outer loop holds one edge
    // per fillet seam plus one per BRIDGED corner: 2N - #welded.
    let n = plan.cap_kind.n_fillets();
    let expected_outer = 2 * n - plan.cap_kind.welded_corner_count();
    if n < 3 || outer_indices.len() != expected_outer || sectors.len() != n {
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
        validate_planner_slots(plan, sector, &sector_name, cap_face)?;
        if sector.face_edge == Some(sector.fillet_end_edge) {
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
        let bridge = sector
            .face_edge
            .map(|face_edge| {
                extract_outer_seam(
                    brep,
                    cap_face,
                    sector,
                    &sector_name,
                    ExactBoundaryOuterSide::FaceEdge,
                    face_edge,
                )
            })
            .transpose()?;

        out.push(ExactBoundaryCapSectorInput {
            sector_index,
            sector_name,
            fillet_index: sector.fillet_index,
            cap_face,
            fillet_end_edge: sector.fillet_end_edge,
            face_edge: sector.face_edge,
            fillet_support_face: fillet.support_face,
            bridge_support_face: bridge.as_ref().map(|seam| seam.support_face),
            fillet_seam: fillet,
            bridge_seam: bridge,
        });
    }
    Ok(out)
}

fn validate_planner_slots(
    plan: &ApexPlan,
    sector: &GregoryCapSector,
    sector_name: &str,
    cap_face: FaceHandle,
) -> Result<(), ExactBoundarySectorDiagnostic> {
    // Resolve the plan edges by ROLE (welded plans are not strictly
    // alternating, so raw 2k/2k+1 slot math is invalid).
    let fillet_edge = plan
        .target_topology
        .fillet_end_edge_index(sector.sector_index)
        .and_then(|index| plan.target_topology.apex_edges.get(index));
    let planned_face_edge = plan
        .target_topology
        .face_edge_edge_index(sector.sector_index)
        .and_then(|index| plan.target_topology.apex_edges.get(index));
    let fillet_role_matches = matches!(
        fillet_edge.map(|edge| edge.role),
        Some(TargetEdgeRole::CapVsCyl { fillet_idx }) if fillet_idx == sector.fillet_index
    );
    // A bridged sector must have a planned CapVsPlane edge; a welded
    // sector must have neither a planned bridge nor a live one.
    let bridge_matches = match sector.face_edge {
        Some(_) => matches!(
            planned_face_edge.map(|edge| edge.role),
            Some(TargetEdgeRole::CapVsPlane { .. })
        ),
        None => planned_face_edge.is_none() && plan.cap_kind.corner_is_welded(sector.sector_index),
    };
    let roles_match = fillet_role_matches && bridge_matches;
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
    let mut cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face).ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingCapHalfEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            None,
        )
    })?;
    let mut cap_half_edge_data = brep.half_edges.get(cap_half_edge).ok_or_else(|| {
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
    let mut support_half_edge_data = brep.half_edges.get(support_half_edge);
    let mut support_face = support_half_edge_data
        .and_then(|data| brep.loops.get(data.loop_))
        .map(|loop_data| loop_data.face)
        .filter(|face| *face != cap_face);
    if support_half_edge_data
        .and_then(|data| brep.loops.get(data.loop_))
        .is_some_and(|loop_data| loop_data.face == cap_face)
    {
        cap_half_edge = support_half_edge;
        cap_half_edge_data = brep.half_edges.get(cap_half_edge).ok_or_else(|| {
            diagnostic(
                ExactBoundarySectorDeclineReason::MissingCapHalfEdge,
                sector,
                sector_name,
                cap_face,
                Some(edge),
                Some(cap_half_edge),
            )
        })?;
        support_half_edge = cap_half_edge_data.twin;
        support_half_edge_data = brep.half_edges.get(support_half_edge);
        support_face = support_half_edge_data
            .and_then(|data| brep.loops.get(data.loop_))
            .map(|loop_data| loop_data.face)
            .filter(|face| *face != cap_face);
    }
    if side == ExactBoundaryOuterSide::FilletEnd
        && (support_half_edge_data.is_none()
            || support_face.is_none_or(|face| {
                brep.faces
                    .get(face)
                    .and_then(|face_data| face_data.surface.as_ref())
                    .is_none()
            }))
        && let Some((semantic_support_half_edge, semantic_support_face)) =
            same_oriented_cylinder_support_side(brep, cap_face, cap_half_edge)
    {
        support_half_edge = semantic_support_half_edge;
        support_half_edge_data = brep.half_edges.get(support_half_edge);
        support_face = Some(semantic_support_face);
    }
    if support_face.is_none_or(|face| face == cap_face)
        && let Some((semantic_support_half_edge, semantic_support_face)) =
            same_oriented_support_side(brep, cap_face, cap_half_edge)
    {
        support_half_edge = semantic_support_half_edge;
        support_half_edge_data = brep.half_edges.get(support_half_edge);
        support_face = Some(semantic_support_face);
    }
    if support_face.is_none_or(|face| face == cap_face)
        && let Some((semantic_support_half_edge, semantic_support_face)) =
            coincident_support_side(brep, cap_face, cap_half_edge)
    {
        support_half_edge = semantic_support_half_edge;
        support_half_edge_data = brep.half_edges.get(support_half_edge);
        support_face = Some(semantic_support_face);
    }
    let support_half_edge_data = support_half_edge_data.ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingSupportHalfEdge,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            Some(support_half_edge),
        )
    })?;
    let support_face = support_face.ok_or_else(|| {
        diagnostic(
            ExactBoundarySectorDeclineReason::MissingSupportFace,
            sector,
            sector_name,
            cap_face,
            Some(edge),
            Some(support_half_edge),
        )
    })?;

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
