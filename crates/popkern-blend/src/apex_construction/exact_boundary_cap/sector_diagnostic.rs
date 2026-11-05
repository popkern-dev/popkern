//! Structured declines of exact-boundary cap sector extraction: stable
//! reason codes and the diagnostic carrying the handles inspected.

use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::apex_construction::error::ApexBuildError;

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
    pub(super) fn new(
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
