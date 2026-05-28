//! Structured diagnostics for exact-boundary cap topology commit.

use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, SharedTrimEdgeFailureReason};

use super::sector::ExactBoundaryCapSectorInput;

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
    /// A support-side face has no surface for pcurve validation/commit.
    SupportSurfaceMissing,
    /// The prepared support surface fits the seam but not the whole current
    /// support loop, so a semantic support-face subdivision is required.
    SupportFaceSubdivisionRequired,
    /// Atomic shared trim validation rejected the committed edge/pcurve triple.
    SharedTrimValidationFailed,
    /// Final BRep validation rejected the committed exact-boundary topology.
    PostCommitValidationFailed,
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
            Self::SupportSurfaceMissing => "exact_boundary_cap.commit.support_surface_missing",
            Self::SupportFaceSubdivisionRequired => {
                "exact_boundary_cap.commit.support_face_subdivision_required"
            }
            Self::SharedTrimValidationFailed => {
                "exact_boundary_cap.commit.shared_trim_validation_failed"
            }
            Self::PostCommitValidationFailed => {
                "exact_boundary_cap.commit.post_commit_validation_failed"
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
    /// Cap, sector, or support face being inspected.
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
    pub(crate) fn new(
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

    pub(crate) fn support_side_missing(
        sector: &ExactBoundaryCapSectorInput,
        side_index: usize,
        cap_half_edge: HalfEdgeHandle,
        support_half_edge: HalfEdgeHandle,
        support_face: Option<FaceHandle>,
    ) -> Self {
        Self {
            reason: ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing,
            code: ExactBoundaryCommitDeclineReason::SupportHalfEdgeMissing.code(),
            sector_name: Some(sector.sector_name.clone()),
            face: support_face.or(Some(sector.cap_face)),
            edge: Some(sector.outer_seams[side_index].edge),
            cap_half_edge: Some(cap_half_edge),
            support_half_edge: Some(support_half_edge),
            shared_trim_failure: None,
            residual: None,
        }
    }

    pub(crate) fn support_surface_missing(
        sector: &ExactBoundaryCapSectorInput,
        side_index: usize,
        cap_half_edge: HalfEdgeHandle,
        support_half_edge: HalfEdgeHandle,
        support_face: FaceHandle,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason: ExactBoundaryCommitDeclineReason::SupportSurfaceMissing,
            code: ExactBoundaryCommitDeclineReason::SupportSurfaceMissing.code(),
            sector_name: Some(sector.sector_name.clone()),
            face: Some(support_face),
            edge: Some(sector.outer_seams[side_index].edge),
            cap_half_edge: Some(cap_half_edge),
            support_half_edge: Some(support_half_edge),
            shared_trim_failure: None,
            residual,
        }
    }

    pub(crate) fn support_face_subdivision_required(
        sector: &ExactBoundaryCapSectorInput,
        side_index: usize,
        cap_half_edge: HalfEdgeHandle,
        support_half_edge: HalfEdgeHandle,
        support_face: FaceHandle,
        residual: f64,
    ) -> Self {
        Self {
            reason: ExactBoundaryCommitDeclineReason::SupportFaceSubdivisionRequired,
            code: ExactBoundaryCommitDeclineReason::SupportFaceSubdivisionRequired.code(),
            sector_name: Some(sector.sector_name.clone()),
            face: Some(support_face),
            edge: Some(sector.outer_seams[side_index].edge),
            cap_half_edge: Some(cap_half_edge),
            support_half_edge: Some(support_half_edge),
            shared_trim_failure: None,
            residual: Some(residual),
        }
    }

    pub(crate) fn shared_trim_validation_failed(
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
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{ExactBoundaryCommitDeclineReason, ExactBoundaryCommitDiagnostic};
    use popkern_core::topology::FaceHandle;

    #[test]
    fn commit_diagnostic_promotes_to_apex_diagnostic() {
        let diagnostic = ExactBoundaryCommitDiagnostic {
            reason: ExactBoundaryCommitDeclineReason::SectorTopologyMissing,
            code: "exact_boundary_cap.commit.sector_topology_missing",
            sector_name: Some("sector_0".to_string()),
            face: Some(FaceHandle::default()),
            edge: None,
            cap_half_edge: None,
            support_half_edge: None,
            shared_trim_failure: None,
            residual: None,
        };

        let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(apex) =
            diagnostic.into_apex_diagnostic()
        else {
            panic!("commit diagnostic must promote through the exact-boundary diagnostic channel");
        };
        assert_eq!(apex.stage, "exact_boundary_cap.commit");
        assert_eq!(
            apex.code,
            "exact_boundary_cap.commit.sector_topology_missing"
        );
        assert_eq!(apex.sector_name, "sector_0");
    }
}
