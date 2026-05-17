use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::topology::{
    BrepRepairActionArtifact, BrepRepairActionRecord, BrepRepairActionSummary, FaceHandle,
    HalfEdgeHandle,
};

/// Structured diagnostics for repair and hardening actions performed during a
/// blend operation.
#[derive(Debug, Clone, Default)]
pub struct BlendRepairReport {
    /// Whole-blend support-face seam relocation diagnostics.
    pub support_seam_relocation: SupportSeamRelocationReport,
}

/// Build a serialization-ready repair-action artifact from blend repair diagnostics.
pub fn build_blend_repair_action_artifact(report: &BlendRepairReport) -> BrepRepairActionArtifact {
    let seam = &report.support_seam_relocation;
    let rolled_back = usize::from(seam.rollback_reason.is_some());
    let mut actions: Vec<_> = seam
        .decisions
        .iter()
        .map(|decision| BrepRepairActionRecord {
            stage: "support_seam_relocation".to_string(),
            outcome: if decision.accepted {
                "authored".to_string()
            } else {
                "declined".to_string()
            },
            subject: format!(
                "cap={:?};half_edge={:?}",
                decision.cap_face, decision.half_edge
            ),
            reason: decision
                .reason
                .map(support_seam_decline_reason_label)
                .map(str::to_string),
            max_motion: decision.max_endpoint_motion,
            max_deviation: decision.max_sampled_disagreement,
            affected_handles: Vec::new(),
            shared_fact_ids: Vec::new(),
        })
        .collect();
    if let Some(reason) = seam.rollback_reason {
        actions.push(BrepRepairActionRecord {
            stage: "support_seam_relocation".to_string(),
            outcome: "rolled_back".to_string(),
            subject: "support_seam_relocation_transaction".to_string(),
            reason: Some(support_seam_rollback_reason_label(reason).to_string()),
            max_motion: seam.max_endpoint_motion,
            max_deviation: seam.max_sampled_disagreement,
            affected_handles: Vec::new(),
            shared_fact_ids: Vec::new(),
        });
    }
    BrepRepairActionArtifact {
        summary: BrepRepairActionSummary {
            attempted_actions: seam.candidates + rolled_back,
            authored_actions: seam.authored,
            declined_actions: seam.declined,
            rolled_back_actions: rolled_back,
            max_motion: seam.max_endpoint_motion,
            max_deviation: seam.max_sampled_disagreement,
        },
        actions,
    }
}

/// Aggregate support-face seam relocation diagnostics.
#[derive(Debug, Clone, Default)]
pub struct SupportSeamRelocationReport {
    /// Number of support seams authored and committed.
    pub authored: usize,
    /// Number of candidate attempts declined before commit.
    pub declined: usize,
    /// Number of candidate attempts inspected.
    pub candidates: usize,
    /// Maximum endpoint movement applied or attempted.
    pub max_endpoint_motion: f64,
    /// Maximum sampled curve/cap/support disagreement observed.
    pub max_sampled_disagreement: f64,
    /// Transaction rollback reason, when the cap-neighborhood job rolled back.
    pub rollback_reason: Option<SupportSeamRollbackReason>,
    /// Per-candidate decisions retained for API and GUI diagnostics.
    pub decisions: Vec<SupportSeamDecision>,
}

impl SupportSeamRelocationReport {
    pub(crate) fn merge(&mut self, other: SupportSeamRelocationReport) {
        self.authored += other.authored;
        self.declined += other.declined;
        self.candidates += other.candidates;
        self.max_endpoint_motion = self.max_endpoint_motion.max(other.max_endpoint_motion);
        self.max_sampled_disagreement = self
            .max_sampled_disagreement
            .max(other.max_sampled_disagreement);
        self.rollback_reason = self.rollback_reason.or(other.rollback_reason);
        self.decisions.extend(other.decisions);
    }

    pub(crate) fn push(&mut self, decision: SupportSeamDecision) {
        self.candidates += 1;
        if decision.accepted {
            self.authored += 1;
        } else {
            self.declined += 1;
        }
        self.max_endpoint_motion = self.max_endpoint_motion.max(decision.max_endpoint_motion);
        self.max_sampled_disagreement = self
            .max_sampled_disagreement
            .max(decision.max_sampled_disagreement);
        self.decisions.push(decision);
    }
}

/// One support-face seam relocation decision.
#[derive(Debug, Clone)]
pub struct SupportSeamDecision {
    /// Candidate half-edge on the cap face.
    pub half_edge: HalfEdgeHandle,
    /// Gregory cap face being repaired.
    pub cap_face: FaceHandle,
    /// Adjacent support face, when the candidate reached support classification.
    pub support_face: Option<FaceHandle>,
    /// Canonical cap boundary side attempted.
    pub boundary_side: Option<BoundarySide>,
    /// Whether this candidate committed.
    pub accepted: bool,
    /// Decline reason for rejected candidates.
    pub reason: Option<SupportSeamDeclineReason>,
    /// Maximum endpoint movement applied or attempted for this candidate.
    pub max_endpoint_motion: f64,
    /// Maximum sampled curve/cap/support disagreement for this candidate.
    pub max_sampled_disagreement: f64,
}

/// Reason a support-face seam relocation candidate declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportSeamDeclineReason {
    /// Candidate was not adjacent to a planar support face.
    NotPlanarSupportCandidate,
    /// Candidate could not be matched to a canonical Gregory boundary.
    NoCanonicalBoundary,
    /// Canonical boundary endpoints were degenerate.
    DegenerateBoundary,
    /// Endpoint relocation exceeded the local motion gate.
    ExcessiveEndpointMotion,
    /// Authored cap boundary did not match relocated edge vertices.
    CapBoundaryMismatch,
    /// Adjacent support face could not be lifted to contain the seam.
    SupportFaceLiftFailed,
    /// Sampled curve/cap/support agreement failed.
    SampledValidationFailed,
    /// Moved vertices were incoherent with incident edge curves.
    MovedVertexCurveIncoherent,
    /// Candidate would remove the local cap-cylinder neighbor evidence.
    CapCylinderNeighborLost,
    /// Candidate would reduce global Gregory cap-cylinder neighbor evidence.
    GlobalCapCylinderNeighborLost,
}

/// Stable serialized label for a support-seam relocation decline.
pub fn support_seam_decline_reason_label(reason: SupportSeamDeclineReason) -> &'static str {
    match reason {
        SupportSeamDeclineReason::NotPlanarSupportCandidate => "not_planar_support_candidate",
        SupportSeamDeclineReason::NoCanonicalBoundary => "no_canonical_boundary",
        SupportSeamDeclineReason::DegenerateBoundary => "degenerate_boundary",
        SupportSeamDeclineReason::ExcessiveEndpointMotion => "excessive_endpoint_motion",
        SupportSeamDeclineReason::CapBoundaryMismatch => "cap_boundary_mismatch",
        SupportSeamDeclineReason::SupportFaceLiftFailed => "support_face_lift_failed",
        SupportSeamDeclineReason::SampledValidationFailed => "sampled_validation_failed",
        SupportSeamDeclineReason::MovedVertexCurveIncoherent => "moved_vertex_curve_incoherent",
        SupportSeamDeclineReason::CapCylinderNeighborLost => "cap_cylinder_neighbor_lost",
        SupportSeamDeclineReason::GlobalCapCylinderNeighborLost => {
            "global_cap_cylinder_neighbor_lost"
        }
    }
}

/// Reason a support-face seam relocation job rolled back after candidate work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportSeamRollbackReason {
    /// Resulting topology failed loop closure validation.
    LoopClosureFailed,
    /// Existing Gregory cap-cylinder boundary hints would have been lost.
    CapCylinderHintsLost,
    /// Existing Gregory cap-cylinder neighbor evidence would have been lost.
    CapCylinderNeighborsLost,
}

/// Stable serialized label for a support-seam relocation rollback.
pub fn support_seam_rollback_reason_label(reason: SupportSeamRollbackReason) -> &'static str {
    match reason {
        SupportSeamRollbackReason::LoopClosureFailed => "loop_closure_failed",
        SupportSeamRollbackReason::CapCylinderHintsLost => "cap_cylinder_hints_lost",
        SupportSeamRollbackReason::CapCylinderNeighborsLost => "cap_cylinder_neighbors_lost",
    }
}
