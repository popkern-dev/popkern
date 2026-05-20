use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::topology::{
    BrepRepairActionArtifact, BrepRepairActionRecord, BrepRepairActionSummary, EdgeHandle,
    FaceHandle, HalfEdgeHandle, VertexHandle,
};

/// Structured diagnostics for repair and hardening actions performed during a
/// blend operation.
#[derive(Debug, Clone, Default)]
pub struct BlendRepairReport {
    /// Whole-blend support-face seam relocation diagnostics.
    pub support_seam_relocation: SupportSeamRelocationReport,
    /// Per-apex construction route decisions for legacy/direct-apex migration tracking.
    pub apex_routes: Vec<ApexConstructionRouteDecision>,
    /// Structured apex construction diagnostics and declined repair candidates.
    pub apex_diagnostics: Vec<ApexConstructionDiagnostic>,
}

/// Construction route used for one smooth multi-fillet apex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApexConstructionRoute {
    /// direct-apex was not requested; legacy `insert_vertex_blend` handled the apex.
    Legacy,
    /// direct-apex planner and builder completed and legacy was skipped.
    DirectApex,
    /// direct-apex was requested but declined; legacy handled the apex.
    DirectApexFallbackToLegacy,
}

/// Per-apex legacy/direct-apex migration diagnostic.
#[derive(Debug, Clone)]
pub struct ApexConstructionRouteDecision {
    /// Original apex vertex being resolved.
    pub vertex: VertexHandle,
    /// Original fillet edges incident to this apex.
    pub incident_edges: Vec<EdgeHandle>,
    /// Route that produced the final topology for this apex.
    pub route: ApexConstructionRoute,
    /// Faces authored by direct-apex, empty for legacy and fallback routes.
    pub built_faces: Vec<FaceHandle>,
    /// Number of corner patch reports produced by this route.
    pub corner_patch_count: usize,
    /// direct-apex decline reason when [`ApexConstructionRoute::DirectApexFallbackToLegacy`] applies.
    pub fallback_reason: Option<String>,
}

/// Structured apex-construction diagnostic emitted by construction planning and
/// repair preflights.
#[derive(Debug, Clone, PartialEq)]
pub enum ApexConstructionDiagnostic {
    /// Target-cylinder trim replacement could not materialize an exact,
    /// ownership-balanced far-side trim.
    CylinderTrimMaterialization(CylinderTrimMaterializationDiagnostic),
    /// Native Gregory cap/support seam continuity did not meet the commercial
    /// G1 gate after construction.
    CapSupportContinuity(CapSupportContinuityDiagnostic),
}

/// Native cap/support continuity diagnostic for a constructed apex cap.
#[derive(Debug, Clone, PartialEq)]
pub struct CapSupportContinuityDiagnostic {
    /// Cap face whose support seam was measured.
    pub cap_face: FaceHandle,
    /// Adjacent support or fillet face across the worst seam.
    pub support_face: FaceHandle,
    /// Topological edge shared by the cap and support face.
    pub edge: EdgeHandle,
    /// Worst sampled unit-normal parallelism deficit.
    pub worst_g1_coplanarity: f64,
    /// Commercial gate used for the measurement.
    pub gate: f64,
    /// Number of cap/support edges measured for the apex.
    pub edges_measured: usize,
}

/// Materialized target-cylinder trim replacement diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub struct CylinderTrimMaterializationDiagnostic {
    /// Fillet index whose target cylinder trim was inspected.
    pub fillet_index: usize,
    /// Target cylinder face being rewritten.
    pub cylinder_face: FaceHandle,
    /// Start vertex of the declined far-side trim.
    pub start: VertexHandle,
    /// End vertex of the declined far-side trim.
    pub end: VertexHandle,
    /// Selected far-side trim strategy.
    pub strategy: CylinderTrimFarSideStrategy,
    /// Support-plane candidates containing the start endpoint.
    pub start_support_faces: Vec<FaceHandle>,
    /// Support-plane candidates containing the end endpoint.
    pub end_support_faces: Vec<FaceHandle>,
    /// Stable decline reason.
    pub reason: CylinderTrimMaterializationDeclineReason,
    /// Whether exact curve authoring was available for the selected span.
    pub exact_authorable: bool,
}

/// Far-side trim direction selected for target-cylinder loop replacement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderTrimFarSideStrategy {
    /// Forward loop path between target-cylinder anchors.
    Forward,
    /// Reverse loop path between target-cylinder anchors.
    Reverse,
}

/// Stable reason a materialized cylinder trim replacement declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderTrimMaterializationDeclineReason {
    /// Endpoints belong to different support planes and no
    /// plane-plane-cylinder transition point exists.
    NoPlanePlaneCylinderTransition,
    /// The preflight declined for a different currently unclassified reason.
    Other,
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
    let mut apex_declined = 0usize;
    let mut apex_max_deviation = 0.0f64;
    for diagnostic in &report.apex_diagnostics {
        match diagnostic {
            ApexConstructionDiagnostic::CylinderTrimMaterialization(_) => {}
            ApexConstructionDiagnostic::CapSupportContinuity(diagnostic) => {
                apex_declined += 1;
                apex_max_deviation = apex_max_deviation.max(diagnostic.worst_g1_coplanarity);
                actions.push(BrepRepairActionRecord {
                    stage: "apex_cap_support_continuity".to_string(),
                    outcome: "declined".to_string(),
                    subject: format!(
                        "cap={:?};support={:?};edge={:?}",
                        diagnostic.cap_face, diagnostic.support_face, diagnostic.edge
                    ),
                    reason: Some("cap_outer_g1_exceeds_gate".to_string()),
                    max_motion: 0.0,
                    max_deviation: diagnostic.worst_g1_coplanarity,
                    affected_handles: vec![
                        format!("face:{:?}", diagnostic.cap_face),
                        format!("face:{:?}", diagnostic.support_face),
                        format!("edge:{:?}", diagnostic.edge),
                    ],
                    shared_fact_ids: Vec::new(),
                });
            }
        }
    }
    BrepRepairActionArtifact {
        summary: BrepRepairActionSummary {
            attempted_actions: seam.candidates + rolled_back + apex_declined,
            authored_actions: seam.authored,
            declined_actions: seam.declined + apex_declined,
            rolled_back_actions: rolled_back,
            max_motion: seam.max_endpoint_motion,
            max_deviation: seam.max_sampled_disagreement.max(apex_max_deviation),
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
