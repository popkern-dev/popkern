use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};
use serde::{Deserialize, Serialize};

use super::ContinuityLevel;
use super::types_repair::ExactBoundaryCandidatePerformanceReport;

const CAP_FAMILY_DEFAULT_G0_TOLERANCE: f64 = 5.0e-3;
const CAP_FAMILY_DEFAULT_G1_NORMAL_DEFICIT_TOLERANCE: f64 = 5.0e-2;
const CAP_FAMILY_DEFAULT_G2_RELATIVE_CURVATURE_TOLERANCE: f64 = 1.0e-2;

/// Semantic cap/support continuity target used by cap-family diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CapContinuityTarget {
    /// Continuity requested by the construction route.
    pub requested: ContinuityLevel,
    /// Lowest continuity level that may be accepted without downgrade.
    pub accepted_floor: ContinuityLevel,
    /// Local positional seam gate.
    pub g0_position_tolerance: f64,
    /// Local tangent-normal parallelism deficit gate.
    pub g1_normal_deficit_tolerance: f64,
    /// Local relative curvature gate for G2 requests.
    pub g2_relative_curvature_tolerance: Option<f64>,
}

impl CapContinuityTarget {
    /// Native cap-family default: require G1 after the cap/support G0 gate.
    pub fn native_cap_g1(g0_position_tolerance: f64, g1_normal_deficit_tolerance: f64) -> Self {
        Self {
            requested: ContinuityLevel::G1,
            accepted_floor: ContinuityLevel::G1,
            g0_position_tolerance,
            g1_normal_deficit_tolerance,
            g2_relative_curvature_tolerance: Some(
                CAP_FAMILY_DEFAULT_G2_RELATIVE_CURVATURE_TOLERANCE,
            ),
        }
    }
}

impl Default for CapContinuityTarget {
    fn default() -> Self {
        Self::native_cap_g1(
            CAP_FAMILY_DEFAULT_G0_TOLERANCE,
            CAP_FAMILY_DEFAULT_G1_NORMAL_DEFICIT_TOLERANCE,
        )
    }
}

/// Measured or pending cap/support continuity evidence for reports.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CapContinuityReport {
    /// Target and local gates used by the route.
    pub target: CapContinuityTarget,
    /// Best continuity level proven by measurement, when measured.
    pub achieved: Option<ContinuityLevel>,
    /// Number of cap/support seams measured.
    pub measured_edges: usize,
    /// Worst sampled positional seam gap.
    pub worst_g0_gap: Option<f64>,
    /// Worst sampled tangent-normal parallelism deficit.
    pub worst_g1_normal_deficit: Option<f64>,
    /// Worst sampled relative curvature mismatch for G2 requests.
    pub worst_g2_relative_curvature_error: Option<f64>,
    /// Stable decline code when the target was not achieved.
    pub decline_code: Option<String>,
}

impl CapContinuityReport {
    /// Report a target that is part of the route contract but has not yet been measured.
    pub fn not_measured(target: CapContinuityTarget) -> Self {
        Self {
            target,
            achieved: None,
            measured_edges: 0,
            worst_g0_gap: None,
            worst_g1_normal_deficit: None,
            worst_g2_relative_curvature_error: None,
            decline_code: None,
        }
    }

    /// Report a native cap/support G1 decline after G0 has been measured.
    pub fn declined_g1(
        target: CapContinuityTarget,
        measured_edges: usize,
        worst_g0_gap: f64,
        worst_g1_normal_deficit: f64,
    ) -> Self {
        Self {
            target,
            achieved: Some(ContinuityLevel::G0),
            measured_edges,
            worst_g0_gap: Some(worst_g0_gap),
            worst_g1_normal_deficit: Some(worst_g1_normal_deficit),
            worst_g2_relative_curvature_error: None,
            decline_code: Some("cap_outer_g1_exceeds_gate".to_string()),
        }
    }

    /// Report a native cap/support decline before positional continuity was met.
    pub fn declined_g0(
        target: CapContinuityTarget,
        measured_edges: usize,
        worst_g0_gap: f64,
        worst_g1_normal_deficit: f64,
    ) -> Self {
        Self {
            target,
            achieved: None,
            measured_edges,
            worst_g0_gap: Some(worst_g0_gap),
            worst_g1_normal_deficit: Some(worst_g1_normal_deficit),
            worst_g2_relative_curvature_error: None,
            decline_code: Some("cap_outer_g0_exceeds_gate".to_string()),
        }
    }
}

impl Default for CapContinuityReport {
    fn default() -> Self {
        Self::not_measured(CapContinuityTarget::default())
    }
}

/// Versioned programmatic and GUI report for an exact-boundary cap candidate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExactBoundaryCapReport {
    /// Schema version for serialized artifacts and downstream UI adapters.
    pub schema_version: u32,
    /// Route that produced or skipped the candidate evidence.
    pub route_label: String,
    /// Candidate outcome in user-facing terms.
    pub candidate_status: ExactBoundaryCandidateOutcome,
    /// Ordered route decisions made while preparing the cap candidate.
    pub decision_chain: Vec<ExactBoundaryRouteDecisionStep>,
    /// Report producer, transaction status, and compatibility metadata.
    pub provenance: ExactBoundaryReportProvenance,
    /// Sector ownership records extracted from live topology.
    pub sectors: Vec<ExactBoundaryReportSector>,
    /// Sector faces committed by an accepted candidate.
    pub committed_faces: Vec<FaceHandle>,
    /// Whether this report has enough source handles for a future preview.
    pub preview_available: bool,
    /// Tolerances and residuals used to explain accept/decline decisions.
    pub tolerance_evidence: ExactBoundaryToleranceEvidence,
    /// Cap/support continuity target and measured status for this route.
    pub continuity: CapContinuityReport,
    /// Candidate-only performance counters and deterministic skip evidence.
    pub performance: ExactBoundaryCandidatePerformanceReport,
    /// Structured diagnostics promoted from candidate stages.
    pub diagnostics: Vec<ExactBoundaryReportDiagnostic>,
    /// Repair fact ids shared by diagnostics in this report.
    pub source_fact_ids: Vec<usize>,
}

impl ExactBoundaryCapReport {
    /// Report for routes where exact-boundary construction was not attempted.
    pub fn not_applicable(route_label: &str) -> Self {
        Self {
            schema_version: 1,
            route_label: route_label.to_string(),
            candidate_status: ExactBoundaryCandidateOutcome::NotApplicable,
            decision_chain: vec![ExactBoundaryRouteDecisionStep::not_attempted(route_label)],
            provenance: ExactBoundaryReportProvenance::not_started(),
            sectors: Vec::new(),
            committed_faces: Vec::new(),
            preview_available: false,
            tolerance_evidence: ExactBoundaryToleranceEvidence::default(),
            continuity: CapContinuityReport::default(),
            performance: ExactBoundaryCandidatePerformanceReport::default(),
            diagnostics: Vec::new(),
            source_fact_ids: Vec::new(),
        }
    }

    /// Return a copy tagged for the final route report.
    pub fn with_route_label(mut self, route_label: &str) -> Self {
        self.route_label = route_label.to_string();
        self
    }

    /// Return a copy with an explicit ordered route decision chain.
    pub fn with_decision_chain(
        mut self,
        decision_chain: Vec<ExactBoundaryRouteDecisionStep>,
    ) -> Self {
        self.decision_chain = decision_chain;
        self
    }
}

/// One ordered route decision in an exact-boundary cap attempt.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactBoundaryRouteDecisionStep {
    /// Stable route label for UI/API filtering.
    pub route_label: String,
    /// Candidate-local outcome for this route.
    pub outcome: ExactBoundaryRouteDecisionOutcome,
    /// First diagnostic code explaining a decline, when any.
    pub diagnostic_code: Option<String>,
}

impl ExactBoundaryRouteDecisionStep {
    /// Build an unattempted route step.
    pub fn not_attempted(route_label: &str) -> Self {
        Self {
            route_label: route_label.to_string(),
            outcome: ExactBoundaryRouteDecisionOutcome::NotAttempted,
            diagnostic_code: None,
        }
    }

    /// Build an accepted route step.
    pub fn accepted(route_label: &str) -> Self {
        Self {
            route_label: route_label.to_string(),
            outcome: ExactBoundaryRouteDecisionOutcome::Accepted,
            diagnostic_code: None,
        }
    }

    /// Build a commit-eligible candidate-only route step.
    pub fn commit_eligible(route_label: &str) -> Self {
        Self {
            route_label: route_label.to_string(),
            outcome: ExactBoundaryRouteDecisionOutcome::CommitEligible,
            diagnostic_code: None,
        }
    }

    /// Build a declined route step with its first deterministic diagnostic.
    pub fn declined(route_label: &str, diagnostic_code: Option<String>) -> Self {
        Self {
            route_label: route_label.to_string(),
            outcome: ExactBoundaryRouteDecisionOutcome::Declined,
            diagnostic_code,
        }
    }
}

/// Candidate-local outcome for an exact-boundary route step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExactBoundaryRouteDecisionOutcome {
    /// Route did not run.
    NotAttempted,
    /// Route passed validation and committed live topology.
    Accepted,
    /// Route passed candidate-only gates and awaits commit preflight.
    CommitEligible,
    /// Route declined with structured diagnostics.
    Declined,
}

/// Versioned producer and transaction metadata for exact-boundary reports.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactBoundaryReportProvenance {
    /// Stable report family name for compatibility checks.
    pub schema_name: String,
    /// Stable producer stage that authored this report.
    pub producer: String,
    /// Transaction outcome for the exact-boundary candidate.
    pub transaction_status: ExactBoundaryReportTransactionStatus,
    /// Whether candidate topology mutation was contained by snapshot rollback.
    pub topology_rollback_performed: bool,
}

impl ExactBoundaryReportProvenance {
    /// Metadata for routes where no exact-boundary transaction was started.
    pub fn not_started() -> Self {
        Self::for_status(ExactBoundaryCandidateOutcome::NotApplicable)
    }

    /// Metadata for a completed exact-boundary candidate outcome.
    pub fn for_status(status: ExactBoundaryCandidateOutcome) -> Self {
        let transaction_status = match status {
            ExactBoundaryCandidateOutcome::NotApplicable => {
                ExactBoundaryReportTransactionStatus::NotStarted
            }
            ExactBoundaryCandidateOutcome::Accepted => {
                ExactBoundaryReportTransactionStatus::Committed
            }
            ExactBoundaryCandidateOutcome::Declined => {
                ExactBoundaryReportTransactionStatus::RolledBack
            }
        };
        Self {
            schema_name: "exact_boundary_cap_report".to_string(),
            producer: "exact_boundary_cap_candidate".to_string(),
            transaction_status,
            topology_rollback_performed: transaction_status
                == ExactBoundaryReportTransactionStatus::RolledBack,
        }
    }
}

/// Transaction outcome for exact-boundary cap construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExactBoundaryReportTransactionStatus {
    /// The route did not start an exact-boundary candidate transaction.
    NotStarted,
    /// Candidate postconditions passed and topology was committed.
    Committed,
    /// Candidate declined and restored the input topology snapshot.
    RolledBack,
}

impl Default for ExactBoundaryCapReport {
    fn default() -> Self {
        Self::not_applicable("not_applicable")
    }
}

/// Exact-boundary candidate outcome for UX and JSON consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExactBoundaryCandidateOutcome {
    /// The route did not attempt an exact-boundary cap candidate.
    NotApplicable,
    /// Candidate passed validation and committed exact-boundary sector faces.
    Accepted,
    /// Candidate declined and restored the input topology.
    Declined,
}

/// Source ownership for one exact-boundary cap sector.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactBoundaryReportSector {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable semantic name for diagnostics and UI highlighting.
    pub sector_name: String,
    /// Fillet index owned by this sector.
    pub fillet_index: usize,
    /// Cap face that owns the cap-side half-edges.
    pub cap_face: FaceHandle,
    /// Fillet-end outer edge.
    pub fillet_end_edge: EdgeHandle,
    /// Plane/support bridge outer edge.
    pub face_edge: EdgeHandle,
    /// Support face across the fillet-end seam.
    pub fillet_support_face: FaceHandle,
    /// Support face across the bridge seam.
    pub bridge_support_face: FaceHandle,
    /// Cap-side half-edges, `[FilletEnd, FaceEdge]`.
    pub cap_half_edges: [HalfEdgeHandle; 2],
    /// Support-side half-edges, `[FilletEnd, FaceEdge]`.
    pub support_half_edges: [HalfEdgeHandle; 2],
    /// Cap-side pcurve availability, `[FilletEnd, FaceEdge]`.
    pub cap_pcurves_present: [bool; 2],
    /// Support-side pcurve availability, `[FilletEnd, FaceEdge]`.
    pub support_pcurves_present: [bool; 2],
}

/// Tolerance and residual evidence explaining candidate accept/decline.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExactBoundaryToleranceEvidence {
    /// BRep model linear tolerance active when the candidate ran.
    pub model_linear_tolerance: Option<f64>,
    /// Candidate-local inner-spoke G0 gate.
    pub inner_spoke_g0_tolerance: f64,
    /// Candidate-local minimum generated spoke length.
    pub inner_spoke_min_length: f64,
    /// Largest residual reported by candidate diagnostics, when any.
    pub max_reported_residual: Option<f64>,
    /// Deterministic decline code used by future skip policies.
    pub deterministic_skip_code: Option<String>,
}

impl Default for ExactBoundaryToleranceEvidence {
    fn default() -> Self {
        Self {
            model_linear_tolerance: None,
            inner_spoke_g0_tolerance: 1.0e-7,
            inner_spoke_min_length: 1.0e-6,
            max_reported_residual: None,
            deterministic_skip_code: None,
        }
    }
}

/// JSON-ready diagnostic row for exact-boundary reports.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExactBoundaryReportDiagnostic {
    /// Diagnostic row schema version for downstream compatibility.
    pub schema_version: u32,
    /// Candidate stage that emitted the diagnostic.
    pub stage: String,
    /// Stable diagnostic code.
    pub code: String,
    /// Optional detailed sub-code for the same diagnostic family.
    pub detail_code: Option<String>,
    /// Stable sector name where available.
    pub sector_name: String,
    /// Source edge being prepared.
    pub edge: Option<EdgeHandle>,
    /// Support half-edge being sampled, when applicable.
    pub support_half_edge: Option<HalfEdgeHandle>,
    /// Support face being sampled, when applicable.
    pub support_face: Option<FaceHandle>,
    /// Measured residual, when available.
    pub residual: Option<f64>,
    /// Sector-local outer side for seam-scoped diagnostics.
    pub outer_side: Option<String>,
    /// Normalized seam sample parameter for sample-scoped diagnostics.
    pub sample_parameter: Option<f64>,
    /// Local tolerance gate used by the diagnostic, when available.
    pub tolerance: Option<f64>,
    /// Shared BRep repair knowledge facts that explain the source evidence.
    pub shared_fact_ids: Vec<usize>,
}

#[cfg(test)]
mod tests {
    use super::{
        CapContinuityReport, CapContinuityTarget, ExactBoundaryCandidateOutcome,
        ExactBoundaryCapReport,
    };
    use crate::types::ContinuityLevel;

    #[test]
    fn not_applicable_report_is_versioned_and_json_ready() {
        let report = ExactBoundaryCapReport::not_applicable("legacy");

        assert_eq!(report.schema_version, 1);
        assert_eq!(report.route_label, "legacy");
        assert_eq!(
            report.candidate_status,
            ExactBoundaryCandidateOutcome::NotApplicable
        );
        assert!(!report.preview_available);
        assert!(report.sectors.is_empty());
        assert_eq!(
            report.provenance.transaction_status,
            super::ExactBoundaryReportTransactionStatus::NotStarted
        );
        assert!(!report.provenance.topology_rollback_performed);
        assert_eq!(report.continuity.target.requested, ContinuityLevel::G1);
        assert_eq!(report.continuity.achieved, None);

        let json = serde_json::to_string(&report).unwrap();
        assert!(json.contains("\"candidate_status\":\"NotApplicable\""));
        assert!(json.contains("\"decision_chain\""));
        assert!(json.contains("\"outcome\":\"NotAttempted\""));
        assert!(json.contains("\"schema_name\":\"exact_boundary_cap_report\""));
        assert!(json.contains("\"requested\":\"G1\""));
    }

    #[test]
    fn cap_continuity_decline_report_records_target_and_residuals() {
        let target = CapContinuityTarget::native_cap_g1(5.0e-3, 5.0e-2);
        let report = CapContinuityReport::declined_g1(target, 3, 1.0e-4, 6.0e-2);

        assert_eq!(report.target.requested, ContinuityLevel::G1);
        assert_eq!(report.target.accepted_floor, ContinuityLevel::G1);
        assert_eq!(report.achieved, Some(ContinuityLevel::G0));
        assert_eq!(report.measured_edges, 3);
        assert_eq!(report.worst_g0_gap, Some(1.0e-4));
        assert_eq!(report.worst_g1_normal_deficit, Some(6.0e-2));
        assert_eq!(
            report.decline_code.as_deref(),
            Some("cap_outer_g1_exceeds_gate")
        );
    }
}
