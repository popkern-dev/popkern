use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};
use serde::{Deserialize, Serialize};

use super::types_repair::ExactBoundaryCandidatePerformanceReport;

/// Versioned programmatic and GUI report for an exact-boundary cap candidate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExactBoundaryCapReport {
    /// Schema version for serialized artifacts and downstream UI adapters.
    pub schema_version: u32,
    /// Route that produced or skipped the candidate evidence.
    pub route_label: String,
    /// Candidate outcome in user-facing terms.
    pub candidate_status: ExactBoundaryCandidateOutcome,
    /// Sector ownership records extracted from live topology.
    pub sectors: Vec<ExactBoundaryReportSector>,
    /// Sector faces committed by an accepted candidate.
    pub committed_faces: Vec<FaceHandle>,
    /// Whether this report has enough source handles for a future preview.
    pub preview_available: bool,
    /// Tolerances and residuals used to explain accept/decline decisions.
    pub tolerance_evidence: ExactBoundaryToleranceEvidence,
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
            sectors: Vec::new(),
            committed_faces: Vec::new(),
            preview_available: false,
            tolerance_evidence: ExactBoundaryToleranceEvidence::default(),
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
    /// Candidate stage that emitted the diagnostic.
    pub stage: String,
    /// Stable diagnostic code.
    pub code: String,
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
    /// Shared BRep repair knowledge facts that explain the source evidence.
    pub shared_fact_ids: Vec<usize>,
}

#[cfg(test)]
mod tests {
    use super::{ExactBoundaryCandidateOutcome, ExactBoundaryCapReport};

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

        let json = serde_json::to_string(&report).unwrap();
        assert!(json.contains("\"candidate_status\":\"NotApplicable\""));
    }
}
