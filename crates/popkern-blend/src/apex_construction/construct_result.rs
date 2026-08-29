//! Build result type for apex construction.

use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

use super::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;
use super::telemetry::PlanTelemetry;
use crate::types::{ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapReport};

#[derive(Debug, Clone)]
/// Topology and telemetry produced by an apex construction attempt.
pub struct ApexBuildResult {
    /// New cap faces created by the builder.
    pub created_faces: Vec<FaceHandle>,
    /// Realized apex vertices, indexed like the plan target vertices.
    pub apex_vertices: Vec<VertexHandle>,
    /// Realized apex edges, indexed like the plan target edges.
    pub apex_edges: Vec<Option<EdgeHandle>>,
    /// Corner-patch reports created by cap surface assignment.
    pub corner_patches: Vec<crate::types::CornerPatchInfo>,
    /// Exact-boundary sector ownership facts extracted from live cap topology.
    pub(crate) exact_boundary_sectors: Vec<ExactBoundaryCapSectorInput>,
    /// Number of exact-boundary curves prepared without approximation.
    pub(crate) exact_boundary_prepared_curve_count: usize,
    /// Number of generated inner spokes prepared without topology mutation.
    pub(crate) exact_boundary_inner_spoke_count: usize,
    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub(crate) exact_boundary_ribbon_count: usize,
    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub(crate) exact_boundary_solved_surface_count: usize,
    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub(crate) exact_boundary_pcurve_pair_count: usize,
    /// Candidate-only performance counters and budget evidence for exact-boundary caps.
    pub(crate) exact_boundary_performance: ExactBoundaryCandidatePerformanceReport,
    /// Versioned exact-boundary report for GUI and programmatic consumers.
    pub(crate) exact_boundary_report: ExactBoundaryCapReport,
    /// Structured diagnostics emitted while planning or declining apex repairs.
    pub diagnostics: Vec<crate::types::ApexConstructionDiagnostic>,
    /// Planner telemetry carried through to the build result.
    pub telemetry: PlanTelemetry,
}
