//! Build result type for apex construction.

use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

use super::telemetry::PlanTelemetry;

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
    /// Planner telemetry carried through to the build result.
    pub telemetry: PlanTelemetry,
}
