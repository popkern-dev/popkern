use popkern_core::topology::{BrepRepairFactId, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Exact-boundary cap preparation diagnostic for candidate construction.
#[derive(Debug, Clone, PartialEq)]
pub struct ExactBoundaryCapDiagnostic {
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
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}
