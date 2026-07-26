//! Candidate status and result types for the exact-boundary cap
//! facade — the record every route (one-patch, multipatch, PnS
//! substitution) returns through, split from `facade.rs`.

use popkern_core::topology::FaceHandle;

use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;
use crate::types::{
    ApexConstructionDiagnostic, ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapReport,
};

/// Status of one exact-boundary cap candidate attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryCapCandidateStatus {
    /// Candidate passed validation and committed exact-boundary sector faces.
    Accepted,
    /// Candidate declined and restored the input BRep snapshot.
    Declined,
}

/// Reusable exact-boundary candidate result for construction, diagnostics, and UX.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryCapCandidateResult {
    /// Candidate status after validation/commit or rollback.
    pub status: ExactBoundaryCapCandidateStatus,
    /// Sector ownership records. On accept these reference committed sector faces.
    pub sectors: Vec<ExactBoundaryCapSectorInput>,
    /// Faces committed by the accepted candidate.
    pub committed_faces: Vec<FaceHandle>,
    /// Candidate-only performance counters and budget evidence.
    pub performance: ExactBoundaryCandidatePerformanceReport,
    /// Versioned report for UI and programmatic consumers.
    pub report: ExactBoundaryCapReport,
    /// Structured diagnostics emitted during decline.
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

impl ExactBoundaryCapCandidateResult {
    /// Number of exact curves prepared without approximation.
    pub(crate) fn prepared_curve_count(&self) -> usize {
        self.performance.prepared_curve_count
    }

    /// Number of generated inner spokes prepared without topology mutation.
    pub(crate) fn inner_spoke_count(&self) -> usize {
        self.performance.inner_spoke_count
    }

    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub(crate) fn ribbon_count(&self) -> usize {
        self.performance.support_ribbon_count
    }

    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub(crate) fn solved_surface_count(&self) -> usize {
        self.performance.solved_surface_count
    }

    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub(crate) fn pcurve_pair_count(&self) -> usize {
        self.performance.pcurve_pair_count
    }
}
