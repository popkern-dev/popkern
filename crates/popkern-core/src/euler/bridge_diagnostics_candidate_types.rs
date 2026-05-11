//! Bridge diagnostic construction-candidate data types.

use super::{DiagnosticOperationId, PostconditionPredictionId, ResourceConflictId};

/// Local graph id for a construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ConstructionCandidateId(pub usize);

/// One construction candidate selected from graph facts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstructionCandidate {
    /// Local construction candidate id.
    pub id: ConstructionCandidateId,
    /// Candidate construction family.
    pub family: ConstructionFamily,
    /// Candidate feasibility.
    pub feasibility: ConstructionFeasibility,
    /// Diagnostic operations participating in the candidate.
    pub operations: Vec<DiagnosticOperationId>,
    /// Resource conflicts participating in the candidate.
    pub conflicts: Vec<ResourceConflictId>,
    /// Requirements that must be satisfied before mutation.
    pub requirements: Vec<ConstructionRequirement>,
    /// Optional postcondition preview for this candidate.
    pub postcondition_preview: Option<PostconditionPredictionId>,
}

/// Candidate construction family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionFamily {
    /// Construction family has not been classified yet.
    Unknown,
    /// Disjoint cross-loop bridge.
    DisjointCrossLoopBridge,
    /// Disjoint same-loop bridge.
    DisjointSameLoopBridge,
    /// Paired same-loop twin exposure.
    PairedSameLoopTwinExposure,
    /// Boundary-side separator bridge.
    BoundarySideSeparatorBridge,
}

/// Candidate feasibility state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionFeasibility {
    /// Feasibility has not been classified yet.
    Unknown,
    /// Existing operator can execute the candidate.
    ExecutableByExistingOperator,
    /// Existing operator can execute after known normalization.
    ExecutableAfterKnownNormalization,
    /// Candidate requires a separator subcomplex.
    RequiresSeparatorSubcomplex,
    /// Candidate requires custom aggregate authoring.
    RequiresCustomAggregateAuthoring,
    /// Candidate is invalid because it would create non-manifold topology.
    InvalidNonManifold,
    /// Candidate family is unsupported or not yet classified.
    UnsupportedUnknownFamily,
}

/// Requirement that must be satisfied before candidate mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionRequirement {
    /// No requirement has been classified yet.
    Unknown,
    /// Revalidate the candidate against the live topology.
    LiveTopologyRevalidation,
    /// Verify postconditions before commit.
    PostconditionVerification,
}
