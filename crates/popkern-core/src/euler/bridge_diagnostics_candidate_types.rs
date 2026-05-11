//! Bridge diagnostic construction-candidate data types.

use super::{DiagnosticOperationId, PostconditionPredictionId, ResourceConflictId};

pub(in crate::euler) struct ConstructionCandidateDraft {
    pub(super) family: ConstructionFamily,
    pub(super) feasibility: ConstructionFeasibility,
    pub(super) feasibility_reason: &'static str,
    pub(super) operations: Vec<DiagnosticOperationId>,
    pub(super) conflicts: Vec<ResourceConflictId>,
    pub(super) requirements: Vec<ConstructionRequirement>,
    pub(super) postcondition_preview: Option<PostconditionPredictionId>,
    pub(super) geometry_facts: Vec<GeometryFeasibilityFact>,
}

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
    /// Stable reason for the feasibility classification.
    pub feasibility_reason: &'static str,
    /// Diagnostic operations participating in the candidate.
    pub operations: Vec<DiagnosticOperationId>,
    /// Resource conflicts participating in the candidate.
    pub conflicts: Vec<ResourceConflictId>,
    /// Requirements that must be satisfied before mutation.
    pub requirements: Vec<ConstructionRequirement>,
    /// Optional postcondition preview for this candidate.
    pub postcondition_preview: Option<PostconditionPredictionId>,
    /// Geometry feasibility facts for this candidate.
    pub geometry_facts: Vec<GeometryFeasibilityFact>,
}

/// Geometry feasibility state for a construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GeometryFeasibility {
    /// Geometry state has not been classified yet.
    Unknown,
    /// No geometry facts are needed for this candidate.
    NotApplicable,
    /// Candidate geometry is ready for construction.
    GeometryReady,
    /// Candidate uses an unsupported geometry family.
    UnsupportedGeometryFamily,
    /// Geometry repair must run before construction.
    RepairRequiredBeforeConstruction,
}

/// One geometry feasibility fact attached to a construction candidate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryFeasibilityFact {
    /// Geometry feasibility state.
    pub feasibility: GeometryFeasibility,
    /// Stable reason label.
    pub reason: &'static str,
    /// Local explicit tolerance used by the fact.
    pub tolerance: DiagnosticTolerance,
}

/// Explicit local tolerance recorded by a diagnostic fact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagnosticTolerance {
    /// Source of the tolerance.
    pub source: DiagnosticToleranceSource,
    /// Magnitude in nanometers when the fact used a numeric tolerance.
    pub nanometers: Option<i64>,
}

impl DiagnosticTolerance {
    /// No numeric tolerance is relevant to this fact.
    pub const NOT_APPLICABLE: Self = Self {
        source: DiagnosticToleranceSource::NotApplicable,
        nanometers: None,
    };
}

/// Source of a diagnostic tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticToleranceSource {
    /// No tolerance is relevant to this fact.
    NotApplicable,
    /// Tolerance came from a local BRep or entity tolerance.
    LocalTopology,
    /// Tolerance came from a local geometry predicate.
    LocalGeometryPredicate,
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
    /// Candidate requires regrouping bridge cycles before mutation.
    RequiresCycleRegrouping,
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
