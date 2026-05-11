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
    pub(super) validity_repair_facts: Vec<ValidityRepairFact>,
    pub(super) operation_dependencies: Vec<OperationDependency>,
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
    /// Validity and repair readiness facts for this candidate.
    pub validity_repair_facts: Vec<ValidityRepairFact>,
    /// Neighboring operation dependencies for this candidate.
    pub operation_dependencies: Vec<OperationDependency>,
}

/// Dependency on a neighboring primitive or modeling operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OperationDependency {
    /// Dependency family.
    pub readiness: OperationDependencyReadiness,
    /// Dependency status.
    pub status: OperationDependencyStatus,
    /// Stable reason label.
    pub reason: &'static str,
    /// Named primitive or operation needed by this dependency.
    pub primitive: &'static str,
}

/// Operation dependency family for a construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperationDependencyReadiness {
    /// Dependency has not been classified.
    Unknown,
    /// No neighboring operation dependency is relevant.
    NotApplicable,
    /// Dependency is already executable.
    Executable,
    /// Dependency is ambiguous under current facts.
    Ambiguous,
    /// Dependency is invalid for this candidate.
    Invalid,
    /// A named lower-level primitive is missing.
    MissingPrimitive,
}

/// Execution status for an operation dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperationDependencyStatus {
    /// Status has not been classified.
    Unknown,
    /// Dependency is satisfied by existing facts/operators.
    Satisfied,
    /// Dependency is not needed.
    NotRequired,
    /// Dependency needs a named primitive before construction.
    RequiredMissingPrimitive,
    /// Dependency is blocked by invalid input.
    BlockedInvalidInput,
    /// Dependency is ambiguous and needs a resolver.
    AmbiguousNeedsResolver,
}

/// Validity and repair readiness for a construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValidityRepairReadiness {
    /// Validity/repair state has not been classified yet.
    Unknown,
    /// No validity/repair fact is relevant.
    NotApplicable,
    /// Candidate input is valid as-is.
    ValidAsInput,
    /// Candidate is valid after non-mutating normalization.
    ValidAfterNonMutatingNormalization,
    /// Candidate requires topology repair.
    RequiresTopologyRepair,
    /// Candidate requires geometry repair.
    RequiresGeometryRepair,
    /// Candidate repair is rejected because it would alter design intent.
    RejectedDesignIntentChangingRepair,
}

/// One validity/repair fact attached to a construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ValidityRepairFact {
    /// Readiness state.
    pub readiness: ValidityRepairReadiness,
    /// Stable reason label.
    pub reason: &'static str,
    /// Semantic effect of any required repair.
    pub status: ValidityRepairStatus,
    /// Whether required repair would change semantic ownership.
    pub changes_semantic_ownership: bool,
    /// Whether required repair would change user-visible design intent.
    pub changes_user_visible_intent: bool,
}

/// Repair status for a validity/repair fact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ValidityRepairStatus {
    /// Repair status has not been classified yet.
    Unknown,
    /// No repair is needed.
    NoRepairRequired,
    /// Non-mutating normalization can satisfy validity.
    NonMutatingNormalizationAvailable,
    /// Mutating repair is required before construction.
    RepairRequired,
    /// Repair is rejected for this candidate.
    Rejected,
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
