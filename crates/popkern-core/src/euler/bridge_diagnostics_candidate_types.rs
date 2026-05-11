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
    pub(super) downstream_facts: Vec<DownstreamReadinessFact>,
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
    /// Downstream readiness facts for export, meshing, naming, and feature history.
    pub downstream_facts: Vec<DownstreamReadinessFact>,
}

impl ConstructionCandidate {
    /// Validate the minimum facets required before planner-driven mutation.
    pub fn validate_mutation_gate(&self) -> Result<(), CandidateCompletenessError> {
        if self.family == ConstructionFamily::Unknown {
            return Err(self.missing_mutation_gate_facet(CandidateMissingFacet::Family));
        }
        if self.feasibility == ConstructionFeasibility::Unknown {
            return Err(self.missing_mutation_gate_facet(CandidateMissingFacet::Feasibility));
        }
        if self.operations.is_empty() {
            return Err(self.missing_mutation_gate_facet(CandidateMissingFacet::TopologyResource));
        }
        if self.postcondition_preview.is_none() {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::PostconditionPreview)
            );
        }
        if self.geometry_facts.is_empty() {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::GeometryFeasibility)
            );
        }
        if self.validity_repair_facts.is_empty() {
            return Err(self.missing_mutation_gate_facet(CandidateMissingFacet::ValidityRepair));
        }
        if self.operation_dependencies.is_empty() {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::OperationDependency)
            );
        }
        if !self
            .operation_dependencies
            .iter()
            .any(|dependency| dependency.primitive.contains("history"))
        {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::SemanticHistoryDependency)
            );
        }
        if self.downstream_facts.is_empty() {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::DownstreamReadiness)
            );
        }
        if !self
            .requirements
            .contains(&ConstructionRequirement::LiveTopologyRevalidation)
        {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::LiveTopologyRevalidation)
            );
        }
        if !self
            .requirements
            .contains(&ConstructionRequirement::PostconditionVerification)
        {
            return Err(
                self.missing_mutation_gate_facet(CandidateMissingFacet::PostconditionVerification)
            );
        }
        Ok(())
    }

    fn missing_mutation_gate_facet(
        &self,
        missing_facet: CandidateMissingFacet,
    ) -> CandidateCompletenessError {
        CandidateCompletenessError {
            candidate: self.id,
            missing_facet,
        }
    }
}

/// Candidate completeness failure for planner-driven mutation gates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateCompletenessError {
    /// Candidate that failed the mutation gate.
    pub candidate: ConstructionCandidateId,
    /// Missing facet or field name.
    pub missing_facet: CandidateMissingFacet,
}

impl std::fmt::Display for CandidateCompletenessError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "candidate {:?} is missing {}",
            self.candidate,
            self.missing_facet.diagnostic_label()
        )
    }
}

impl std::error::Error for CandidateCompletenessError {}

/// Missing candidate facet required by the planner mutation gate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CandidateMissingFacet {
    /// Candidate family was not classified.
    Family,
    /// Candidate feasibility was not classified.
    Feasibility,
    /// Candidate has no topology/resource operation membership.
    TopologyResource,
    /// Candidate lacks a postcondition preview.
    PostconditionPreview,
    /// Candidate lacks geometry feasibility facts.
    GeometryFeasibility,
    /// Candidate lacks validity/repair facts.
    ValidityRepair,
    /// Candidate lacks neighboring operation dependencies.
    OperationDependency,
    /// Candidate lacks semantic history dependencies.
    SemanticHistoryDependency,
    /// Candidate lacks downstream readiness facts.
    DownstreamReadiness,
    /// Candidate lacks live topology revalidation.
    LiveTopologyRevalidation,
    /// Candidate lacks postcondition verification.
    PostconditionVerification,
}

impl CandidateMissingFacet {
    /// Stable diagnostic label for this missing facet.
    pub fn diagnostic_label(self) -> &'static str {
        match self {
            Self::Family => "candidate-family",
            Self::Feasibility => "candidate-feasibility",
            Self::TopologyResource => "topology-resource-facet",
            Self::PostconditionPreview => "postcondition-preview",
            Self::GeometryFeasibility => "geometry-feasibility",
            Self::ValidityRepair => "validity-repair",
            Self::OperationDependency => "operation-dependency",
            Self::SemanticHistoryDependency => "semantic-history-dependency",
            Self::DownstreamReadiness => "downstream-readiness",
            Self::LiveTopologyRevalidation => "live-topology-revalidation",
            Self::PostconditionVerification => "postcondition-verification",
        }
    }
}

/// Downstream readiness fact for generated topology and geometry consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DownstreamReadinessFact {
    /// Downstream readiness family.
    pub readiness: DownstreamReadiness,
    /// Downstream readiness status.
    pub status: DownstreamReadinessStatus,
    /// Stable reason label.
    pub reason: &'static str,
}

/// Downstream readiness for export, meshing, naming, and feature history.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DownstreamReadiness {
    /// Downstream readiness has not been classified.
    Unknown,
    /// Downstream readiness is not applicable.
    NotApplicable,
    /// Downstream consumers can use the candidate output as-is.
    Ready,
    /// Downstream risk is known and bounded by candidate requirements.
    RiskManaged,
    /// Candidate is blocked by a missing primitive.
    BlockedByMissingPrimitive,
    /// Candidate is blocked until repair is performed.
    BlockedByRepairRequired,
    /// Candidate output would be invalid for export or meshing.
    InvalidForExportOrMesh,
}

/// Downstream readiness status for generated topology and geometry consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DownstreamReadinessStatus {
    /// Downstream status has not been classified.
    Unknown,
    /// Downstream consumers are ready.
    Ready,
    /// Downstream consumers are blocked.
    Blocked,
    /// Downstream risk is known and bounded.
    RiskManaged,
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
