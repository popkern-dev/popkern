//! Euler operators — the only way to modify BRep topology.
//!
//! Each operator maintains the Euler-Poincaré invariant by construction.

mod boundary_side_separator_bridge;
mod boundary_side_separator_bridge_apply;
mod boundary_side_separator_bridge_apply_types;
mod bridge_diagnostics;
mod construction_residual_retire;
mod contract;
mod corner_apex_diagnostics;
mod cylinder_mediated_support;
mod cylinder_mediated_support_apply;
mod cylinder_mediated_support_partition;
mod face_boundary_fragment;
mod face_boundary_retire;
mod face_family_bridge;
mod face_merge;
mod face_merge_split;
mod kemr;
mod kfmrh;
mod loop_bridge;
mod loop_bridge_apply;
mod mef;
mod mev;
mod mvfs;
mod same_loop_boundary_exposure;
mod same_loop_bridge_dissolve;
mod same_loop_face_bridge;
mod sweep;
pub use boundary_side_separator_bridge::{
    BoundarySideSeparatorBridgePlan, BoundarySideSeparatorBridgeRequest,
    BoundarySideSeparatorBridgeTopologyDelta,
};
pub use boundary_side_separator_bridge_apply::BoundarySideSeparatorBridgeResult;
pub use bridge_diagnostics::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticReport, BridgeDiagnosticRequest,
    BridgeDiagnosticSummary, CandidateCompletenessError, CandidateMissingFacet,
    ConstructionCandidate, ConstructionCandidateId, ConstructionFamily, ConstructionFeasibility,
    ConstructionRequirement, DiagnosticHandle, DiagnosticOperation, DiagnosticOperationId,
    DiagnosticOperationKind, DiagnosticOperationPurpose, DiagnosticOperatorFamily, DiagnosticOwner,
    DiagnosticResource, DiagnosticResourceId, DiagnosticResourceKind, DiagnosticResourceStatus,
    DiagnosticTolerance, DiagnosticToleranceSource, DiagnosticTopologyDelta, DownstreamReadiness,
    DownstreamReadinessFact, DownstreamReadinessStatus, GeometryFeasibility,
    GeometryFeasibilityFact, OperationDependency, OperationDependencyReadiness,
    OperationDependencyStatus, OperatorReadiness, OperatorReadinessStatus, PostconditionPrediction,
    PostconditionPredictionId, ResourceConflict, ResourceConflictId, ResourceConflictKind,
    ResourceConflictStatus, TopologyDiagnosticGraph, TopologyGraphFingerprint, ValidityRepairFact,
    ValidityRepairReadiness, ValidityRepairStatus,
};
pub use construction_residual_retire::ConstructionResidualRetireResult;
pub use corner_apex_diagnostics::{
    ApexFeasibilityScoreFact, ApexFeasibilityScoreIngredient, ApexFeasibilityScoreMissing,
    ApexSizingFact, ApexSizingInput, ApexSizingMode, ApexSizingState, BlendKind,
    BlendTypeCompatibility, BlendValence, CapRunoutMode, CenterLocusFeasibilityFact,
    CenterLocusFeasibilityState, CliffTerminationMode, ConstructionRiskAction,
    ConstructionRiskFact, ConstructionRiskKind, ConstructionRiskState, CornerApexConditionCategory,
    CornerApexConditionId, CornerApexCoverageFacetPayload, CornerApexCoverageRow,
    CornerApexCoverageValidationError, CornerApexDiagnosticReport, CornerApexDiagnosticRequest,
    CornerApexFactContract, CornerApexStrategyFamily, CornerApexStrategyReadiness,
    CornerGeometryAvailability, CornerGeometryFact, CornerSupportAdequacyFact,
    CornerSupportAdequacyState, CornerTopologyRecord, CornerTopologyState, DiagnosticArtifactKind,
    EdgeConvexity, GeometryProbeDependency, IntegrationValidationAction, IntegrationValidationFact,
    IntegrationValidationKind, IntegrationValidationState, LifecycleAction, LifecycleFact,
    LifecycleState, LocalBlendEdgeSelection, LocalCornerScope, MutatorReadinessRow,
    OverflowHandlingMode, PatchCandidateApplicability, PatchCandidateDowngradePolicy,
    PatchCandidateFact, PatchCandidateFeasibility, PatchCandidateQuality, StrategyMetadataFact,
    StrategyOptionReadiness, StrategyReadinessFact, SupportFootprintFeasibilityFact,
    SupportFootprintFeasibilityState, SurfaceCompatibility, SurfaceRepresentation,
    TopologicalConfiguration, TrialInteractionFact, TrialInteractionState,
    VertexClassificationFact, VertexConvexityClass, classify_apex_sizing,
    classify_blend_type_compatibility, classify_blend_valence, classify_construction_risk,
    classify_convexity, classify_integration_validation, classify_lifecycle,
    classify_strategy_readiness, classify_support_adequacy, classify_topological_configuration,
    classify_trial_interaction, classify_vertex_classification, corner_apex_condition_ids,
    corner_apex_fact_contracts, validate_corner_apex_coverage_row, validate_mutator_readiness_row,
};
pub use cylinder_mediated_support::{
    CylinderMediatedOutsideIntervalInteraction, CylinderMediatedSupportContinuationAggregateDelta,
    CylinderMediatedSupportContinuationBatchPlan, CylinderMediatedSupportContinuationDelta,
    CylinderMediatedSupportContinuationPlan, CylinderMediatedSupportContinuationRequest,
};
pub use cylinder_mediated_support_apply::{
    CylinderMediatedContextOwnerFace, CylinderMediatedOwnerRepairSpan,
    CylinderMediatedSupportContinuationBatchResult,
    CylinderMediatedSupportContinuationPostconditions,
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement, CylinderMediatedUnsupportedOwnerFace,
};
pub use cylinder_mediated_support_partition::{
    CylinderMediatedOutsideIntervalPartition, CylinderMediatedOutsideIntervalPartitionSegment,
    CylinderMediatedOutsideLoopPartition,
};
pub use face_boundary_fragment::BoundaryFragmentRetireResult;
pub use face_boundary_retire::BoundaryFaceRetireResult;
pub use face_family_bridge::{
    FaceFamilyBridgeCycleRequest, FaceFamilyBridgeCycleResult, FaceFamilyBridgePlan,
    FaceFamilyBridgeResult,
};
pub use face_merge::FaceMergeResult;
pub use face_merge_split::{FaceMergeSplitRequest, FaceMergeSplitResult};
pub use kemr::{KemrResult, MekrResult};
pub use kfmrh::MfkrhResult;
pub use loop_bridge::{TwoLoopBridgePlan, TwoLoopBridgeRequest};
pub use loop_bridge_apply::TwoLoopBridgeResult;
pub use mef::MefResult;
pub use mev::MevResult;
pub use mvfs::MvfsResult;
pub use same_loop_boundary_exposure::SameLoopBoundaryExposureResult;
pub use same_loop_face_bridge::{
    SameLoopFaceBridgePlan, SameLoopFaceBridgeRequest, SameLoopFaceBridgeResult,
};
pub use sweep::ExtrudeResult;
pub use sweep::RevolveResult;
