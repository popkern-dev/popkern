//! Read-only corner/apex diagnostic facts for fillet and blend planning.

#[path = "corner_apex_diagnostics_classifiers.rs"]
mod classifiers;
#[path = "corner_apex_diagnostics_live.rs"]
mod live;
#[path = "corner_apex_diagnostics_live_helpers.rs"]
mod live_helpers;
#[path = "corner_apex_diagnostics_live_primitives.rs"]
mod live_primitives;
#[path = "corner_apex_diagnostics_manifest.rs"]
mod manifest;
#[path = "corner_apex_diagnostics_patch.rs"]
mod patch;
#[path = "corner_apex_diagnostics_patch_types.rs"]
mod patch_types;
#[path = "corner_apex_diagnostics_score_types.rs"]
mod score_types;
#[path = "corner_apex_diagnostics_selection.rs"]
mod selection;
#[path = "corner_apex_diagnostics_selection_types.rs"]
mod selection_types;
#[path = "corner_apex_diagnostics_strategy.rs"]
mod strategy;
#[path = "corner_apex_diagnostics_strategy_types.rs"]
mod strategy_types;
#[path = "corner_apex_diagnostics_support_types.rs"]
mod support_types;
#[path = "corner_apex_diagnostics_three_face_tangency.rs"]
mod three_face_tangency;
#[path = "corner_apex_diagnostics_topology_neighborhood.rs"]
mod topology_neighborhood;
#[path = "corner_apex_diagnostics_topology_validity.rs"]
mod topology_validity;
#[path = "corner_apex_diagnostics_types.rs"]
mod types;

pub use classifiers::{
    classify_apex_sizing, classify_blend_type_compatibility, classify_blend_valence,
    classify_construction_risk, classify_convexity, classify_integration_validation,
    classify_lifecycle, classify_strategy_readiness, classify_support_adequacy,
    classify_topological_configuration, classify_trial_interaction, classify_vertex_classification,
};
pub use manifest::{
    CornerApexConditionCategory, CornerApexConditionId, CornerApexFactContract,
    corner_apex_condition_ids, corner_apex_fact_contracts,
};
#[allow(unused_imports)]
pub use patch_types::{
    PatchCandidateApplicability, PatchCandidateDowngradePolicy, PatchCandidateFact,
    PatchCandidateFeasibility, PatchCandidateQuality,
};
#[allow(unused_imports)]
pub use score_types::{
    ApexFeasibilityScoreFact, ApexFeasibilityScoreIngredient, ApexFeasibilityScoreMissing,
};
#[allow(unused_imports)]
pub use selection_types::{
    CornerSelectionIntentFact, CornerSelectionOptionState, CornerSelectionPropagation,
    CornerSelectionSource,
};
#[allow(unused_imports)]
pub use strategy_types::{
    AdvancedStrategyMode, CapRunoutMode, CliffTerminationMode, OverflowHandlingMode,
    StrategyMetadataFact, StrategyOptionReadiness,
};
#[allow(unused_imports)]
pub use support_types::{
    CenterLocusFeasibilityFact, CenterLocusFeasibilityState, SupportFootprintFeasibilityFact,
    SupportFootprintFeasibilityState,
};
pub use types::{
    ApexSizingFact, ApexSizingInput, ApexSizingMode, ApexSizingState, BlendKind,
    BlendTypeCompatibility, BlendValence, ConstructionRiskFact, ConstructionRiskKind,
    ConstructionRiskState, CornerApexDiagnosticReport, CornerApexDiagnosticRequest,
    CornerApexStrategyFamily, CornerApexStrategyReadiness, CornerGeometryAvailability,
    CornerGeometryFact, CornerSupportAdequacyFact, CornerSupportAdequacyState,
    CornerTopologyRecord, CornerTopologyState, DiagnosticArtifactKind, EdgeConvexity,
    GeometryProbeDependency, IntegrationValidationFact, IntegrationValidationKind,
    IntegrationValidationState, LifecycleFact, LifecycleState, LocalBlendEdgeSelection,
    LocalCornerScope, StrategyReadinessFact, SurfaceCompatibility, SurfaceRepresentation,
    TopologicalConfiguration, TrialInteractionFact, TrialInteractionState,
    VertexClassificationFact, VertexConvexityClass,
};

#[cfg(test)]
#[path = "corner_apex_diagnostics_high_valence_tests.rs"]
mod high_valence_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_live_primitive_tests.rs"]
mod live_primitive_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_live_tests.rs"]
mod live_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_selection_tests.rs"]
mod selection_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_tests.rs"]
mod tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_three_face_live_tests.rs"]
mod three_face_live_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_three_face_negative_tests.rs"]
mod three_face_negative_tests;
#[cfg(test)]
#[path = "corner_apex_diagnostics_topology_validity_tests.rs"]
mod topology_validity_tests;
