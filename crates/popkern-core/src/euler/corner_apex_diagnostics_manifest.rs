//! Corner/apex diagnostic coverage manifest and fact contracts.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerApexConditionCategory {
    TopologyExtraction,
    SurfaceCompatibility,
    GeometryFactAvailability,
    ValenceConfiguration,
    Convexity,
    BlendTypeCompatibility,
    ConservativeExposure,
    SizingRangeClearance,
    TrialInteraction,
    InsideTight,
    SupportAdequacy,
    StrategyReadiness,
    ConstructionRisk,
    IntegrationValidation,
    ContinuityVerification,
    SemanticSurgeryHistory,
    DownstreamReadiness,
    LifecycleFreshness,
    DiagnosticArtifactCleanup,
}

impl CornerApexConditionCategory {
    pub fn diagnostic_label(self) -> &'static str {
        match self {
            Self::TopologyExtraction => "topology-extraction",
            Self::SurfaceCompatibility => "surface-compatibility",
            Self::GeometryFactAvailability => "geometry-fact-availability",
            Self::ValenceConfiguration => "valence-configuration",
            Self::Convexity => "convexity",
            Self::BlendTypeCompatibility => "blend-type-compatibility",
            Self::ConservativeExposure => "conservative-exposure",
            Self::SizingRangeClearance => "sizing-range-clearance",
            Self::TrialInteraction => "trial-interaction",
            Self::InsideTight => "inside-tight",
            Self::SupportAdequacy => "support-adequacy",
            Self::StrategyReadiness => "strategy-readiness",
            Self::ConstructionRisk => "construction-risk",
            Self::IntegrationValidation => "integration-validation",
            Self::ContinuityVerification => "continuity-verification",
            Self::SemanticSurgeryHistory => "semantic-surgery-history",
            Self::DownstreamReadiness => "downstream-readiness",
            Self::LifecycleFreshness => "lifecycle-freshness",
            Self::DiagnosticArtifactCleanup => "diagnostic-artifact-cleanup",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CornerApexConditionId {
    pub category: CornerApexConditionCategory,
    pub id: &'static str,
}

impl CornerApexConditionId {
    pub const fn new(category: CornerApexConditionCategory, id: &'static str) -> Self {
        Self { category, id }
    }

    pub fn diagnostic_label(self) -> String {
        format!("corner_{}_{}", self.category.diagnostic_label(), self.id)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CornerApexFactContract {
    pub condition: CornerApexConditionId,
    pub source_inputs: &'static [&'static str],
    pub output_states: &'static [&'static str],
    pub tolerance_source: &'static str,
    pub not_applicable: &'static str,
    pub unknown: &'static str,
    pub compact_label_fragment: &'static str,
    pub matrix_assertion_shape: &'static str,
}

macro_rules! condition {
    ($category:ident, $id:literal) => {
        CornerApexConditionId::new(CornerApexConditionCategory::$category, $id)
    };
}

pub const CORNER_APEX_CONDITION_IDS: [CornerApexConditionId; 45] = [
    condition!(TopologyExtraction, "local-neighborhood"),
    condition!(TopologyExtraction, "stale-or-missing-handle"),
    condition!(SurfaceCompatibility, "analytic-support"),
    condition!(SurfaceCompatibility, "conditional-support"),
    condition!(GeometryFactAvailability, "tangent-normal-dihedral"),
    condition!(GeometryFactAvailability, "blocked-geometry-dependency"),
    condition!(ValenceConfiguration, "free-end"),
    condition!(ValenceConfiguration, "through-vertex"),
    condition!(ValenceConfiguration, "standard-corner"),
    condition!(ValenceConfiguration, "high-valence"),
    condition!(ValenceConfiguration, "invalid-topology"),
    condition!(Convexity, "all-convex"),
    condition!(Convexity, "all-concave"),
    condition!(Convexity, "mixed"),
    condition!(Convexity, "near-tangent"),
    condition!(Convexity, "variable-stop"),
    condition!(BlendTypeCompatibility, "homogeneous"),
    condition!(BlendTypeCompatibility, "parameter-varying"),
    condition!(BlendTypeCompatibility, "heterogeneous"),
    condition!(BlendTypeCompatibility, "incompatible"),
    condition!(ConservativeExposure, "overflow-cliff-collision"),
    condition!(SizingRangeClearance, "range-ready"),
    condition!(SizingRangeClearance, "range-too-large"),
    condition!(SizingRangeClearance, "adjacent-collision"),
    condition!(SizingRangeClearance, "suggest-smaller-radius"),
    condition!(TrialInteraction, "mutual-trim-conditioning"),
    condition!(TrialInteraction, "missing-intersection"),
    condition!(InsideTight, "inside-tight"),
    condition!(InsideTight, "missing-center-locus"),
    condition!(SupportAdequacy, "footprint-fits"),
    condition!(SupportAdequacy, "overflow-cliff-hole-trim"),
    condition!(StrategyReadiness, "core-strategy-family"),
    condition!(StrategyReadiness, "cap-deferred-fallback-quality"),
    condition!(ConstructionRisk, "intersection-tangency"),
    condition!(ConstructionRisk, "patch-boundary-continuity"),
    condition!(ConstructionRisk, "approximation-mode-retry"),
    condition!(ConstructionRisk, "numerical-solver-failure"),
    condition!(IntegrationValidation, "watertightness-intersection"),
    condition!(IntegrationValidation, "tolerance-exceeded"),
    condition!(ContinuityVerification, "g0-g1-g2"),
    condition!(SemanticSurgeryHistory, "topology-surgery-preview"),
    condition!(SemanticSurgeryHistory, "merge-readiness"),
    condition!(DownstreamReadiness, "mesh-export-risk"),
    condition!(LifecycleFreshness, "stale-reuse-refusal"),
    condition!(DiagnosticArtifactCleanup, "artifact-cleanup"),
];

pub fn corner_apex_condition_ids() -> &'static [CornerApexConditionId] {
    &CORNER_APEX_CONDITION_IDS
}

pub const CORNER_APEX_FACT_CONTRACTS: [CornerApexFactContract; 9] = [
    CornerApexFactContract {
        condition: condition!(TopologyExtraction, "local-neighborhood"),
        source_inputs: &["corner-request", "brep-topology", "semantic-roles"],
        output_states: &["ready", "stale-handle", "missing-handle", "non-manifold"],
        tolerance_source: "not-applicable",
        not_applicable: "request-has-no-corner-scope",
        unknown: "unknown-input-fact",
        compact_label_fragment: "corner_topology",
        matrix_assertion_shape: "condition-id+six-facet-status",
    },
    CornerApexFactContract {
        condition: condition!(GeometryFactAvailability, "tangent-normal-dihedral"),
        source_inputs: &["edge-curve", "support-surface", "local-tolerance"],
        output_states: &["available", "missing-input-fact", "blocked-by-dependency"],
        tolerance_source: "local-geometry-predicate",
        not_applicable: "topology-only-row",
        unknown: "unknown-input-fact",
        compact_label_fragment: "corner_geometry_facts",
        matrix_assertion_shape: "geometry-feasibility+operation-dependency",
    },
    CornerApexFactContract {
        condition: condition!(Convexity, "mixed"),
        source_inputs: &["dihedral-facts", "selected-blend-edges"],
        output_states: &["all-convex", "all-concave", "mixed", "near-tangent"],
        tolerance_source: "local-geometry-predicate",
        not_applicable: "not-blend-vertex",
        unknown: "unknown-input-fact",
        compact_label_fragment: "corner_convexity",
        matrix_assertion_shape: "classification+strategy-readiness",
    },
    CornerApexFactContract {
        condition: condition!(SizingRangeClearance, "range-ready"),
        source_inputs: &["radius-law", "dihedral", "edge-length", "face-margin"],
        output_states: &["ready", "range-too-large", "adjacent-collision"],
        tolerance_source: "local-topology-or-geometry",
        not_applicable: "no-sizing-requested",
        unknown: "unknown-tolerance-model",
        compact_label_fragment: "apex_sizing",
        matrix_assertion_shape: "geometry-feasibility+tolerance-source",
    },
    CornerApexFactContract {
        condition: condition!(TrialInteraction, "mutual-trim-conditioning"),
        source_inputs: &["candidate-pipes", "intersection-query", "tolerance"],
        output_states: &["well-conditioned", "near-tangent", "ambiguous"],
        tolerance_source: "local-geometry-predicate",
        not_applicable: "not-valence-two",
        unknown: "unknown-input-fact",
        compact_label_fragment: "trial_interaction",
        matrix_assertion_shape: "operation-dependency+geometry-feasibility",
    },
    CornerApexFactContract {
        condition: condition!(SupportAdequacy, "overflow-cliff-hole-trim"),
        source_inputs: &["support-footprint", "face-trim-domain", "inner-loops"],
        output_states: &["fits", "overflow", "cliff", "trim-domain-exit"],
        tolerance_source: "local-geometry-predicate",
        not_applicable: "no-support-face",
        unknown: "unknown-input-fact",
        compact_label_fragment: "support_adequacy",
        matrix_assertion_shape: "geometry-feasibility+downstream-readiness",
    },
    CornerApexFactContract {
        condition: condition!(StrategyReadiness, "core-strategy-family"),
        source_inputs: &["classification", "sizing", "support-adequacy"],
        output_states: &["ready", "managed-risk", "missing-primitive", "invalid"],
        tolerance_source: "derived-from-input-facts",
        not_applicable: "no-construction-candidate",
        unknown: "unknown-input-fact",
        compact_label_fragment: "strategy_readiness",
        matrix_assertion_shape: "candidate-family+expected-delta",
    },
    CornerApexFactContract {
        condition: condition!(ConstructionRisk, "patch-boundary-continuity"),
        source_inputs: &["strategy", "boundary-curves", "continuity-target"],
        output_states: &["no-risk", "managed-retry", "missing-primitive", "invalid"],
        tolerance_source: "local-geometry-predicate",
        not_applicable: "non-patch-strategy",
        unknown: "unknown-input-fact",
        compact_label_fragment: "construction_risk",
        matrix_assertion_shape: "operation-dependency+downstream-readiness",
    },
    CornerApexFactContract {
        condition: condition!(IntegrationValidation, "watertightness-intersection"),
        source_inputs: &["postcondition-preview", "continuity-target", "tolerance"],
        output_states: &["ready", "repair-required", "invalid", "unknown-input-fact"],
        tolerance_source: "postcondition-preview",
        not_applicable: "diagnostic-only-refusal",
        unknown: "unknown-input-fact",
        compact_label_fragment: "integration_validation",
        matrix_assertion_shape: "postcondition-preview+downstream-readiness",
    },
];

pub fn corner_apex_fact_contracts() -> &'static [CornerApexFactContract] {
    &CORNER_APEX_FACT_CONTRACTS
}
