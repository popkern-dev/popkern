use super::*;
use crate::euler::{DiagnosticHandle, DiagnosticTolerance, DiagnosticToleranceSource};

fn tolerance() -> DiagnosticTolerance {
    DiagnosticTolerance {
        source: DiagnosticToleranceSource::LocalGeometryPredicate,
        nanometers: Some(10),
    }
}

fn edge(blend_kind: BlendKind, radius_nanometers: Option<i64>) -> LocalBlendEdgeSelection {
    LocalBlendEdgeSelection {
        edge: DiagnosticHandle::None,
        blend_kind,
        radius_nanometers,
    }
}

fn topology(incident_blend_edges: usize, state: CornerTopologyState) -> CornerTopologyRecord {
    CornerTopologyRecord {
        vertex: DiagnosticHandle::None,
        incident_blend_edges,
        incident_non_blend_edges: 0,
        surrounding_faces: 3,
        face_adjacency_edges: 3,
        neighborhood: Default::default(),
        validity: Default::default(),
        state,
        geometry: CornerGeometryFact {
            tangent: CornerGeometryAvailability::Available,
            normal: CornerGeometryAvailability::Available,
            dihedral: CornerGeometryAvailability::Available,
            tolerance: tolerance(),
        },
        surface_representation: SurfaceRepresentation::Plane,
        surface_compatibility: SurfaceCompatibility::Full,
    }
}

#[test]
fn corner_apex_condition_manifest_names_every_required_group() {
    let ids = corner_apex_condition_ids();
    assert!(ids.len() >= 40);
    for id in ids {
        let label = id.diagnostic_label();
        assert!(
            label.starts_with("corner_"),
            "condition id must be stable and scoped: {label}"
        );
        assert!(!id.category.diagnostic_label().is_empty());
        assert!(!id.id.is_empty());
    }
    for expected in [
        CornerApexConditionCategory::TopologyExtraction,
        CornerApexConditionCategory::SurfaceCompatibility,
        CornerApexConditionCategory::GeometryFactAvailability,
        CornerApexConditionCategory::ValenceConfiguration,
        CornerApexConditionCategory::Convexity,
        CornerApexConditionCategory::BlendTypeCompatibility,
        CornerApexConditionCategory::ConservativeExposure,
        CornerApexConditionCategory::SizingRangeClearance,
        CornerApexConditionCategory::TrialInteraction,
        CornerApexConditionCategory::InsideTight,
        CornerApexConditionCategory::SupportAdequacy,
        CornerApexConditionCategory::StrategyReadiness,
        CornerApexConditionCategory::ConstructionRisk,
        CornerApexConditionCategory::IntegrationValidation,
        CornerApexConditionCategory::ContinuityVerification,
        CornerApexConditionCategory::SemanticSurgeryHistory,
        CornerApexConditionCategory::DownstreamReadiness,
        CornerApexConditionCategory::LifecycleFreshness,
        CornerApexConditionCategory::DiagnosticArtifactCleanup,
    ] {
        assert!(
            ids.iter().any(|id| id.category == expected),
            "missing condition category {}",
            expected.diagnostic_label()
        );
    }
}

#[test]
fn corner_apex_fact_contracts_are_implementation_ready() {
    for contract in corner_apex_fact_contracts() {
        assert!(!contract.source_inputs.is_empty());
        assert!(!contract.output_states.is_empty());
        assert!(!contract.tolerance_source.is_empty());
        assert!(!contract.not_applicable.is_empty());
        assert!(!contract.unknown.is_empty());
        assert!(!contract.compact_label_fragment.is_empty());
        assert!(!contract.matrix_assertion_shape.is_empty());
        assert!(
            corner_apex_condition_ids().contains(&contract.condition),
            "contract references condition not present in manifest: {}",
            contract.condition.diagnostic_label()
        );
    }
}

#[test]
fn valence_topology_convexity_and_blend_type_classifiers_are_named() {
    assert_eq!(classify_blend_valence(0), BlendValence::NotBlendVertex);
    assert_eq!(classify_blend_valence(1), BlendValence::FreeEnd);
    assert_eq!(classify_blend_valence(2), BlendValence::ThroughVertex);
    assert_eq!(classify_blend_valence(3), BlendValence::StandardCorner);
    assert_eq!(classify_blend_valence(4), BlendValence::HighValenceCorner);

    assert_eq!(
        classify_convexity(&[EdgeConvexity::Convex, EdgeConvexity::Concave]),
        VertexConvexityClass::Mixed
    );
    assert_eq!(
        classify_convexity(&[EdgeConvexity::Convex, EdgeConvexity::Convex]),
        VertexConvexityClass::AllConvex
    );
    assert_eq!(
        classify_convexity(&[EdgeConvexity::Concave, EdgeConvexity::Concave]),
        VertexConvexityClass::AllConcave
    );
    assert_eq!(
        classify_convexity(&[EdgeConvexity::VariableNearVertex]),
        VertexConvexityClass::VariableConvexityStop
    );

    let homogeneous = [
        edge(BlendKind::FilletConstantRadius, Some(1)),
        edge(BlendKind::FilletConstantRadius, Some(1)),
    ];
    assert_eq!(
        classify_blend_type_compatibility(&homogeneous, BlendValence::ThroughVertex),
        BlendTypeCompatibility::Homogeneous
    );
    let selected = [
        edge(BlendKind::FilletConstantRadius, Some(1)),
        edge(BlendKind::FilletConstantRadius, Some(2)),
    ];
    assert_eq!(
        classify_blend_type_compatibility(&selected, BlendValence::ThroughVertex),
        BlendTypeCompatibility::ParameterVarying
    );
    let heterogeneous = [
        edge(BlendKind::FilletConstantRadius, Some(1)),
        edge(BlendKind::Chamfer, None),
    ];
    assert_eq!(
        classify_blend_type_compatibility(&heterogeneous, BlendValence::ThroughVertex),
        BlendTypeCompatibility::Heterogeneous
    );
    let incompatible = [
        edge(BlendKind::G2Blend, Some(1)),
        edge(BlendKind::Chamfer, None),
    ];
    assert_eq!(
        classify_blend_type_compatibility(&incompatible, BlendValence::ThroughVertex),
        BlendTypeCompatibility::Incompatible("g2-chamfer-corner")
    );
    assert_eq!(
        classify_topological_configuration(&topology(1, CornerTopologyState::Ready)),
        TopologicalConfiguration::FreeEnd
    );
    assert_eq!(
        classify_topological_configuration(&topology(2, CornerTopologyState::Ready)),
        TopologicalConfiguration::ThroughVertex
    );
    assert_eq!(
        classify_topological_configuration(&topology(3, CornerTopologyState::Ready)),
        TopologicalConfiguration::SimpleCorner
    );
    assert_eq!(
        classify_topological_configuration(&topology(4, CornerTopologyState::Ready)),
        TopologicalConfiguration::HighValence
    );
}

#[test]
fn sizing_trial_support_and_strategy_classifiers_expose_missing_primitives() {
    let sizing = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::ConstantRadius,
        requested_extent_nm: Some(20),
        available_extent_nm: Some(10),
        adjacent_clearance_nm: Some(100),
        curvature_radius_nm: Some(100),
        tolerance: tolerance(),
    });
    assert_eq!(sizing.state, ApexSizingState::RadiusTooLarge);
    assert_eq!(sizing.max_feasible_radius_nm, Some(10));

    let trial = classify_trial_interaction(BlendValence::ThroughVertex, false, None, tolerance());
    assert_eq!(
        trial.state,
        TrialInteractionState::MissingIntersectionPrimitive
    );
    assert_eq!(
        trial.dependency,
        Some(GeometryProbeDependency::Intersection)
    );

    let support = classify_support_adequacy(
        CornerSupportAdequacyState::Cliff,
        Some(GeometryProbeDependency::SupportFootprint),
        tolerance(),
    );
    let classification = classify_vertex_classification(
        &topology(2, CornerTopologyState::Ready),
        &[edge(BlendKind::FilletConstantRadius, Some(1))],
        &[EdgeConvexity::Convex],
        support.state,
    );
    let strategy = classify_strategy_readiness(classification, sizing, support);
    assert_eq!(strategy.family, CornerApexStrategyFamily::Setback);
    assert_eq!(
        strategy.readiness,
        CornerApexStrategyReadiness::RequiresMissingPrimitive("radius-or-setback-adjustment")
    );
}

#[test]
fn construction_integration_and_lifecycle_classifiers_are_typed() {
    let risk = classify_construction_risk(
        ConstructionRiskKind::PatchBoundaryIncompatible,
        false,
        tolerance(),
    );
    assert_eq!(
        risk.state,
        ConstructionRiskState::RequiresMissingPrimitive("patch-boundary-repair")
    );
    let numerical = classify_construction_risk(
        ConstructionRiskKind::NumericalSolverFailure,
        false,
        tolerance(),
    );
    assert_eq!(
        numerical.state,
        ConstructionRiskState::RequiresMissingPrimitive("numerical-solver-failure")
    );

    let integration = classify_integration_validation(
        IntegrationValidationKind::Continuity,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        integration.state,
        IntegrationValidationState::RepairRequired("continuity-repair-or-downgrade")
    );

    let lifecycle = classify_lifecycle(false, false, DiagnosticArtifactKind::DiagnosticArtifact);
    assert_eq!(lifecycle.state, LifecycleState::TopologyFingerprintMismatch);
    assert_eq!(
        lifecycle.artifact,
        DiagnosticArtifactKind::DiagnosticArtifact
    );
}
