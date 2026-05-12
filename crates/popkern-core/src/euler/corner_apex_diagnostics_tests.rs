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
        dihedral_margin_nm: None,
        face_area_margin_nm: None,
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
fn sizing_modes_remain_distinct_from_radius_scalar() {
    for mode in [
        ApexSizingMode::RangeRange,
        ApexSizingMode::RangeAngle,
        ApexSizingMode::AngleRange,
        ApexSizingMode::VariableRadiusEndpoint,
    ] {
        let sizing = classify_apex_sizing(ApexSizingInput {
            mode,
            requested_extent_nm: Some(10),
            available_extent_nm: Some(20),
            adjacent_clearance_nm: Some(20),
            curvature_radius_nm: None,
            dihedral_margin_nm: None,
            face_area_margin_nm: None,
            tolerance: tolerance(),
        });
        assert_eq!(sizing.mode, mode);
        assert_eq!(sizing.state, ApexSizingState::Ready);
    }

    let missing = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::RangeAngle,
        requested_extent_nm: None,
        available_extent_nm: Some(20),
        adjacent_clearance_nm: Some(20),
        curvature_radius_nm: None,
        dihedral_margin_nm: None,
        face_area_margin_nm: None,
        tolerance: tolerance(),
    });
    assert_eq!(missing.mode, ApexSizingMode::RangeAngle);
    assert_eq!(missing.state, ApexSizingState::MissingSolverFact);
}

#[test]
fn sizing_clearance_reports_radius_suggestions_and_apex_collisions() {
    let oversized = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::ConstantRadius,
        requested_extent_nm: Some(30),
        available_extent_nm: Some(20),
        adjacent_clearance_nm: Some(100),
        curvature_radius_nm: None,
        dihedral_margin_nm: None,
        face_area_margin_nm: None,
        tolerance: tolerance(),
    });
    assert_eq!(oversized.state, ApexSizingState::RadiusTooLarge);
    assert_eq!(oversized.max_feasible_radius_nm, Some(20));
    assert_eq!(oversized.margin_nm, Some(-10));

    let close_apex = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::ConstantRadius,
        requested_extent_nm: Some(30),
        available_extent_nm: Some(100),
        adjacent_clearance_nm: Some(20),
        curvature_radius_nm: None,
        dihedral_margin_nm: None,
        face_area_margin_nm: None,
        tolerance: tolerance(),
    });
    assert_eq!(close_apex.state, ApexSizingState::AdjacentApexCollision);
    assert_eq!(close_apex.max_feasible_radius_nm, Some(20));
    assert_eq!(close_apex.margin_nm, Some(70));
}

#[test]
fn sizing_feasibility_score_reports_complete_and_missing_ingredients() {
    let complete = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::ConstantRadius,
        requested_extent_nm: Some(10),
        available_extent_nm: Some(20),
        adjacent_clearance_nm: Some(20),
        curvature_radius_nm: Some(50),
        dihedral_margin_nm: Some(5),
        face_area_margin_nm: Some(7),
        tolerance: tolerance(),
    });
    assert_eq!(complete.state, ApexSizingState::Ready);
    assert_eq!(complete.score.apex_range_margin_nm, Some(10));
    assert_eq!(complete.score.curvature_radius_margin_nm, Some(30));
    assert_eq!(complete.score.dihedral_margin_nm, Some(5));
    assert_eq!(complete.score.face_area_margin_nm, Some(7));
    assert_eq!(complete.score.missing_ingredients(), Vec::new());

    let missing = classify_apex_sizing(ApexSizingInput {
        mode: ApexSizingMode::ConstantRadius,
        requested_extent_nm: Some(10),
        available_extent_nm: None,
        adjacent_clearance_nm: None,
        curvature_radius_nm: None,
        dihedral_margin_nm: None,
        face_area_margin_nm: None,
        tolerance: tolerance(),
    });
    assert_eq!(missing.state, ApexSizingState::MissingSolverFact);
    assert_eq!(
        missing.score.missing_ingredients(),
        vec![
            ApexFeasibilityScoreIngredient::ApexRangeMargin,
            ApexFeasibilityScoreIngredient::DihedralMargin,
            ApexFeasibilityScoreIngredient::CurvatureRadiusMargin,
            ApexFeasibilityScoreIngredient::FaceAreaMargin,
        ]
    );
}

#[test]
fn support_adequacy_separates_inside_tight_from_missing_center_locus() {
    let inside_tight =
        classify_support_adequacy(CornerSupportAdequacyState::InsideTight, None, tolerance());
    assert_eq!(inside_tight.state, CornerSupportAdequacyState::InsideTight);
    assert_eq!(
        inside_tight.center_locus.state,
        CenterLocusFeasibilityState::Penetrates
    );
    assert_eq!(inside_tight.dependency, None);

    let missing_locus = classify_support_adequacy(
        CornerSupportAdequacyState::UnknownInputFact,
        Some(GeometryProbeDependency::CenterLocus),
        tolerance(),
    );
    assert_eq!(
        missing_locus.center_locus.state,
        CenterLocusFeasibilityState::MissingPrimitive("rolling-ball-center-locus")
    );
    assert_eq!(
        missing_locus.dependency,
        Some(GeometryProbeDependency::CenterLocus)
    );
}

#[test]
fn support_footprint_facts_do_not_mask_each_other() {
    let overflow =
        classify_support_adequacy(CornerSupportAdequacyState::Overflow, None, tolerance());
    assert_eq!(
        overflow.footprint.overflow,
        SupportFootprintFeasibilityState::Overflow
    );
    assert_eq!(
        overflow.footprint.cliff,
        SupportFootprintFeasibilityState::Clear
    );

    let cliff = classify_support_adequacy(CornerSupportAdequacyState::Cliff, None, tolerance());
    assert_eq!(
        cliff.footprint.cliff,
        SupportFootprintFeasibilityState::Cliff
    );
    assert_eq!(
        cliff.footprint.overflow,
        SupportFootprintFeasibilityState::Clear
    );
    assert_eq!(
        cliff.footprint.boundary_reachability,
        SupportFootprintFeasibilityState::BoundaryUnreachable
    );

    let collision = classify_support_adequacy(
        CornerSupportAdequacyState::AdjacentBlendCollision,
        None,
        tolerance(),
    );
    assert_eq!(
        collision.footprint.blend_collision,
        SupportFootprintFeasibilityState::BlendTerritoryCollision
    );

    let trim_exit = classify_support_adequacy(
        CornerSupportAdequacyState::TrimDomainExit,
        None,
        tolerance(),
    );
    assert_eq!(
        trim_exit.footprint.boundary_reachability,
        SupportFootprintFeasibilityState::BoundaryUnreachable
    );
}

#[test]
fn support_footprint_reports_missing_projection_primitive() {
    let support = classify_support_adequacy(
        CornerSupportAdequacyState::UnknownInputFact,
        Some(GeometryProbeDependency::Projection),
        tolerance(),
    );
    assert_eq!(
        support.footprint.overflow,
        SupportFootprintFeasibilityState::MissingPrimitive("support-projection")
    );
    assert_eq!(
        support.footprint.boundary_reachability,
        SupportFootprintFeasibilityState::MissingPrimitive("support-projection")
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
    assert_eq!(
        risk.action,
        ConstructionRiskAction::MissingPrimitive("patch-boundary-repair")
    );
    let degeneracy = classify_construction_risk(
        ConstructionRiskKind::IntersectionDegenerate,
        false,
        tolerance(),
    );
    assert_eq!(
        degeneracy.state,
        ConstructionRiskState::ManagedRetry("intersection-refinement")
    );
    assert_eq!(
        degeneracy.action,
        ConstructionRiskAction::Retry("intersection-refinement")
    );
    let curvature =
        classify_construction_risk(ConstructionRiskKind::CurvatureExplosion, false, tolerance());
    assert_eq!(
        curvature.action,
        ConstructionRiskAction::Reclassify("continuity-relaxation")
    );
    let self_intersection =
        classify_construction_risk(ConstructionRiskKind::SelfIntersection, false, tolerance());
    assert_eq!(
        self_intersection.action,
        ConstructionRiskAction::Invalid("self-intersection")
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
    assert_eq!(
        numerical.action,
        ConstructionRiskAction::RollbackRequired("transactional-construction-rollback")
    );

    let lifecycle = classify_lifecycle(false, false, DiagnosticArtifactKind::DiagnosticArtifact);
    assert_eq!(lifecycle.state, LifecycleState::TopologyFingerprintMismatch);
    assert_eq!(
        lifecycle.artifact,
        DiagnosticArtifactKind::DiagnosticArtifact
    );
}
