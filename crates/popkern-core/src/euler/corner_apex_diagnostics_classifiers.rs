//! Pure classifiers for corner/apex diagnostic facts.
#![allow(missing_docs)]

use super::score_types::{ApexFeasibilityScoreFact, ApexFeasibilityScoreMissing};
use super::support_types::CenterLocusFeasibilityFact;
use super::types::*;
use crate::euler::DiagnosticTolerance;

pub fn classify_blend_valence(incident_blend_edges: usize) -> BlendValence {
    match incident_blend_edges {
        0 => BlendValence::NotBlendVertex,
        1 => BlendValence::FreeEnd,
        2 => BlendValence::ThroughVertex,
        3 => BlendValence::StandardCorner,
        _ => BlendValence::HighValenceCorner,
    }
}

pub fn classify_convexity(edges: &[EdgeConvexity]) -> VertexConvexityClass {
    if edges.is_empty() {
        return VertexConvexityClass::NotApplicable;
    }
    if edges.contains(&EdgeConvexity::Unknown) {
        return VertexConvexityClass::UnknownInputFact;
    }
    if edges.contains(&EdgeConvexity::VariableNearVertex) {
        return VertexConvexityClass::VariableConvexityStop;
    }
    if edges.contains(&EdgeConvexity::NearTangent) {
        return VertexConvexityClass::NearTangent;
    }
    let convex = edges.iter().all(|edge| *edge == EdgeConvexity::Convex);
    let concave = edges.iter().all(|edge| *edge == EdgeConvexity::Concave);
    if convex {
        VertexConvexityClass::AllConvex
    } else if concave {
        VertexConvexityClass::AllConcave
    } else {
        VertexConvexityClass::Mixed
    }
}

pub fn classify_blend_type_compatibility(
    selected_edges: &[LocalBlendEdgeSelection],
    blend_valence: BlendValence,
) -> BlendTypeCompatibility {
    if selected_edges.is_empty() {
        return BlendTypeCompatibility::NotApplicable;
    }
    if selected_edges
        .iter()
        .any(|edge| edge.blend_kind == BlendKind::Unknown)
    {
        return BlendTypeCompatibility::UnknownInputFact;
    }
    let first_kind = selected_edges[0].blend_kind;
    let same_kind = selected_edges
        .iter()
        .all(|edge| edge.blend_kind == first_kind);
    let radii = selected_edges
        .iter()
        .map(|edge| edge.radius_nanometers)
        .collect::<std::collections::BTreeSet<_>>();
    if same_kind && radii.len() <= 1 {
        BlendTypeCompatibility::Homogeneous
    } else if same_kind {
        BlendTypeCompatibility::ParameterVarying
    } else if blend_valence == BlendValence::HighValenceCorner
        && selected_edges
            .iter()
            .any(|edge| edge.blend_kind == BlendKind::Setback)
    {
        BlendTypeCompatibility::Incompatible("high-valence-mixed-sizing")
    } else if selected_edges
        .iter()
        .any(|edge| edge.blend_kind == BlendKind::G2Blend)
        && selected_edges
            .iter()
            .any(|edge| edge.blend_kind == BlendKind::Chamfer)
    {
        BlendTypeCompatibility::Incompatible("g2-chamfer-corner")
    } else {
        BlendTypeCompatibility::Heterogeneous
    }
}

pub fn classify_topological_configuration(
    record: &CornerTopologyRecord,
) -> TopologicalConfiguration {
    match record.state {
        CornerTopologyState::Ready if record.incident_blend_edges == 1 => {
            TopologicalConfiguration::FreeEnd
        }
        CornerTopologyState::Ready if record.incident_blend_edges == 2 => {
            TopologicalConfiguration::ThroughVertex
        }
        CornerTopologyState::Ready if record.incident_blend_edges > 3 => {
            TopologicalConfiguration::HighValence
        }
        CornerTopologyState::Ready if record.face_adjacency_edges >= 2 => {
            TopologicalConfiguration::SimpleCorner
        }
        CornerTopologyState::Ready => TopologicalConfiguration::UnknownInputFact,
        CornerTopologyState::NonManifold => TopologicalConfiguration::NonManifold,
        CornerTopologyState::Boundary => TopologicalConfiguration::BoundaryVertex,
        CornerTopologyState::OpenShell => TopologicalConfiguration::OpenShellOrSheet,
        CornerTopologyState::WireOrSheet => TopologicalConfiguration::WireLike,
        CornerTopologyState::StaleHandle | CornerTopologyState::MissingHandle => {
            TopologicalConfiguration::UnknownInputFact
        }
    }
}

pub fn classify_vertex_classification(
    record: &CornerTopologyRecord,
    selected_edges: &[LocalBlendEdgeSelection],
    convexities: &[EdgeConvexity],
    exposure: CornerSupportAdequacyState,
) -> VertexClassificationFact {
    let valence = classify_blend_valence(record.incident_blend_edges);
    VertexClassificationFact {
        valence,
        convexity: classify_convexity(convexities),
        blend_types: classify_blend_type_compatibility(selected_edges, valence),
        topology: classify_topological_configuration(record),
        exposure,
    }
}

pub fn classify_apex_sizing(input: ApexSizingInput) -> ApexSizingFact {
    if input.mode == ApexSizingMode::NotApplicable {
        return sizing_fact(input, ApexSizingState::NotApplicable, None, None);
    }
    let Some(requested) = input.requested_extent_nm else {
        return sizing_fact(input, ApexSizingState::MissingSolverFact, None, None);
    };
    if requested <= 0 {
        return sizing_fact(input, ApexSizingState::ZeroOrNegativeRange, None, None);
    }
    let Some(available) = input.available_extent_nm else {
        return sizing_fact(input, ApexSizingState::MissingSolverFact, None, None);
    };
    let margin = available - requested;
    if margin < 0 {
        return sizing_fact(
            input,
            ApexSizingState::RadiusTooLarge,
            Some(available),
            Some(margin),
        );
    }
    if input
        .adjacent_clearance_nm
        .is_some_and(|clearance| requested > clearance)
    {
        return sizing_fact(
            input,
            ApexSizingState::AdjacentApexCollision,
            input.adjacent_clearance_nm,
            Some(margin),
        );
    }
    if input
        .curvature_radius_nm
        .is_some_and(|radius| radius > 0 && requested * 2 > radius)
    {
        return sizing_fact(
            input,
            ApexSizingState::HighCurvatureRisk,
            Some(requested),
            Some(margin),
        );
    }
    sizing_fact(input, ApexSizingState::Ready, Some(requested), Some(margin))
}

fn sizing_fact(
    input: ApexSizingInput,
    state: ApexSizingState,
    max_feasible_radius_nm: Option<i64>,
    margin_nm: Option<i64>,
) -> ApexSizingFact {
    let score = apex_feasibility_score(&input);
    ApexSizingFact {
        mode: input.mode,
        state,
        max_feasible_radius_nm,
        margin_nm,
        score,
        tolerance: input.tolerance,
    }
}

fn apex_feasibility_score(input: &ApexSizingInput) -> ApexFeasibilityScoreFact {
    let apex_range_margin_nm = input
        .requested_extent_nm
        .zip(input.available_extent_nm)
        .map(|(requested, available)| available - requested);
    let curvature_radius_margin_nm = input
        .requested_extent_nm
        .zip(input.curvature_radius_nm)
        .map(|(requested, curvature_radius)| curvature_radius - requested.saturating_mul(2));
    ApexFeasibilityScoreFact {
        apex_range_margin_nm,
        dihedral_margin_nm: input.dihedral_margin_nm,
        curvature_radius_margin_nm,
        face_area_margin_nm: input.face_area_margin_nm,
        missing: ApexFeasibilityScoreMissing {
            apex_range_margin: apex_range_margin_nm.is_none(),
            dihedral_margin: input.dihedral_margin_nm.is_none(),
            curvature_radius_margin: curvature_radius_margin_nm.is_none(),
            face_area_margin: input.face_area_margin_nm.is_none(),
        },
    }
}

pub fn classify_trial_interaction(
    valence: BlendValence,
    has_intersection_primitive: bool,
    residual_nm: Option<i64>,
    tolerance: DiagnosticTolerance,
) -> TrialInteractionFact {
    if valence != BlendValence::ThroughVertex {
        return TrialInteractionFact {
            state: TrialInteractionState::NotApplicable,
            dependency: None,
            tolerance,
        };
    }
    if !has_intersection_primitive {
        return TrialInteractionFact {
            state: TrialInteractionState::MissingIntersectionPrimitive,
            dependency: Some(GeometryProbeDependency::Intersection),
            tolerance,
        };
    }
    match residual_nm {
        None => TrialInteractionFact {
            state: TrialInteractionState::Ambiguous,
            dependency: Some(GeometryProbeDependency::Intersection),
            tolerance,
        },
        Some(residual) if tolerance.nanometers.is_some_and(|limit| residual <= limit) => {
            TrialInteractionFact {
                state: TrialInteractionState::WellConditioned,
                dependency: None,
                tolerance,
            }
        }
        Some(_) => TrialInteractionFact {
            state: TrialInteractionState::NearTangentUpgrade,
            dependency: Some(GeometryProbeDependency::Intersection),
            tolerance,
        },
    }
}

pub fn classify_support_adequacy(
    footprint: CornerSupportAdequacyState,
    dependency: Option<GeometryProbeDependency>,
    tolerance: DiagnosticTolerance,
) -> CornerSupportAdequacyFact {
    let center_locus = center_locus_feasibility(footprint, dependency);
    CornerSupportAdequacyFact {
        state: footprint,
        center_locus,
        dependency,
        tolerance,
    }
}

fn center_locus_feasibility(
    footprint: CornerSupportAdequacyState,
    dependency: Option<GeometryProbeDependency>,
) -> CenterLocusFeasibilityFact {
    if dependency == Some(GeometryProbeDependency::CenterLocus) {
        CenterLocusFeasibilityFact::missing_primitive("rolling-ball-center-locus")
    } else if footprint == CornerSupportAdequacyState::InsideTight {
        CenterLocusFeasibilityFact::penetrates()
    } else {
        CenterLocusFeasibilityFact::clear()
    }
}

pub fn classify_strategy_readiness(
    classification: VertexClassificationFact,
    sizing: ApexSizingFact,
    support: CornerSupportAdequacyFact,
) -> StrategyReadinessFact {
    let (family, readiness) = if let BlendTypeCompatibility::Incompatible(reason) =
        classification.blend_types
    {
        (
            CornerApexStrategyFamily::VertexPatch,
            CornerApexStrategyReadiness::Invalid(reason),
        )
    } else if classification.blend_types == BlendTypeCompatibility::Heterogeneous {
        (
            CornerApexStrategyFamily::VertexPatch,
            CornerApexStrategyReadiness::RequiresMissingPrimitive("mixed-blend-family-meshing"),
        )
    } else if classification.topology == TopologicalConfiguration::NonManifold {
        (
            CornerApexStrategyFamily::VertexPatch,
            CornerApexStrategyReadiness::Invalid("non-manifold-corner"),
        )
    } else if matches!(
        sizing.state,
        ApexSizingState::RadiusTooLarge
            | ApexSizingState::AdjacentApexCollision
            | ApexSizingState::ZeroOrNegativeRange
    ) {
        (
            CornerApexStrategyFamily::Setback,
            CornerApexStrategyReadiness::RequiresMissingPrimitive("radius-or-setback-adjustment"),
        )
    } else if support.state == CornerSupportAdequacyState::Cliff {
        (
            CornerApexStrategyFamily::CliffTermination,
            CornerApexStrategyReadiness::ReadyWithManagedRisk,
        )
    } else {
        match classification.valence {
            BlendValence::FreeEnd => (
                CornerApexStrategyFamily::CapRunout,
                CornerApexStrategyReadiness::Ready,
            ),
            BlendValence::ThroughVertex => (
                CornerApexStrategyFamily::MutualTrim,
                CornerApexStrategyReadiness::Ready,
            ),
            BlendValence::StandardCorner => (
                CornerApexStrategyFamily::ThreeFaceBlend,
                CornerApexStrategyReadiness::Ready,
            ),
            BlendValence::HighValenceCorner => (
                CornerApexStrategyFamily::HighValencePatch,
                CornerApexStrategyReadiness::RequiresMissingPrimitive("n-sided-corner-patch"),
            ),
            BlendValence::NotBlendVertex => (
                CornerApexStrategyFamily::VertexPatch,
                CornerApexStrategyReadiness::UnknownInputFact,
            ),
        }
    };
    StrategyReadinessFact {
        family,
        readiness,
        continuity_target: "G1",
        deferred: family == CornerApexStrategyFamily::DeferredVertexPatch,
        fallback_eligible: matches!(
            family,
            CornerApexStrategyFamily::VertexPatch | CornerApexStrategyFamily::HighValencePatch
        ),
    }
}

pub fn classify_construction_risk(
    kind: ConstructionRiskKind,
    primitive_available: bool,
    tolerance: DiagnosticTolerance,
) -> ConstructionRiskFact {
    let state = if primitive_available {
        ConstructionRiskState::NoRisk
    } else {
        ConstructionRiskState::RequiresMissingPrimitive(match kind {
            ConstructionRiskKind::IntersectionDegenerate => "intersection-refinement",
            ConstructionRiskKind::TangencyUnsolved => "three-face-tangency-solver",
            ConstructionRiskKind::BoundaryCurveMissing => "boundary-curve-extraction",
            ConstructionRiskKind::SetbackSplitInfeasible => "setback-split-construction",
            ConstructionRiskKind::BoundaryConditionsIncomplete => "boundary-condition-assembly",
            ConstructionRiskKind::TwistCompatibility => "twist-compatibility-solver",
            ConstructionRiskKind::PatchBoundaryIncompatible => "patch-boundary-repair",
            ConstructionRiskKind::CurvatureExplosion => "continuity-relaxation",
            ConstructionRiskKind::SelfIntersection => "self-intersection-repair",
            ConstructionRiskKind::VariableRadiusCollapse => "variable-radius-collapse-handler",
            ConstructionRiskKind::ChamferRoundMeshing => "mixed-chamfer-round-meshing",
            ConstructionRiskKind::FacetedSmoothMismatch => "geometry-mode-conversion",
            ConstructionRiskKind::ApproximationResidual => "approximation-refinement",
            ConstructionRiskKind::NumericalSolverFailure => "numerical-solver-failure",
        })
    };
    ConstructionRiskFact {
        kind,
        state,
        tolerance,
    }
}

pub fn classify_integration_validation(
    kind: IntegrationValidationKind,
    passes: Option<bool>,
    tolerance: DiagnosticTolerance,
) -> IntegrationValidationFact {
    let state = match passes {
        Some(true) => IntegrationValidationState::Ready,
        Some(false) => IntegrationValidationState::RepairRequired(match kind {
            IntegrationValidationKind::Continuity => "continuity-repair-or-downgrade",
            IntegrationValidationKind::MergeReadiness => "merge-refusal-preserves-history",
            _ => "integration-repair-required",
        }),
        None => IntegrationValidationState::UnknownInputFact,
    };
    IntegrationValidationFact {
        kind,
        state,
        tolerance,
    }
}

pub fn classify_lifecycle(
    fingerprint_matches: bool,
    reused_after_mutation: bool,
    artifact: DiagnosticArtifactKind,
) -> LifecycleFact {
    let state = if reused_after_mutation {
        LifecycleState::ReusedAfterMutation
    } else if fingerprint_matches {
        LifecycleState::ScopedDerivedDiscardable
    } else {
        LifecycleState::TopologyFingerprintMismatch
    };
    LifecycleFact { state, artifact }
}
