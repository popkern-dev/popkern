//! Live BRep extraction for corner/apex diagnostics.

use super::classifiers::{
    classify_apex_sizing, classify_construction_risk, classify_integration_validation,
    classify_lifecycle, classify_strategy_readiness, classify_support_adequacy,
    classify_vertex_classification,
};
use super::live_helpers::*;
use super::live_primitives::assess_live_primitives;
use super::topology_neighborhood::{extract_topology_neighborhood, stale_topology_neighborhood};
use super::types::*;
use crate::euler::{DiagnosticHandle, DiagnosticTolerance};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, TopologyRole, VertexHandle};

impl BRep {
    /// Build a read-only corner/apex diagnostic report from live BRep topology.
    pub fn diagnose_corner_apex(
        &self,
        request: CornerApexDiagnosticRequest,
    ) -> CornerApexDiagnosticReport {
        let tolerance = local_tolerance(self, request.scope.vertex);
        let topology = extract_topology(self, &request, tolerance);
        let support = classify_support_adequacy(
            support_state(self, &topology),
            support_dependency(&topology),
            tolerance,
        );
        let convexities = edge_convexities(self, &request, topology.vertex, tolerance);
        let classification = classify_vertex_classification(
            &topology,
            &request.selected_edges,
            &convexities,
            support.state,
        );
        let sizing = classify_apex_sizing(sizing_input(self, &request, topology.vertex, tolerance));
        let primitive_assessment = assess_live_primitives(
            self,
            &request,
            classification,
            classify_strategy_readiness(classification, sizing, support),
            tolerance,
        );
        let trial_interaction = primitive_assessment.trial_interaction;
        let strategy = primitive_assessment.strategy;
        let construction_risk = construction_risk(&classification, &strategy, tolerance);
        let integration_validation = integration_validation(&strategy, tolerance);
        let lifecycle = classify_lifecycle(
            topology.state == CornerTopologyState::Ready,
            false,
            DiagnosticArtifactKind::AcceptedGraphFact,
        );

        CornerApexDiagnosticReport {
            request,
            topology,
            classification,
            sizing,
            trial_interaction,
            support_adequacy: support,
            strategy,
            construction_risk,
            integration_validation,
            lifecycle,
        }
    }
}

fn extract_topology(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    tolerance: DiagnosticTolerance,
) -> CornerTopologyRecord {
    let DiagnosticHandle::Vertex(vertex) = request.scope.vertex else {
        return missing_topology(
            request.scope.vertex,
            CornerTopologyState::MissingHandle,
            tolerance,
        );
    };
    if !brep.vertices.contains_key(vertex) {
        return missing_topology(
            request.scope.vertex,
            CornerTopologyState::StaleHandle,
            tolerance,
        );
    }

    let selected_edge_stale = request
        .selected_edges
        .iter()
        .any(|selection| live_edge(selection.edge, brep).is_none());
    let incident_edges = incident_edges_at_vertex(brep, vertex);
    let incident_faces = incident_faces_at_vertex(brep, vertex);
    let incident_blend_edges = request
        .selected_edges
        .iter()
        .filter_map(|selection| live_edge(selection.edge, brep))
        .filter(|edge| incident_edges.contains(edge))
        .count();
    let incident_non_blend_edges = incident_edges.len().saturating_sub(incident_blend_edges);
    let state = if selected_edge_stale {
        CornerTopologyState::StaleHandle
    } else {
        topology_state(brep, vertex, &incident_edges)
    };
    let geometry = geometry_fact(brep, &incident_edges, &incident_faces, tolerance);
    let surface_representation = surface_representation(brep, &incident_faces);
    let surface_compatibility = surface_compatibility(brep, &incident_faces);
    let mut neighborhood =
        extract_topology_neighborhood(brep, vertex, &incident_edges, &request.selected_edges);
    if selected_edge_stale {
        neighborhood
            .stale_sources
            .push(CornerTopologyFactSource::SelectedEdge);
    }

    CornerTopologyRecord {
        vertex: DiagnosticHandle::Vertex(vertex),
        incident_blend_edges,
        incident_non_blend_edges,
        surrounding_faces: incident_faces.len(),
        face_adjacency_edges: face_adjacency_edges(brep, &incident_edges),
        neighborhood,
        state,
        geometry,
        surface_representation,
        surface_compatibility,
    }
}

fn missing_topology(
    vertex: DiagnosticHandle,
    state: CornerTopologyState,
    tolerance: DiagnosticTolerance,
) -> CornerTopologyRecord {
    CornerTopologyRecord {
        vertex,
        incident_blend_edges: 0,
        incident_non_blend_edges: 0,
        surrounding_faces: 0,
        face_adjacency_edges: 0,
        neighborhood: stale_topology_neighborhood(CornerTopologyFactSource::Vertex),
        state,
        geometry: CornerGeometryFact {
            tangent: CornerGeometryAvailability::MissingInputFact,
            normal: CornerGeometryAvailability::MissingInputFact,
            dihedral: CornerGeometryAvailability::MissingInputFact,
            tolerance,
        },
        surface_representation: SurfaceRepresentation::UnknownImported,
        surface_compatibility: SurfaceCompatibility::Unknown("missing-corner-topology"),
    }
}

fn topology_state(
    brep: &BRep,
    vertex: VertexHandle,
    incident_edges: &[EdgeHandle],
) -> CornerTopologyState {
    if incident_edges.is_empty() {
        return CornerTopologyState::WireOrSheet;
    }
    if incident_edges.iter().any(|edge| {
        let half_edge = brep.edges[*edge].half_edge;
        brep.nonmanifold_radial.contains_key(&half_edge)
            || brep
                .half_edges
                .get(half_edge)
                .is_some_and(|he| brep.nonmanifold_radial.contains_key(&he.twin))
    }) {
        return CornerTopologyState::NonManifold;
    }
    if incident_edges
        .iter()
        .any(|edge| brep.boundary_edges.contains(edge))
    {
        return CornerTopologyState::Boundary;
    }
    if brep.vertices[vertex].outgoing == HalfEdgeHandle::default() {
        return CornerTopologyState::WireOrSheet;
    }
    CornerTopologyState::Ready
}

fn geometry_fact(
    brep: &BRep,
    incident_edges: &[EdgeHandle],
    incident_faces: &[FaceHandle],
    tolerance: DiagnosticTolerance,
) -> CornerGeometryFact {
    let tangent = if incident_edges.is_empty() {
        CornerGeometryAvailability::NotApplicable
    } else if incident_edges
        .iter()
        .all(|edge| brep.edges[*edge].curve.is_some())
    {
        CornerGeometryAvailability::Available
    } else {
        CornerGeometryAvailability::BlockedByDependency(GeometryProbeDependency::CurveTangent)
    };
    let normal = if incident_faces.is_empty() {
        CornerGeometryAvailability::NotApplicable
    } else if incident_faces
        .iter()
        .all(|face| brep.faces[*face].surface.is_some())
    {
        CornerGeometryAvailability::Available
    } else {
        CornerGeometryAvailability::BlockedByDependency(GeometryProbeDependency::SurfaceNormal)
    };
    let dihedral = if incident_faces.len() < 2 {
        CornerGeometryAvailability::NotApplicable
    } else if tangent == CornerGeometryAvailability::Available
        && normal == CornerGeometryAvailability::Available
    {
        CornerGeometryAvailability::Available
    } else {
        CornerGeometryAvailability::BlockedByDependency(GeometryProbeDependency::DihedralAngle)
    };
    CornerGeometryFact {
        tangent,
        normal,
        dihedral,
        tolerance,
    }
}

fn surface_representation(brep: &BRep, faces: &[FaceHandle]) -> SurfaceRepresentation {
    let mut representation = None;
    for face in faces {
        let Some(surface) = &brep.faces[*face].surface else {
            return SurfaceRepresentation::UnknownImported;
        };
        let current = match surface.geometry {
            SurfaceGeometry::Plane { .. } => SurfaceRepresentation::Plane,
            SurfaceGeometry::Cylinder { .. } => SurfaceRepresentation::Cylinder,
            SurfaceGeometry::Cone { .. } => SurfaceRepresentation::Cone,
            SurfaceGeometry::Sphere { .. } => SurfaceRepresentation::Sphere,
            SurfaceGeometry::Torus { .. } => SurfaceRepresentation::Torus,
            SurfaceGeometry::Nurbs(_) => SurfaceRepresentation::Spline,
            SurfaceGeometry::Gregory(_) => SurfaceRepresentation::ExistingBlend,
        };
        if representation.is_some_and(|previous| previous != current) {
            return SurfaceRepresentation::SweptOrRevolved;
        }
        representation = Some(current);
    }
    representation.unwrap_or(SurfaceRepresentation::UnknownImported)
}

fn surface_compatibility(brep: &BRep, faces: &[FaceHandle]) -> SurfaceCompatibility {
    if faces.is_empty() {
        return SurfaceCompatibility::Unknown("no-incident-support-faces");
    }
    let mut has_existing_blend = false;
    for face in faces {
        let Some(surface) = &brep.faces[*face].surface else {
            return SurfaceCompatibility::Unknown("missing-support-surface");
        };
        if brep.faces[*face].role == TopologyRole::ConstructionResidual {
            return SurfaceCompatibility::Incompatible("construction-residual-support");
        }
        has_existing_blend |= matches!(surface.geometry, SurfaceGeometry::Gregory(_));
    }
    if has_existing_blend {
        SurfaceCompatibility::Conditional("nested-blend-support")
    } else {
        SurfaceCompatibility::Full
    }
}

fn sizing_input(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    vertex: DiagnosticHandle,
    tolerance: DiagnosticTolerance,
) -> ApexSizingInput {
    let requested = request
        .selected_edges
        .iter()
        .filter_map(|selection| selection.radius_nanometers)
        .max();
    let available = live_vertex(vertex, brep).and_then(|v| {
        request
            .selected_edges
            .iter()
            .filter_map(|selection| live_edge(selection.edge, brep))
            .filter(|edge| edge_contains_vertex(brep, *edge, v))
            .filter_map(|edge| edge_extent_nm(brep, edge, v))
            .min()
    });
    let mode = sizing_mode(&request.selected_edges);
    ApexSizingInput {
        mode,
        requested_extent_nm: requested,
        available_extent_nm: available,
        adjacent_clearance_nm: available.map(|extent| extent / 2),
        curvature_radius_nm: None,
        tolerance,
    }
}

fn sizing_mode(selected_edges: &[LocalBlendEdgeSelection]) -> ApexSizingMode {
    if selected_edges.is_empty() {
        return ApexSizingMode::NotApplicable;
    }
    if selected_edges
        .iter()
        .any(|edge| edge.blend_kind == BlendKind::FilletVariableRadius)
    {
        ApexSizingMode::VariableRadiusEndpoint
    } else {
        ApexSizingMode::ConstantRadius
    }
}

fn support_state(brep: &BRep, topology: &CornerTopologyRecord) -> CornerSupportAdequacyState {
    if topology.state == CornerTopologyState::Boundary {
        return CornerSupportAdequacyState::Cliff;
    }
    if topology.state != CornerTopologyState::Ready {
        return CornerSupportAdequacyState::UnknownInputFact;
    }
    if topology.geometry.normal != CornerGeometryAvailability::Available {
        return CornerSupportAdequacyState::UnknownInputFact;
    }
    let DiagnosticHandle::Vertex(vertex) = topology.vertex else {
        return CornerSupportAdequacyState::UnknownInputFact;
    };
    if incident_faces_at_vertex(brep, vertex)
        .iter()
        .any(|face| !brep.faces[*face].inner_loops.is_empty())
    {
        return CornerSupportAdequacyState::HoleOrInnerLoopCollision;
    }
    CornerSupportAdequacyState::Fits
}

fn support_dependency(topology: &CornerTopologyRecord) -> Option<GeometryProbeDependency> {
    if topology.geometry.normal == CornerGeometryAvailability::Available {
        None
    } else {
        Some(GeometryProbeDependency::SupportFootprint)
    }
}

fn edge_convexities(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    vertex: DiagnosticHandle,
    tolerance: DiagnosticTolerance,
) -> Vec<EdgeConvexity> {
    let Some(vertex) = live_vertex(vertex, brep) else {
        return vec![EdgeConvexity::Unknown];
    };
    request
        .selected_edges
        .iter()
        .filter_map(|selection| live_edge(selection.edge, brep))
        .filter(|edge| edge_contains_vertex(brep, *edge, vertex))
        .map(|edge| edge_convexity(brep, edge, vertex, tolerance))
        .collect()
}

fn edge_convexity(
    brep: &BRep,
    edge: EdgeHandle,
    vertex: VertexHandle,
    tolerance: DiagnosticTolerance,
) -> EdgeConvexity {
    let faces = edge_faces(brep, edge);
    if faces.len() < 2 {
        return EdgeConvexity::Unknown;
    }
    let Some(first) = face_normal_at_vertex(brep, faces[0], vertex) else {
        return EdgeConvexity::Unknown;
    };
    let Some(second) = face_normal_at_vertex(brep, faces[1], vertex) else {
        return EdgeConvexity::Unknown;
    };
    let Some(first) = first.normalized() else {
        return EdgeConvexity::Unknown;
    };
    let Some(second) = second.normalized() else {
        return EdgeConvexity::Unknown;
    };
    let near_tangent_limit = tolerance.nanometers.unwrap_or(10).max(1) as f64 * 1e-9;
    if (1.0 - first.dot(&second).abs()) <= near_tangent_limit.max(1e-9) {
        EdgeConvexity::NearTangent
    } else {
        EdgeConvexity::Convex
    }
}

fn face_normal_at_vertex(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Option<crate::geometry::Vec3> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let position = brep.vertices.get(vertex)?.position;
    let (u, v) = project_point_to_uv(surface, &position)?;
    Some(surface.normal_at(u, v))
}

fn construction_risk(
    classification: &VertexClassificationFact,
    strategy: &StrategyReadinessFact,
    tolerance: DiagnosticTolerance,
) -> ConstructionRiskFact {
    let kind = match strategy.family {
        CornerApexStrategyFamily::MutualTrim => ConstructionRiskKind::IntersectionDegenerate,
        CornerApexStrategyFamily::ThreeFaceBlend => ConstructionRiskKind::TangencyUnsolved,
        CornerApexStrategyFamily::HighValencePatch | CornerApexStrategyFamily::VertexPatch => {
            ConstructionRiskKind::BoundaryConditionsIncomplete
        }
        CornerApexStrategyFamily::Setback => ConstructionRiskKind::SetbackSplitInfeasible,
        CornerApexStrategyFamily::CliffTermination
        | CornerApexStrategyFamily::OverflowSmoothing => {
            ConstructionRiskKind::ApproximationResidual
        }
        CornerApexStrategyFamily::CapRunout | CornerApexStrategyFamily::DeferredVertexPatch => {
            ConstructionRiskKind::BoundaryCurveMissing
        }
    };
    let primitive_available = matches!(
        strategy.readiness,
        CornerApexStrategyReadiness::Ready | CornerApexStrategyReadiness::ReadyWithManagedRisk
    ) && !matches!(
        classification.blend_types,
        BlendTypeCompatibility::Heterogeneous | BlendTypeCompatibility::Incompatible(_)
    );
    classify_construction_risk(kind, primitive_available, tolerance)
}

fn integration_validation(
    strategy: &StrategyReadinessFact,
    tolerance: DiagnosticTolerance,
) -> IntegrationValidationFact {
    let passes = matches!(
        strategy.readiness,
        CornerApexStrategyReadiness::Ready | CornerApexStrategyReadiness::ReadyWithManagedRisk
    );
    classify_integration_validation(
        IntegrationValidationKind::Watertightness,
        Some(passes),
        tolerance,
    )
}
