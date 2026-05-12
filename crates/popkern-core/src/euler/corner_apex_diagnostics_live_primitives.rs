//! Live read-only geometry primitives for corner/apex diagnostics.

use super::classifiers::classify_trial_interaction;
use super::live_helpers::*;
use super::three_face_tangency::{
    solve_three_face_tangency, tangency_supports, three_face_tangency_support_missing_primitive,
};
use super::types::*;
use crate::euler::{DiagnosticHandle, DiagnosticTolerance};
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

pub(super) struct LivePrimitiveAssessment {
    pub trial_interaction: TrialInteractionFact,
    pub strategy: StrategyReadinessFact,
}

pub(super) fn assess_live_primitives(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    classification: VertexClassificationFact,
    base_strategy: StrategyReadinessFact,
    tolerance: DiagnosticTolerance,
) -> LivePrimitiveAssessment {
    match base_strategy.family {
        CornerApexStrategyFamily::MutualTrim => {
            assess_mutual_trim_probe(brep, request, base_strategy, tolerance)
        }
        CornerApexStrategyFamily::ThreeFaceBlend => {
            assess_three_face_tangency(brep, request, base_strategy, tolerance)
        }
        _ => LivePrimitiveAssessment {
            trial_interaction: classify_trial_interaction(
                classification.valence,
                false,
                None,
                tolerance,
            ),
            strategy: base_strategy,
        },
    }
}

fn assess_mutual_trim_probe(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    mut strategy: StrategyReadinessFact,
    tolerance: DiagnosticTolerance,
) -> LivePrimitiveAssessment {
    let selected = live_selected_edges(brep, request);
    let missing = mutual_trim_missing_reason(brep, request, &selected);
    let trial_interaction = if missing.is_none() && selected.len() == 2 {
        let angle = edge_tangent_angle_at_vertex(brep, request.scope.vertex, selected);
        match angle {
            Some(angle) if angle.sin().abs() > angle_tolerance(tolerance) => {
                classify_trial_interaction(BlendValence::ThroughVertex, true, Some(0), tolerance)
            }
            Some(_) => TrialInteractionFact {
                state: TrialInteractionState::NearTangentUpgrade,
                dependency: Some(GeometryProbeDependency::Intersection),
                tolerance,
            },
            None => classify_trial_interaction(BlendValence::ThroughVertex, false, None, tolerance),
        }
    } else {
        classify_trial_interaction(BlendValence::ThroughVertex, false, None, tolerance)
    };
    strategy.readiness = match (missing, trial_interaction.state) {
        (Some(reason), _) => CornerApexStrategyReadiness::RequiresMissingPrimitive(reason),
        (_, TrialInteractionState::WellConditioned) => CornerApexStrategyReadiness::Ready,
        (_, TrialInteractionState::NearTangentUpgrade) => {
            CornerApexStrategyReadiness::RequiresMissingPrimitive("near-tangent-mutual-trim")
        }
        (_, TrialInteractionState::SelfIntersecting) => {
            CornerApexStrategyReadiness::Invalid("coincident-mutual-trim-intersection")
        }
        (_, TrialInteractionState::Ambiguous) => {
            CornerApexStrategyReadiness::RequiresMissingPrimitive("ambiguous-mutual-trim-branch")
        }
        _ => CornerApexStrategyReadiness::RequiresMissingPrimitive(
            "missing-mutual-trim-intersection-refinement",
        ),
    };
    LivePrimitiveAssessment {
        trial_interaction,
        strategy,
    }
}

fn assess_three_face_tangency(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    mut strategy: StrategyReadinessFact,
    tolerance: DiagnosticTolerance,
) -> LivePrimitiveAssessment {
    let selected = live_selected_edges(brep, request);
    strategy.readiness = match three_face_tangency_missing_reason(brep, request, &selected) {
        Some(reason) => CornerApexStrategyReadiness::RequiresMissingPrimitive(reason),
        None if solve_live_three_face_tangency(brep, request, tolerance).is_some() => {
            CornerApexStrategyReadiness::Ready
        }
        None if tolerance.nanometers == Some(0) => {
            CornerApexStrategyReadiness::Invalid("three-face-tangency-residual-limit-failure")
        }
        None => {
            CornerApexStrategyReadiness::RequiresMissingPrimitive("three-face-tangency-no-solution")
        }
    };
    LivePrimitiveAssessment {
        trial_interaction: classify_trial_interaction(
            BlendValence::StandardCorner,
            false,
            None,
            tolerance,
        ),
        strategy,
    }
}

fn mutual_trim_missing_reason(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    selected: &[EdgeHandle],
) -> Option<&'static str> {
    if selected.len() != 2 {
        return Some("blend-pipe-intersection-selection");
    }
    if !selected_radii_equal(request) {
        return Some("variable-radius-mutual-trim-probe");
    }
    if selected.iter().any(|edge| !is_line_edge(brep, *edge)) {
        return Some("nonlinear-spine-mutual-trim-probe");
    }
    if selected
        .iter()
        .any(|edge| !edge_has_two_planar_supports(brep, *edge))
    {
        return Some("nonplanar-support-mutual-trim-probe");
    }
    None
}

fn three_face_tangency_missing_reason(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    selected: &[EdgeHandle],
) -> Option<&'static str> {
    if selected.len() != 3 {
        return Some("three-face-tangency-selection");
    }
    if !selected_radii_equal(request) {
        return Some("variable-radius-three-face-tangency-solver");
    }
    let Some(vertex) = live_vertex(request.scope.vertex, brep) else {
        return Some("three-face-tangency-stale-vertex");
    };
    let faces = incident_faces_at_vertex(brep, vertex);
    if faces.len() < 3 {
        return Some("three-face-tangency-support-count");
    }
    if let Some(reason) = faces
        .iter()
        .filter_map(|face| brep.faces[*face].surface.as_ref())
        .map(three_face_tangency_support_missing_primitive)
        .find(|reason| *reason != "three-face-tangency-analytic-support-ready")
    {
        return Some(reason);
    }
    let vertex_point = brep.vertices[vertex].position;
    if faces
        .iter()
        .any(|face| has_singular_tangency_normal(brep, *face, vertex_point))
    {
        return Some("three-face-tangency-singular-normal");
    }
    if has_multiple_solution_ambiguity(brep, &faces) {
        return Some("three-face-tangency-multiple-solution-ambiguity");
    }
    if let Err(reason) = tangency_supports(brep, vertex) {
        return Some(reason);
    }
    None
}

fn has_singular_tangency_normal(brep: &BRep, face: FaceHandle, vertex: Point3) -> bool {
    let Some(surface) = brep.faces[face].surface.as_ref() else {
        return false;
    };
    let near_zero = |vector: Vec3| vector.length() <= f64::EPSILON;
    match surface.geometry {
        SurfaceGeometry::Sphere { center, .. } => near_zero(vertex - center),
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let Some(axis) = axis.normalized() else {
                return true;
            };
            let from_origin = vertex - origin;
            near_zero(from_origin - axis * from_origin.dot(&axis))
        }
        SurfaceGeometry::Cone { apex, .. } => near_zero(vertex - apex),
        SurfaceGeometry::Plane { normal, .. } => normal.normalized().is_none(),
        _ => false,
    }
}

fn has_multiple_solution_ambiguity(brep: &BRep, faces: &[FaceHandle]) -> bool {
    faces.iter().take(3).all(|face| {
        brep.faces[*face]
            .surface
            .as_ref()
            .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Sphere { .. }))
    })
}

fn solve_live_three_face_tangency(
    brep: &BRep,
    request: &CornerApexDiagnosticRequest,
    tolerance: DiagnosticTolerance,
) -> Option<Point3> {
    let vertex = live_vertex(request.scope.vertex, brep)?;
    let radius = selected_radius_meters(request)?;
    let supports = tangency_supports(brep, vertex).ok()?;
    solve_three_face_tangency(&supports, brep.vertices[vertex].position, radius, tolerance)
}

fn live_selected_edges(brep: &BRep, request: &CornerApexDiagnosticRequest) -> Vec<EdgeHandle> {
    let Some(vertex) = live_vertex(request.scope.vertex, brep) else {
        return Vec::new();
    };
    request
        .selected_edges
        .iter()
        .filter_map(|selection| live_edge(selection.edge, brep))
        .filter(|edge| edge_contains_vertex(brep, *edge, vertex))
        .collect()
}

fn selected_radii_equal(request: &CornerApexDiagnosticRequest) -> bool {
    let mut radii = request
        .selected_edges
        .iter()
        .map(|selection| selection.radius_nanometers);
    let Some(first) = radii.next() else {
        return false;
    };
    first.is_some() && radii.all(|radius| radius == first)
}

fn selected_radius_meters(request: &CornerApexDiagnosticRequest) -> Option<f64> {
    request
        .selected_edges
        .first()
        .and_then(|selection| selection.radius_nanometers)
        .map(|radius| radius as f64 * 1e-9)
        .filter(|radius| radius.is_finite() && *radius > 0.0)
}

fn is_line_edge(brep: &BRep, edge: EdgeHandle) -> bool {
    brep.edges
        .get(edge)
        .and_then(|edge| edge.curve.as_ref())
        .is_some_and(|curve| matches!(curve.geometry, CurveGeometry::Line { .. }))
}

fn edge_has_two_planar_supports(brep: &BRep, edge: EdgeHandle) -> bool {
    let faces = edge_faces(brep, edge);
    faces.len() >= 2 && faces.iter().all(|face| plane_record(brep, *face).is_some())
}

fn edge_tangent_angle_at_vertex(
    brep: &BRep,
    vertex: DiagnosticHandle,
    edges: Vec<EdgeHandle>,
) -> Option<f64> {
    let vertex = live_vertex(vertex, brep)?;
    let first = edge_tangent_at_vertex(brep, edges[0], vertex)?;
    let second = edge_tangent_at_vertex(brep, edges[1], vertex)?;
    Some(first.angle_to(&second))
}

fn edge_tangent_at_vertex(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> Option<Vec3> {
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let parameter = curve
        .geometry
        .closest_point(&brep.vertices.get(vertex)?.position)?;
    curve.geometry.tangent_at(parameter).normalized()
}

#[derive(Clone, Copy)]
struct PlaneRecord;

fn plane_record(brep: &BRep, face: FaceHandle) -> Option<PlaneRecord> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    match surface {
        SurfaceGeometry::Plane { normal, .. } => normal.normalized().map(|_| PlaneRecord),
        _ => None,
    }
}

fn angle_tolerance(tolerance: DiagnosticTolerance) -> f64 {
    tolerance_meters(tolerance).max(1e-9)
}

fn tolerance_meters(tolerance: DiagnosticTolerance) -> f64 {
    tolerance.nanometers.unwrap_or(10).max(1) as f64 * 1e-9
}
