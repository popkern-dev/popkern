//! Trim-endpoint candidate generation for edge/surface vertex-star repair.

use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, SurfaceGeometry};
use crate::topology::brep_edge_surface_compatibility::BrepEdgeSurfaceCompatibilityOptions;
use crate::topology::brep_edge_surface_compatibility_vertex_candidate::BrepEdgeSurfaceVertexCandidateOutcome;
use crate::topology::{BRep, HalfEdgeHandle};

pub(super) struct TrimEndpointCandidate {
    pub point: Option<Point3>,
    pub outcome: BrepEdgeSurfaceVertexCandidateOutcome,
    pub anchor_motion: f64,
}

pub(super) fn trim_endpoint_candidates_for_side(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    current: Point3,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> Vec<TrimEndpointCandidate> {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return vec![trim_decline(
            BrepEdgeSurfaceVertexCandidateOutcome::MissingTopology,
        )];
    };
    let Some(next) = brep.half_edges.get(half_edge_data.next) else {
        return vec![trim_decline(
            BrepEdgeSurfaceVertexCandidateOutcome::MissingTopology,
        )];
    };
    let Some(pcurve) = half_edge_data.pcurve.as_ref() else {
        return Vec::new();
    };
    let Some(origin) = brep.vertices.get(half_edge_data.origin) else {
        return vec![trim_decline(
            BrepEdgeSurfaceVertexCandidateOutcome::MissingTopology,
        )];
    };
    let Some(next_origin) = brep.vertices.get(next.origin) else {
        return vec![trim_decline(
            BrepEdgeSurfaceVertexCandidateOutcome::MissingTopology,
        )];
    };
    let expected_start = half_edge_data.origin == next.origin
        || current.distance_to(&origin.position) <= current.distance_to(&next_origin.position);
    let Some(endpoints) = trim_endpoint_candidates(pcurve, surface, current, expected_start) else {
        return vec![trim_decline(
            BrepEdgeSurfaceVertexCandidateOutcome::EvaluationFailed,
        )];
    };
    endpoints
        .into_iter()
        .map(|mut candidate| {
            if candidate.anchor_motion > options.max_anchor_motion {
                candidate.point = None;
                candidate.outcome = BrepEdgeSurfaceVertexCandidateOutcome::AnchorMotionExceeded;
            }
            candidate
        })
        .filter(|candidate| candidate.point.is_some())
        .collect()
}

pub(super) fn trim_decline(
    outcome: BrepEdgeSurfaceVertexCandidateOutcome,
) -> TrimEndpointCandidate {
    TrimEndpointCandidate {
        point: None,
        outcome,
        anchor_motion: 0.0,
    }
}

fn trim_endpoint_candidates(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    current: Point3,
    expected_start: bool,
) -> Option<[TrimEndpointCandidate; 2]> {
    let (t_start, t_end) = pcurve.parameter_range();
    Some([
        trim_endpoint_at(pcurve, surface, current, t_start, expected_start)?,
        trim_endpoint_at(pcurve, surface, current, t_end, !expected_start)?,
    ])
}

fn trim_endpoint_at(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    current: Point3,
    parameter: f64,
    expected: bool,
) -> Option<TrimEndpointCandidate> {
    let uv = pcurve.point_at(parameter);
    let point = surface.point_at(uv.x, uv.y);
    point_is_finite(point).then(|| TrimEndpointCandidate {
        point: Some(point),
        outcome: if expected {
            BrepEdgeSurfaceVertexCandidateOutcome::Accepted
        } else {
            BrepEdgeSurfaceVertexCandidateOutcome::ReconciledEndpoint
        },
        anchor_motion: point.distance_to(&current),
    })
}

fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}
