//! Single-surface edge repair for edge/surface compatibility.

use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::brep_edge_surface_compatibility_helpers::{
    edge_vertices, move_vertex_with_tolerance, projected_endpoint_uvs,
    proposed_vertex_incident_endpoint_gap,
};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, edge_segment};
use popkern_math::PortableMath;

pub(super) fn commit_single_surface_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    options: BrepEdgeSurfaceCompatibilityOptions,
    mut evidence: BrepEdgeSurfaceRepairEvidence,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    let Some((start, end)) = edge_segment(brep, edge) else {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
            evidence,
        );
    };
    let (repair_start, repair_end) =
        projected_surface_endpoints(surface, start, end).unwrap_or((start, end));
    let endpoint_gap = start
        .distance_to(&repair_start)
        .max(end.distance_to(&repair_end));
    if endpoint_gap > options.max_anchor_motion {
        evidence.anchor_motion = endpoint_gap;
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::AnchorMotionExceeded,
            evidence,
        );
    }
    let Some(curve) = single_surface_curve(surface, repair_start, repair_end, options.tolerance)
    else {
        evidence.boundary_failure = Some("unsupported_single_surface_trim".to_string());
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingAdjacentSurfaces,
            evidence,
        );
    };
    let residual = crate::healing::deviation::edge_surface_deviation(
        &curve.geometry,
        (curve.t_start, curve.t_end),
        surface,
        options.samples,
    )
    .max_deviation;
    evidence.residual_deviation_a = residual;
    if residual > options.tolerance {
        return BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::ResidualDeviation,
            ),
            endpoint_gap,
            residual_deviation: residual,
            evidence,
            shared_fact_ids: Vec::new(),
        };
    }
    if endpoint_gap > options.tolerance {
        let Some((start_vertex, end_vertex)) = edge_vertices(brep, edge) else {
            return declined(
                edge,
                BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
                evidence,
            );
        };
        let incident_gap = proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            start_vertex,
            repair_start,
            options.tolerance,
        )
        .max(proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            end_vertex,
            repair_end,
            options.tolerance,
        ));
        if incident_gap > options.max_anchor_motion {
            evidence.max_incident_endpoint_gap = incident_gap;
            return BrepEdgeSurfaceCompatibilityRepairDecision {
                edge,
                outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                    BrepEdgeSurfaceCompatibilityDeclineReason::IncidentEndpointConflict,
                ),
                endpoint_gap,
                residual_deviation: residual,
                evidence,
                shared_fact_ids: Vec::new(),
            };
        }
        evidence.anchor_motion = endpoint_gap;
        move_vertex_with_tolerance(brep, start_vertex, repair_start, endpoint_gap);
        move_vertex_with_tolerance(brep, end_vertex, repair_end, endpoint_gap);
    }
    let Some((he_start, he_end)) = half_edge_segment(brep, half_edge) else {
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::MissingTopology,
            evidence,
        );
    };
    let Some(pcurve) = single_surface_pcurve(surface, he_start, he_end) else {
        evidence.boundary_failure = Some("single_surface_pcurve_projection_failed".to_string());
        return declined(
            edge,
            BrepEdgeSurfaceCompatibilityDeclineReason::SurfaceBoundaryReconstructionFailed,
            evidence,
        );
    };
    brep.edges[edge].curve = Some(curve);
    brep.half_edges[half_edge].pcurve = Some(pcurve);
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary,
        endpoint_gap,
        residual_deviation: residual,
        evidence,
        shared_fact_ids: Vec::new(),
    }
}

fn projected_surface_endpoints(
    surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Option<(Point3, Point3)> {
    let (start_uv, end_uv) = projected_endpoint_uvs(surface, start, end)?;
    Some((
        surface.point_at(start_uv.0, start_uv.1),
        surface.point_at(end_uv.0, end_uv.1),
    ))
}

fn single_surface_curve(
    surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<CurveOnEdge> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let start_distance = (start - *origin).dot(normal).abs();
            let end_distance = (end - *origin).dot(normal).abs();
            (start_distance.max(end_distance) <= tolerance.max(1.0e-9))
                .then(|| CurveOnEdge::line_segment(start, end, tolerance))
                .flatten()
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => cylinder_surface_curve(*origin, *axis, *radius, start, end, tolerance),
        _ => None,
    }
}

fn cylinder_surface_curve(
    origin: Point3,
    axis: Vec3,
    radius: f64,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<CurveOnEdge> {
    let axis = axis.normalized()?;
    let start_axis = (start - origin).dot(&axis);
    let end_axis = (end - origin).dot(&axis);
    let start_center = origin + axis * start_axis;
    let end_center = origin + axis * end_axis;
    let start_radial = start - start_center;
    let end_radial = end - end_center;
    let start_radius = start_radial.length();
    let end_radius = end_radial.length();
    if (start_radius - radius)
        .abs()
        .max((end_radius - radius).abs())
        > tolerance.max(1.0e-8)
    {
        return None;
    }
    let start_dir = start_radial.normalized()?;
    let end_dir = end_radial.normalized()?;
    if start_dir.cross(&end_dir).length() <= tolerance.max(1.0e-8) {
        return CurveOnEdge::line_segment(start, end, tolerance);
    }
    if (start_axis - end_axis).abs() > tolerance.max(1.0e-8) {
        return None;
    }
    let center = start_center;
    let (u_axis, v_axis) = circle_basis_for_normal(axis);
    let start_angle = start_radial.dot(&v_axis).patan2(start_radial.dot(&u_axis));
    let mut end_angle = end_radial.dot(&v_axis).patan2(end_radial.dot(&u_axis));
    while end_angle - start_angle > std::f64::consts::PI {
        end_angle -= std::f64::consts::TAU;
    }
    while end_angle - start_angle < -std::f64::consts::PI {
        end_angle += std::f64::consts::TAU;
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal: axis,
            radius,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    })
}

fn circle_basis_for_normal(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

fn single_surface_pcurve(surface: &SurfaceGeometry, start: Point3, end: Point3) -> Option<Pcurve> {
    let (start_uv, end_uv) = projected_endpoint_uvs(surface, start, end)?;
    Some(Pcurve::linear(start_uv, end_uv, 0.0, 1.0))
}

fn half_edge_segment(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<(Point3, Point3)> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let start = brep.vertices.get(half_edge_data.origin)?.position;
    let end_vertex = brep.half_edges.get(half_edge_data.next)?.origin;
    let end = brep.vertices.get(end_vertex)?.position;
    Some((start, end))
}

fn declined(
    edge: EdgeHandle,
    reason: BrepEdgeSurfaceCompatibilityDeclineReason,
    evidence: BrepEdgeSurfaceRepairEvidence,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(reason),
        endpoint_gap: 0.0,
        residual_deviation: 0.0,
        evidence,
        shared_fact_ids: Vec::new(),
    }
}
