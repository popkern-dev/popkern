//! Pcurve-trace repair candidate for edge/surface compatibility.

use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::topology::brep_curve_trim_reconstruction_helpers::project_edge_curve_to_surface_pcurve;
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::brep_edge_surface_compatibility_helpers::{
    edge_vertices, endpoint_gap, move_vertex_with_tolerance, proposed_vertex_incident_endpoint_gap,
    residuals_for_candidate,
};
use crate::topology::{
    BRep, EdgeHandle, HalfEdgeHandle, SharedTrimEdgeInput, commit_shared_trim_edge_transactionally,
};
use popkern_math::PortableMath;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_commit_pcurve_trace_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    options: BrepEdgeSurfaceCompatibilityOptions,
    mut evidence: BrepEdgeSurfaceRepairEvidence,
) -> Option<BrepEdgeSurfaceCompatibilityRepairDecision> {
    if !supports_surface_trim_authority(surf_a) && !supports_surface_trim_authority(surf_b) {
        return None;
    }
    let candidates = [
        cylinder_supported_surface_candidate_curve(brep, edge, surf_a, surf_b),
        supports_surface_trim_authority(surf_a)
            .then(|| pcurve_trace_curve(brep, edge, he_a, surf_a))
            .flatten(),
        supports_surface_trim_authority(surf_b)
            .then(|| pcurve_trace_curve(brep, edge, he_b, surf_b))
            .flatten(),
    ];
    let (curve, residual_a, residual_b, endpoint_gap) = candidates
        .into_iter()
        .flatten()
        .filter_map(|curve| {
            let (residual_a, residual_b) = residuals_for_candidate(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                surf_a,
                surf_b,
                options.samples,
            );
            let (start, end) = crate::topology::edge_segment(brep, edge)?;
            let endpoint_gap =
                endpoint_gap(&curve.geometry, curve.t_start, curve.t_end, start, end);
            (residual_a.max(residual_b) <= options.tolerance
                && endpoint_gap <= options.max_anchor_motion)
                .then_some((curve, residual_a, residual_b, endpoint_gap))
        })
        .min_by(|left, right| {
            left.1
                .max(left.2)
                .max(left.3)
                .total_cmp(&right.1.max(right.2).max(right.3))
        })?;
    evidence.residual_deviation_a = residual_a;
    evidence.residual_deviation_b = residual_b;
    evidence.anchor_motion = endpoint_gap;
    let rollback_snapshot = brep.clone();
    if endpoint_gap > options.tolerance {
        let (start_vertex, end_vertex) = edge_vertices(brep, edge)?;
        let repaired_start = curve.geometry.point_at(curve.t_start);
        let repaired_end = curve.geometry.point_at(curve.t_end);
        let incident_gap = proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            start_vertex,
            repaired_start,
            options.tolerance,
        )
        .max(proposed_vertex_incident_endpoint_gap(
            brep,
            edge,
            end_vertex,
            repaired_end,
            options.tolerance,
        ));
        if incident_gap > options.max_anchor_motion {
            evidence.max_incident_endpoint_gap = incident_gap;
            return Some(BrepEdgeSurfaceCompatibilityRepairDecision {
                edge,
                outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                    BrepEdgeSurfaceCompatibilityDeclineReason::IncidentEndpointConflict,
                ),
                endpoint_gap,
                residual_deviation: residual_a.max(residual_b),
                evidence,
                shared_fact_ids: Vec::new(),
            });
        }
        move_vertex_with_tolerance(brep, start_vertex, repaired_start, endpoint_gap);
        move_vertex_with_tolerance(brep, end_vertex, repaired_end, endpoint_gap);
    }
    let mut pcurve_source = brep.clone();
    pcurve_source.edges[edge].curve = Some(curve.clone());
    let Some(pcurve_a) =
        project_edge_curve_to_surface_pcurve(&pcurve_source, he_a, options.samples)
    else {
        *brep = rollback_snapshot;
        return None;
    };
    let Some(pcurve_b) =
        project_edge_curve_to_surface_pcurve(&pcurve_source, he_b, options.samples)
    else {
        *brep = rollback_snapshot;
        return None;
    };
    if let Err(failure) = commit_shared_trim_edge_transactionally(
        brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a: he_a,
            half_edge_b: he_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: options.tolerance,
            validate_brep: options.validate_output,
            provenance: "edge_surface_compatibility_pcurve_trace_repair",
        },
    ) {
        *brep = rollback_snapshot;
        evidence.shared_trim_failure = Some(failure.reason);
        evidence.shared_trim_curve_surface_deviation =
            failure.validation.max_curve_surface_deviation;
        evidence.shared_trim_side_agreement_deviation =
            failure.validation.max_side_agreement_deviation;
        evidence.shared_trim_endpoint_gap = failure.validation.max_endpoint_gap;
        return Some(BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::SharedTrimValidationFailed,
            ),
            endpoint_gap: 0.0,
            residual_deviation: residual_a.max(residual_b),
            evidence,
            shared_fact_ids: Vec::new(),
        });
    }
    Some(BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary,
        endpoint_gap,
        residual_deviation: residual_a.max(residual_b),
        evidence,
        shared_fact_ids: Vec::new(),
    })
}

fn supports_surface_trim_authority(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
    )
}

fn cylinder_supported_surface_candidate_curve(
    brep: &BRep,
    edge: EdgeHandle,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let cylinder = match (surf_a, surf_b) {
        (SurfaceGeometry::Cylinder { .. }, SurfaceGeometry::Gregory(_))
        | (SurfaceGeometry::Cylinder { .. }, SurfaceGeometry::Nurbs(_)) => surf_a,
        (SurfaceGeometry::Gregory(_), SurfaceGeometry::Cylinder { .. })
        | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Cylinder { .. }) => surf_b,
        _ => return None,
    };
    let (start, end) = crate::topology::edge_segment(brep, edge)?;
    cylinder_curve_between(cylinder, start, end)
}

fn cylinder_curve_between(
    cylinder: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = cylinder
    else {
        return None;
    };
    let axis_unit = axis.normalized()?;
    let start_axis = (start - *origin).dot(&axis_unit);
    let end_axis = (end - *origin).dot(&axis_unit);
    let (u_axis, v_axis) = circle_basis_for_normal(axis_unit);
    let start_center = *origin + axis_unit * start_axis;
    let end_center = *origin + axis_unit * end_axis;
    let start_radial = start - start_center;
    let end_radial = end - end_center;
    let start_angle = start_radial.dot(&v_axis).patan2(start_radial.dot(&u_axis));
    let mut end_angle = end_radial.dot(&v_axis).patan2(end_radial.dot(&u_axis));
    while end_angle - start_angle > std::f64::consts::PI {
        end_angle -= std::f64::consts::TAU;
    }
    while end_angle - start_angle < -std::f64::consts::PI {
        end_angle += std::f64::consts::TAU;
    }
    if (start_axis - end_axis).abs() <= 1.0e-9 {
        return Some(CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center: *origin + axis_unit * start_axis,
                normal: axis_unit,
                radius: *radius,
                start_angle,
                end_angle,
            },
            t_start: start_angle,
            t_end: end_angle,
        });
    }
    let samples = 33usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let angle = start_angle + (end_angle - start_angle) * fraction;
        let station = start_axis + (end_axis - start_axis) * fraction;
        let center = *origin + axis_unit * station;
        points.push(center + u_axis * (*radius * angle.pcos()) + v_axis * (*radius * angle.psin()));
        parameters.push(fraction);
    }
    let degree = 5.min(points.len().saturating_sub(1));
    let geometry = crate::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)?;
    Some(CurveOnEdge {
        geometry,
        t_start: 0.0,
        t_end: 1.0,
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

fn pcurve_trace_curve(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let pcurve = brep.half_edges.get(half_edge)?.pcurve.as_ref()?;
    let edge_data = brep.edges.get(edge)?;
    let forward = half_edge == edge_data.half_edge;
    let (t_start, t_end) = pcurve.parameter_range();
    let samples = 33usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let pcurve_fraction = if forward { fraction } else { 1.0 - fraction };
        let parameter = t_start + (t_end - t_start) * pcurve_fraction;
        let uv = pcurve.point_at(parameter);
        points.push(surface.point_at(uv.x, uv.y));
        parameters.push(fraction);
    }
    let degree = 5.min(points.len().saturating_sub(1));
    let geometry = crate::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)?;
    Some(CurveOnEdge {
        geometry,
        t_start: 0.0,
        t_end: 1.0,
    })
}
