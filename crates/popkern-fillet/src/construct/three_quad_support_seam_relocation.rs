//! Relocation-backed support seam authoring helpers.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::fitting;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "three_quad_support_seam_validation.rs"]
mod validation;
use crate::types::SupportSeamDeclineReason;

const SAMPLE_COUNT: usize = 17;
const SUPPORT_SEAM_TOLERANCE: f64 = 2.0e-4;

pub(super) struct SupportSeamAttempt<'a> {
    pub(super) brep: &'a mut BRep,
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) cap_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) edge: EdgeHandle,
    pub(super) cap_surface: &'a SurfaceGeometry,
    pub(super) boundary_side: BoundarySide,
    pub(super) min_cylinder_neighbors: usize,
}

pub(super) struct SupportSeamAttemptOutcome {
    pub(super) accepted: bool,
    pub(super) reason: Option<SupportSeamDeclineReason>,
    pub(super) max_endpoint_motion: f64,
    pub(super) max_sampled_disagreement: f64,
}

impl SupportSeamAttemptOutcome {
    fn accepted(max_endpoint_motion: f64, max_sampled_disagreement: f64) -> Self {
        Self {
            accepted: true,
            reason: None,
            max_endpoint_motion,
            max_sampled_disagreement,
        }
    }

    fn declined(reason: SupportSeamDeclineReason) -> Self {
        Self {
            accepted: false,
            reason: Some(reason),
            max_endpoint_motion: 0.0,
            max_sampled_disagreement: 0.0,
        }
    }

    fn declined_with_metrics(
        reason: SupportSeamDeclineReason,
        max_endpoint_motion: f64,
        max_sampled_disagreement: f64,
    ) -> Self {
        Self {
            accepted: false,
            reason: Some(reason),
            max_endpoint_motion,
            max_sampled_disagreement,
        }
    }
}

pub(super) fn try_author_one(attempt: SupportSeamAttempt<'_>) -> SupportSeamAttemptOutcome {
    let brep = attempt.brep;
    let before = brep.clone();
    let max_endpoint_motion = match relocate_support_vertices_to_boundary(
        brep,
        attempt.half_edge,
        attempt.cap_surface,
        attempt.boundary_side,
    ) {
        Ok(motion) => motion,
        Err(reason) => {
            super::trace_decline(
                attempt.half_edge,
                "support boundary vertex relocation failed",
            );
            return SupportSeamAttemptOutcome::declined(reason);
        }
    };
    let Some((curve, cap_pcurve)) = build_cap_boundary_curve(
        brep,
        attempt.half_edge,
        attempt.cap_surface,
        attempt.boundary_side,
    ) else {
        super::trace_decline(
            attempt.half_edge,
            "cap boundary endpoints do not match edge vertices",
        );
        return SupportSeamAttemptOutcome::declined_with_metrics(
            SupportSeamDeclineReason::CapBoundaryMismatch,
            max_endpoint_motion,
            f64::INFINITY,
        );
    };
    brep.edges[attempt.edge].curve = Some(curve);
    brep.half_edges[attempt.half_edge].pcurve = Some(cap_pcurve);
    if lift_support_face_to_edge_curve(brep, attempt.support_face, attempt.edge).is_err() {
        super::trace_decline(attempt.half_edge, "support face lift failed");
        *brep = before;
        return SupportSeamAttemptOutcome::declined_with_metrics(
            SupportSeamDeclineReason::SupportFaceLiftFailed,
            max_endpoint_motion,
            f64::INFINITY,
        );
    }
    let twin = brep.half_edges[attempt.half_edge].twin;
    write_lifted_support_pcurve(brep, twin);
    let seam_disagreement = validation::validate_support_seam(
        brep,
        attempt.half_edge,
        twin,
        attempt.cap_surface,
        SUPPORT_SEAM_TOLERANCE,
        SAMPLE_COUNT,
    );
    let vertices_valid = validation::validate_moved_vertex_curve_coherence(
        brep,
        attempt.half_edge,
        SUPPORT_SEAM_TOLERANCE,
    );
    let cap_neighbor_valid = super::cap_face_preserves_cylinder_neighbor(brep, attempt.cap_face);
    let global_neighbors_valid =
        super::gregory_cap_cylinder_neighbor_count(brep) >= attempt.min_cylinder_neighbors;
    let max_sampled_disagreement = match seam_disagreement {
        Ok(value) | Err(value) => value,
    };
    if seam_disagreement.is_ok() && vertices_valid && cap_neighbor_valid && global_neighbors_valid {
        return SupportSeamAttemptOutcome::accepted(max_endpoint_motion, max_sampled_disagreement);
    }
    let reason = if seam_disagreement.is_err() {
        super::trace_decline(attempt.half_edge, "support seam sampled validation failed");
        SupportSeamDeclineReason::SampledValidationFailed
    } else if !vertices_valid {
        super::trace_decline(
            attempt.half_edge,
            "moved vertex incident curve coherence failed",
        );
        SupportSeamDeclineReason::MovedVertexCurveIncoherent
    } else if !cap_neighbor_valid {
        super::trace_decline(
            attempt.half_edge,
            "cap-cylinder neighbor preservation failed",
        );
        SupportSeamDeclineReason::CapCylinderNeighborLost
    } else {
        super::trace_decline(
            attempt.half_edge,
            "global cap-cylinder neighbor preservation failed",
        );
        SupportSeamDeclineReason::GlobalCapCylinderNeighborLost
    };
    *brep = before;
    SupportSeamAttemptOutcome::declined_with_metrics(
        reason,
        max_endpoint_motion,
        max_sampled_disagreement,
    )
}

fn lift_support_face_to_edge_curve(
    brep: &mut BRep,
    support_face: FaceHandle,
    edge: EdgeHandle,
) -> Result<(), ()> {
    let support_surface = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .ok_or(())?;
    let SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    } = support_surface.geometry
    else {
        return Err(());
    };
    let normal = plane_normal.normalized().ok_or(())?;
    let curve = brep
        .edges
        .get(edge)
        .and_then(|edge| edge.curve.as_ref())
        .ok_or(())?;
    let curve_nurbs = curve
        .geometry
        .to_nurbs_in_range((curve.t_start, curve.t_end));
    let mut control_points = Vec::with_capacity(curve_nurbs.control_points.len() * 2);
    control_points.extend(curve_nurbs.control_points.iter().copied());
    control_points.extend(curve_nurbs.control_points.iter().map(|point| {
        let offset = *point - plane_origin;
        *point - normal * offset.dot(&normal)
    }));
    let mut weights = Vec::with_capacity(curve_nurbs.weights.len() * 2);
    weights.extend(curve_nurbs.weights.iter().copied());
    weights.extend(curve_nurbs.weights.iter().copied());
    let surface = NurbsSurface {
        degree_u: 1,
        degree_v: curve_nurbs.degree,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: curve_nurbs.knots,
        control_points,
        num_u: 2,
        num_v: curve_nurbs.control_points.len(),
        weights,
    };
    let mut lifted = support_surface.clone();
    lifted.geometry = SurfaceGeometry::Nurbs(surface);
    brep.faces[support_face].surface = Some(lifted);
    Ok(())
}

fn write_lifted_support_pcurve(brep: &mut BRep, support_half_edge: HalfEdgeHandle) {
    let Some(edge) = brep.half_edges.get(support_half_edge).map(|he| he.edge) else {
        return;
    };
    let Some(curve) = brep.edges.get(edge).and_then(|edge| edge.curve.as_ref()) else {
        return;
    };
    let t_start = curve.t_start;
    let t_end = curve.t_end;
    let Some(surface) = brep
        .half_edges
        .get(support_half_edge)
        .and_then(|he| brep.loops.get(he.loop_))
        .and_then(|loop_| brep.faces.get(loop_.face))
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return;
    };
    let u_boundary = best_lifted_support_u_boundary(&surface, curve);
    let uv_start = (u_boundary, t_start);
    let uv_end = (u_boundary, t_end);
    brep.half_edges[support_half_edge].pcurve =
        Some(Pcurve::linear(uv_start, uv_end, t_start, t_end));
}

fn best_lifted_support_u_boundary(surface: &SurfaceGeometry, curve: &CurveOnEdge) -> f64 {
    let mut best = (0.0, f64::INFINITY);
    for u in [0.0, 1.0] {
        let mut worst = 0.0_f64;
        for i in 0..SAMPLE_COUNT {
            let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
            let t = curve.t_start + s * (curve.t_end - curve.t_start);
            let p_curve = curve.geometry.point_at(t);
            let p_surface = surface.point_at(u, t);
            worst = worst.max(p_curve.distance_to(&p_surface));
        }
        if worst < best.1 {
            best = (u, worst);
        }
    }
    if crate::a2_trace::cap_g0_enabled() {
        eprintln!(
            "[a2-trace cap-g0] lifted support best boundary u={:.1} worst={:.6e}",
            best.0, best.1
        );
    }
    best.0
}

fn relocate_support_vertices_to_boundary(
    brep: &mut BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    boundary_side: BoundarySide,
) -> Result<f64, SupportSeamDeclineReason> {
    let Some((target_start, target_end)) =
        oriented_boundary_endpoints_for_half_edge(brep, half_edge, cap_surface, boundary_side)
    else {
        return Err(SupportSeamDeclineReason::NoCanonicalBoundary);
    };
    let Some(twin) = brep.half_edges.get(half_edge).map(|he| he.twin) else {
        return Err(SupportSeamDeclineReason::NoCanonicalBoundary);
    };
    let Some(start_vertex) = brep.half_edges.get(half_edge).map(|he| he.origin) else {
        return Err(SupportSeamDeclineReason::NoCanonicalBoundary);
    };
    let Some(end_vertex) = brep.half_edges.get(twin).map(|he| he.origin) else {
        return Err(SupportSeamDeclineReason::NoCanonicalBoundary);
    };
    let old_start = brep.vertices[start_vertex].position;
    let old_end = brep.vertices[end_vertex].position;
    if target_start.distance_to(&target_end) <= SUPPORT_SEAM_TOLERANCE {
        return Err(SupportSeamDeclineReason::DegenerateBoundary);
    }
    let current_span = old_start.distance_to(&old_end).max(SUPPORT_SEAM_TOLERANCE);
    let allowed_move = current_span.max(target_start.distance_to(&target_end)) * 2.0;
    if old_start
        .distance_to(&target_start)
        .max(old_end.distance_to(&target_end))
        > allowed_move
    {
        return Err(SupportSeamDeclineReason::ExcessiveEndpointMotion);
    }
    let max_motion = old_start
        .distance_to(&target_start)
        .max(old_end.distance_to(&target_end));
    brep.vertices[start_vertex].position = target_start;
    brep.vertices[end_vertex].position = target_end;
    snap_incident_edge_curve_endpoints(brep, start_vertex, old_start, target_start);
    snap_incident_edge_curve_endpoints(brep, end_vertex, old_end, target_end);
    Ok(max_motion)
}

fn oriented_boundary_endpoints_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    boundary_side: BoundarySide,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    let start = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let end = brep
        .vertices
        .get(brep.half_edges.get(twin)?.origin)?
        .position;
    let (uv0, uv1) = boundary_uv_range(boundary_side)?;
    let p0 = cap_surface.point_at(uv0.0, uv0.1);
    let p1 = cap_surface.point_at(uv1.0, uv1.1);
    if start.distance_to(&p0) + end.distance_to(&p1)
        <= start.distance_to(&p1) + end.distance_to(&p0)
    {
        Some((p0, p1))
    } else {
        Some((p1, p0))
    }
}

fn snap_incident_edge_curve_endpoints(
    brep: &mut BRep,
    vertex: VertexHandle,
    old_position: popkern_core::geometry::Point3,
    new_position: popkern_core::geometry::Point3,
) {
    let incident_edges: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter_map(|(edge_handle, edge)| {
            let primary = edge.half_edge;
            let twin = brep.half_edges.get(primary)?.twin;
            let a = brep.half_edges.get(primary)?.origin;
            let b = brep.half_edges.get(twin)?.origin;
            (a == vertex || b == vertex).then_some(edge_handle)
        })
        .collect();
    for edge_handle in incident_edges {
        let Some(curve_on_edge) = brep
            .edges
            .get(edge_handle)
            .and_then(|edge| edge.curve.clone())
        else {
            continue;
        };
        let mut nurbs = curve_on_edge
            .geometry
            .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));
        let first_distance = nurbs
            .control_points
            .first()
            .map(|point| point.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        let last_distance = nurbs
            .control_points
            .last()
            .map(|point| point.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        if first_distance <= last_distance {
            if let Some(first) = nurbs.control_points.first_mut() {
                *first = new_position;
            }
        } else if let Some(last) = nurbs.control_points.last_mut() {
            *last = new_position;
        }
        let (t_start, t_end) = nurbs.parameter_range();
        brep.edges[edge_handle].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start,
            t_end,
        });
    }
}

fn build_cap_boundary_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    boundary_side: BoundarySide,
) -> Option<(CurveOnEdge, Pcurve)> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    let start = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let end = brep
        .vertices
        .get(brep.half_edges.get(twin)?.origin)?
        .position;
    let (uv0, uv1) = boundary_uv_range(boundary_side)?;
    let p0 = cap_surface.point_at(uv0.0, uv0.1);
    let p1 = cap_surface.point_at(uv1.0, uv1.1);
    let forward = start.distance_to(&p0) + end.distance_to(&p1)
        <= start.distance_to(&p1) + end.distance_to(&p0);
    let (uv_start, uv_end) = if forward { (uv0, uv1) } else { (uv1, uv0) };
    let start_gap = start.distance_to(&cap_surface.point_at(uv_start.0, uv_start.1));
    let end_gap = end.distance_to(&cap_surface.point_at(uv_end.0, uv_end.1));
    if start_gap.max(end_gap) > SUPPORT_SEAM_TOLERANCE {
        return None;
    }

    let mut points = Vec::with_capacity(SAMPLE_COUNT);
    let mut params = Vec::with_capacity(SAMPLE_COUNT);
    for i in 0..SAMPLE_COUNT {
        let s = i as f64 / (SAMPLE_COUNT - 1) as f64;
        let u = uv_start.0 + s * (uv_end.0 - uv_start.0);
        let v = uv_start.1 + s * (uv_end.1 - uv_start.1);
        points.push(cap_surface.point_at(u, v));
        params.push(s);
    }
    let curve = fitting::interpolate_curve(&points, &params, 5).ok()?;
    Some((
        CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve),
            t_start: 0.0,
            t_end: 1.0,
        },
        Pcurve::linear(uv_start, uv_end, 0.0, 1.0),
    ))
}

fn boundary_uv_range(boundary_side: BoundarySide) -> Option<((f64, f64), (f64, f64))> {
    match boundary_side {
        BoundarySide::UMin => Some(((0.0, 0.0), (0.0, 1.0))),
        BoundarySide::UMax => Some(((1.0, 0.0), (1.0, 1.0))),
        BoundarySide::VMin => Some(((0.0, 0.0), (1.0, 0.0))),
        BoundarySide::VMax => Some(((0.0, 1.0), (1.0, 1.0))),
        BoundarySide::Indexed(_) => None,
    }
}
