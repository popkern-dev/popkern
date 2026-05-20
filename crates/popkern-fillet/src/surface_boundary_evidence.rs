//! Final construction audit for bounded-surface boundary evidence.

#[path = "surface_boundary_audit_events.rs"]
mod surface_boundary_audit_events;

use self::surface_boundary_audit_events::{
    BoundaryAuditRecordInput, record_boundary_audit, surface_boundary_shared_trim_failure_reason,
    surface_boundary_shared_trim_validation,
};
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::edge_triple::{
    BoundarySide, EdgeAnchors, EdgeHint, EdgeTripleOptions, SurfaceSide, compute_edge_geometry,
};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{
    BRep, FaceHandle, HalfEdgeHandle, SharedTrimEdgeInput, SurfaceBoundaryAuditOutcome,
    SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason,
    commit_shared_trim_edge_transactionally,
};

pub(crate) fn audit_gregory_analytic_boundary_evidence(brep: &mut BRep) {
    let mut half_edges = Vec::new();
    for (_, face) in &brep.faces {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Gregory(_)) {
            continue;
        }
        let start = brep.loops[face.outer_loop].half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            half_edges.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }
    for half_edge in half_edges {
        audit_half_edge(brep, half_edge);
    }
}

fn audit_half_edge(brep: &mut BRep, half_edge: HalfEdgeHandle) {
    if let Some(boundary) = brep.surface_boundary_hint(half_edge) {
        if authored_boundary_hint_is_fresh(brep, half_edge, boundary) {
            record_boundary_audit(
                brep,
                BoundaryAuditRecordInput {
                    half_edge,
                    outcome: SurfaceBoundaryAuditOutcome::FreshHint,
                    boundary: Some(boundary),
                    ..BoundaryAuditRecordInput::default()
                },
            );
            return;
        }
        record_boundary_audit(
            brep,
            BoundaryAuditRecordInput {
                half_edge,
                outcome: SurfaceBoundaryAuditOutcome::StaleHintCleared,
                boundary: Some(boundary),
                ..BoundaryAuditRecordInput::default()
            },
        );
        brep.clear_surface_boundary_hint(half_edge);
    }
    if brep.surface_boundary_decline(half_edge).is_some() {
        record_boundary_audit(
            brep,
            BoundaryAuditRecordInput {
                half_edge,
                outcome: SurfaceBoundaryAuditOutcome::ReauditDecline,
                ..BoundaryAuditRecordInput::default()
            },
        );
    }
    if !has_analytic_neighbor(brep, half_edge) {
        return;
    }
    if try_author_exact_boundary_owned_triple(brep, half_edge) {
        return;
    }
    let Some(pcurve) = brep.half_edges[half_edge].pcurve.clone() else {
        brep.set_surface_boundary_decline(
            half_edge,
            SurfaceBoundaryHintDecline {
                reason: SurfaceBoundaryHintDeclineReason::MissingPcurve,
                best_boundary_deviation: f64::INFINITY,
                best_boundary: None,
                endpoint_gap: 0.0,
            },
        );
        return;
    };
    let classification = classify_rectangular_boundary_pcurve(&pcurve);
    if let Some(boundary) = classification.boundary {
        brep.set_surface_boundary_hint(half_edge, boundary);
        record_boundary_audit(
            brep,
            BoundaryAuditRecordInput {
                half_edge,
                outcome: SurfaceBoundaryAuditOutcome::ClassifiedHint,
                boundary: Some(boundary),
                tolerance: Some(brep.tolerance.linear().max(1.0e-7)),
                ..BoundaryAuditRecordInput::default()
            },
        );
    } else {
        brep.set_surface_boundary_decline(
            half_edge,
            SurfaceBoundaryHintDecline {
                reason: SurfaceBoundaryHintDeclineReason::NotOnCanonicalBoundary,
                best_boundary_deviation: classification.best_deviation,
                best_boundary: classification.boundary,
                endpoint_gap: 0.0,
            },
        );
        record_boundary_audit(
            brep,
            BoundaryAuditRecordInput {
                half_edge,
                outcome: SurfaceBoundaryAuditOutcome::ClassifiedDecline,
                boundary: classification.boundary,
                decline_reason: Some(SurfaceBoundaryHintDeclineReason::NotOnCanonicalBoundary),
                tolerance: Some(brep.tolerance.linear().max(1.0e-7)),
                ..BoundaryAuditRecordInput::default()
            },
        );
    }
}

fn authored_boundary_hint_is_fresh(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    boundary: BoundarySide,
) -> bool {
    let Some((cap_face, _, _)) = adjacent_faces(brep, half_edge) else {
        return false;
    };
    let Some(surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return false;
    };
    let Some(pcurve) = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| half_edge.pcurve.as_ref())
    else {
        return false;
    };
    let classification = classify_rectangular_boundary_pcurve(pcurve);
    let Some(best_boundary) = classification.boundary else {
        return false;
    };
    if best_boundary != boundary {
        return false;
    }
    if classification.best_deviation > brep.tolerance.linear().max(1.0e-7) {
        return false;
    }
    boundary_pcurve_endpoint_gap(brep, half_edge, surface, pcurve)
        <= brep.tolerance.linear().max(1.0e-7)
}

fn try_author_exact_boundary_owned_triple(brep: &mut BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some((cap_face, other_face, twin)) = adjacent_faces(brep, half_edge) else {
        return false;
    };
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return false;
    };
    if !matches!(cap_surface, SurfaceGeometry::Gregory(_)) {
        return false;
    }
    let Some(other_surface) = brep
        .faces
        .get(other_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| surface.geometry.clone())
    else {
        return false;
    };
    if !is_supported_analytic_surface(&other_surface) {
        return false;
    }
    let Some(start) = brep
        .half_edges
        .get(half_edge)
        .and_then(|data| brep.vertices.get(data.origin).map(|vertex| vertex.position))
    else {
        return false;
    };
    let Some(end) = brep
        .half_edges
        .get(twin)
        .and_then(|data| brep.vertices.get(data.origin).map(|vertex| vertex.position))
    else {
        return false;
    };
    let boundaries = [
        BoundarySide::UMin,
        BoundarySide::UMax,
        BoundarySide::VMin,
        BoundarySide::VMax,
    ];
    for boundary in boundaries {
        let hint = EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary,
        };
        let Ok(triple) = compute_edge_geometry(
            &cap_surface,
            &other_surface,
            EdgeAnchors::Open { start, end },
            Some(hint),
            &brep.tolerance,
            &EdgeTripleOptions::default(),
        ) else {
            continue;
        };
        let Some(edge) = brep.half_edges.get(half_edge).map(|data| data.edge) else {
            return false;
        };
        let curve = CurveOnEdge {
            geometry: triple.curve_3d,
            t_start: triple.t_start,
            t_end: triple.t_end,
        };
        let curve = orient_curve_to_canonical_edge(brep, edge, curve);
        let pcurve_a = orient_pcurve_to_half_edge(brep, half_edge, &cap_surface, triple.pcurve_a);
        let pcurve_b = orient_pcurve_to_half_edge(brep, twin, &other_surface, triple.pcurve_b);
        let shared_trim_result = commit_shared_trim_edge_transactionally(
            brep,
            SharedTrimEdgeInput {
                edge,
                half_edge_a: half_edge,
                half_edge_b: twin,
                curve,
                pcurve_a,
                pcurve_b,
                tolerance: brep.tolerance.linear().max(1.0e-7),
                validate_brep: false,
                provenance: "gregory_analytic_boundary_audit",
            },
        );
        if shared_trim_result.is_ok() {
            record_boundary_audit(
                brep,
                BoundaryAuditRecordInput {
                    half_edge,
                    outcome: SurfaceBoundaryAuditOutcome::SharedTrimAccepted,
                    boundary: Some(boundary),
                    tolerance: Some(brep.tolerance.linear().max(1.0e-7)),
                    ..BoundaryAuditRecordInput::default()
                },
            );
            return brep.set_surface_boundary_hint(half_edge, boundary);
        }
        if let Err(failure) = shared_trim_result {
            record_boundary_audit(
                brep,
                BoundaryAuditRecordInput {
                    half_edge,
                    outcome: SurfaceBoundaryAuditOutcome::SharedTrimDeclined,
                    boundary: Some(boundary),
                    shared_trim_failure: Some(surface_boundary_shared_trim_failure_reason(
                        failure.reason,
                    )),
                    shared_trim_validation: Some(surface_boundary_shared_trim_validation(
                        failure.validation,
                    )),
                    tolerance: Some(brep.tolerance.linear().max(1.0e-7)),
                    ..BoundaryAuditRecordInput::default()
                },
            );
        }
    }
    false
}

fn orient_pcurve_to_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: Pcurve,
) -> Pcurve {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return pcurve;
    };
    let Some(next) = brep.half_edges.get(half_edge_data.next) else {
        return pcurve;
    };
    let Some(start) = brep
        .vertices
        .get(half_edge_data.origin)
        .map(|vertex| vertex.position)
    else {
        return pcurve;
    };
    let Some(end) = brep.vertices.get(next.origin).map(|vertex| vertex.position) else {
        return pcurve;
    };
    let p_start = surface_point_at_pcurve_endpoint(surface, &pcurve, pcurve.t_start);
    let p_end = surface_point_at_pcurve_endpoint(surface, &pcurve, pcurve.t_end);
    let forward_gap = p_start.distance_to(&start).max(p_end.distance_to(&end));
    let reversed_gap = p_start.distance_to(&end).max(p_end.distance_to(&start));
    if reversed_gap < forward_gap {
        pcurve.reversed()
    } else {
        pcurve
    }
}

fn orient_curve_to_canonical_edge(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    mut curve: CurveOnEdge,
) -> CurveOnEdge {
    let Some(edge_data) = brep.edges.get(edge) else {
        return curve;
    };
    let Some(half_edge_data) = brep.half_edges.get(edge_data.half_edge) else {
        return curve;
    };
    let Some(next) = brep.half_edges.get(half_edge_data.next) else {
        return curve;
    };
    let Some(start) = brep
        .vertices
        .get(half_edge_data.origin)
        .map(|vertex| vertex.position)
    else {
        return curve;
    };
    let Some(end) = brep.vertices.get(next.origin).map(|vertex| vertex.position) else {
        return curve;
    };
    let p_start = curve.geometry.point_at(curve.t_start);
    let p_end = curve.geometry.point_at(curve.t_end);
    let forward_gap = p_start.distance_to(&start).max(p_end.distance_to(&end));
    let reversed_gap = p_start.distance_to(&end).max(p_end.distance_to(&start));
    if reversed_gap < forward_gap {
        std::mem::swap(&mut curve.t_start, &mut curve.t_end);
    }
    curve
}

fn surface_point_at_pcurve_endpoint(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pcurve.point_at(parameter);
    surface.point_at(uv.x, uv.y)
}

fn boundary_pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> f64 {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return f64::INFINITY;
    };
    let Some(twin) = brep.half_edges.get(half_edge_data.twin) else {
        return f64::INFINITY;
    };
    let Some(start) = brep
        .vertices
        .get(half_edge_data.origin)
        .map(|vertex| vertex.position)
    else {
        return f64::INFINITY;
    };
    let Some(end) = brep.vertices.get(twin.origin).map(|vertex| vertex.position) else {
        return f64::INFINITY;
    };
    let p_start = surface_point_at_pcurve_endpoint(surface, pcurve, pcurve.t_start);
    let p_end = surface_point_at_pcurve_endpoint(surface, pcurve, pcurve.t_end);
    start.distance_to(&p_start).max(end.distance_to(&p_end))
}

fn adjacent_faces(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, FaceHandle, HalfEdgeHandle)> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let cap_face = brep.loops.get(half_edge_data.loop_)?.face;
    let twin = half_edge_data.twin;
    let other_face = brep.loops.get(brep.half_edges.get(twin)?.loop_)?.face;
    (cap_face != other_face).then_some((cap_face, other_face, twin))
}

fn has_analytic_neighbor(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(half_edge_data.twin) else {
        return false;
    };
    let Some(twin_face) = brep.loops.get(twin.loop_).map(|loop_data| loop_data.face) else {
        return false;
    };
    brep.faces
        .get(twin_face)
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| is_supported_analytic_surface(&surface.geometry))
}

fn is_supported_analytic_surface(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    )
}

struct BoundaryPcurveClassification {
    boundary: Option<BoundarySide>,
    best_deviation: f64,
}

fn classify_rectangular_boundary_pcurve(pcurve: &Pcurve) -> BoundaryPcurveClassification {
    let (t0, t1) = pcurve.parameter_range();
    let mut max_u_min = 0.0_f64;
    let mut max_u_max = 0.0_f64;
    let mut max_v_min = 0.0_f64;
    let mut max_v_max = 0.0_f64;
    for i in 0..9 {
        let s = i as f64 / 8.0;
        let uv = pcurve.point_at(t0 + s * (t1 - t0));
        max_u_min = max_u_min.max(uv.x.abs());
        max_u_max = max_u_max.max((1.0 - uv.x).abs());
        max_v_min = max_v_min.max(uv.y.abs());
        max_v_max = max_v_max.max((1.0 - uv.y).abs());
    }
    let candidates = [
        (BoundarySide::UMin, max_u_min),
        (BoundarySide::UMax, max_u_max),
        (BoundarySide::VMin, max_v_min),
        (BoundarySide::VMax, max_v_max),
    ];
    let Some(best) = candidates.into_iter().min_by(|a, b| a.1.total_cmp(&b.1)) else {
        return BoundaryPcurveClassification {
            boundary: None,
            best_deviation: f64::INFINITY,
        };
    };
    BoundaryPcurveClassification {
        boundary: (best.1 <= 1.0e-6).then_some(best.0),
        best_deviation: best.1,
    }
}

#[cfg(test)]
#[path = "surface_boundary_evidence_tests.rs"]
mod tests;
