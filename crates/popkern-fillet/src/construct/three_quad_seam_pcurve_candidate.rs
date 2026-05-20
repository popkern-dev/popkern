//! Boundary-owned seam pcurve candidate validation.

use super::*;
use crate::construct::helpers;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::trim_domain::{
    TrimDomainCandidateKind, classify_trim_domain_candidate,
};

const BOUNDARY_AUTHORING_SAMPLE_COUNT: usize = 16;

pub(super) enum BoundaryOwnedCandidate {
    AcceptedBoundary {
        half_edge: popkern_core::topology::HalfEdgeHandle,
        boundary_side: BoundarySide,
        curve: CurveOnEdge,
        cap_pcurve: Pcurve,
        other_pcurve: Pcurve,
    },
    AcceptedInterior {
        half_edge: popkern_core::topology::HalfEdgeHandle,
        curve: CurveOnEdge,
        cap_pcurve: Pcurve,
        other_pcurve: Pcurve,
    },
    Declined {
        half_edge: popkern_core::topology::HalfEdgeHandle,
        max_disagreement: f64,
        boundary_side: BoundarySide,
        endpoint_gap: f64,
    },
    Unavailable,
}

pub(super) fn boundary_owned_candidate(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &SurfaceOnFace,
    segment: &BoundarySegment,
    segment_offset: usize,
) -> BoundaryOwnedCandidate {
    let Some(twin) = brep.half_edges.get(he).map(|half_edge| half_edge.twin) else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    let Some(twin_face) = brep
        .half_edges
        .get(twin)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_| loop_.face)
    else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    let Some(twin_surface) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref()) else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    let Some(reference_curve) = boundary_reference_curve_for_candidate(brep, he, segment) else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    let other_pcurve =
        project_reference_curve_to_half_edge_pcurve(brep, twin, twin_surface, &reference_curve);
    let Some(other_pcurve) = other_pcurve else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    let tolerance = tolerance_for_boundary_authoring(brep);
    let mut best: Option<(BoundarySide, Pcurve, f64, f64)> = None;
    for boundary_side in boundary_sides_for_segment_offset(segment_offset) {
        let cap_pcurve =
            canonical_gregory_boundary_pcurve_for_side(brep, he, cap_surface, *boundary_side);
        let endpoint_gap = boundary_pcurve_endpoint_gap(brep, he, cap_surface, &cap_pcurve);
        let max_disagreement = boundary_owned_seam_disagreement(BoundarySeamEvidence {
            brep,
            he,
            twin,
            cap_surface,
            twin_surface,
            reference_curve: &reference_curve,
            cap_pcurve: &cap_pcurve,
            other_pcurve: &other_pcurve,
        });
        match best {
            Some((_, _, best_disagreement, best_endpoint_gap))
                if max_disagreement.max(endpoint_gap)
                    >= best_disagreement.max(best_endpoint_gap) => {}
            _ => best = Some((*boundary_side, cap_pcurve, max_disagreement, endpoint_gap)),
        }
    }
    let Some((boundary_side, cap_pcurve, max_disagreement, endpoint_gap)) = best else {
        return BoundaryOwnedCandidate::Unavailable;
    };
    if max_disagreement <= tolerance && endpoint_gap <= tolerance {
        return BoundaryOwnedCandidate::AcceptedBoundary {
            half_edge: he,
            boundary_side,
            curve: reference_curve.curve_on_edge(),
            cap_pcurve,
            other_pcurve,
        };
    }
    if let Some(accepted) = projected_trim_domain_candidate(ProjectedTrimCandidateInput {
        brep,
        he,
        twin,
        cap_surface,
        twin_surface,
        reference_curve: &reference_curve,
        other_pcurve: &other_pcurve,
        tolerance,
    }) {
        return accepted;
    }
    if crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!(
            "[direct-apex-trace cap-g0] boundary-owned pcurve: declined he={he:?} \
             disagreement={max_disagreement:.6e} endpoint_gap={endpoint_gap:.6e}"
        );
    }
    BoundaryOwnedCandidate::Declined {
        half_edge: he,
        max_disagreement,
        boundary_side,
        endpoint_gap,
    }
}

struct ProjectedTrimCandidateInput<'a> {
    brep: &'a BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    twin: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &'a SurfaceOnFace,
    twin_surface: &'a SurfaceOnFace,
    reference_curve: &'a BoundaryReferenceCurve<'a>,
    other_pcurve: &'a Pcurve,
    tolerance: f64,
}

fn projected_trim_domain_candidate(
    input: ProjectedTrimCandidateInput<'_>,
) -> Option<BoundaryOwnedCandidate> {
    let ProjectedTrimCandidateInput {
        brep,
        he,
        twin,
        cap_surface,
        twin_surface,
        reference_curve,
        other_pcurve,
        tolerance,
    } = input;
    let cap_pcurve =
        project_reference_curve_to_half_edge_pcurve(brep, he, cap_surface, reference_curve)?;
    let anchors = half_edge_topology_endpoints(brep, he)?;
    let classification = classify_trim_domain_candidate(
        &cap_surface.geometry,
        &cap_pcurve,
        Some(anchors),
        tolerance,
        BOUNDARY_AUTHORING_SAMPLE_COUNT,
    )
    .ok()?;
    let max_disagreement = boundary_owned_seam_disagreement(BoundarySeamEvidence {
        brep,
        he,
        twin,
        cap_surface,
        twin_surface,
        reference_curve,
        cap_pcurve: &cap_pcurve,
        other_pcurve,
    });
    if max_disagreement > tolerance {
        return None;
    }
    match classification.kind {
        TrimDomainCandidateKind::BoundarySubspan { side } => {
            Some(BoundaryOwnedCandidate::AcceptedBoundary {
                half_edge: he,
                boundary_side: side,
                curve: reference_curve.curve_on_edge(),
                cap_pcurve,
                other_pcurve: other_pcurve.clone(),
            })
        }
        TrimDomainCandidateKind::Interior => Some(BoundaryOwnedCandidate::AcceptedInterior {
            half_edge: he,
            curve: reference_curve.curve_on_edge(),
            cap_pcurve,
            other_pcurve: other_pcurve.clone(),
        }),
    }
}

pub(super) fn tolerance_for_boundary_authoring(brep: &BRep) -> f64 {
    brep.tolerance.linear().max(2.0e-4)
}

fn boundary_sides_for_segment_offset(segment_offset: usize) -> &'static [BoundarySide] {
    match segment_offset {
        0 => &[BoundarySide::UMin],
        _ => &[BoundarySide::UMax, BoundarySide::VMin, BoundarySide::VMax],
    }
}

struct BoundarySeamEvidence<'a> {
    brep: &'a BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    twin: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &'a SurfaceOnFace,
    twin_surface: &'a SurfaceOnFace,
    reference_curve: &'a BoundaryReferenceCurve<'a>,
    cap_pcurve: &'a Pcurve,
    other_pcurve: &'a Pcurve,
}

struct BoundaryReferenceCurve<'a> {
    curve: &'a dyn popkern_core::geometry::curve::CurveEval,
    geometry: CurveGeometry,
    t_start: f64,
    t_end: f64,
}

impl BoundaryReferenceCurve<'_> {
    fn curve_on_edge(&self) -> CurveOnEdge {
        CurveOnEdge {
            geometry: self.geometry.clone(),
            t_start: self.t_start,
            t_end: self.t_end,
        }
    }
}

fn boundary_reference_curve_for_candidate<'a>(
    brep: &'a BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    segment: &'a BoundarySegment,
) -> Option<BoundaryReferenceCurve<'a>> {
    use popkern_core::geometry::curve::CurveEval;

    if matches!(segment, BoundarySegment::FaceEdge { .. })
        && let Some(curve_on_edge) = brep
            .half_edges
            .get(he)
            .and_then(|half_edge| brep.edges.get(half_edge.edge))
            .and_then(|edge| edge.curve.as_ref())
    {
        return Some(BoundaryReferenceCurve {
            curve: &curve_on_edge.geometry,
            geometry: curve_on_edge.geometry.clone(),
            t_start: curve_on_edge.t_start,
            t_end: curve_on_edge.t_end,
        });
    }

    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t_start, t_end) = curve.parameter_range();
    Some(BoundaryReferenceCurve {
        curve,
        geometry: CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    })
}

fn boundary_owned_seam_disagreement(evidence: BoundarySeamEvidence<'_>) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;

    let t0 = evidence.reference_curve.t_start;
    let t1 = evidence.reference_curve.t_end;
    let p0 = evidence.reference_curve.curve.point_at(t0);
    let p1 = evidence.reference_curve.curve.point_at(t1);
    let he_origin = evidence.brep.vertices[evidence.brep.half_edges[evidence.he].origin].position;
    let twin_origin =
        evidence.brep.vertices[evidence.brep.half_edges[evidence.twin].origin].position;
    let cap_forward = he_origin.distance_to(&p0) <= he_origin.distance_to(&p1);
    let twin_forward = twin_origin.distance_to(&p0) <= twin_origin.distance_to(&p1);

    let mut worst = 0.0_f64;
    for i in 0..BOUNDARY_AUTHORING_SAMPLE_COUNT {
        let s = i as f64 / (BOUNDARY_AUTHORING_SAMPLE_COUNT - 1) as f64;
        let cap_curve_s = if cap_forward { s } else { 1.0 - s };
        let twin_curve_s = if twin_forward { s } else { 1.0 - s };
        let p_cap_ref = evidence
            .reference_curve
            .curve
            .point_at(t0 + cap_curve_s * (t1 - t0));
        let p_twin_ref = evidence
            .reference_curve
            .curve
            .point_at(t0 + twin_curve_s * (t1 - t0));
        let cap_uv = evidence.cap_pcurve.point_at(
            evidence.cap_pcurve.t_start
                + s * (evidence.cap_pcurve.t_end - evidence.cap_pcurve.t_start),
        );
        let other_uv = evidence.other_pcurve.point_at(
            evidence.other_pcurve.t_start
                + s * (evidence.other_pcurve.t_end - evidence.other_pcurve.t_start),
        );
        let p_cap = evidence.cap_surface.geometry.point_at(cap_uv.x, cap_uv.y);
        let p_other = evidence
            .twin_surface
            .geometry
            .point_at(other_uv.x, other_uv.y);
        worst = worst
            .max(p_cap.distance_to(&p_cap_ref))
            .max(p_other.distance_to(&p_twin_ref));
    }
    worst
}

fn canonical_gregory_boundary_pcurve_for_side(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &SurfaceOnFace,
    boundary_side: BoundarySide,
) -> Pcurve {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::surface::SurfaceEval;

    let (start_uv, end_uv) = match boundary_side {
        BoundarySide::UMin => (Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        BoundarySide::UMax => (Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        BoundarySide::VMin => (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        BoundarySide::VMax => (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        BoundarySide::Indexed(_) => return Pcurve::linear((0.0, 0.0), (0.0, 0.0), 0.0, 1.0),
    };
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let p_start = cap_surface.geometry.point_at(start_uv.x, start_uv.y);
    let p_end = cap_surface.geometry.point_at(end_uv.x, end_uv.y);
    let (start, end) = if he_origin.distance_to(&p_start) <= he_origin.distance_to(&p_end) {
        (start_uv, end_uv)
    } else {
        (end_uv, start_uv)
    };
    Pcurve::linear((start.x, start.y), (end.x, end.y), 0.0, 1.0)
}

fn boundary_pcurve_endpoint_gap(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &SurfaceOnFace,
    cap_pcurve: &Pcurve,
) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;

    let Some(half_edge) = brep.half_edges.get(he) else {
        return f64::INFINITY;
    };
    let Some(twin) = brep.half_edges.get(half_edge.twin) else {
        return f64::INFINITY;
    };
    let start = brep.vertices[half_edge.origin].position;
    let end = brep.vertices[twin.origin].position;
    let uv_start = cap_pcurve.point_at(cap_pcurve.t_start);
    let uv_end = cap_pcurve.point_at(cap_pcurve.t_end);
    let boundary_start = cap_surface.geometry.point_at(uv_start.x, uv_start.y);
    let boundary_end = cap_surface.geometry.point_at(uv_end.x, uv_end.y);
    start
        .distance_to(&boundary_start)
        .max(end.distance_to(&boundary_end))
}

fn half_edge_topology_endpoints(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    let half_edge = brep.half_edges.get(he)?;
    let next = brep.half_edges.get(half_edge.next)?;
    Some((
        brep.vertices.get(half_edge.origin)?.position,
        brep.vertices.get(next.origin)?.position,
    ))
}

fn project_reference_curve_to_half_edge_pcurve(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    surface: &SurfaceOnFace,
    reference_curve: &BoundaryReferenceCurve<'_>,
) -> Option<Pcurve> {
    use popkern_core::geometry::Point3;

    let t0 = reference_curve.t_start;
    let t1 = reference_curve.t_end;
    let p0 = reference_curve.curve.point_at(t0);
    let p1 = reference_curve.curve.point_at(t1);
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let forward = he_origin.distance_to(&p0) <= he_origin.distance_to(&p1);

    let n_samples = 12usize;
    let mut uv_points = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let p = reference_curve.curve.point_at(t0 + curve_s * (t1 - t0));
        let (u, v) = helpers::project_to_uv(&surface.geometry, &p);
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &params, 5)
        .ok()
        .map(Into::into)
}
