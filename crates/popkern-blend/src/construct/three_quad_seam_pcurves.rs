//! Boundary-owned seam pcurve authoring for three-quad cap faces.

#[path = "three_quad_seam_pcurve_candidate.rs"]
mod three_quad_seam_pcurve_candidate;
#[path = "three_quad_support_seams.rs"]
mod three_quad_support_seams;
pub(crate) use crate::types::SupportSeamRelocationReport;

use super::three_quad_face_matching::{
    face_has_cylinder_surface, face_surface_matches, half_edge_endpoints_match_curve,
};
use super::*;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::topology::FaceHandle;
use popkern_core::topology::{
    SharedTrimEdgeFailureReason, SharedTrimEdgeInput, SharedTrimEdgeValidation,
    SurfaceBoundaryAuditEvent, SurfaceBoundaryAuditOutcome, SurfaceBoundarySharedTrimFailureReason,
    SurfaceBoundarySharedTrimValidation, commit_shared_trim_edge_transactionally,
};
use popkern_core::topology::{SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason};
use three_quad_seam_pcurve_candidate::{
    BoundaryOwnedCandidate, boundary_owned_candidate, tolerance_for_boundary_authoring,
};

pub(super) fn author_boundary_owned_seam_pcurves(
    brep: &mut BRep,
    sub_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    subface_index: usize,
    semantic_cap_index: Option<usize>,
) -> Vec<crate::types::ApexConstructionDiagnostic> {
    let mut diagnostics = Vec::new();
    let segment_cap_index = semantic_cap_index.unwrap_or(subface_index);
    if boundary.segments.len() < 2 * (segment_cap_index + 1) {
        return diagnostics;
    }
    for (segment_offset, segment_index) in [2 * segment_cap_index, 2 * segment_cap_index + 1]
        .into_iter()
        .enumerate()
    {
        let Some(segment) = boundary.segments.get(segment_index) else {
            continue;
        };
        let candidates = subface_half_edge_candidates_for_boundary_segment(
            brep,
            sub_face,
            segment,
            fillet_surfaces,
        );
        let Some(cap_surface) = brep.faces.get(sub_face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let mut best_decline = None;
        let Some((he, boundary_side, curve, cap_pcurve, other_pcurve)) =
            candidates.into_iter().find_map(|he| {
                match boundary_owned_candidate(brep, he, cap_surface, segment, segment_offset) {
                    BoundaryOwnedCandidate::AcceptedBoundary {
                        half_edge,
                        boundary_side,
                        curve,
                        cap_pcurve,
                        other_pcurve,
                    } => Some((
                        half_edge,
                        Some(boundary_side),
                        curve,
                        cap_pcurve,
                        other_pcurve,
                    )),
                    BoundaryOwnedCandidate::AcceptedInterior {
                        half_edge,
                        curve,
                        cap_pcurve,
                        other_pcurve,
                    } => Some((half_edge, None, curve, cap_pcurve, other_pcurve)),
                    BoundaryOwnedCandidate::Declined {
                        half_edge,
                        max_disagreement,
                        boundary_side,
                        endpoint_gap,
                    } => {
                        if brep.surface_boundary_hint(half_edge).is_none() {
                            match best_decline {
                                Some((_, best_disagreement, _, _))
                                    if max_disagreement >= best_disagreement => {}
                                _ => {
                                    best_decline = Some((
                                        half_edge,
                                        max_disagreement,
                                        boundary_side,
                                        endpoint_gap,
                                    ))
                                }
                            }
                        }
                        None
                    }
                    BoundaryOwnedCandidate::Unavailable => None,
                }
            })
        else {
            if let Some((half_edge, max_disagreement, boundary_side, endpoint_gap)) = best_decline
                && brep.surface_boundary_hint(half_edge).is_none()
            {
                diagnostics.push(
                    crate::construct::three_quad_boundary_pcurve_diagnostics::boundary_owned_pcurve_decline_diagnostic(
                        brep,
                        half_edge,
                        segment,
                        segment_index,
                        max_disagreement,
                        endpoint_gap,
                        tolerance_for_boundary_authoring(brep),
                    ),
                );
                brep.set_surface_boundary_decline(
                    half_edge,
                    SurfaceBoundaryHintDecline {
                        reason: if endpoint_gap > tolerance_for_boundary_authoring(brep) {
                            SurfaceBoundaryHintDeclineReason::AuthoredHintEndpointMismatch
                        } else {
                            SurfaceBoundaryHintDeclineReason::NotOnCanonicalBoundary
                        },
                        best_boundary_deviation: max_disagreement,
                        best_boundary: Some(boundary_side),
                        endpoint_gap,
                    },
                );
            }
            if crate::direct_apex_trace::cap_g0_enabled() {
                eprintln!(
                    "[direct-apex-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                     segment[{segment_index}] no validated half-edge"
                );
            }
            continue;
        };
        let twin = brep.half_edges[he].twin;
        let edge = brep.half_edges[he].edge;
        let curve = orient_curve_to_canonical_edge(brep, edge, curve);
        let commit_result = commit_shared_trim_edge_transactionally(
            brep,
            SharedTrimEdgeInput {
                edge,
                half_edge_a: he,
                half_edge_b: twin,
                curve,
                pcurve_a: cap_pcurve,
                pcurve_b: other_pcurve,
                tolerance: tolerance_for_boundary_authoring(brep),
                validate_brep: false,
                provenance: "boundary_owned_seam_pcurves",
            },
        );
        match commit_result {
            Ok(_) => record_boundary_owned_audit(
                brep,
                BoundaryOwnedAuditInput {
                    half_edge: he,
                    outcome: SurfaceBoundaryAuditOutcome::SharedTrimAccepted,
                    boundary: boundary_side,
                    tolerance: Some(tolerance_for_boundary_authoring(brep)),
                    ..BoundaryOwnedAuditInput::default()
                },
            ),
            Err(failure) => {
                record_boundary_owned_audit(
                    brep,
                    BoundaryOwnedAuditInput {
                        half_edge: he,
                        outcome: SurfaceBoundaryAuditOutcome::SharedTrimDeclined,
                        boundary: boundary_side,
                        shared_trim_failure: Some(surface_boundary_shared_trim_failure_reason(
                            failure.reason,
                        )),
                        shared_trim_validation: Some(surface_boundary_shared_trim_validation(
                            failure.validation,
                        )),
                        tolerance: Some(tolerance_for_boundary_authoring(brep)),
                    },
                );
                continue;
            }
        }
        if !brep.half_edges.contains_key(he) {
            continue;
        }
        if let Some(boundary_side) = boundary_side {
            brep.set_surface_boundary_hint(he, boundary_side);
        }
        if crate::direct_apex_trace::cap_g0_enabled() {
            eprintln!(
                "[direct-apex-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                 segment[{segment_index}] wrote he={he:?} twin={twin:?} edge={:?}",
                brep.half_edges[he].edge
            );
        }
    }
    diagnostics
}

struct BoundaryOwnedAuditInput {
    half_edge: popkern_core::topology::HalfEdgeHandle,
    outcome: SurfaceBoundaryAuditOutcome,
    boundary: Option<popkern_core::geometry::edge_triple::BoundarySide>,
    shared_trim_failure: Option<SurfaceBoundarySharedTrimFailureReason>,
    shared_trim_validation: Option<SurfaceBoundarySharedTrimValidation>,
    tolerance: Option<f64>,
}

impl Default for BoundaryOwnedAuditInput {
    fn default() -> Self {
        Self {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            outcome: SurfaceBoundaryAuditOutcome::ClassifiedDecline,
            boundary: None,
            shared_trim_failure: None,
            shared_trim_validation: None,
            tolerance: None,
        }
    }
}

fn record_boundary_owned_audit(brep: &mut BRep, input: BoundaryOwnedAuditInput) {
    let Some(half_edge_data) = brep.half_edges.get(input.half_edge) else {
        return;
    };
    let edge = half_edge_data.edge;
    let face = brep
        .loops
        .get(half_edge_data.loop_)
        .map(|loop_data| loop_data.face);
    let twin_face = brep
        .half_edges
        .get(half_edge_data.twin)
        .and_then(|twin| brep.loops.get(twin.loop_))
        .map(|loop_data| loop_data.face);
    brep.record_surface_boundary_audit_event(SurfaceBoundaryAuditEvent {
        half_edge: input.half_edge,
        edge: Some(edge),
        face,
        twin_face,
        outcome: input.outcome,
        boundary: input.boundary,
        decline_reason: None,
        shared_trim_failure: input.shared_trim_failure,
        shared_trim_validation: input.shared_trim_validation,
        tolerance: input.tolerance,
        provenance: "boundary_owned_seam_pcurves".to_string(),
    });
}

fn surface_boundary_shared_trim_failure_reason(
    reason: SharedTrimEdgeFailureReason,
) -> SurfaceBoundarySharedTrimFailureReason {
    match reason {
        SharedTrimEdgeFailureReason::MissingTopology => {
            SurfaceBoundarySharedTrimFailureReason::MissingTopology
        }
        SharedTrimEdgeFailureReason::HalfEdgeEdgeMismatch => {
            SurfaceBoundarySharedTrimFailureReason::HalfEdgeEdgeMismatch
        }
        SharedTrimEdgeFailureReason::HalfEdgesNotTwins => {
            SurfaceBoundarySharedTrimFailureReason::HalfEdgesNotTwins
        }
        SharedTrimEdgeFailureReason::EndpointMismatch => {
            SurfaceBoundarySharedTrimFailureReason::EndpointMismatch
        }
        SharedTrimEdgeFailureReason::CurveSurfaceMismatch => {
            SurfaceBoundarySharedTrimFailureReason::CurveSurfaceMismatch
        }
        SharedTrimEdgeFailureReason::SideAgreementMismatch => {
            SurfaceBoundarySharedTrimFailureReason::SideAgreementMismatch
        }
        SharedTrimEdgeFailureReason::ValidationFailed => {
            SurfaceBoundarySharedTrimFailureReason::ValidationFailed
        }
    }
}

fn surface_boundary_shared_trim_validation(
    validation: SharedTrimEdgeValidation,
) -> SurfaceBoundarySharedTrimValidation {
    SurfaceBoundarySharedTrimValidation {
        max_curve_surface_deviation: validation.max_curve_surface_deviation,
        max_side_agreement_deviation: validation.max_side_agreement_deviation,
        max_endpoint_gap: validation.max_endpoint_gap,
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

pub(super) fn author_post_assignment_support_face_boundary_seams(
    brep: &mut BRep,
    sub_faces: &[FaceHandle],
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    semantic_cap_indices: &[Option<usize>],
) -> SupportSeamRelocationReport {
    let before = brep.clone();
    let mut report = SupportSeamRelocationReport::default();
    // Already-authored no-op (PK-234): when EVERY bridge segment is
    // already carried by some sub-face's ring — true by construction
    // after the saddle boundary-ring carve — relocation can only drag
    // exact junction vertices onto mis-mapped targets (measured: the
    // per-sub-face semantic mapping lands rotated on the carved ring
    // and all six ring vertices moved 0.066–0.222 toward neighboring
    // junctions). The mitred classes this pass serves sit ~4e-2 off
    // their segments and still author.
    let all_bridges_carried = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FaceEdge { .. }))
        .all(|segment| {
            sub_faces.iter().any(|sub_face| {
                three_quad_support_seams::segment_already_authored(brep, *sub_face, segment)
            })
        });
    if all_bridges_carried {
        return report;
    }
    let required_hints = gregory_cap_cylinder_boundary_hint_count(brep);
    let required_neighbors = three_quad_support_seams::gregory_cap_cylinder_neighbor_count(brep);
    for (subface_index, sub_face) in sub_faces.iter().enumerate() {
        let segment_cap_index = semantic_cap_indices
            .get(subface_index)
            .and_then(|cap_index| *cap_index)
            .unwrap_or(subface_index);
        for (segment_offset, segment_index) in [2 * segment_cap_index, 2 * segment_cap_index + 1]
            .into_iter()
            .enumerate()
        {
            let Some(segment @ BoundarySegment::FaceEdge { .. }) =
                boundary.segments.get(segment_index)
            else {
                continue;
            };
            let candidates = subface_half_edge_candidates_for_boundary_segment(
                brep,
                *sub_face,
                segment,
                fillet_surfaces,
            );
            let result = three_quad_support_seams::try_author_support_face_boundary_seam(
                brep,
                *sub_face,
                segment_offset,
                &candidates,
            );
            if let Some(authored_seam) = result.authored {
                brep.set_surface_boundary_hint(
                    authored_seam.half_edge,
                    authored_seam.boundary_side,
                );
            }
            report.merge(result.report);
        }
    }
    let hints_preserved = gregory_cap_cylinder_boundary_hint_count(brep) >= required_hints;
    let neighbors_preserved =
        three_quad_support_seams::gregory_cap_cylinder_neighbor_count(brep) >= required_neighbors;
    let loops_closed = brep.validate_loop_closure().is_ok();
    if !hints_preserved || !neighbors_preserved || !loops_closed {
        *brep = before;
        report.rollback_reason = Some(if !loops_closed {
            crate::types::SupportSeamRollbackReason::LoopClosureFailed
        } else if !hints_preserved {
            crate::types::SupportSeamRollbackReason::CapCylinderHintsLost
        } else {
            crate::types::SupportSeamRollbackReason::CapCylinderNeighborsLost
        });
        report.authored = 0;
    }
    report
}

fn gregory_cap_cylinder_boundary_hint_count(brep: &BRep) -> usize {
    let mut hinted = 0usize;
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
            let half_edge = &brep.half_edges[cursor];
            let twin_face = brep.loops[brep.half_edges[half_edge.twin].loop_].face;
            if face_has_cylinder_surface(brep, twin_face)
                && brep.surface_boundary_hint(cursor).is_some()
            {
                hinted += 1;
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
    }
    hinted
}

fn subface_half_edge_candidates_for_boundary_segment(
    brep: &BRep,
    sub_face: FaceHandle,
    segment: &BoundarySegment,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Vec<popkern_core::topology::HalfEdgeHandle> {
    let mut primary = Vec::new();
    let mut fallback = Vec::new();
    let Some(face_data) = brep.faces.get(sub_face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..(brep.half_edges.len() + 1) {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let Some(twin_half_edge) = brep.half_edges.get(he.twin) else {
            break;
        };
        let Some(twin_loop) = brep.loops.get(twin_half_edge.loop_) else {
            break;
        };
        let twin_face = twin_loop.face;
        let matches = match segment {
            BoundarySegment::FaceEdge { face, curve } => {
                twin_face == *face || half_edge_endpoints_match_curve(brep, cursor, curve, 1.0e-5)
            }
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => {
                half_edge_endpoints_match_curve(brep, cursor, curve, 1.0e-5)
                    || fillet_surfaces.get(*fillet_index).is_some_and(|fillet| {
                        face_surface_matches(brep, twin_face, &fillet.surface)
                    })
                    || face_has_cylinder_surface(brep, twin_face)
            }
        };
        if matches {
            primary.push(cursor);
        } else if matches!(segment, BoundarySegment::FaceEdge { .. }) {
            fallback.push(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    primary.extend(fallback);
    primary
}
