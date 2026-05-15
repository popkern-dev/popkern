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
use popkern_core::topology::FaceHandle;
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
) {
    if boundary.segments.len() < 2 * (subface_index + 1) {
        return;
    }
    for (segment_offset, segment_index) in [2 * subface_index, 2 * subface_index + 1]
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
        let Some((he, boundary_side, cap_pcurve, other_pcurve)) =
            candidates.into_iter().find_map(|he| {
                match boundary_owned_candidate(brep, he, cap_surface, segment, segment_offset) {
                    BoundaryOwnedCandidate::AcceptedBoundary {
                        half_edge,
                        boundary_side,
                        cap_pcurve,
                        other_pcurve,
                    } => Some((half_edge, Some(boundary_side), cap_pcurve, other_pcurve)),
                    BoundaryOwnedCandidate::AcceptedInterior {
                        half_edge,
                        cap_pcurve,
                        other_pcurve,
                    } => Some((half_edge, None, cap_pcurve, other_pcurve)),
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
            if crate::a2_trace::cap_g0_enabled() {
                eprintln!(
                    "[a2-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                     segment[{segment_index}] no validated half-edge"
                );
            }
            continue;
        };
        let twin = brep.half_edges[he].twin;
        if let Some(boundary_side) = boundary_side {
            brep.set_surface_boundary_hint(he, boundary_side);
        }
        brep.half_edges[he].pcurve = Some(cap_pcurve);
        brep.half_edges[twin].pcurve = Some(other_pcurve);
        if crate::a2_trace::cap_g0_enabled() {
            eprintln!(
                "[a2-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                 segment[{segment_index}] wrote he={he:?} twin={twin:?} edge={:?}",
                brep.half_edges[he].edge
            );
        }
    }
}

pub(super) fn author_post_assignment_support_face_boundary_seams(
    brep: &mut BRep,
    sub_faces: &[FaceHandle],
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> SupportSeamRelocationReport {
    let before = brep.clone();
    let mut report = SupportSeamRelocationReport::default();
    let required_hints = gregory_cap_cylinder_boundary_hint_count(brep);
    let required_neighbors = three_quad_support_seams::gregory_cap_cylinder_neighbor_count(brep);
    for (subface_index, sub_face) in sub_faces.iter().enumerate() {
        for (segment_offset, segment_index) in [2 * subface_index, 2 * subface_index + 1]
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
