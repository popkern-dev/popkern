//! Transactional support-face seam authoring for three-quad Gregory caps.

use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::types::{SupportSeamDecision, SupportSeamDeclineReason, SupportSeamRelocationReport};

#[path = "three_quad_support_seam_relocation.rs"]
mod relocation;

pub(super) struct AuthoredSupportSeam {
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) boundary_side: BoundarySide,
}

pub(super) struct SupportSeamAuthoringResult {
    pub(super) authored: Option<AuthoredSupportSeam>,
    pub(super) report: SupportSeamRelocationReport,
}

pub(super) fn try_author_support_face_boundary_seam(
    brep: &mut BRep,
    cap_face: FaceHandle,
    segment_offset: usize,
    candidates: &[HalfEdgeHandle],
) -> SupportSeamAuthoringResult {
    let mut report = SupportSeamRelocationReport::default();
    let Some(cap_surface) = brep
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    else {
        return SupportSeamAuthoringResult {
            authored: None,
            report,
        };
    };
    if !matches!(cap_surface, SurfaceGeometry::Gregory(_)) {
        return SupportSeamAuthoringResult {
            authored: None,
            report,
        };
    }
    try_author_support_face_boundary_seam_with_relocation(
        brep,
        cap_face,
        segment_offset,
        candidates,
        cap_surface,
        &mut report,
    )
}

fn try_author_support_face_boundary_seam_with_relocation(
    brep: &mut BRep,
    cap_face: FaceHandle,
    segment_offset: usize,
    candidates: &[HalfEdgeHandle],
    cap_surface: SurfaceGeometry,
    report: &mut SupportSeamRelocationReport,
) -> SupportSeamAuthoringResult {
    let min_cylinder_neighbors = gregory_cap_cylinder_neighbor_count(brep);
    for &half_edge in candidates {
        let Some((support_face, edge)) = support_plane_face_and_edge(brep, cap_face, half_edge)
        else {
            trace_decline(half_edge, "not a planar support candidate");
            report.push(SupportSeamDecision {
                half_edge,
                cap_face,
                support_face: None,
                boundary_side: None,
                accepted: false,
                reason: Some(SupportSeamDeclineReason::NotPlanarSupportCandidate),
                max_endpoint_motion: 0.0,
                max_sampled_disagreement: 0.0,
            });
            continue;
        };
        for boundary_side in support_boundary_sides(segment_offset) {
            let outcome = relocation::try_author_one(relocation::SupportSeamAttempt {
                brep,
                half_edge,
                cap_face,
                support_face,
                edge,
                cap_surface: &cap_surface,
                boundary_side: *boundary_side,
                min_cylinder_neighbors,
            });
            report.push(SupportSeamDecision {
                half_edge,
                cap_face,
                support_face: Some(support_face),
                boundary_side: Some(*boundary_side),
                accepted: outcome.accepted,
                reason: outcome.reason,
                max_endpoint_motion: outcome.max_endpoint_motion,
                max_sampled_disagreement: outcome.max_sampled_disagreement,
            });
            if outcome.accepted {
                return SupportSeamAuthoringResult {
                    authored: Some(AuthoredSupportSeam {
                        half_edge,
                        boundary_side: *boundary_side,
                    }),
                    report: std::mem::take(report),
                };
            }
        }
    }
    SupportSeamAuthoringResult {
        authored: None,
        report: std::mem::take(report),
    }
}

fn support_plane_face_and_edge(
    brep: &BRep,
    cap_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, EdgeHandle)> {
    let he = brep.half_edges.get(half_edge)?;
    let twin = he.twin;
    let support_face = brep.loops.get(brep.half_edges.get(twin)?.loop_)?.face;
    if support_face == cap_face {
        return None;
    }
    let support_surface = brep.faces.get(support_face)?.surface.as_ref()?;
    if !matches!(support_surface.geometry, SurfaceGeometry::Plane { .. }) {
        return None;
    }
    Some((support_face, he.edge))
}

fn support_boundary_sides(segment_offset: usize) -> &'static [BoundarySide] {
    match segment_offset {
        0 => &[],
        _ => &[BoundarySide::UMax, BoundarySide::VMin, BoundarySide::VMax],
    }
}

pub(super) fn gregory_cap_cylinder_neighbor_count(brep: &BRep) -> usize {
    let mut count = 0usize;
    for (_, face) in &brep.faces {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Gregory(_)) {
            continue;
        }
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(half_edge) = brep.half_edges.get(cursor) else {
                break;
            };
            let twin = half_edge.twin;
            let Some(twin_face) = brep
                .half_edges
                .get(twin)
                .and_then(|he| brep.loops.get(he.loop_))
                .map(|loop_| loop_.face)
            else {
                break;
            };
            if brep
                .faces
                .get(twin_face)
                .and_then(|face| face.surface.as_ref())
                .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
            {
                count += 1;
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
    }
    count
}

fn cap_face_preserves_cylinder_neighbor(brep: &BRep, cap_face: FaceHandle) -> bool {
    let Some(face) = brep.faces.get(cap_face) else {
        return false;
    };
    let Some(sof) = face.surface.as_ref() else {
        return false;
    };
    if !matches!(sof.geometry, SurfaceGeometry::Gregory(_)) {
        return false;
    }
    let Some(loop_data) = brep.loops.get(face.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            return false;
        };
        let twin = half_edge.twin;
        let Some(twin_face) = brep
            .half_edges
            .get(twin)
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|loop_| loop_.face)
        else {
            return false;
        };
        if brep
            .faces
            .get(twin_face)
            .and_then(|face| face.surface.as_ref())
            .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
        {
            return true;
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    false
}

fn trace_decline(half_edge: HalfEdgeHandle, reason: &str) {
    if crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!("[direct-apex-trace cap-g0] support seam declined he={half_edge:?}: {reason}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn support_seam_report_merges_counts_and_extrema() {
        let mut report = SupportSeamRelocationReport::default();
        report.push(SupportSeamDecision {
            half_edge: HalfEdgeHandle::default(),
            cap_face: FaceHandle::default(),
            support_face: None,
            boundary_side: None,
            accepted: false,
            reason: Some(SupportSeamDeclineReason::NotPlanarSupportCandidate),
            max_endpoint_motion: 0.01,
            max_sampled_disagreement: 0.02,
        });

        let mut other = SupportSeamRelocationReport::default();
        other.push(SupportSeamDecision {
            half_edge: HalfEdgeHandle::default(),
            cap_face: FaceHandle::default(),
            support_face: Some(FaceHandle::default()),
            boundary_side: Some(BoundarySide::VMax),
            accepted: true,
            reason: None,
            max_endpoint_motion: 0.03,
            max_sampled_disagreement: 0.004,
        });

        report.merge(other);

        assert_eq!(report.authored, 1);
        assert_eq!(report.declined, 1);
        assert_eq!(report.candidates, 2);
        assert_eq!(report.decisions.len(), 2);
        assert_eq!(report.max_endpoint_motion, 0.03);
        assert_eq!(report.max_sampled_disagreement, 0.02);
    }
}

/// PK-234 idempotence gate for the post-assignment support-seam
/// relocation: `true` when `sub_face`'s outer ring already carries an
/// edge lying ON the boundary segment's curve within the
/// boundary-authoring tolerance (2e-4 — the seam re-anchor's band).
/// True by construction after the saddle boundary-ring carve; the
/// mitred classes the relocation serves sit ~4e-2 off and still
/// author.
pub(super) fn segment_already_authored(
    brep: &popkern_core::topology::BRep,
    sub_face: popkern_core::topology::FaceHandle,
    segment: &crate::blend::boundary::BoundarySegment,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    let curve = match segment {
        crate::blend::boundary::BoundarySegment::FaceEdge { curve, .. }
        | crate::blend::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
    };
    let (t0, t1) = curve.parameter_range();
    let Some(face) = brep.faces.get(sub_face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..64 {
        let Some(hd) = brep.half_edges.get(cursor) else {
            return false;
        };
        if let Some(edge_curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
            let eval: &dyn CurveEval = &edge_curve.geometry;
            let (e0, e1) = (edge_curve.t_start, edge_curve.t_end);
            // Worst distance from 9 segment samples to the edge
            // curve (65-sample scan + windowed refinement).
            let mut worst = 0.0_f64;
            for i in 0..=8 {
                let p = curve.point_at(t0 + (t1 - t0) * i as f64 / 8.0);
                let mut d_min = f64::INFINITY;
                let mut t_min = e0;
                for j in 0..=64 {
                    let t = e0 + (e1 - e0) * j as f64 / 64.0;
                    let d = eval.point_at(t).distance_to(&p);
                    if d < d_min {
                        d_min = d;
                        t_min = t;
                    }
                }
                let mut window = (e1 - e0).abs() / 64.0;
                for _ in 0..16 {
                    let (lo, hi) = (
                        (t_min - window).max(e0.min(e1)),
                        (t_min + window).min(e0.max(e1)),
                    );
                    for j in 0..=8 {
                        let t = lo + (hi - lo) * j as f64 / 8.0;
                        let d = eval.point_at(t).distance_to(&p);
                        if d < d_min {
                            d_min = d;
                            t_min = t;
                        }
                    }
                    window *= 0.4;
                }
                worst = worst.max(d_min);
            }
            if worst <= 2.0e-4 {
                return true;
            }
        }
        cursor = hd.next;
        if cursor == start {
            break;
        }
    }
    false
}
