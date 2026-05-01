//! Phase 6c-A2.4 Gregory cap surface assembly.
//!
//! Bridges from A2's [`ApexPlan`] + post-stitch BRep state to the
//! legacy [`apply_three_quad_decomposition_to_blend`] entry point
//! (the K-P promotion ladder).  The bridge:
//!
//! 1. Walks `plan.target_topology.outer_loop_edge_indices` and emits
//!    one [`BoundarySegment`] per apex_edge — `FilletEnd` for
//!    [`TargetEdgeRole::CapVsCyl`], `FaceEdge` for
//!    [`TargetEdgeRole::CapVsPlane`].
//! 2. Resolves the `face` field of each `FaceEdge` by walking the
//!    BRep edge's twin half-edge.
//! 3. Calls the legacy K-P promotion ladder via
//!    `apply_three_quad_decomposition_to_blend`.
//!
//! Plus topological helpers used by `construct_apex_inner` to
//! identify the cap face after stage 3 and merge any unmerged
//! corner sub-faces left by the V-spoke kef sequence on Gregory
//! caps with mixed-host mef sites.
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md`
//! §A2.4.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::{TargetCapFace, TargetEdgeRole};
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::construct::three_quad_apply::{
    ThreeQuadFailure, apply_three_quad_decomposition_to_blend,
};
use crate::vertex_classifier::VertexBlendClassification;

/// **A2.4 Gregory cap surface assembly** — bridge from A2's plan +
/// post-stitch BRep state to legacy `three_quad_apply::
/// apply_three_quad_decomposition_to_blend`.
///
/// Builds a `BlendBoundary` from `plan.target_topology.apex_edges`
/// (each apex_edge becomes a `BoundarySegment::FilletEnd` or
/// `FaceEdge` based on its role + curve), then invokes the
/// existing K-P promotion ladder.  The ladder produces N or 4N
/// sub-quad Gregory patches whose surfaces are assigned to the
/// poked sub-faces.
pub(super) fn assemble_gregory_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
    apex_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    if fillet_surfaces.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: "Gregory cap surface assembly needs fillet_surfaces; caller \
                     used construct_apex (legacy entry) instead of \
                     construct_apex_with_surfaces"
                .into(),
        });
    }

    let outer_indices = match &plan.target_topology.cap_face {
        TargetCapFace::Gregory {
            outer_loop_edge_indices,
        } => outer_loop_edge_indices,
        other => {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "assemble_gregory_cap_surface called with non-Gregory cap_face: {other:?}"
                ),
            });
        }
    };
    if outer_indices.len() < 6 || !outer_indices.len().is_multiple_of(2) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory cap outer loop must have 2N >= 6 edges, got {}",
                outer_indices.len()
            ),
        });
    }
    let n = outer_indices.len() / 2;

    // Walk the planned 2N outer-loop edge indices and emit one
    // BoundarySegment per edge, alternating FilletEnd / FaceEdge.
    let mut segments: Vec<BoundarySegment> = Vec::with_capacity(2 * n);
    for (k, &edge_idx) in outer_indices.iter().enumerate() {
        let target_edge = plan
            .target_topology
            .apex_edges
            .get(edge_idx)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "outer_loop_edge_indices[{k}]={edge_idx} out of range for \
                         apex_edges (len {})",
                    plan.target_topology.apex_edges.len(),
                ),
            })?;
        let segment = match target_edge.role {
            TargetEdgeRole::CapVsCyl { fillet_idx } => {
                if k % 2 != 0 {
                    return Err(ApexBuildError::PlanMismatch {
                        reason: format!(
                            "Gregory cap loop expects FilletEnd at even index {k}, \
                                 found CapVsCyl at odd position"
                        ),
                    });
                }
                BoundarySegment::FilletEnd {
                    fillet_index: fillet_idx,
                    curve: target_edge.curve.clone(),
                }
            }
            TargetEdgeRole::CapVsPlane { .. } => {
                if k % 2 != 1 {
                    return Err(ApexBuildError::PlanMismatch {
                        reason: format!(
                            "Gregory cap loop expects FaceEdge at odd index {k}, \
                                 found CapVsPlane at even position"
                        ),
                    });
                }
                let edge = *apex_edge_handles.get(edge_idx).ok_or_else(|| {
                    ApexBuildError::PlanMismatch {
                        reason: format!(
                            "apex_edge_handles missing index {edge_idx} \
                                 (len {}) for FaceEdge bridge",
                            apex_edge_handles.len(),
                        ),
                    }
                })?;
                let face = find_face_across_edge(brep, edge, cap_face).ok_or_else(|| {
                    let edge_alive = brep.edges.contains_key(edge);
                    let adj = if edge_alive {
                        let primary = brep.edges[edge].half_edge;
                        let pa = brep
                            .half_edges
                            .get(primary)
                            .and_then(|he| brep.loops.get(he.loop_))
                            .map(|l| l.face);
                        let twin = brep.half_edges.get(primary).map(|he| he.twin);
                        let pb = twin
                            .and_then(|t| brep.half_edges.get(t))
                            .and_then(|he| brep.loops.get(he.loop_))
                            .map(|l| l.face);
                        format!("primary_face={pa:?}, twin_face={pb:?}, cap_face={cap_face:?}")
                    } else {
                        "edge no longer exists in BRep".to_string()
                    };
                    ApexBuildError::PlanMismatch {
                        reason: format!(
                            "FaceEdge[{k}] (apex_edge[{edge_idx}], BRep edge {edge:?}) \
                                 has no twin face on the non-cap side ({adj})"
                        ),
                    }
                })?;
                BoundarySegment::FaceEdge {
                    face,
                    curve: target_edge.curve.clone(),
                }
            }
            TargetEdgeRole::InnerSpoke => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory cap outer loop should not contain InnerSpoke at index {k}"
                    ),
                });
            }
        };
        segments.push(segment);
    }

    // Vertex position: centroid of the 2N apex vertices.  Used by
    // some BlendBoundary consumers (e.g. sphere-fit substitution)
    // but not by `apply_three_quad_decomposition_to_blend`'s core
    // path; centroid is a safe geometric default.
    let denom = apex_vertex_handles.len() as f64;
    let vertex_position = if denom > 0.0 {
        let mut sx = 0.0;
        let mut sy = 0.0;
        let mut sz = 0.0;
        for &vh in apex_vertex_handles {
            let p = brep.vertices[vh].position;
            sx += p.x;
            sy += p.y;
            sz += p.z;
        }
        Point3::new(sx / denom, sy / denom, sz / denom)
    } else {
        Point3::new(0.0, 0.0, 0.0)
    };

    let boundary = BlendBoundary {
        vertex: plan.apex_vertex,
        vertex_position,
        segments,
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };

    // Classification — A2 fires only at SameSideConvex apices today
    // (cube/pyramid mixed-radius corners).  A2.5 will lift this to
    // honor the planner's own classification once OppositeSide /
    // saddle handling lands.
    let classification = VertexBlendClassification::SameSideConvex;

    // Drive the legacy K-P promotion ladder with all rungs enabled.
    // The ladder picks the highest-quality rung that converges; for
    // cube_mixed N=3 the typical landing is biquintic-full or two-
    // stripe with K-P binonic promotion above.  These flags mirror
    // the legacy default-on dispatch when the ladder is fully
    // exercised; subsequent passes can plumb FilletOptions through
    // for user-driven knobs.
    let result = apply_three_quad_decomposition_to_blend(
        brep,
        cap_face,
        &boundary,
        fillet_surfaces,
        plan.apex_vertex,
        n,
        classification,
        true, // use_deg5
        true, // use_biquintic_full
        true, // use_biseptic
        true, // use_two_stripe
        true, // use_kp_binonic
        true, // use_kp_polar
    );

    match result {
        Ok(_patches) => Ok(()),
        Err(ThreeQuadFailure::Overlap) => Err(ApexBuildError::PlanMismatch {
            reason: "Gregory three-quad decomposition self-intersects (overlap remediation)".into(),
        }),
        Err(ThreeQuadFailure::Other(s)) => Err(ApexBuildError::PlanMismatch {
            reason: format!("Gregory three-quad decomposition failed: {s}"),
        }),
    }
}

/// Sweep across the BRep and `kef` residual edges that separate
/// cap-side corner sub-faces left over by stage 3's V-spoke kef
/// sequence.  Two passes:
///
/// 1. **Cross-cap-accumulator pass** — kef any edge whose two
///    adjacent faces both contain ≥2 apex_edges in their outer
///    loop (cap-side accumulators).  Trivially correct since
///    multi-apex-edge faces are by construction interior to the
///    cap region.
/// 2. **Bisected-spoke pass** — kef any edge that originates from
///    an apex_vertex (i.e., a vertex appearing on `apex_edge_
///    handles`'s endpoints) AND whose adjacent faces both
///    contain ≥1 apex_edge in their outer loop AND whose adjacent
///    faces are different.  These are the "leftover" V-spoke
///    halves that `resolve_apex_vertices`'s `split_edge` left
///    alive (V→apex_vertex was kef'd by stage 3, but apex_vertex→
///    perp_anchor remains, separating two corner sub-faces).
///
/// Iterates the combined sweep until no candidate boundary
/// remains.  Cap at `apex_edges.len() * 2 + 8` iterations.
pub(super) fn merge_unmerged_apex_subfaces(
    brep: &mut BRep,
    apex_edge_handles: &[EdgeHandle],
) -> Result<(), String> {
    let cap_iters = apex_edge_handles.len() * 2 + 8;
    let apex_set: std::collections::HashSet<EdgeHandle> =
        apex_edge_handles.iter().copied().collect();
    let apex_vertex_set = apex_vertex_set(brep, &apex_set);
    for _ in 0..cap_iters {
        // Pass 1: cross-cap-accumulator merge.
        let cap_side_faces = faces_with_multi_apex_edges(brep, &apex_set);
        let mut victim: Option<HalfEdgeHandle> = None;
        if cap_side_faces.len() >= 2 {
            victim = find_kef_candidate_between(brep, &cap_side_faces, &apex_set);
        }
        // Pass 2: bisected-spoke kef.  Walks every face containing
        // any apex_edge and looks for an outgoing non-apex edge from
        // an apex_vertex whose twin face is a different apex-
        // containing face.
        if victim.is_none() {
            victim = find_bisected_spoke_kef_candidate(brep, &apex_set, &apex_vertex_set);
        }
        let Some(he) = victim else {
            return Ok(());
        };
        brep.kef(he)
            .map_err(|e| format!("kef on cap-merging edge: {e:?}"))?;
    }
    Err("merge_unmerged_apex_subfaces: iteration cap hit".to_string())
}

/// Collect the set of vertices that appear as endpoints of any
/// apex_edge in `apex_set`.  These are the apex_vertices.
fn apex_vertex_set(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> std::collections::HashSet<VertexHandle> {
    let mut out = std::collections::HashSet::new();
    for &eh in apex_set {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let primary = edge.half_edge;
        if let Some(he) = brep.half_edges.get(primary) {
            out.insert(he.origin);
            if let Some(twin) = brep.half_edges.get(he.twin) {
                out.insert(twin.origin);
            }
        }
    }
    out
}

/// Pass 1 helper: find any non-apex edge between two distinct
/// cap-side accumulator faces.
fn find_kef_candidate_between(
    brep: &BRep,
    cap_side_faces: &[FaceHandle],
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> Option<HalfEdgeHandle> {
    for fh in cap_side_faces {
        let face = brep.faces.get(*fh)?;
        let loop_data = brep.loops.get(face.outer_loop)?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if !apex_set.contains(&he.edge) {
                let twin_loop = brep.half_edges.get(he.twin).map(|t| t.loop_);
                let twin_face = twin_loop.and_then(|tl| brep.loops.get(tl).map(|l| l.face));
                if let Some(tf) = twin_face
                    && tf != *fh
                    && cap_side_faces.contains(&tf)
                {
                    return Some(cursor);
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    None
}

/// Pass 2 helper: find a bisected-spoke edge.  Walks every face
/// containing an apex_edge and looks for a non-apex edge whose
/// origin is an apex_vertex AND whose twin face also contains at
/// least one apex_edge.  This catches the case where a corner
/// sub-face has only one apex_edge (so it's not a "cap-side
/// accumulator") but is still on the cap perimeter and must merge.
fn find_bisected_spoke_kef_candidate(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
    apex_vertex_set: &std::collections::HashSet<VertexHandle>,
) -> Option<HalfEdgeHandle> {
    let containing_faces: Vec<FaceHandle> = brep
        .faces
        .iter()
        .filter_map(|(fh, face)| {
            let loop_data = brep.loops.get(face.outer_loop)?;
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut iters = 0;
            let max_iters = brep.half_edges.len() + 1;
            while let Some(he) = brep.half_edges.get(cursor) {
                if apex_set.contains(&he.edge) {
                    return Some(fh);
                }
                cursor = he.next;
                iters += 1;
                if cursor == start || iters >= max_iters {
                    return None;
                }
            }
            None
        })
        .collect();
    for fh in &containing_faces {
        let face = brep.faces.get(*fh)?;
        let loop_data = brep.loops.get(face.outer_loop)?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            if !apex_set.contains(&he.edge) && apex_vertex_set.contains(&he.origin) {
                let twin = he.twin;
                let twin_origin = brep.half_edges.get(twin).map(|t| t.origin);
                let twin_loop = brep.half_edges.get(twin).map(|t| t.loop_);
                let twin_face = twin_loop.and_then(|tl| brep.loops.get(tl).map(|l| l.face));
                if let (Some(to), Some(tf)) = (twin_origin, twin_face)
                    && apex_vertex_set.contains(&to)
                    && tf != *fh
                    && containing_faces.contains(&tf)
                {
                    return Some(cursor);
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    None
}

fn faces_with_multi_apex_edges(
    brep: &BRep,
    apex_set: &std::collections::HashSet<EdgeHandle>,
) -> Vec<FaceHandle> {
    let mut out = Vec::new();
    for (fh, face) in &brep.faces {
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        let mut count = 0usize;
        while let Some(he) = brep.half_edges.get(cursor) {
            if apex_set.contains(&he.edge) {
                count += 1;
                if count >= 2 {
                    out.push(fh);
                    break;
                }
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
    out
}

/// Find the face whose outer loop walks every edge in
/// `apex_edge_handles`.  Returns `None` when no such face exists.
///
/// The cap face is uniquely determined by the property that all
/// 2N (or N) plan apex_edges are on its outer loop in the planned
/// cap-loop order.  After stage 3 (kef + kev), exactly one such
/// face exists in a well-formed BRep.  This helper identifies it
/// independently of which side `kef_cross_face_v_edges` happened
/// to keep — the legacy helper's "blend_face" return chases the
/// last kef's `he_face`, which on Gregory caps with mixed-host
/// mef sites can be a non-cap remnant (e.g. a far-side plane sub-
/// face whose interior was unrelated to the apex region).
pub(super) fn find_face_bounded_by_all_apex_edges(
    brep: &BRep,
    apex_edge_handles: &[EdgeHandle],
) -> Option<FaceHandle> {
    if apex_edge_handles.is_empty() {
        return None;
    }
    for (fh, face) in &brep.faces {
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut loop_edges: Vec<EdgeHandle> = Vec::new();
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            loop_edges.push(he.edge);
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if apex_edge_handles.iter().all(|eh| loop_edges.contains(eh)) {
            return Some(fh);
        }
    }
    None
}

/// Given an edge handle whose two adjacent faces are `cap_face` and
/// some other face, return the other face.  Returns `None` when
/// the edge isn't manifold-2 or `cap_face` doesn't bound it.
fn find_face_across_edge(
    brep: &BRep,
    edge: EdgeHandle,
    cap_face: FaceHandle,
) -> Option<FaceHandle> {
    let edge_data = brep.edges.get(edge)?;
    let primary = edge_data.half_edge;
    let primary_he = brep.half_edges.get(primary)?;
    let twin_he = brep.half_edges.get(primary_he.twin)?;
    let face_a = brep.loops.get(primary_he.loop_)?.face;
    let face_b = brep.loops.get(twin_he.loop_)?.face;
    if face_a == cap_face {
        Some(face_b)
    } else if face_b == cap_face {
        Some(face_a)
    } else {
        None
    }
}
