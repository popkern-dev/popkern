//! Retired repeated-apex residual topology handling.

use super::*;
use popkern_core::topology::{Face, TopologyRole};

pub(super) struct RetiredBaseResidualLoopPlan {
    arc_index: usize,
    apex: VertexHandle,
    sides: Vec<RetiredBaseResidualSide>,
}

struct RetiredBaseResidualSide {
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
}

impl RetiredBaseResidualLoopPlan {
    pub(super) fn author_loop(
        &self,
        brep: &mut BRep,
    ) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
        let mut loop_half_edges = Vec::with_capacity(self.sides.len() + 2);
        for (side_index, side) in self.sides.iter().enumerate() {
            let (edge, half_edge) = insert_boundary_half_edge(brep, side.start);
            brep.edges[edge].role = TopologyRole::ConstructionResidual;
            brep.edges[edge].curve = Some(side.curve.clone());
            if support_split_half_edge_end_vertex(brep, half_edge) != Some(side.end) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition durable base residual arc[{}] side[{}] \
                         authored {:?}->{:?} but did not resolve that endpoint",
                        self.arc_index, side_index, side.start, side.end
                    ),
                });
            }
            loop_half_edges.push(half_edge);
        }

        let (seam_edge, seam_a, seam_b) = insert_unwired_edge_pair(brep, self.apex, self.apex);
        brep.edges[seam_edge].role = TopologyRole::ConstructionResidual;
        loop_half_edges.push(seam_a);
        loop_half_edges.push(seam_b);
        validate_half_edge_cycle(
            brep,
            &loop_half_edges,
            "durable support repeated-apex base residual loop",
            FaceHandle::default(),
        )?;
        Ok(loop_half_edges)
    }
}

pub(super) fn plan_retired_base_residual_loop(
    brep: &BRep,
    span: &RetiredApexSpan,
) -> Result<RetiredBaseResidualLoopPlan, ApexBuildError> {
    let first = span
        .half_edges
        .first()
        .copied()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition base span arc[{}] has no residual sides",
                span.arc_index
            ),
        })?;
    let apex = brep
        .half_edges
        .get(first)
        .map(|data| data.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition base span arc[{}] starts with stale side {first:?}",
                span.arc_index
            ),
        })?;
    let mut sides = Vec::with_capacity(span.half_edges.len());
    for (side_index, half_edge) in span.half_edges.iter().copied().enumerate() {
        let half_edge_data =
            brep.half_edges
                .get(half_edge)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition base span arc[{}] side[{side_index}] \
                         references stale half-edge {half_edge:?}",
                        span.arc_index
                    ),
                })?;
        let end = support_split_half_edge_end_vertex(brep, half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition base span arc[{}] side[{side_index}] \
                     cannot resolve an endpoint for {half_edge:?}",
                    span.arc_index
                ),
            }
        })?;
        let curve = brep
            .edges
            .get(half_edge_data.edge)
            .and_then(|edge| edge.curve.clone())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition base span arc[{}] side[{side_index}] \
                     cannot preserve residual authoring data for {half_edge:?}: source edge has \
                     no curve",
                    span.arc_index
                ),
            })?;
        sides.push(RetiredBaseResidualSide {
            start: half_edge_data.origin,
            end,
            curve,
        });
    }
    Ok(RetiredBaseResidualLoopPlan {
        arc_index: span.arc_index,
        apex,
        sides,
    })
}

pub(super) fn residual_strip_loop_with_duplicated_seam(
    brep: &mut BRep,
    ownership: &FilletEndOwnership,
    shared_seam_override: Option<HalfEdgeHandle>,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    if let Some(shared) = shared_seam_override {
        return Ok(vec![
            shared,
            ownership.retired_apex_cell.t_a_to_apex,
            ownership.retired_apex_cell.apex_to_t_b,
        ]);
    }
    let source = shared_seam_override.unwrap_or(ownership.cylinder_candidate_half_edge);
    let shared = duplicate_boundary_half_edge_for_residual_strip(
        brep,
        source,
        ownership.fillet_index,
        false,
    )?;
    Ok(vec![
        shared,
        ownership.retired_apex_cell.t_a_to_apex,
        ownership.retired_apex_cell.apex_to_t_b,
    ])
}

pub(super) fn mark_loop_edges_construction_residual(
    brep: &mut BRep,
    half_edges: &[HalfEdgeHandle],
) {
    for half_edge in half_edges {
        if let Some(edge) = brep.half_edges.get(*half_edge).map(|data| data.edge)
            && let Some(edge_data) = brep.edges.get_mut(edge)
        {
            edge_data.role = TopologyRole::ConstructionResidual;
        }
    }
}

pub(super) fn demote_retired_apex_residual_faces(
    brep: &mut BRep,
    cap_face: FaceHandle,
    apex_vertex: VertexHandle,
    retired_spans: &[RetiredApexSpan],
    fillet_ends: &[FilletEndOwnership],
) {
    let mut residual_faces = std::collections::BTreeSet::new();
    for span in retired_spans {
        for half_edge in &span.half_edges {
            if let Some(face) = face_for_half_edge(brep, *half_edge) {
                residual_faces.insert(face);
            }
        }
    }
    for ownership in fillet_ends {
        for half_edge in [
            ownership.retired_apex_cell.t_a_to_apex,
            ownership.retired_apex_cell.apex_to_t_b,
        ] {
            if let Some(face) = face_for_half_edge(brep, half_edge) {
                residual_faces.insert(face);
            }
        }
    }
    let target_cylinder_faces = fillet_ends
        .iter()
        .map(|ownership| ownership.target_cylinder_face)
        .collect::<std::collections::BTreeSet<_>>();
    residual_faces.extend(brep.faces.iter().filter_map(|(face, face_data)| {
        if face_data.surface.is_none() || face == cap_face || target_cylinder_faces.contains(&face)
        {
            return None;
        }
        face_contains_vertex(brep, face, apex_vertex).then_some(face)
    }));
    residual_faces.extend(brep.faces.iter().filter_map(|(face, face_data)| {
        if face == cap_face || target_cylinder_faces.contains(&face) {
            return None;
        }
        (face_data.surface.is_none() && face_has_construction_boundary_edge(brep, face_data))
            .then_some(face)
    }));
    for &face in &residual_faces {
        if face == cap_face {
            continue;
        }
        if (face_contains_vertex(brep, face, apex_vertex)
            || brep
                .faces
                .get(face)
                .is_some_and(|face_data| face_has_construction_boundary_edge(brep, face_data)))
            && let Some(face_data) = brep.faces.get_mut(face)
        {
            face_data.role = TopologyRole::ConstructionResidual;
            face_data.surface = None;
        }
    }
}

pub(super) fn validate_existing_target_cylinder_strip(
    brep: &BRep,
    ownership: &FilletEndOwnership,
    apex_vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    let face = brep
        .faces
        .get(ownership.target_cylinder_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] selected stale target cylinder \
                 face {:?} as the authoritative strip",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        })?;
    let cycle = collect_loop_half_edges(brep, face.outer_loop)?;
    validate_half_edge_cycle(
        brep,
        &cycle,
        "authoritative target-cylinder strip loop",
        ownership.target_cylinder_face,
    )?;
    if cycle.iter().any(|half_edge| {
        brep.half_edges
            .get(*half_edge)
            .is_some_and(|data| data.origin == apex_vertex)
    }) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} still \
                 retains apex vertex {apex_vertex:?} after the cylinder trim rewrite",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        });
    }
    let expected = fillet_surfaces.get(ownership.fillet_index).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] has no fillet surface for \
                 target-cylinder strip validation",
                ownership.fillet_index
            ),
        }
    })?;
    let actual = face
        .surface
        .as_ref()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} has no \
                 production surface after the cylinder trim rewrite",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        })?;
    if std::mem::discriminant(&actual.geometry) != std::mem::discriminant(&expected.surface) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target face {:?} is not the \
                 expected cylinder strip surface",
                ownership.fillet_index, ownership.target_cylinder_face
            ),
        });
    }
    Ok(())
}

fn duplicate_boundary_half_edge_for_residual_strip(
    brep: &mut BRep,
    source: HalfEdgeHandle,
    fillet_index: usize,
    reverse: bool,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let source_data =
        brep.half_edges
            .get(source)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot duplicate \
                     stale residual strip seam {source:?}"
                ),
            })?;
    let curve = brep
        .edges
        .get(source_data.edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot duplicate \
                 residual strip seam {source:?}: source edge has no curve"
            ),
        })?;
    let start = if reverse {
        support_split_half_edge_end_vertex(brep, source).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot duplicate \
                     residual strip seam {source:?}: source endpoint is not recoverable"
                ),
            }
        })?
    } else {
        source_data.origin
    };
    let (edge, half_edge) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].role = TopologyRole::ConstructionResidual;
    brep.edges[edge].curve = Some(curve);
    Ok(half_edge)
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let loop_handle = brep.half_edges.get(half_edge)?.loop_;
    brep.loops.get(loop_handle).map(|loop_data| loop_data.face)
}

fn face_contains_vertex(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    brep.faces
        .get(face)
        .and_then(|face_data| collect_loop_half_edges(brep, face_data.outer_loop).ok())
        .is_some_and(|cycle| {
            cycle.iter().any(|half_edge| {
                brep.half_edges
                    .get(*half_edge)
                    .is_some_and(|data| data.origin == vertex)
            })
        })
}

fn face_has_construction_boundary_edge(brep: &BRep, face: &Face) -> bool {
    collect_loop_half_edges(brep, face.outer_loop)
        .ok()
        .is_some_and(|cycle| {
            cycle.iter().any(|half_edge| {
                brep.half_edges.get(*half_edge).is_some_and(|data| {
                    brep.boundary_edges.contains(&data.edge)
                        || brep
                            .half_edges
                            .get(data.twin)
                            .and_then(|twin| brep.loops.get(twin.loop_))
                            .is_none()
                })
            })
        })
}

/// Role coherence for retired subcomplexes: an edge whose EVERY
/// incident face is `ConstructionResidual` is itself residual. The
/// loop-wiring helper stamps edge roles from the owning face's role
/// DURING face creation, when the face still carries the default
/// (`Production`) role — the caller demotes the face immediately
/// afterwards, leaving residual-interior edges (e.g. the retired base
/// span's apex seam self-loop) mislabeled `Production` and visible to
/// production diagnostics (measured: pyramid4_zoo_mixed shipped a
/// zero-length "Production" edge at the apex).
pub(super) fn demote_residual_interior_edges(brep: &mut BRep) -> usize {
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    let mut demoted = 0usize;
    for eh in edge_handles {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        if edge.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let face_of = |h| {
            brep.half_edges
                .get(h)
                .and_then(|d| brep.loops.get(d.loop_))
                .map(|l| l.face)
        };
        let (Some(f_a), Some(f_b)) = (face_of(he), face_of(twin)) else {
            continue;
        };
        let residual = |f| {
            brep.faces
                .get(f)
                .is_some_and(|d| d.role == TopologyRole::ConstructionResidual)
        };
        if residual(f_a)
            && residual(f_b)
            && let Some(edge) = brep.edges.get_mut(eh)
        {
            edge.role = TopologyRole::ConstructionResidual;
            demoted += 1;
        }
    }
    demoted
}
