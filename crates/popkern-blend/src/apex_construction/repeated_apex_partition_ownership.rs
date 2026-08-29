use super::super::constructive_topology;
use super::*;

#[derive(Clone)]
pub(super) struct FaceEdgeOwnership {
    pub(super) fillet_index: usize,
    pub(super) edge: EdgeHandle,
    pub(super) apex_to_t_b: HalfEdgeHandle,
    pub(super) face_edge_cap_half_edge: HalfEdgeHandle,
    pub(super) support_half_edge: HalfEdgeHandle,
    pub(super) support_face: FaceHandle,
    pub(super) t_a_next_to_apex: HalfEdgeHandle,
}

impl FaceEdgeOwnership {
    pub(super) fn from_retired_span(
        brep: &BRep,
        span: &RetiredApexSpan,
        cap_face: FaceHandle,
        fillet_index: usize,
        n_fillets: usize,
    ) -> Result<Self, ApexBuildError> {
        if span.half_edges.len() != 3 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition expected FaceEdge arc[{fillet_index}] to \
                     have 3 half-edges, got {}",
                    span.half_edges.len()
                ),
            });
        }
        let face_edge_cap_half_edge = span.half_edges[1];
        let face = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .and_then(|half_edge| {
                brep.loops
                    .get(half_edge.loop_)
                    .map(|loop_data| loop_data.face)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     middle half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        if face != cap_face {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] middle \
                     half-edge is on {face:?}, expected cap face {cap_face:?}"
                ),
            });
        }
        let edge = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     cap half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        let support_half_edge = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     cap half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        let support_face = brep
            .half_edges
            .get(support_half_edge)
            .and_then(|half_edge| {
                brep.loops
                    .get(half_edge.loop_)
                    .map(|loop_data| loop_data.face)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     support half-edge {support_half_edge:?}"
                ),
            })?;
        if support_face == cap_face {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] support side \
                     is still on cap face {cap_face:?}"
                ),
            });
        }
        if fillet_index >= n_fillets {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc fillet index {fillet_index} \
                     out of range for {n_fillets} fillets"
                ),
            });
        }
        Ok(Self {
            fillet_index,
            edge,
            apex_to_t_b: span.half_edges[0],
            face_edge_cap_half_edge,
            support_half_edge,
            support_face,
            t_a_next_to_apex: span.half_edges[2],
        })
    }
}

pub(super) struct FilletEndOwnership {
    pub(super) fillet_index: usize,
    pub(super) edge: EdgeHandle,
    pub(super) cap_half_edge: HalfEdgeHandle,
    pub(super) cylinder_candidate_half_edge: HalfEdgeHandle,
    pub(super) target_cylinder_face: FaceHandle,
    pub(super) retired_apex_cell: RetiredApexCell,
}

pub(super) struct RetiredApexCell {
    pub(super) previous_face_edge_index: usize,
    pub(super) next_face_edge_index: usize,
    pub(super) t_a_to_apex: HalfEdgeHandle,
    pub(super) apex_to_t_b: HalfEdgeHandle,
}

pub(super) fn validate_partition_ownership_records(
    brep: &BRep,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
    n_corners: usize,
) -> Result<(), ApexBuildError> {
    for face_edge in face_edges {
        validate_face_edge_ownership_record(brep, face_edge)?;
    }
    for fillet_end in fillet_ends {
        validate_fillet_end_ownership_record(brep, fillet_end, n_corners)?;
    }
    Ok(())
}

fn validate_face_edge_ownership_record(
    brep: &BRep,
    face_edge: &FaceEdgeOwnership,
) -> Result<(), ApexBuildError> {
    let cap_twin = brep
        .half_edges
        .get(face_edge.face_edge_cap_half_edge)
        .map(|half_edge| half_edge.twin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] cap half-edge {:?} is stale",
                face_edge.fillet_index, face_edge.face_edge_cap_half_edge
            ),
        })?;
    if cap_twin != face_edge.support_half_edge {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] support half-edge mismatch: \
                 recorded {:?}, actual {:?}",
                face_edge.fillet_index, face_edge.support_half_edge, cap_twin
            ),
        });
    }
    let support_face = brep
        .half_edges
        .get(face_edge.support_half_edge)
        .and_then(|half_edge| {
            brep.loops
                .get(half_edge.loop_)
                .map(|loop_data| loop_data.face)
        })
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] support half-edge {:?} is stale",
                face_edge.fillet_index, face_edge.support_half_edge
            ),
        })?;
    if support_face != face_edge.support_face {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] support face changed from \
                 {:?} to {:?}",
                face_edge.fillet_index, face_edge.support_face, support_face
            ),
        });
    }
    let cap_edge = brep
        .half_edges
        .get(face_edge.face_edge_cap_half_edge)
        .map(|half_edge| half_edge.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] cap half-edge {:?} is stale",
                face_edge.fillet_index, face_edge.face_edge_cap_half_edge
            ),
        })?;
    if cap_edge != face_edge.edge {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FaceEdge[{}] edge changed from {:?} to {:?}",
                face_edge.fillet_index, face_edge.edge, cap_edge
            ),
        });
    }
    Ok(())
}

fn validate_fillet_end_ownership_record(
    brep: &BRep,
    fillet_end: &FilletEndOwnership,
    n_corners: usize,
) -> Result<(), ApexBuildError> {
    let cap_start = brep
        .half_edges
        .get(fillet_end.cap_half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                fillet_end.fillet_index, fillet_end.cap_half_edge
            ),
        })?;
    let cap_end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has stale twin",
                fillet_end.fillet_index, fillet_end.cap_half_edge
            ),
        }
    })?;
    let cap_twin = brep
        .half_edges
        .get(fillet_end.cap_half_edge)
        .map(|half_edge| half_edge.twin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                fillet_end.fillet_index, fillet_end.cap_half_edge
            ),
        })?;
    if cap_twin != fillet_end.cylinder_candidate_half_edge {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cylinder candidate mismatch: \
                 recorded {:?}, actual {:?}",
                fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge, cap_twin
            ),
        });
    }
    let cap_edge = brep
        .half_edges
        .get(fillet_end.cap_half_edge)
        .map(|half_edge| half_edge.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                fillet_end.fillet_index, fillet_end.cap_half_edge
            ),
        })?;
    if cap_edge != fillet_end.edge {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] edge changed from {:?} to {:?}",
                fillet_end.fillet_index, fillet_end.edge, cap_edge
            ),
        });
    }
    if !brep
        .faces
        .get(fillet_end.target_cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} \
                 is not a live cylinder face",
                fillet_end.fillet_index, fillet_end.target_cylinder_face
            ),
        });
    }
    let cylinder_start = brep
        .half_edges
        .get(fillet_end.cylinder_candidate_half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cylinder candidate {:?} is stale",
                fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge
            ),
        })?;
    let cylinder_end = half_edge_end_vertex(brep, fillet_end.cylinder_candidate_half_edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cylinder candidate {:?} has stale twin",
                fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge
            ),
        })?;
    if cylinder_start != cap_end || cylinder_end != cap_start {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cylinder candidate is not \
                 the reverse of the cap half-edge",
                fillet_end.fillet_index
            ),
        });
    }
    validate_retired_apex_cell_shape(brep, fillet_end, cap_start, cap_end)?;
    // Cell indices are CORNER indices (welded corners have no
    // FaceEdge but are still corners — PK-219).
    if fillet_end.retired_apex_cell.previous_face_edge_index >= n_corners
        || fillet_end.retired_apex_cell.next_face_edge_index >= n_corners
    {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] has invalid retired apex \
                 cell corner indices {} and {} for {n_corners} corners",
                fillet_end.fillet_index,
                fillet_end.retired_apex_cell.previous_face_edge_index,
                fillet_end.retired_apex_cell.next_face_edge_index,
            ),
        });
    }
    Ok(())
}

fn validate_retired_apex_cell_shape(
    brep: &BRep,
    fillet_end: &FilletEndOwnership,
    cap_start: VertexHandle,
    cap_end: VertexHandle,
) -> Result<(), ApexBuildError> {
    let retired_a_start = brep
        .half_edges
        .get(fillet_end.retired_apex_cell.t_a_to_apex)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] retired T_a->apex half-edge \
                 {:?} is stale",
                fillet_end.fillet_index, fillet_end.retired_apex_cell.t_a_to_apex
            ),
        })?;
    let retired_apex = half_edge_end_vertex(brep, fillet_end.retired_apex_cell.t_a_to_apex)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] retired T_a->apex half-edge \
                 {:?} has stale twin",
                fillet_end.fillet_index, fillet_end.retired_apex_cell.t_a_to_apex
            ),
        })?;
    let retired_b_start = brep
        .half_edges
        .get(fillet_end.retired_apex_cell.apex_to_t_b)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] retired apex->T_b half-edge \
                 {:?} is stale",
                fillet_end.fillet_index, fillet_end.retired_apex_cell.apex_to_t_b
            ),
        })?;
    let retired_b_end = half_edge_end_vertex(brep, fillet_end.retired_apex_cell.apex_to_t_b)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] retired apex->T_b half-edge \
                 {:?} has stale twin",
                fillet_end.fillet_index, fillet_end.retired_apex_cell.apex_to_t_b
            ),
        })?;
    if retired_a_start != cap_start || retired_apex != retired_b_start || retired_b_end != cap_end {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] retired apex cell is not \
                 T_a -> apex -> T_b for cap endpoints {:?}->{:?}",
                fillet_end.fillet_index, cap_start, cap_end
            ),
        });
    }
    Ok(())
}

pub(super) fn validate_target_cylinder_loop_ownership_ready(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
        let cap_start = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        let cap_end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has \
                     stale twin",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            }
        })?;
        let Some(vertices) = constructive_topology::collect_outer_loop_vertices(
            brep,
            fillet_end.target_cylinder_face,
        ) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} \
                     has no walkable outer loop",
                    fillet_end.fillet_index, fillet_end.target_cylinder_face
                ),
            });
        };
        if !vertices.contains(&cap_start) || !vertices.contains(&cap_end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot yet author shared \
                     cap/cylinder ownership: target cylinder face {:?} outer loop does not \
                     contain both cap endpoints {:?}->{:?}",
                    fillet_end.fillet_index, fillet_end.target_cylinder_face, cap_start, cap_end
                ),
            });
        }
    }
    Ok(())
}
