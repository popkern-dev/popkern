use super::*;

pub(super) fn materialized_shared_fillet_end_cylinder_half_edge(
    brep: &mut BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    topology_plan: &CylinderTrimRewriteTopologyPlan,
    start: VertexHandle,
    end: VertexHandle,
    reuse_shared_cap_cylinder_seam: bool,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let existing = existing_shared_fillet_end_half_edge_for_segment(brep, rebuild_plan, start, end);
    if existing.is_none()
        || (!reuse_shared_cap_cylinder_seam
            && matches!(
                topology_plan.seam_rewrite,
                FilletEndSeamRewrite::ReuseCurrentTwinPair { .. }
            ))
    {
        let source =
            existing.unwrap_or_else(|| seam_rewrite_current_cylinder_half_edge(topology_plan));
        let original_edge = seam_rewrite_edge(topology_plan);
        let curve = brep
            .edges
            .get(original_edge)
            .and_then(|edge| edge.curve.clone())
            .map(Ok)
            .unwrap_or_else(|| {
                exact_cylinder_arc_trim_curve_between(
                    brep,
                    rebuild_plan.cylinder_face,
                    start,
                    end,
                    "duplicated shared FilletEnd seam",
                    rebuild_plan.fillet_index,
                )
            })?;
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(curve);
        if support_split_half_edge_end_vertex(brep, half_edge) != Some(end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] duplicated shared cylinder \
                     seam has wrong orientation: expected {start:?}->{end:?}, got {half_edge:?}",
                    rebuild_plan.fillet_index
                ),
            });
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/cylinder-loop-replace] FilletEnd[{}] duplicated shared \
                 cylinder seam {:?}->{:?}: original {:?}, cylinder boundary {:?}",
                rebuild_plan.fillet_index, start, end, source, half_edge
            );
        }
        return Ok(half_edge);
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/cylinder-loop-replace] FilletEnd[{}] reuses shared \
             cap/cylinder seam {:?}->{:?}: {:?} ({})",
            rebuild_plan.fillet_index,
            start,
            end,
            existing,
            topology_plan.seam_rewrite.diagnostic_label()
        );
    }
    existing.ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{}] expected reusable shared \
             cap/cylinder seam {:?}->{:?}, but no semantic half-edge was found",
            rebuild_plan.fillet_index, start, end
        ),
    })
}

fn existing_shared_fillet_end_half_edge_for_segment(
    brep: &BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let (_edge, cap_half_edge, cylinder_half_edge) =
        rebuild_plan.shared_fillet_end_half_edges().ok()?;
    if let Some(recorded) = [cap_half_edge, cylinder_half_edge]
        .into_iter()
        .find(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .is_some_and(|data| data.origin == start)
                && half_edge_end_vertex(brep, *half_edge)
                    .is_some_and(|actual_end| actual_end == end)
        })
    {
        return Some(recorded);
    }
    semantic_shared_fillet_end_cylinder_half_edge(brep, rebuild_plan, start, end)
}

fn seam_rewrite_current_cylinder_half_edge(
    topology_plan: &CylinderTrimRewriteTopologyPlan,
) -> HalfEdgeHandle {
    match topology_plan.seam_rewrite {
        FilletEndSeamRewrite::ReuseCurrentTwinPair {
            current_cylinder_half_edge,
            ..
        }
        | FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
            current_cylinder_half_edge,
            ..
        } => current_cylinder_half_edge,
    }
}

fn seam_rewrite_edge(topology_plan: &CylinderTrimRewriteTopologyPlan) -> EdgeHandle {
    match topology_plan.seam_rewrite {
        FilletEndSeamRewrite::ReuseCurrentTwinPair { edge, .. }
        | FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide { edge, .. } => edge,
    }
}

fn semantic_shared_fillet_end_cylinder_half_edge(
    brep: &BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let face = brep.faces.get(rebuild_plan.cylinder_face)?;
    collect_loop_half_edges(brep, face.outer_loop)
        .ok()?
        .into_iter()
        .find(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .is_some_and(|data| data.origin == start)
                && half_edge_end_vertex(brep, *half_edge)
                    .is_some_and(|actual_end| actual_end == end)
        })
}
