use super::*;

pub(super) fn materialized_shared_fillet_end_cylinder_half_edge(
    brep: &mut BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    topology_plan: &CylinderTrimRewriteTopologyPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let existing =
        existing_shared_fillet_end_half_edge_for_segment(brep, rebuild_plan, start, end)?;
    if !matches!(
        topology_plan.seam_rewrite,
        FilletEndSeamRewrite::ReuseCurrentTwinPair { .. }
    ) {
        return Ok(existing);
    }
    let original_edge = brep
        .half_edges
        .get(existing)
        .map(|half_edge| half_edge.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot duplicate stale shared \
                 cylinder seam {existing:?}",
                rebuild_plan.fillet_index
            ),
        })?;
    let curve = brep
        .edges
        .get(original_edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot duplicate shared cylinder \
                 seam {existing:?}: original edge {original_edge:?} has no curve",
                rebuild_plan.fillet_index
            ),
        })?;
    let (edge, half_edge) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    if support_split_half_edge_end_vertex(brep, half_edge) != Some(end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] duplicated shared cylinder seam \
                 has wrong orientation: expected {start:?}->{end:?}, got {half_edge:?}",
                rebuild_plan.fillet_index
            ),
        });
    }
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace §14.9/cylinder-loop-replace] FilletEnd[{}] duplicated shared cylinder \
             seam {:?}->{:?}: original {:?}, cylinder boundary {:?}",
            rebuild_plan.fillet_index, start, end, existing, half_edge
        );
    }
    Ok(half_edge)
}

fn existing_shared_fillet_end_half_edge_for_segment(
    brep: &BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let (_edge, cap_half_edge, cylinder_half_edge) = rebuild_plan.shared_fillet_end_half_edges()?;
    [cap_half_edge, cylinder_half_edge]
        .into_iter()
        .find(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .is_some_and(|data| data.origin == start)
                && half_edge_end_vertex(brep, *half_edge)
                    .is_some_and(|actual_end| actual_end == end)
        })
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot find existing shared \
                 FilletEnd half-edge {:?}->{:?} while replacing target-cylinder loop",
                rebuild_plan.fillet_index, start, end
            ),
        })
}
