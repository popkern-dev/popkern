//! Curve refresh helpers for repeated-apex partition endpoint adoption.

use super::*;

pub(super) fn refresh_fillet_end_curves_after_endpoint_adoption(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_ends: &[FilletEndOwnership],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), ApexBuildError> {
    for ownership in fillet_ends {
        let sector_index = plan
            .target_topology
            .apex_edges
            .chunks(2)
            .position(|chunk| {
                chunk.first().is_some_and(|edge| {
                    matches!(
                        edge.role,
                        crate::apex_construction::target_topology::TargetEdgeRole::CapVsCyl { fillet_idx }
                            if fillet_idx == ownership.fillet_index
                    )
                })
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot refresh FilletEnd[{}] curve after \
                     endpoint adoption: no matching target edge",
                    ownership.fillet_index
                ),
            })?;
        let target_edge_index = 2 * sector_index;
        let curve = &plan.target_topology.apex_edges[target_edge_index].curve;
        attach_fillet_end_curve_to_edge(
            brep,
            ownership.edge,
            curve,
            ownership.fillet_index,
            fillet_surfaces,
        )?;
    }
    Ok(())
}
