//! Curve refresh for cap-boundary rewrites that reauthor shared FilletEnd seams.

use super::*;

pub(super) fn refresh_reauthored_fillet_end_curves_for_plan(
    required: bool,
    brep: &mut BRep,
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> Result<(), ApexBuildError> {
    if !required {
        return Ok(());
    }
    for plan in topology_plans {
        let FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
            edge,
            current_cylinder_half_edge,
            ..
        } = plan.seam_rewrite
        else {
            continue;
        };
        let start = brep
            .half_edges
            .get(current_cylinder_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot refresh reauthored FilletEnd[{}] \
                     curve: stale cap-side half-edge {current_cylinder_half_edge:?}",
                    plan.fillet_index
                ),
            })?;
        let end = half_edge_end_vertex(brep, current_cylinder_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot refresh reauthored FilletEnd[{}] \
                     curve: half-edge {current_cylinder_half_edge:?} has stale twin",
                    plan.fillet_index
                ),
            }
        })?;
        let curve = exact_cylinder_arc_trim_curve_between(
            brep,
            plan.cylinder_face,
            start,
            end,
            "reauthored shared FilletEnd seam",
            plan.fillet_index,
        )
        .or_else(|_| {
            exact_cylinder_trim_curve_between(
                brep,
                plan.cylinder_face,
                start,
                end,
                "reauthored shared FilletEnd seam",
                plan.fillet_index,
            )
        })?;
        brep.edges[edge].curve = Some(curve);
    }
    Ok(())
}
