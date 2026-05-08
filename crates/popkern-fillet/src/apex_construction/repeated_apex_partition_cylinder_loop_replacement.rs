//! Materialized target-cylinder loop replacement for repeated-apex rewrites.

use super::*;

pub(super) fn apply_materialized_cylinder_loop_replacements(
    brep: &mut BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    subcomplex_plan: &TrimRewriteSubcomplexPlan,
    pre_materialized: &[PreMaterializedTrimSide],
) -> Result<(), ApexBuildError> {
    let mut all_support_owned_trim_sides = Vec::new();
    let mut pending_cylinder_loop_replacements = Vec::new();
    for (rebuild_plan, topology_plan) in rebuild_plans.iter().zip(topology_plans) {
        let allowed_support_faces = target_cylinder_support_faces_with_semantic(
            brep,
            rebuild_plan.cylinder_face,
            cap_face,
            face_edges,
        );
        let Some(fragment) =
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)
        else {
            continue;
        };
        if !fragment.materialized_anchor_replacement {
            continue;
        }
        let mut cycle = Vec::with_capacity(topology_plan.loop_candidate.segments.len());
        let mut support_owned_trim_sides = Vec::new();
        for segment in &topology_plan.loop_candidate.segments {
            if segment.label.contains("shared FilletEnd") {
                let half_edge = materialized_shared_fillet_end_cylinder_half_edge(
                    brep,
                    rebuild_plan,
                    topology_plan,
                    segment.start,
                    segment.end,
                )?;
                cycle.push(half_edge);
            } else if segment.label.contains("far-side") {
                let trim_sides = if let Some(trim_sides) = pre_materialized_trim_sides(
                    pre_materialized,
                    rebuild_plan.cylinder_face,
                    rebuild_plan.fillet_index,
                    segment.start,
                    segment.end,
                ) {
                    trim_sides
                } else {
                    materialized_far_side_trim_sides(
                        brep,
                        rebuild_plan.cylinder_face,
                        segment.start,
                        segment.end,
                        segment.label,
                        rebuild_plan.fillet_index,
                        &allowed_support_faces,
                    )?
                };
                support_owned_trim_sides.extend(trim_sides.iter().cloned());
                cycle.extend(trim_sides.iter().map(|trim_side| trim_side.half_edge));
            } else {
                let planned_support_face = subcomplex_plan.planned_generator_support_face(
                    topology_plan.cylinder_face,
                    topology_plan.fillet_index,
                    segment.start,
                    segment.end,
                    segment.label,
                )?;
                let allowed_planned_support =
                    std::collections::BTreeSet::from([planned_support_face]);
                let support_face = support_face_for_generator_trim_side(
                    brep,
                    GeneratorTrimSupportOwnerQuery {
                        cylinder_face: rebuild_plan.cylinder_face,
                        cap_face,
                        start: segment.start,
                        end: segment.end,
                        label: segment.label,
                        fillet_index: rebuild_plan.fillet_index,
                        allowed_support_faces: &allowed_planned_support,
                    },
                )?;
                let (half_edge, actual_support_face, support_side_owned, support_source) =
                    materialized_generator_trim_side(
                        brep,
                        rebuild_plan.cylinder_face,
                        support_face,
                        segment.start,
                        segment.end,
                        segment.label,
                        rebuild_plan.fillet_index,
                    )?;
                cycle.push(half_edge);
                support_owned_trim_sides.push(SupportOwnedTrimSide {
                    cylinder_face: topology_plan.cylinder_face,
                    support_face: actual_support_face,
                    fillet_index: topology_plan.fillet_index,
                    half_edge,
                    support_side_owned,
                    support_source,
                    requires_support_ownership: true,
                });
            }
        }
        validate_half_edge_cycle(
            brep,
            &cycle,
            "materialized target-cylinder replacement loop",
            topology_plan.cylinder_face,
        )?;
        all_support_owned_trim_sides.extend(support_owned_trim_sides);
        pending_cylinder_loop_replacements.push(PendingCylinderLoopReplacement {
            fillet_index: topology_plan.fillet_index,
            cylinder_face: topology_plan.cylinder_face,
            kind: topology_plan.loop_candidate.kind,
            cycle,
        });
    }
    for replacement in pending_cylinder_loop_replacements {
        let loop_handle = brep
            .faces
            .get(replacement.cylinder_face)
            .map(|face| face.outer_loop)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot replace stale target \
                     cylinder face {:?}",
                    replacement.fillet_index, replacement.cylinder_face
                ),
            })?;
        validate_half_edge_cycle(
            brep,
            &replacement.cycle,
            "materialized target-cylinder replacement loop",
            replacement.cylinder_face,
        )?;
        wire_existing_loop(
            brep,
            loop_handle,
            replacement.cylinder_face,
            &replacement.cycle,
        )?;
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/cylinder-loop-replace] FilletEnd[{}] {:?}/{:?} -> {:?}",
                replacement.fillet_index,
                replacement.cylinder_face,
                replacement.kind,
                replacement.cycle
            );
        }
    }
    materialize_closed_support_trim_cycles(brep, face_edges, &all_support_owned_trim_sides)?;
    for trim_side in all_support_owned_trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        validate_replacement_loop_support_side_ownership(
            brep,
            trim_side.cylinder_face,
            trim_side.fillet_index,
            &[trim_side.half_edge],
        )?;
    }
    Ok(())
}
