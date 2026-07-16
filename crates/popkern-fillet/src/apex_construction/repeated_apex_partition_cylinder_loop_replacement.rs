//! Materialized target-cylinder loop replacement for repeated-apex rewrites.

use super::*;

pub(super) struct MaterializedCylinderLoopReplacementRequest<'a> {
    pub(super) cap_face: FaceHandle,
    pub(super) face_edges: &'a [FaceEdgeOwnership],
    pub(super) rebuild_plans: &'a [SemanticCylinderTrimRebuildPlan],
    pub(super) topology_plans: &'a [CylinderTrimRewriteTopologyPlan],
    pub(super) subcomplex_plan: &'a TrimRewriteSubcomplexPlan,
    pub(super) pre_materialized: &'a [PreMaterializedTrimSide],
    pub(super) side_support_bridge_loops: &'a [LoopHandle],
    pub(super) include_selected_nonmaterialized: bool,
    pub(super) reuse_shared_cap_cylinder_seam: bool,
    pub(super) support_ownership_already_materialized: bool,
}

pub(super) fn apply_materialized_cylinder_loop_replacements(
    brep: &mut BRep,
    request: MaterializedCylinderLoopReplacementRequest<'_>,
) -> Result<Vec<PendingCylinderLoopReplacement>, ApexBuildError> {
    let cap_face = request.cap_face;
    let face_edges = request.face_edges;
    let rebuild_plans = request.rebuild_plans;
    let topology_plans = request.topology_plans;
    let subcomplex_plan = request.subcomplex_plan;
    let pre_materialized = request.pre_materialized;
    let side_support_bridge_loops = request.side_support_bridge_loops;
    let include_selected_nonmaterialized = request.include_selected_nonmaterialized;
    let reuse_shared_cap_cylinder_seam = request.reuse_shared_cap_cylinder_seam;
    let support_ownership_already_materialized = request.support_ownership_already_materialized;

    let mut all_support_owned_trim_sides = Vec::new();
    let mut residual_boundary_spans = Vec::new();
    let mut existing_residual_sides = Vec::new();
    let mut pending_cylinder_loop_replacements = Vec::new();
    if crate::direct_apex_trace::enabled() && !side_support_bridge_loops.is_empty() {
        eprintln!(
            "[direct-apex-trace §14.9/cylinder-loop-replace] side-support continuation bridge loops {:?}",
            side_support_bridge_loops
        );
    }
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
        if !fragment.materialized_anchor_replacement && !include_selected_nonmaterialized {
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
                    reuse_shared_cap_cylinder_seam,
                )?;
                push_if_boundary_span(brep, half_edge, &mut residual_boundary_spans);
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
                for trim_side in &trim_sides {
                    if !trim_side.requires_support_ownership {
                        push_if_boundary_span(
                            brep,
                            trim_side.half_edge,
                            &mut residual_boundary_spans,
                        );
                        push_if_unwired_residual_twin(
                            brep,
                            trim_side.half_edge,
                            &mut existing_residual_sides,
                        );
                    }
                    residual_boundary_spans
                        .extend(trim_side.authored_boundary_spans.iter().copied());
                }
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
                let materialized = materialized_generator_trim_side(
                    brep,
                    GeneratorTrimMaterializationRequest {
                        cylinder_face: rebuild_plan.cylinder_face,
                        support_face,
                        start: segment.start,
                        end: segment.end,
                        label: segment.label,
                        fillet_index: rebuild_plan.fillet_index,
                        side_support_bridge_loops,
                        skip_interior_patch: reuse_shared_cap_cylinder_seam,
                    },
                )?;
                residual_boundary_spans
                    .extend(materialized.authored_boundary_spans.iter().copied());
                cycle.push(materialized.half_edge);
                support_owned_trim_sides.push(SupportOwnedTrimSide {
                    cylinder_face: topology_plan.cylinder_face,
                    support_face: materialized.support_face,
                    fillet_index: topology_plan.fillet_index,
                    half_edge: materialized.half_edge,
                    support_side_owned: materialized.support_side_owned,
                    support_source: materialized.support_source,
                    authored_boundary_spans: materialized.authored_boundary_spans,
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
    if support_ownership_already_materialized {
        let _ =
            repair_support_source_loops_after_cylinder_rewire(brep, &all_support_owned_trim_sides)?;
    }
    for replacement in &pending_cylinder_loop_replacements {
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
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/cylinder-loop-replace] FilletEnd[{}] {:?}/{:?} -> {:?}",
                replacement.fillet_index,
                replacement.cylinder_face,
                replacement.kind,
                replacement.cycle
            );
        }
    }
    if !support_ownership_already_materialized {
        materialize_closed_support_trim_cycles(
            brep,
            face_edges,
            &all_support_owned_trim_sides,
            reuse_shared_cap_cylinder_seam,
        )?;
    }
    materialize_authored_support_boundary_residual_cycles(
        brep,
        &residual_boundary_spans,
        &existing_residual_sides,
    )?;
    repair_stale_shell_anchors(brep)?;
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
    Ok(pending_cylinder_loop_replacements)
}

fn push_if_boundary_span(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    boundary_spans: &mut Vec<HalfEdgeHandle>,
) {
    let Some(edge) = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.edge)
    else {
        return;
    };
    if brep.boundary_edges.contains(&edge) {
        boundary_spans.push(half_edge);
    }
}

fn push_if_unwired_residual_twin(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    residual_sides: &mut Vec<HalfEdgeHandle>,
) {
    let Some(twin) = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.twin)
    else {
        return;
    };
    if brep
        .half_edges
        .get(twin)
        .is_some_and(|twin_data| twin_data.loop_ == LoopHandle::default())
    {
        residual_sides.push(twin);
    }
}
