//! Transactional target-cylinder trim rewrite orchestration.

use super::*;

pub(super) fn validate_cylinder_trim_rebuild_plans(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    let rebuild_plans = plan_semantic_cylinder_trim_rebuilds(brep, fillet_ends)?;
    for rebuild_plan in rebuild_plans {
        if rebuild_plan.current_cylinder_order_is_constructive() {
            continue;
        }
        let reverse_status = rebuild_plan.reverse_path_status();
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cylinder trim rebuild plan \
                 found generator anchors {:?}->{:?}, but the cylinder half-edge orientation \
                 would keep off-cylinder path vertices {:?}; {reverse_status}; semantic \
                 rebuild segments required: {}",
                rebuild_plan.fillet_index,
                rebuild_plan.start_anchor,
                rebuild_plan.end_anchor,
                rebuild_plan.forward_off_surface_vertices(),
                rebuild_plan.segment_summary()
            ),
        });
    }
    Ok(())
}

pub(super) fn validate_cylinder_trim_rewrite_readiness(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    let rebuild_plans = plan_semantic_cylinder_trim_rebuilds(brep, fillet_ends)?;
    for rebuild_plan in rebuild_plans {
        rebuild_plan.validate_rewrite_readiness()?;
    }
    Ok(())
}

pub(super) fn apply_cylinder_trim_rewrites_transactionally(
    brep: &mut BRep,
    cap_face: FaceHandle,
    cap_loop: LoopHandle,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
    require_v_free_repeated_apex: bool,
) -> Result<CylinderTrimRewriteResult, ApexBuildError> {
    let rollback_snapshot = brep.clone();
    let result = (|| {
        materialize_missing_generator_anchors(brep, cap_face, fillet_ends, face_edges)?;
        let rebuild_plans = plan_semantic_cylinder_trim_rebuilds(brep, fillet_ends)?;
        for rebuild_plan in &rebuild_plans {
            rebuild_plan.validate_rewrite_readiness()?;
            rebuild_plan.validate_replacement_curve_authoring_ready(brep)?;
        }
        let topology_plans = rebuild_plans
            .iter()
            .map(|plan| plan.plan_cylinder_trim_rewrite_topology(brep))
            .collect::<Result<Vec<_>, _>>()?;
        let (include_selected_nonmaterialized, diagnostics) =
            selected_far_side_trim_materialization_preflight(
                brep,
                cap_face,
                face_edges,
                &rebuild_plans,
                &topology_plans,
            );
        let subcomplex_plan = plan_trim_subcomplex_for_cylinder_rewrite(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            include_selected_nonmaterialized,
        )?;
        validate_cylinder_trim_rewrite_topology_preflight(CylinderTrimRewritePreflight {
            brep,
            cap_face,
            face_edges,
            rebuild_plans: &rebuild_plans,
            topology_plans: &topology_plans,
            subcomplex_plan: &subcomplex_plan,
            include_selected_nonmaterialized,
            require_v_free_repeated_apex,
        })?;
        let pre_materialized = pre_materialize_far_side_trim_sides(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            include_selected_nonmaterialized,
        )?;
        let side_support_result = apply_selected_side_support_continuation(brep, &subcomplex_plan)?;
        let support_ownership_already_materialized = matches!(
            side_support_result,
            SideSupportRegionContinuationOutcome::Applied(_)
        );
        let side_support_bridge_loops = side_support_result.bridge_loops();
        let cylinder_loop_replacements = apply_materialized_cylinder_loop_replacements(
            brep,
            MaterializedCylinderLoopReplacementRequest {
                cap_face,
                face_edges,
                rebuild_plans: &rebuild_plans,
                topology_plans: &topology_plans,
                subcomplex_plan: &subcomplex_plan,
                pre_materialized: &pre_materialized,
                side_support_bridge_loops,
                include_selected_nonmaterialized,
                reuse_shared_cap_cylinder_seam: require_v_free_repeated_apex,
                support_ownership_already_materialized,
            },
        )?;
        let cap_boundary_plan = plan_cap_boundary_rewrite_for_cylinder_trim(
            brep,
            cap_face,
            cap_loop,
            face_edges,
            fillet_ends,
            &topology_plans,
        )?;
        cap_boundary_plan.validate_rewrite_postconditions(brep, face_edges, &topology_plans)?;
        cap_boundary_plan.validate_adjacent_loop_rewrite_readiness()?;
        let mut boundary_outcome = cap_boundary_plan.boundary_outcome(&topology_plans);
        if require_v_free_repeated_apex {
            for replacement in &cylinder_loop_replacements {
                boundary_outcome.mark_existing_target_cylinder_strip(replacement.fillet_index);
            }
        }
        if !require_v_free_repeated_apex || cylinder_loop_replacements.len() != fillet_ends.len() {
            boundary_outcome.plan_strip_loop_overrides(brep, fillet_ends)?;
        }
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/cylinder-trim-rewrite] topology_plans=[{}]; side_support_continuation={}; {}; {}",
                topology_plans
                    .iter()
                    .map(CylinderTrimRewriteTopologyPlan::diagnostic_label)
                    .collect::<Vec<_>>()
                    .join("; "),
                side_support_result.diagnostic_label(),
                cap_boundary_plan.diagnostic_label(),
                boundary_outcome.diagnostic_label()
            );
        }
        cap_boundary_plan.apply_adjacent_loop_rewrites(brep)?;
        Ok(CylinderTrimRewriteResult {
            boundary_outcome,
            diagnostics,
        })
    })();
    if result.is_err() {
        *brep = rollback_snapshot;
    }
    result
}

pub(super) struct CylinderTrimRewriteResult {
    pub(super) boundary_outcome: PartitionBoundaryOutcome,
    pub(super) diagnostics: Vec<crate::types::ApexConstructionDiagnostic>,
}

struct CylinderTrimRewritePreflight<'a> {
    brep: &'a BRep,
    cap_face: FaceHandle,
    face_edges: &'a [FaceEdgeOwnership],
    rebuild_plans: &'a [SemanticCylinderTrimRebuildPlan],
    topology_plans: &'a [CylinderTrimRewriteTopologyPlan],
    subcomplex_plan: &'a TrimRewriteSubcomplexPlan,
    include_selected_nonmaterialized: bool,
    require_v_free_repeated_apex: bool,
}

fn validate_cylinder_trim_rewrite_topology_preflight(
    preflight: CylinderTrimRewritePreflight<'_>,
) -> Result<(), ApexBuildError> {
    let brep = preflight.brep;
    let cap_face = preflight.cap_face;
    let face_edges = preflight.face_edges;
    let rebuild_plans = preflight.rebuild_plans;
    let topology_plans = preflight.topology_plans;
    let subcomplex_plan = preflight.subcomplex_plan;
    let include_selected_nonmaterialized = preflight.include_selected_nonmaterialized;
    let require_v_free_repeated_apex = preflight.require_v_free_repeated_apex;

    let before = topology_count_summary(brep);
    let mut candidate = brep.clone();
    let pre_materialized = pre_materialize_far_side_trim_sides(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
        include_selected_nonmaterialized,
    )?;
    let side_support_diagnostics =
        SideSupportRegionDiagnosticSnapshot::build_read_only(&candidate, subcomplex_plan);
    debug_assert!(side_support_diagnostics.validate_fresh(&candidate));
    let side_support_result =
        apply_selected_side_support_continuation(&mut candidate, subcomplex_plan)?;
    let support_ownership_already_materialized = matches!(
        side_support_result,
        SideSupportRegionContinuationOutcome::Applied(_)
    );
    let side_support_bridge_loops = side_support_result.bridge_loops();
    match apply_materialized_cylinder_loop_replacements(
        &mut candidate,
        MaterializedCylinderLoopReplacementRequest {
            cap_face,
            face_edges,
            rebuild_plans,
            topology_plans,
            subcomplex_plan,
            pre_materialized: &pre_materialized,
            side_support_bridge_loops,
            include_selected_nonmaterialized,
            reuse_shared_cap_cylinder_seam: require_v_free_repeated_apex,
            support_ownership_already_materialized,
        },
    ) {
        Ok(_) => {
            let after = topology_count_summary(&candidate);
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace §14.9/topology-preflight] accepted before {}; after {}; delta {}; side_support_continuation={}; {}; {}",
                    before.diagnostic_label(),
                    after.diagnostic_label(),
                    TopologyCountDelta::between(before, after).diagnostic_label(),
                    side_support_result.diagnostic_label(),
                    subcomplex_plan.diagnostic_label(),
                    side_support_diagnostics.diagnostic_label()
                );
            }
            Ok(())
        }
        Err(error) => {
            let after = topology_count_summary(&candidate);
            Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition target-cylinder trim rewrite preflight \
                     refused to commit topology mutation: {error}. Probe counts before {}; after \
                     {}; delta {}. Required balance: {}. Planned trim subcomplex: {}. Side-support diagnostics: {}. This is a construction blocker, not a validation blocker: the \
                     implementation must plan the Euler-balanced replacement subcomplex before \
                     changing the live BRep.",
                    before.diagnostic_label(),
                    after.diagnostic_label(),
                    TopologyCountDelta::between(before, after).diagnostic_label(),
                    TopologyBalanceRequirement::from_deficit_increase(
                        after.deficit_to_genus_zero() - before.deficit_to_genus_zero()
                    )
                    .diagnostic_label(),
                    subcomplex_plan.diagnostic_label(),
                    side_support_diagnostics.diagnostic_label()
                ),
            })
        }
    }
}

fn plan_trim_subcomplex_for_cylinder_rewrite(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    include_selected_nonmaterialized: bool,
) -> Result<TrimRewriteSubcomplexPlan, ApexBuildError> {
    let mut scratch = brep.clone();
    let mut builder = TrimRewriteSubcomplexBuilder::new();
    for (rebuild_plan, topology_plan) in rebuild_plans.iter().zip(topology_plans) {
        let allowed_support_faces =
            target_cylinder_support_faces(brep, rebuild_plan.cylinder_face, cap_face, face_edges);
        let Some(fragment) =
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)
        else {
            continue;
        };
        if !fragment.materialized_anchor_replacement && !include_selected_nonmaterialized {
            continue;
        }
        for segment in &topology_plan.loop_candidate.segments {
            if segment.label.contains("shared FilletEnd") {
                builder.push_existing_seam(
                    topology_plan.fillet_index,
                    topology_plan.cylinder_face,
                    cap_face,
                    segment.start,
                    segment.end,
                );
            } else if segment.label.contains("far-side") {
                let before = topology_count_summary(&scratch);
                let trim_sides = materialized_far_side_trim_sides(
                    &mut scratch,
                    rebuild_plan.cylinder_face,
                    segment.start,
                    segment.end,
                    segment.label,
                    rebuild_plan.fillet_index,
                    &allowed_support_faces,
                )?;
                let after = topology_count_summary(&scratch);
                builder.add_topology_delta(TopologyCountDelta::between(before, after));
                for trim_side in trim_sides {
                    if !trim_side.requires_support_ownership {
                        continue;
                    }
                    let start = scratch
                        .half_edges
                        .get(trim_side.half_edge)
                        .map(|half_edge| half_edge.origin)
                        .ok_or_else(|| ApexBuildError::PlanMismatch {
                            reason: format!(
                                "planned repeated-apex partition FilletEnd[{}] trim \
                                     subcomplex planner produced stale far-side half-edge {:?}",
                                topology_plan.fillet_index, trim_side.half_edge
                            ),
                        })?;
                    let end =
                        half_edge_end_vertex(&scratch, trim_side.half_edge).ok_or_else(|| {
                            ApexBuildError::PlanMismatch {
                                reason: format!(
                                    "planned repeated-apex partition FilletEnd[{}] trim subcomplex \
                                 planner produced far-side half-edge {:?} with stale twin",
                                    topology_plan.fillet_index, trim_side.half_edge
                                ),
                            }
                        })?;
                    builder.push_new_support_cylinder_trim(
                        topology_plan.fillet_index,
                        trim_side.cylinder_face,
                        trim_side.support_face,
                        start,
                        end,
                    );
                }
            } else {
                let support_face = support_face_for_generator_trim_side(
                    brep,
                    GeneratorTrimSupportOwnerQuery {
                        cylinder_face: rebuild_plan.cylinder_face,
                        cap_face,
                        start: segment.start,
                        end: segment.end,
                        label: segment.label,
                        fillet_index: rebuild_plan.fillet_index,
                        allowed_support_faces: &allowed_support_faces,
                    },
                )?;
                builder.push_new_generator_trim(
                    topology_plan.fillet_index,
                    topology_plan.cylinder_face,
                    support_face,
                    segment.start,
                    segment.end,
                );
            }
        }
    }
    Ok(builder.finish())
}

fn selected_far_side_trim_materialization_preflight(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> (bool, Vec<crate::types::ApexConstructionDiagnostic>) {
    let mut candidate = brep.clone();
    match pre_materialize_far_side_trim_sides(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
        true,
    ) {
        Ok(_) => (true, Vec::new()),
        Err(error) => {
            let diagnostics = materialization_preflight_diagnostics(
                brep,
                cap_face,
                face_edges,
                rebuild_plans,
                topology_plans,
                &error.to_string(),
            );
            (false, diagnostics)
        }
    }
}

fn materialization_preflight_diagnostics(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    error: &str,
) -> Vec<crate::types::ApexConstructionDiagnostic> {
    rebuild_plans
        .iter()
        .zip(topology_plans)
        .filter_map(|(rebuild_plan, topology_plan)| {
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)?;
            let far_side = topology_plan
                .loop_candidate
                .segments
                .iter()
                .find(|segment| segment.label.contains("far-side"))?;
            let allowed_support_faces = target_cylinder_support_faces(
                brep,
                rebuild_plan.cylinder_face,
                cap_face,
                face_edges,
            );
            let start_support_faces =
                support_faces_containing_vertex(brep, far_side.start, &allowed_support_faces);
            let end_support_faces =
                support_faces_containing_vertex(brep, far_side.end, &allowed_support_faces);
            let reason = if error.contains("no plane-plane-cylinder transition point") {
                crate::types::CylinderTrimMaterializationDeclineReason::NoPlanePlaneCylinderTransition
            } else {
                crate::types::CylinderTrimMaterializationDeclineReason::Other
            };
            let strategy = match topology_plan.loop_candidate.far_side_strategy {
                FarSideTrimStrategy::Forward => crate::types::CylinderTrimFarSideStrategy::Forward,
                FarSideTrimStrategy::Reverse => crate::types::CylinderTrimFarSideStrategy::Reverse,
            };
            Some(crate::types::ApexConstructionDiagnostic::CylinderTrimMaterialization(
                crate::types::CylinderTrimMaterializationDiagnostic {
                    fillet_index: topology_plan.fillet_index,
                    cylinder_face: topology_plan.cylinder_face,
                    start: far_side.start,
                    end: far_side.end,
                    strategy,
                    start_support_faces,
                    end_support_faces,
                    reason,
                    exact_authorable: false,
                },
            ))
        })
        .collect()
}

fn support_faces_containing_vertex(
    brep: &BRep,
    vertex: VertexHandle,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<FaceHandle> {
    let Some(point) = brep.vertices.get(vertex).map(|vertex| vertex.position) else {
        return Vec::new();
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    let mut faces = allowed_support_faces
        .iter()
        .copied()
        .filter(|face| {
            let Some(surface) = brep.faces.get(*face).and_then(|face| face.surface.as_ref()) else {
                return false;
            };
            let SurfaceGeometry::Plane { origin, normal, .. } = surface.geometry else {
                return false;
            };
            (point - origin).dot(&normal).abs() <= tolerance
        })
        .collect::<Vec<_>>();
    faces.sort();
    faces.dedup();
    faces
}

fn plan_semantic_cylinder_trim_rebuilds(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<Vec<SemanticCylinderTrimRebuildPlan>, ApexBuildError> {
    fillet_ends
        .iter()
        .map(|fillet_end| SemanticCylinderTrimRebuildPlan::build(brep, fillet_end))
        .collect()
}
