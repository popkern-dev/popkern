//! Transactional target-cylinder trim rewrite orchestration.

#[path = "repeated_apex_partition_cylinder_trim_transaction_diagnostics.rs"]
mod transaction_diagnostics;
#[path = "repeated_apex_partition_cylinder_trim_subcomplex.rs"]
mod trim_subcomplex;
use trim_subcomplex::{
    plan_trim_subcomplex_for_cylinder_rewrite, twinless_boundary_spans_of_loops,
};

use super::*;
use transaction_diagnostics::{
    plan_semantic_cylinder_trim_rebuilds, selected_far_side_trim_materialization_preflight,
};

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
        let (include_selected_nonmaterialized, mut diagnostics) = if require_v_free_repeated_apex {
            selected_far_side_trim_materialization_preflight(
                brep,
                cap_face,
                face_edges,
                &rebuild_plans,
                &topology_plans,
            )
        } else {
            (false, Vec::new())
        };
        let subcomplex_plan = plan_trim_subcomplex_for_cylinder_rewrite(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            include_selected_nonmaterialized,
        )?;
        let default_preflight =
            validate_cylinder_trim_rewrite_topology_preflight(CylinderTrimRewritePreflight {
                brep,
                cap_face,
                face_edges,
                rebuild_plans: &rebuild_plans,
                topology_plans: &topology_plans,
                subcomplex_plan: &subcomplex_plan,
                include_selected_nonmaterialized,
                require_v_free_repeated_apex,
            });
        // Side-support bridge connectors twin the CURRENT cylinder
        // winding, so a fillet end with ANY opposing reserved bridge
        // span must be planned with that winding BEFORE anything
        // applies — under the fully-reversed Reauthored loop the
        // reserved span co-directs with the planned generator and
        // cannot twin (measured: cube_mixed strict seam census 3/6
        // when no end flipped; cube_extreme support 13v1 chain
        // 28->23 vs 28->24 when FilletEnd[2] kept Reauthored while
        // its (22,25) connector was reserved). Spans without a
        // reserved connector author through the mef/patch path,
        // which works in either winding. Ends with NO reserved
        // spans keep the default Reauthored path. Probe the
        // decision on a clone: the continuation's apply path is
        // its own readiness authority.
        let (probe_bridge_spans, probe_preserves_current_winding): (
            Vec<(VertexHandle, VertexHandle)>,
            bool,
        ) = {
            let mut probe = brep.clone();
            match apply_selected_side_support_continuation(&mut probe, &subcomplex_plan) {
                Ok(outcome) => (
                    twinless_boundary_spans_of_loops(&probe, outcome.bridge_loops()),
                    matches!(
                        &outcome,
                        SideSupportRegionContinuationOutcome::Applied(result)
                            if result.preserves_current_winding()
                    ),
                ),
                Err(_) => (Vec::new(), false),
            }
        };
        // The PK-167 in-place continuation preserves the CURRENT
        // shell winding by construction (face identity kept, no
        // reserved twinless spans exist to signal it) — every fillet
        // end must plan with the current winding, or the Reauthored
        // reauthor steals the cap sides and hands each FilletEnd
        // seam twin to the retired apex-cell residual (measured:
        // cube_extreme cap x unsurfaced-residual seams, V-free
        // acceptance declined with 3 missing support pcurves).
        let prefer_current_winding: Vec<bool> = rebuild_plans
            .iter()
            .map(|plan| {
                probe_preserves_current_winding
                    || plan
                        .current_winding_generator_spans()
                        .iter()
                        .any(|(start, end)| probe_bridge_spans.contains(&(*end, *start)))
            })
            .collect();
        // The default (Reauthored-winding) preflight is advisory when a
        // replan is available: a first-preflight failure with reserved
        // connector spans is exactly the mixed-winding defect the
        // replan repairs, so the transaction must not abort before
        // trying it. Without a replan trigger the default preflight
        // failure is final.
        let (topology_plans, subcomplex_plan, include_selected_nonmaterialized) =
            if !prefer_current_winding.iter().any(|prefer| *prefer) {
                default_preflight?;
                (
                    topology_plans,
                    subcomplex_plan,
                    include_selected_nonmaterialized,
                )
            } else {
                let replanned = (|| -> Result<_, ApexBuildError> {
                    let topology_plans = rebuild_plans
                        .iter()
                        .zip(&prefer_current_winding)
                        .map(|(plan, prefer)| {
                            plan.plan_cylinder_trim_rewrite_topology_with_winding(brep, *prefer)
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    let (include_selected_nonmaterialized, replanned_diagnostics) =
                        if require_v_free_repeated_apex {
                            selected_far_side_trim_materialization_preflight(
                                brep,
                                cap_face,
                                face_edges,
                                &rebuild_plans,
                                &topology_plans,
                            )
                        } else {
                            (false, Vec::new())
                        };
                    let subcomplex_plan = plan_trim_subcomplex_for_cylinder_rewrite(
                        brep,
                        cap_face,
                        face_edges,
                        &rebuild_plans,
                        &topology_plans,
                        include_selected_nonmaterialized,
                    )?;
                    validate_cylinder_trim_rewrite_topology_preflight(
                        CylinderTrimRewritePreflight {
                            brep,
                            cap_face,
                            face_edges,
                            rebuild_plans: &rebuild_plans,
                            topology_plans: &topology_plans,
                            subcomplex_plan: &subcomplex_plan,
                            include_selected_nonmaterialized,
                            require_v_free_repeated_apex,
                        },
                    )?;
                    Ok((
                        topology_plans,
                        subcomplex_plan,
                        include_selected_nonmaterialized,
                        replanned_diagnostics,
                    ))
                })();
                match (replanned, default_preflight) {
                    (Ok((plans, subcomplex, include, replanned_diagnostics)), _) => {
                        diagnostics.extend(replanned_diagnostics);
                        (plans, subcomplex, include)
                    }
                    // Replan failed but the default plans preflighted
                    // clean — keep them rather than aborting a viable
                    // transaction.
                    (Err(_), Ok(())) => (
                        topology_plans,
                        subcomplex_plan,
                        include_selected_nonmaterialized,
                    ),
                    // Both shapes failed: report the replanned wall —
                    // it is the shape the reserved connectors demand.
                    (Err(replanned_error), Err(_)) => return Err(replanned_error),
                }
            };
        let side_support_result = apply_selected_side_support_continuation(brep, &subcomplex_plan)?;
        if std::env::var_os("PK_SLIR_TRACE").is_some() {
            eprintln!(
                "[slir] transaction FINAL apply: continuation outcome {}",
                side_support_result.diagnostic_label()
            );
        }
        if let Some(diagnostic) = side_support_result.applied_diagnostic(cap_face) {
            diagnostics.push(diagnostic);
        }
        // Re-resolve plan handles through the PK-167 continuation's
        // coincident-duplicate vertex remap (plans authored pre-merge).
        let (topology_plans, subcomplex_plan) =
            if let SideSupportRegionContinuationOutcome::Applied(result) = &side_support_result {
                (
                    topology_plans
                        .iter()
                        .map(|plan| plan.with_canonical_vertices(|v| result.canonical_vertex(v)))
                        .collect::<Vec<_>>(),
                    subcomplex_plan.with_canonical_vertices(|v| result.canonical_vertex(v)),
                )
            } else {
                (topology_plans, subcomplex_plan)
            };
        let support_ownership_already_materialized = matches!(
            side_support_result,
            SideSupportRegionContinuationOutcome::Applied(_)
        );
        let side_support_bridge_loops = side_support_result.bridge_loops();
        let pre_materialized = pre_materialize_far_side_trim_sides(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            include_selected_nonmaterialized,
        )?;
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
        let mut boundary_outcome = if require_v_free_repeated_apex {
            materialized_cap_boundary_from_target_cylinder_replacements(
                brep,
                cap_face,
                face_edges,
                fillet_ends,
                &cylinder_loop_replacements,
                require_v_free_repeated_apex,
            )?
            .unwrap_or_else(|| PartitionBoundaryOutcome::from_current_cap_cycle(Vec::new()))
        } else {
            PartitionBoundaryOutcome::from_current_cap_cycle(Vec::new())
        };
        let cap_boundary_plan = if boundary_outcome.cap_cycle.is_empty() {
            let plan_result = plan_cap_boundary_rewrite_for_cylinder_trim(
                brep,
                cap_face,
                cap_loop,
                face_edges,
                fillet_ends,
                &topology_plans,
            );
            match plan_result {
                Ok(plan) => {
                    plan.validate_rewrite_postconditions(brep, face_edges, &topology_plans)?;
                    plan.refresh_reauthored_fillet_end_curves(brep, &topology_plans)?;
                    plan.validate_adjacent_loop_rewrite_readiness()?;
                    boundary_outcome = plan.boundary_outcome(&topology_plans);
                    Some(plan)
                }
                Err(plan_error) => {
                    let Some(materialized) =
                        materialized_cap_boundary_from_target_cylinder_replacements(
                            brep,
                            cap_face,
                            face_edges,
                            fillet_ends,
                            &cylinder_loop_replacements,
                            require_v_free_repeated_apex,
                        )?
                    else {
                        return Err(plan_error);
                    };
                    boundary_outcome = materialized;
                    None
                }
            }
        } else {
            None
        };
        if require_v_free_repeated_apex {
            for replacement in &cylinder_loop_replacements {
                boundary_outcome.mark_existing_target_cylinder_strip(replacement.fillet_index);
            }
        }
        if !require_v_free_repeated_apex || cylinder_loop_replacements.len() != fillet_ends.len() {
            boundary_outcome.plan_strip_loop_overrides(brep, fillet_ends)?;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/cylinder-trim-rewrite] topology_plans=[{}]; side_support_continuation={}; {}; {}",
                topology_plans
                    .iter()
                    .map(CylinderTrimRewriteTopologyPlan::diagnostic_label)
                    .collect::<Vec<_>>()
                    .join("; "),
                side_support_result.diagnostic_label(),
                cap_boundary_plan
                    .as_ref()
                    .map(|plan| plan.diagnostic_label())
                    .unwrap_or_else(|| {
                        "materialized target-cylinder cap boundary from selected seam twins"
                            .to_string()
                    }),
                boundary_outcome.diagnostic_label()
            );
        }
        if let Some(cap_boundary_plan) = &cap_boundary_plan {
            cap_boundary_plan.apply_adjacent_loop_rewrites(brep)?;
            author_materialized_shared_fillet_end_trims_for_plan(brep, &topology_plans)?;
        }
        // PK-074: the support-cell surface solve runs as the LATE pocket
        // completion in direct_apex_dispatch (acceptance preparation) —
        // mid-transaction the model still carries transiently broken
        // loops that fail mef's global postconditions, and later stages
        // keep re-trimming the pocket complex. See docs/decisions/
        // 2026-07-11-side-support-pocket-floor-reconstruction.md.
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
    let require_v_free_repeated_apex = preflight.require_v_free_repeated_apex;
    let include_selected_nonmaterialized =
        preflight.include_selected_nonmaterialized && require_v_free_repeated_apex;

    let before = topology_count_summary(brep);
    let mut candidate = brep.clone();
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
    // Re-resolve plan handles through the PK-167 merge remap.
    let remapped_plans: Vec<CylinderTrimRewriteTopologyPlan>;
    let remapped_subcomplex: TrimRewriteSubcomplexPlan;
    let (topology_plans, subcomplex_plan) =
        if let SideSupportRegionContinuationOutcome::Applied(result) = &side_support_result {
            remapped_plans = topology_plans
                .iter()
                .map(|plan| plan.with_canonical_vertices(|v| result.canonical_vertex(v)))
                .collect();
            remapped_subcomplex =
                subcomplex_plan.with_canonical_vertices(|v| result.canonical_vertex(v));
            (&remapped_plans[..], &remapped_subcomplex)
        } else {
            (topology_plans, subcomplex_plan)
        };
    let pre_materialized = pre_materialize_far_side_trim_sides(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
        include_selected_nonmaterialized,
    )?;
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
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/topology-preflight] accepted before {}; after {}; delta {}; side_support_continuation={}; {}; {}",
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
