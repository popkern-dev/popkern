#[path = "repeated_apex_partition_cap_boundary_adjacent_loop_repair.rs"]
mod adjacent_loop_repair;

use super::*;
use adjacent_loop_repair::author_boundary_preserving_adjacent_loop_repair;

pub(super) fn plan_cap_boundary_rewrite_for_cylinder_trim(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_loop: popkern_core::topology::LoopHandle,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> Result<CapBoundaryRewritePlan, ApexBuildError> {
    let cap_boundary_rewrite_required = topology_plans.iter().any(|plan| {
        matches!(
            plan.seam_rewrite,
            FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide { .. }
        )
    });
    if cap_boundary_rewrite_required && face_edges.len() != topology_plans.len() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition target-cylinder trim rewrite needs a whole cap \
                 boundary rewrite, but has {} FaceEdge ownership record(s) for {} FilletEnd \
                 topology plan(s)",
                face_edges.len(),
                topology_plans.len()
            ),
        });
    }
    if face_edges.len() != fillet_ends.len() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition target-cylinder trim rewrite cannot plan cap \
                 boundary: {} FaceEdge ownership record(s) for {} FilletEnd ownership record(s)",
                face_edges.len(),
                fillet_ends.len()
            ),
        });
    }
    let current_cap_cycle = intended_cap_boundary_cycle(face_edges, fillet_ends);
    validate_half_edge_cycle(
        brep,
        &current_cap_cycle,
        "current intended cap boundary",
        cap_face,
    )?;
    let reauthored_cap_cycle = if cap_boundary_rewrite_required {
        let cycle = reversed_twin_cycle(brep, &current_cap_cycle, cap_face)?;
        validate_half_edge_cycle(brep, &cycle, "reauthored cap boundary", cap_face)?;
        cycle
    } else {
        current_cap_cycle.clone()
    };
    let adjacent_loop_rewrites = if cap_boundary_rewrite_required {
        plan_adjacent_loop_rewrites_for_reauthored_cap_cycle(brep, cap_face, &reauthored_cap_cycle)?
    } else {
        Vec::new()
    };
    Ok(CapBoundaryRewritePlan {
        cap_face,
        cap_loop,
        required: cap_boundary_rewrite_required,
        fillet_end_count: topology_plans.len(),
        face_edge_count: face_edges.len(),
        current_cap_cycle,
        reauthored_cap_cycle,
        adjacent_loop_rewrites,
    })
}

fn plan_adjacent_loop_rewrites_for_reauthored_cap_cycle(
    brep: &BRep,
    cap_face: FaceHandle,
    reauthored_cap_cycle: &[HalfEdgeHandle],
) -> Result<Vec<AdjacentLoopRewriteRequirement>, ApexBuildError> {
    let mut rewrites = Vec::new();
    for stolen_cap_side in reauthored_cap_cycle {
        let Some(stolen) = brep.half_edges.get(*stolen_cap_side) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite references stale \
                     reauthored cap half-edge {stolen_cap_side:?}"
                ),
            });
        };
        let replacement_neighbor_side = stolen.twin;
        let neighbor_loop = brep
            .loops
            .contains_key(stolen.loop_)
            .then_some(stolen.loop_);
        let neighbor_face =
            neighbor_loop.and_then(|loop_handle| brep.loops.get(loop_handle).map(|data| data.face));
        if neighbor_face == Some(cap_face) {
            continue;
        }
        let current_neighbor_cycle = if let Some(loop_handle) = neighbor_loop {
            collect_loop_half_edges(brep, loop_handle)?
        } else {
            Vec::new()
        };
        let has_boundary_neighbor_side = current_neighbor_cycle.iter().any(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .is_some_and(|data| !brep.half_edges.contains_key(data.twin))
        });
        let rewrite_kind = if neighbor_loop.is_none() {
            AdjacentLoopRewriteKind::AttachTwinToLaterStripLoop
        } else if has_boundary_neighbor_side {
            AdjacentLoopRewriteKind::BoundaryPreservingAdjacentLoopRepairRequired
        } else {
            AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle
        };
        rewrites.push(AdjacentLoopRewriteRequirement {
            stolen_cap_side: *stolen_cap_side,
            replacement_neighbor_side,
            neighbor_loop,
            neighbor_face,
            current_neighbor_cycle: current_neighbor_cycle.clone(),
            reauthored_neighbor_cycle: if rewrite_kind
                == AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle
            {
                let reauthored = reversed_twin_cycle(brep, &current_neighbor_cycle, cap_face)?;
                validate_half_edge_cycle(
                    brep,
                    &reauthored,
                    "reauthored adjacent twin boundary",
                    neighbor_face.unwrap_or(cap_face),
                )?;
                reauthored
            } else {
                Vec::new()
            },
            rewrite_kind,
        });
    }
    Ok(rewrites)
}

pub(super) struct CapBoundaryRewritePlan {
    cap_face: FaceHandle,
    cap_loop: popkern_core::topology::LoopHandle,
    required: bool,
    fillet_end_count: usize,
    face_edge_count: usize,
    current_cap_cycle: Vec<HalfEdgeHandle>,
    reauthored_cap_cycle: Vec<HalfEdgeHandle>,
    adjacent_loop_rewrites: Vec<AdjacentLoopRewriteRequirement>,
}

struct AdjacentLoopRewriteRequirement {
    stolen_cap_side: HalfEdgeHandle,
    replacement_neighbor_side: HalfEdgeHandle,
    neighbor_loop: Option<popkern_core::topology::LoopHandle>,
    neighbor_face: Option<FaceHandle>,
    current_neighbor_cycle: Vec<HalfEdgeHandle>,
    reauthored_neighbor_cycle: Vec<HalfEdgeHandle>,
    rewrite_kind: AdjacentLoopRewriteKind,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AdjacentLoopRewriteKind {
    ReverseAdjacentTwinCycle,
    BoundaryPreservingAdjacentLoopRepairRequired,
    AttachTwinToLaterStripLoop,
}

impl CapBoundaryRewritePlan {
    pub(super) fn apply_adjacent_loop_rewrites(
        &self,
        brep: &mut BRep,
    ) -> Result<(), ApexBuildError> {
        if !self.required {
            return Ok(());
        }
        for rewrite in &self.adjacent_loop_rewrites {
            let Some(loop_handle) = rewrite.neighbor_loop else {
                continue;
            };
            let Some(face) = rewrite.neighbor_face else {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition adjacent loop rewrite has loop \
                         {loop_handle:?} but no owning face: {}",
                        rewrite.diagnostic_label()
                    ),
                });
            };
            match rewrite.rewrite_kind {
                AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle => {
                    if rewrite.reauthored_neighbor_cycle.is_empty() {
                        return Err(ApexBuildError::PlanMismatch {
                            reason: format!(
                                "planned repeated-apex partition adjacent loop rewrite has no \
                                 reauthored cycle: {}",
                                rewrite.diagnostic_label()
                            ),
                        });
                    }
                    wire_existing_loop(
                        brep,
                        loop_handle,
                        face,
                        &rewrite.reauthored_neighbor_cycle,
                    )?;
                }
                AdjacentLoopRewriteKind::BoundaryPreservingAdjacentLoopRepairRequired => {
                    author_boundary_preserving_adjacent_loop_repair(
                        brep,
                        rewrite,
                        loop_handle,
                        face,
                    )?;
                }
                AdjacentLoopRewriteKind::AttachTwinToLaterStripLoop => {}
            }
        }
        Ok(())
    }

    pub(super) fn boundary_outcome(
        &self,
        topology_plans: &[CylinderTrimRewriteTopologyPlan],
    ) -> PartitionBoundaryOutcome {
        if !self.required {
            return PartitionBoundaryOutcome::from_current_cap_cycle(
                self.current_cap_cycle.clone(),
            );
        }
        let mut strip_side_overrides = std::collections::BTreeMap::new();
        for plan in topology_plans {
            let cap_half_edge = match plan.seam_rewrite {
                FilletEndSeamRewrite::ReuseCurrentTwinPair { cap_half_edge, .. }
                | FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                    cap_half_edge, ..
                } => cap_half_edge,
            };
            strip_side_overrides.insert(plan.fillet_index, cap_half_edge);
        }
        PartitionBoundaryOutcome::from_reauthored_cap_cycle(
            self.reauthored_cap_cycle.clone(),
            strip_side_overrides,
        )
    }

    pub(super) fn validate_adjacent_loop_rewrite_readiness(&self) -> Result<(), ApexBuildError> {
        if !self.required {
            return Ok(());
        }
        let missing_adjacent_cycles = self
            .adjacent_loop_rewrites
            .iter()
            .filter(|rewrite| {
                rewrite.rewrite_kind == AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle
                    && rewrite.reauthored_neighbor_cycle.is_empty()
            })
            .map(AdjacentLoopRewriteRequirement::diagnostic_label)
            .collect::<Vec<_>>();
        if !missing_adjacent_cycles.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite cannot yet mutate: \
                     adjacent support/cylinder loop rewrite plans are incomplete before stealing \
                     live half-edges for the cap boundary: {}",
                    missing_adjacent_cycles.join("; ")
                ),
            });
        }
        Ok(())
    }

    pub(super) fn validate_rewrite_postconditions(
        &self,
        brep: &BRep,
        face_edges: &[FaceEdgeOwnership],
        topology_plans: &[CylinderTrimRewriteTopologyPlan],
    ) -> Result<(), ApexBuildError> {
        if !self.required {
            return Ok(());
        }
        for plan in topology_plans {
            self.validate_fillet_end_rewrite_postcondition(brep, plan)?;
        }
        for face_edge in face_edges {
            self.validate_face_edge_rewrite_postcondition(brep, face_edge)?;
        }
        Ok(())
    }

    fn validate_fillet_end_rewrite_postcondition(
        &self,
        brep: &BRep,
        plan: &CylinderTrimRewriteTopologyPlan,
    ) -> Result<(), ApexBuildError> {
        let FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
            cap_half_edge,
            current_cylinder_half_edge,
            ..
        } = plan.seam_rewrite
        else {
            return Ok(());
        };
        if !self
            .reauthored_cap_cycle
            .contains(&current_cylinder_half_edge)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     does not include the reauthored cap-side half-edge \
                     {current_cylinder_half_edge:?}",
                    plan.fillet_index
                ),
            });
        }
        let twin = brep
            .half_edges
            .get(current_cylinder_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     references stale reauthored cap-side half-edge \
                     {current_cylinder_half_edge:?}",
                    plan.fillet_index
                ),
            })?;
        if twin != cap_half_edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     expected reauthored cap-side half-edge {current_cylinder_half_edge:?} \
                     to be opposite selected cylinder-side half-edge {cap_half_edge:?}, got \
                     twin {twin:?}",
                    plan.fillet_index
                ),
            });
        }
        let Some(shared_segment) = plan.loop_candidate.segments.first() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     has no selected cylinder shared-seam segment",
                    plan.fillet_index
                ),
            });
        };
        let selected_start = brep
            .half_edges
            .get(cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     references stale selected cylinder-side half-edge {cap_half_edge:?}",
                    plan.fillet_index
                ),
            })?;
        let selected_end = half_edge_end_vertex(brep, cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for \
                         FilletEnd[{}] selected cylinder-side half-edge {cap_half_edge:?} \
                         has stale twin",
                    plan.fillet_index
                ),
            }
        })?;
        if shared_segment.start != selected_start || shared_segment.end != selected_end {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FilletEnd[{}] \
                     selected cylinder side {cap_half_edge:?} is {:?}->{:?}, but the selected \
                     cylinder loop shared segment is {:?}->{:?}",
                    plan.fillet_index,
                    selected_start,
                    selected_end,
                    shared_segment.start,
                    shared_segment.end
                ),
            });
        }
        Ok(())
    }

    fn validate_face_edge_rewrite_postcondition(
        &self,
        brep: &BRep,
        face_edge: &FaceEdgeOwnership,
    ) -> Result<(), ApexBuildError> {
        if !self
            .reauthored_cap_cycle
            .contains(&face_edge.support_half_edge)
        {
            return Ok(());
        }
        let twin = brep
            .half_edges
            .get(face_edge.support_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FaceEdge[{}] \
                     references stale support half-edge {:?}",
                    face_edge.fillet_index, face_edge.support_half_edge
                ),
            })?;
        if twin != face_edge.face_edge_cap_half_edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cap-boundary rewrite for FaceEdge[{}] \
                     expected reauthored cap side {:?} to be opposite previous cap side {:?}, \
                     got twin {twin:?}",
                    face_edge.fillet_index,
                    face_edge.support_half_edge,
                    face_edge.face_edge_cap_half_edge
                ),
            });
        }
        Ok(())
    }

    pub(super) fn diagnostic_label(&self) -> String {
        if self.required {
            format!(
                "reauthor the whole cap boundary on {:?}/{:?} across {} FilletEnd seam(s) and \
                 {} FaceEdge seam(s) so every cap-side half-edge remains opposite its \
                 support/cylinder-side twin; current_cap_cycle={:?}; \
                 reauthored_cap_cycle={:?}; adjacent_loop_rewrites=[{}]",
                self.cap_face,
                self.cap_loop,
                self.fillet_end_count,
                self.face_edge_count,
                self.current_cap_cycle,
                self.reauthored_cap_cycle,
                self.adjacent_loop_rewrites
                    .iter()
                    .map(AdjacentLoopRewriteRequirement::diagnostic_label)
                    .collect::<Vec<_>>()
                    .join("; ")
            )
        } else {
            format!(
                "reuse current cap boundary on {:?}/{:?} with {} FilletEnd seam(s) and {} \
                 FaceEdge seam(s); cap_cycle={:?}",
                self.cap_face,
                self.cap_loop,
                self.fillet_end_count,
                self.face_edge_count,
                self.current_cap_cycle
            )
        }
    }
}

impl AdjacentLoopRewriteRequirement {
    fn diagnostic_label(&self) -> String {
        format!(
            "stolen_cap_side={:?}, replacement_neighbor_side={:?}, neighbor_loop={:?}, \
             neighbor_face={:?}, rewrite_kind={}, current_neighbor_cycle={:?}, \
             reauthored_neighbor_cycle={:?}: {}",
            self.stolen_cap_side,
            self.replacement_neighbor_side,
            self.neighbor_loop,
            self.neighbor_face,
            self.rewrite_kind.diagnostic_label(),
            self.current_neighbor_cycle,
            self.reauthored_neighbor_cycle,
            self.rewrite_kind.reason()
        )
    }
}
