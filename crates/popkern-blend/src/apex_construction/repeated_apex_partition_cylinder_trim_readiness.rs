//! Readiness and topology selection for semantic cylinder trim rebuilds.

use super::*;

impl SemanticCylinderTrimRebuildPlan {
    pub(super) fn current_cylinder_order_is_constructive(&self) -> bool {
        self.retired_fragments.first().is_some_and(|fragment| {
            !fragment.materialized_anchor_replacement && fragment.off_surface_vertices.is_empty()
        })
    }

    pub(super) fn validate_rewrite_readiness(&self) -> Result<(), ApexBuildError> {
        let mut blocking_fragments = Vec::new();
        for (fragment_index, fragment) in self.retired_fragments.iter().enumerate() {
            if fragment.materialized_anchor_replacement {
                continue;
            }
            if fragment.path.is_empty() || fragment.adjacent_face_runs.is_empty() {
                blocking_fragments.push(format!(
                    "fragment[{fragment_index}] vertices {:?} has no adjacent-face run coverage",
                    fragment.vertices
                ));
            }
            let covered_half_edges = fragment
                .adjacent_face_runs
                .iter()
                .map(|run| run.path.len())
                .sum::<usize>();
            if covered_half_edges != fragment.path.len() {
                blocking_fragments.push(format!(
                    "fragment[{fragment_index}] vertices {:?} has {} run half-edge(s) for {} \
                     path half-edge(s)",
                    fragment.vertices,
                    covered_half_edges,
                    fragment.path.len()
                ));
            }
        }
        if !blocking_fragments.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot yet rewrite target \
                     cylinder face {:?}: {}. The rewrite must have complete adjacent-face run \
                     coverage so every new cylinder trim segment has a unique owning neighbor.",
                    self.fillet_index,
                    self.cylinder_face,
                    blocking_fragments.join("; ")
                ),
            });
        }
        Ok(())
    }

    pub(super) fn validate_replacement_curve_authoring_ready(
        &self,
        brep: &BRep,
    ) -> Result<(), ApexBuildError> {
        let candidate = self.selected_cylinder_loop_orientation_candidate(brep)?;
        let fragments = candidate
            .as_ref()
            .and_then(|candidate| self.fragment_for_strategy(candidate.far_side_strategy))
            .map(std::slice::from_ref)
            .unwrap_or(self.retired_fragments.as_slice());
        let mut missing = Vec::new();
        let mut needs_reownership = Vec::new();
        for fragment in fragments {
            for run in &fragment.adjacent_face_runs {
                if fragment.materialized_anchor_replacement {
                    continue;
                }
                if !run.requires_replacement() {
                    continue;
                }
                if replacement_run_exact_curve_kind(brep, self.cylinder_face, run).is_some() {
                    continue;
                }
                if let Some(plan) = replacement_run_reownership_plan(brep, self.cylinder_face, run)
                {
                    needs_reownership.push(format!(
                        "{:?} {:?}->{:?}: {}",
                        run.adjacent_face,
                        run.start,
                        run.end,
                        plan.diagnostic_label()
                    ));
                    continue;
                }
                missing.push(format!(
                    "{:?} {:?}->{:?} support surface {}; {}; {}",
                    run.adjacent_face,
                    run.start,
                    run.end,
                    face_surface_kind(brep, run.adjacent_face),
                    replacement_run_curve_blocker(brep, self.cylinder_face, run),
                    replacement_run_reownership_blocker(brep, self.cylinder_face, run)
                ));
            }
        }
        if !needs_reownership.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot apply target-cylinder \
                     trim rewrite yet: replacement support-face run(s) must be split and re-owned \
                     before curve authoring: {}. The mutation must author those ownership splits \
                     first, then attach exact support/cylinder curves to each replacement span.",
                    self.fillet_index,
                    needs_reownership.join("; ")
                ),
            });
        }
        if !missing.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot apply target-cylinder \
                     trim rewrite yet: exact replacement curve authoring is missing for \
                     support-face run(s): {}. The rewrite must derive the cylinder/support \
                     intersection curves for these runs before mutating topology; using straight \
                     chords would be a non-commercial doctored trim.",
                    self.fillet_index,
                    missing.join("; ")
                ),
            });
        }
        Ok(())
    }

    pub(super) fn fragment_for_strategy(
        &self,
        strategy: FarSideTrimStrategy,
    ) -> Option<&RetiredCylinderTrimFragment> {
        match strategy {
            FarSideTrimStrategy::Forward => self.retired_fragments.first(),
            FarSideTrimStrategy::Reverse => self.retired_fragments.get(1),
        }
    }

    pub(super) fn plan_cylinder_trim_rewrite_topology(
        &self,
        brep: &BRep,
    ) -> Result<CylinderTrimRewriteTopologyPlan, ApexBuildError> {
        self.plan_cylinder_trim_rewrite_topology_with_winding(brep, false)
    }

    /// The two generator spans of this fillet end as traversed by the
    /// CURRENT-winding cylinder loop. A side-support bridge whose
    /// reserved boundary spans oppose BOTH pairs can twin this loop
    /// directly (the §14.9 trim-subcomplex claims succeed).
    pub(super) fn current_winding_generator_spans(&self) -> [(VertexHandle, VertexHandle); 2] {
        [
            (self.cap_start, self.start_anchor),
            (self.end_anchor, self.cap_end),
        ]
    }

    /// Re-plan variant for after the side-support continuation
    /// applied bridge loops: prefer the Current cylinder winding so
    /// the generator spans twin the bridge connectors.
    pub(super) fn plan_cylinder_trim_rewrite_topology_with_winding(
        &self,
        brep: &BRep,
        prefer_current_winding: bool,
    ) -> Result<CylinderTrimRewriteTopologyPlan, ApexBuildError> {
        let loop_candidate = self
            .selected_cylinder_loop_orientation_candidate_with_winding(
                brep,
                prefer_current_winding,
            )?
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot plan target-cylinder \
                     trim topology rewrite: no exact-authorable far-side loop strategy exists",
                    self.fillet_index
                ),
            })?;
        let seam_rewrite = if loop_candidate.requires_fillet_end_reorientation() {
            let (edge, cap_half_edge, current_cylinder_half_edge) =
                self.shared_fillet_end_half_edges()?;
            self.validate_reauthored_seam_requires_cap_boundary_reorientation(
                brep,
                &loop_candidate,
                cap_half_edge,
            )?;
            FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            }
        } else {
            let (edge, cap_half_edge, current_cylinder_half_edge) =
                self.shared_fillet_end_half_edges()?;
            FilletEndSeamRewrite::ReuseCurrentTwinPair {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            }
        };
        Ok(CylinderTrimRewriteTopologyPlan {
            fillet_index: self.fillet_index,
            cylinder_face: self.cylinder_face,
            loop_candidate,
            seam_rewrite,
        })
    }

    pub(super) fn shared_fillet_end_half_edges(
        &self,
    ) -> Result<(EdgeHandle, HalfEdgeHandle, HalfEdgeHandle), ApexBuildError> {
        self.segments
            .iter()
            .find_map(|segment| match segment {
                SemanticCylinderTrimSegment::SharedFilletEnd {
                    edge,
                    cap_half_edge,
                    cylinder_half_edge,
                    ..
                } => Some((*edge, cap_half_edge, cylinder_half_edge)),
                _ => None,
            })
            .map(|(edge, cap_half_edge, cylinder_half_edge)| {
                (edge, *cap_half_edge, *cylinder_half_edge)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] has no shared FilletEnd seam \
                     in the target-cylinder trim plan",
                    self.fillet_index
                ),
            })
    }

    fn validate_reauthored_seam_requires_cap_boundary_reorientation(
        &self,
        brep: &BRep,
        loop_candidate: &CylinderLoopOrientationCandidate,
        cap_half_edge: HalfEdgeHandle,
    ) -> Result<(), ApexBuildError> {
        let cap_start = brep
            .half_edges
            .get(cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge \
                     {cap_half_edge:?} is stale while planning seam reorientation",
                    self.fillet_index
                ),
            })?;
        let cap_end = half_edge_end_vertex(brep, cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge \
                     {cap_half_edge:?} has stale twin while planning seam reorientation",
                    self.fillet_index
                ),
            }
        })?;
        let Some(shared_segment) = loop_candidate.segments.first() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] selected cylinder loop has no \
                     shared seam segment",
                    self.fillet_index
                ),
            });
        };
        if shared_segment.start == cap_start && shared_segment.end == cap_end {
            return Ok(());
        }
        if shared_segment.start == cap_end && shared_segment.end == cap_start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] selected reauthored cylinder \
                     loop is already opposite the current cap half-edge {:?}->{:?}; this should \
                     have reused the current cap boundary orientation instead of requesting cap \
                     boundary reorientation",
                    self.fillet_index, cap_start, cap_end
                ),
            });
        }
        Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] selected reauthored cylinder loop \
                 seam {:?}->{:?} does not share the current cap seam endpoints {:?}->{:?}",
                self.fillet_index, shared_segment.start, shared_segment.end, cap_start, cap_end
            ),
        })
    }

    fn selected_cylinder_loop_orientation_candidate(
        &self,
        brep: &BRep,
    ) -> Result<Option<CylinderLoopOrientationCandidate>, ApexBuildError> {
        self.selected_cylinder_loop_orientation_candidate_with_winding(brep, false)
    }

    pub(super) fn selected_cylinder_loop_orientation_candidate_with_winding(
        &self,
        brep: &BRep,
        prefer_current_winding: bool,
    ) -> Result<Option<CylinderLoopOrientationCandidate>, ApexBuildError> {
        let Some(strategy) = self.selected_far_side_strategy(brep) else {
            return Ok(None);
        };
        let candidates =
            self.cylinder_loop_orientation_candidates(brep, strategy, prefer_current_winding)?;
        let mut rejected = Vec::new();
        for candidate in candidates {
            match candidate.validate_contiguous(self.fillet_index) {
                Ok(()) => return Ok(Some(candidate)),
                Err(reason) => rejected.push(format!("{}: {reason}", candidate.diagnostic_label())),
            }
        }
        Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] selected {:?} far-side trim \
                 ownership as exact-authorable, but no cylinder loop orientation candidate is \
                 contiguous: {}",
                self.fillet_index,
                strategy,
                rejected.join("; ")
            ),
        })
    }

    fn cylinder_loop_orientation_candidates(
        &self,
        brep: &BRep,
        strategy: FarSideTrimStrategy,
        prefer_current_winding: bool,
    ) -> Result<Vec<CylinderLoopOrientationCandidate>, ApexBuildError> {
        let mut candidates = Vec::new();
        candidates.push(CylinderLoopOrientationCandidate {
            kind: CylinderLoopOrientationKind::CurrentSharedFilletEnd,
            far_side_strategy: strategy,
            segments: self.oriented_cylinder_cycle_for_strategy(brep, strategy)?,
        });
        if strategy == FarSideTrimStrategy::Reverse {
            if prefer_current_winding {
                // Side-support bridge loops exist: their connector
                // directions are FORCED by twinning the source
                // support loops, so they pair only with the CURRENT
                // cylinder winding. `Reverse` names which retired
                // fragment supplies the authorable far-side geometry
                // — not the traversal direction of the far-side
                // segment — so traverse it start->end to keep the
                // Current cycle contiguous (measured on cube_mixed
                // strict: the fully-reversed loop co-directs with
                // all six bridge connectors and the cap seam, nothing
                // twins, and the seam census drops to 3/6).
                candidates.push(CylinderLoopOrientationCandidate {
                    kind: CylinderLoopOrientationKind::CurrentSharedFilletEnd,
                    far_side_strategy: strategy,
                    segments: self
                        .oriented_cylinder_cycle_for_strategy(brep, FarSideTrimStrategy::Forward)?,
                });
            }
            candidates.push(CylinderLoopOrientationCandidate {
                kind: CylinderLoopOrientationKind::ReauthoredSharedFilletEnd,
                far_side_strategy: strategy,
                segments: vec![
                    OrientedCylinderTrimSegment {
                        label: "reauthored shared FilletEnd",
                        start: self.cap_start,
                        end: self.cap_end,
                    },
                    OrientedCylinderTrimSegment {
                        label: "generator to end anchor",
                        start: self.cap_end,
                        end: self.end_anchor,
                    },
                    OrientedCylinderTrimSegment {
                        label: "reverse far-side support-run trim",
                        start: self.end_anchor,
                        end: self.start_anchor,
                    },
                    OrientedCylinderTrimSegment {
                        label: "generator to FilletEnd start",
                        start: self.start_anchor,
                        end: self.cap_start,
                    },
                ],
            });
        }
        Ok(candidates)
    }

    fn selected_far_side_strategy(&self, brep: &BRep) -> Option<FarSideTrimStrategy> {
        let [forward, reverse] = self.retired_fragments.as_slice() else {
            return None;
        };
        let forward_status =
            fragment_replacement_authoring_status(brep, self.cylinder_face, forward);
        let reverse_status =
            fragment_replacement_authoring_status(brep, self.cylinder_face, reverse);
        let _diagnostic_context = (&forward_status.diagnostic, &reverse_status.diagnostic);
        if forward_status.authorable {
            return Some(FarSideTrimStrategy::Forward);
        }
        if !forward_status.authorable && reverse_status.authorable {
            return Some(FarSideTrimStrategy::Reverse);
        }
        None
    }

    fn oriented_cylinder_cycle_for_strategy(
        &self,
        brep: &BRep,
        strategy: FarSideTrimStrategy,
    ) -> Result<Vec<OrientedCylinderTrimSegment>, ApexBuildError> {
        let mut cycle = self.oriented_cylinder_cycle(brep)?;
        let far_side = cycle
            .iter_mut()
            .find(|segment| segment.label == "far-side support-run trim")
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] has no far-side cycle segment",
                    self.fillet_index
                ),
            })?;
        match strategy {
            FarSideTrimStrategy::Forward => {
                far_side.start = self.start_anchor;
                far_side.end = self.end_anchor;
            }
            FarSideTrimStrategy::Reverse => {
                far_side.start = self.end_anchor;
                far_side.end = self.start_anchor;
            }
        }
        Ok(cycle)
    }

    pub(super) fn forward_off_surface_vertices(&self) -> &[VertexHandle] {
        self.retired_fragments
            .first()
            .map(|fragment| fragment.off_surface_vertices.as_slice())
            .unwrap_or(&[])
    }

    pub(super) fn reverse_path_status(&self) -> String {
        let Some(reverse_fragment) = self.retired_fragments.get(1) else {
            return "the reverse loop path was not captured".to_string();
        };
        if reverse_fragment.off_surface_vertices.is_empty() {
            "the reverse loop path is V-free/on-cylinder, so the rewrite must reverse ownership \
             across existing trim edges rather than reuse the current cylinder half-edge order"
                .to_string()
        } else {
            format!(
                "the reverse loop path also has off-cylinder vertices {:?}",
                reverse_fragment.off_surface_vertices
            )
        }
    }

    pub(super) fn segment_summary(&self) -> String {
        let segments = self
            .segments
            .iter()
            .map(SemanticCylinderTrimSegment::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ");
        let retired = self
            .retired_fragments
            .iter()
            .map(RetiredCylinderTrimFragment::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "FilletEnd[{}] on {:?}: segments=[{}], retire=[{}], cap={:?}->{:?}, \
             anchors={:?}->{:?}",
            self.fillet_index,
            self.cylinder_face,
            segments,
            retired,
            self.cap_start,
            self.cap_end,
            self.start_anchor,
            self.end_anchor
        )
    }
}
