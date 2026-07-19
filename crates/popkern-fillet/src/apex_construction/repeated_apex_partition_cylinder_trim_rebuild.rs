//! Semantic rebuild model for target-cylinder trims in repeated-apex rewrites.

use super::*;

#[derive(Clone)]
pub(super) struct SemanticCylinderTrimRebuildPlan {
    pub(super) fillet_index: usize,
    pub(super) cylinder_face: FaceHandle,
    pub(super) start_anchor: VertexHandle,
    pub(super) cap_start: VertexHandle,
    pub(super) cap_end: VertexHandle,
    pub(super) end_anchor: VertexHandle,
    pub(super) segments: Vec<SemanticCylinderTrimSegment>,
    pub(super) retired_fragments: Vec<RetiredCylinderTrimFragment>,
}

#[derive(Clone)]
pub(super) enum SemanticCylinderTrimSegment {
    Generator {
        start: VertexHandle,
        end: VertexHandle,
    },
    SharedFilletEnd {
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        cylinder_half_edge: HalfEdgeHandle,
        start: VertexHandle,
        end: VertexHandle,
    },
    FarSideCylinderTrim {
        start: VertexHandle,
        end: VertexHandle,
        replacement_for: Vec<RetiredCylinderTrimFragment>,
    },
}

impl SemanticCylinderTrimRebuildPlan {
    pub(super) fn build(
        brep: &BRep,
        fillet_end: &FilletEndOwnership,
    ) -> Result<Self, ApexBuildError> {
        let cap_start = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        let cap_end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has \
                     stale twin",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            }
        })?;
        let cap_start_point = brep.vertices[cap_start].position;
        let cap_end_point = brep.vertices[cap_end].position;
        let start_anchor = generator_anchor_on_cylinder_loop(
            brep,
            fillet_end.target_cylinder_face,
            cap_start_point,
        )
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot plan cylinder trim \
                 rebuild: no generator anchor on target cylinder face {:?} for start endpoint \
                 {:?}: {}",
                fillet_end.fillet_index,
                fillet_end.target_cylinder_face,
                cap_start,
                generator_anchor_blocker(brep, fillet_end.target_cylinder_face, cap_start_point)
            ),
        })?;
        let end_anchor =
            generator_anchor_on_cylinder_loop(brep, fillet_end.target_cylinder_face, cap_end_point)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{}] cannot plan cylinder trim \
                         rebuild: no generator anchor on target cylinder face {:?} for end \
                         endpoint {:?}: {}",
                        fillet_end.fillet_index,
                        fillet_end.target_cylinder_face,
                        cap_end,
                        generator_anchor_blocker(
                            brep,
                            fillet_end.target_cylinder_face,
                            cap_end_point
                        )
                    ),
                })?;

        let forward_fragment = if let Some(forward_path) = outer_loop_path_between_vertices(
            brep,
            fillet_end.target_cylinder_face,
            start_anchor,
            end_anchor,
        ) {
            RetiredCylinderTrimFragment::from_path(
                brep,
                fillet_end.target_cylinder_face,
                forward_path,
            )
        } else if vertices_are_materialized_generator_anchors(
            brep,
            fillet_end.target_cylinder_face,
            start_anchor,
            end_anchor,
        ) {
            RetiredCylinderTrimFragment::materialized_anchor_replacement(start_anchor, end_anchor)
        } else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot plan cylinder trim \
                     rebuild: no loop path {:?}->{:?} on target cylinder face {:?}",
                    fillet_end.fillet_index,
                    start_anchor,
                    end_anchor,
                    fillet_end.target_cylinder_face
                ),
            });
        };
        let reverse_fragment = if let Some(reverse_path) = outer_loop_path_between_vertices(
            brep,
            fillet_end.target_cylinder_face,
            end_anchor,
            start_anchor,
        ) {
            RetiredCylinderTrimFragment::from_path(
                brep,
                fillet_end.target_cylinder_face,
                reverse_path,
            )
        } else if vertices_are_materialized_generator_anchors(
            brep,
            fillet_end.target_cylinder_face,
            start_anchor,
            end_anchor,
        ) {
            RetiredCylinderTrimFragment::materialized_anchor_replacement(end_anchor, start_anchor)
        } else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot plan cylinder trim \
                     rebuild: no loop path {:?}->{:?} on target cylinder face {:?}",
                    fillet_end.fillet_index,
                    end_anchor,
                    start_anchor,
                    fillet_end.target_cylinder_face
                ),
            });
        };
        let retired_fragments = vec![forward_fragment.clone(), reverse_fragment.clone()];
        let segments = vec![
            SemanticCylinderTrimSegment::Generator {
                start: start_anchor,
                end: cap_start,
            },
            SemanticCylinderTrimSegment::SharedFilletEnd {
                edge: fillet_end.edge,
                cap_half_edge: fillet_end.cap_half_edge,
                cylinder_half_edge: fillet_end.cylinder_candidate_half_edge,
                start: cap_end,
                end: cap_start,
            },
            SemanticCylinderTrimSegment::Generator {
                start: cap_end,
                end: end_anchor,
            },
            SemanticCylinderTrimSegment::FarSideCylinderTrim {
                start: end_anchor,
                end: start_anchor,
                replacement_for: retired_fragments.clone(),
            },
        ];
        let plan = Self {
            fillet_index: fillet_end.fillet_index,
            cylinder_face: fillet_end.target_cylinder_face,
            start_anchor,
            cap_start,
            cap_end,
            end_anchor,
            segments,
            retired_fragments,
        };
        plan.validate_segment_endpoints(brep)?;
        plan.validate_segment_geometry(brep)?;
        plan.validate_oriented_cylinder_cycle(brep)?;
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/rebuild-fragments] FilletEnd[{}] cylinder {:?} anchors \
                 {:?}->{:?}: forward path {} he / vertices {:?} (off-surface {:?}); reverse path \
                 {} he / vertices {:?} (off-surface {:?})",
                plan.fillet_index,
                plan.cylinder_face,
                plan.start_anchor,
                plan.end_anchor,
                plan.retired_fragments[0].path.len(),
                plan.retired_fragments[0].vertices,
                plan.retired_fragments[0].off_surface_vertices,
                plan.retired_fragments[1].path.len(),
                plan.retired_fragments[1].vertices,
                plan.retired_fragments[1].off_surface_vertices,
            );
        }
        Ok(plan)
    }

    fn validate_segment_endpoints(&self, brep: &BRep) -> Result<(), ApexBuildError> {
        for segment in &self.segments {
            let (start, end) = segment.endpoints();
            let start_point = brep
                .vertices
                .get(start)
                .map(|vertex| vertex.position)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{}] cylinder trim segment \
                         references stale start vertex {start:?}",
                        self.fillet_index
                    ),
                })?;
            let end_point = brep
                .vertices
                .get(end)
                .map(|vertex| vertex.position)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{}] cylinder trim segment \
                         references stale end vertex {end:?}",
                        self.fillet_index
                    ),
                })?;
            if !point_lies_on_cylinder_face(brep, self.cylinder_face, start_point)
                || !point_lies_on_cylinder_face(brep, self.cylinder_face, end_point)
            {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{}] cylinder trim segment \
                         {start:?}->{end:?} is not bounded by points on target cylinder face {:?}",
                        self.fillet_index, self.cylinder_face
                    ),
                });
            }
        }
        Ok(())
    }

    fn validate_segment_geometry(&self, brep: &BRep) -> Result<(), ApexBuildError> {
        let cylinder = CylinderTrimSurface::from_face(brep, self.cylinder_face, self.fillet_index)?;
        for segment in &self.segments {
            match segment {
                SemanticCylinderTrimSegment::Generator { start, end } => {
                    cylinder.validate_generator_segment(brep, self.fillet_index, *start, *end)?;
                }
                SemanticCylinderTrimSegment::FarSideCylinderTrim {
                    replacement_for, ..
                } => {
                    self.validate_far_side_replacement_runs(replacement_for)?;
                }
                SemanticCylinderTrimSegment::SharedFilletEnd {
                    edge,
                    cap_half_edge: _,
                    cylinder_half_edge,
                    ..
                } => {
                    if brep
                        .edges
                        .get(*edge)
                        .is_none_or(|edge_data| edge_data.curve.is_none())
                    {
                        return Err(ApexBuildError::PlanMismatch {
                            reason: format!(
                                "planned repeated-apex partition FilletEnd[{}] shared cylinder \
                                 trim edge {edge:?}/{cylinder_half_edge:?} has no curve geometry",
                                self.fillet_index
                            ),
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn validate_far_side_replacement_runs(
        &self,
        fragments: &[RetiredCylinderTrimFragment],
    ) -> Result<(), ApexBuildError> {
        if fragments
            .iter()
            .any(|fragment| fragment.materialized_anchor_replacement)
        {
            return Ok(());
        }
        let run_count = fragments
            .iter()
            .map(|fragment| fragment.adjacent_face_runs.len())
            .sum::<usize>();
        if run_count == 0 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] far-side cylinder trim has \
                     no support-face run decomposition",
                    self.fillet_index
                ),
            });
        }
        let replacement_runs = fragments
            .iter()
            .flat_map(|fragment| fragment.adjacent_face_runs.iter())
            .filter(|run| run.requires_replacement())
            .count();
        if replacement_runs == 0 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] far-side cylinder trim has \
                     no off-cylinder support-face run to replace",
                    self.fillet_index
                ),
            });
        }
        Ok(())
    }

    fn validate_oriented_cylinder_cycle(&self, brep: &BRep) -> Result<(), ApexBuildError> {
        let cycle = self.oriented_cylinder_cycle(brep)?;
        for i in 0..cycle.len() {
            let current = &cycle[i];
            let next = &cycle[(i + 1) % cycle.len()];
            if current.end != next.start {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{}] target cylinder cycle is \
                         not contiguous: segment {} ends at {:?}, but segment {} starts at {:?}",
                        self.fillet_index, current.label, current.end, next.label, next.start
                    ),
                });
            }
        }
        Ok(())
    }

    pub(super) fn oriented_cylinder_cycle(
        &self,
        brep: &BRep,
    ) -> Result<Vec<OrientedCylinderTrimSegment>, ApexBuildError> {
        let (shared_half_edge, recorded_start, recorded_end) = self
            .segments
            .iter()
            .find_map(|segment| match segment {
                SemanticCylinderTrimSegment::SharedFilletEnd {
                    cylinder_half_edge,
                    start,
                    end,
                    ..
                } => Some((*cylinder_half_edge, *start, *end)),
                _ => None,
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] has no shared cylinder \
                     FilletEnd segment",
                    self.fillet_index
                ),
            })?;
        let shared_start = brep
            .half_edges
            .get(shared_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] shared cylinder half-edge \
                     {shared_half_edge:?} is stale",
                    self.fillet_index
                ),
            })?;
        let shared_end = half_edge_end_vertex(brep, shared_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] shared cylinder half-edge \
                     {shared_half_edge:?} has stale twin",
                    self.fillet_index
                ),
            }
        })?;
        if shared_start != recorded_start || shared_end != recorded_end {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] shared cylinder half-edge \
                     orientation changed: recorded {:?}->{:?}, actual {:?}->{:?}",
                    self.fillet_index, recorded_start, recorded_end, shared_start, shared_end
                ),
            });
        }
        Ok(vec![
            OrientedCylinderTrimSegment {
                label: "shared FilletEnd",
                start: shared_start,
                end: shared_end,
            },
            OrientedCylinderTrimSegment {
                label: "generator to start anchor",
                start: self.cap_start,
                end: self.start_anchor,
            },
            OrientedCylinderTrimSegment {
                label: "far-side support-run trim",
                start: self.start_anchor,
                end: self.end_anchor,
            },
            OrientedCylinderTrimSegment {
                label: "generator to FilletEnd end",
                start: self.end_anchor,
                end: self.cap_end,
            },
        ])
    }
}

impl SemanticCylinderTrimSegment {
    pub(super) fn endpoints(&self) -> (VertexHandle, VertexHandle) {
        match self {
            SemanticCylinderTrimSegment::Generator { start, end }
            | SemanticCylinderTrimSegment::SharedFilletEnd { start, end, .. }
            | SemanticCylinderTrimSegment::FarSideCylinderTrim { start, end, .. } => (*start, *end),
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SemanticCylinderTrimSegment::Generator { start, end } => {
                format!("generator {start:?}->{end:?}")
            }
            SemanticCylinderTrimSegment::SharedFilletEnd {
                edge,
                cap_half_edge: _,
                cylinder_half_edge,
                start,
                end,
            } => format!(
                "shared FilletEnd edge={edge:?} cylinder_he={cylinder_half_edge:?} \
                 {start:?}->{end:?}"
            ),
            SemanticCylinderTrimSegment::FarSideCylinderTrim {
                start,
                end,
                replacement_for,
            } => format!(
                "far-side cylinder trim {start:?}->{end:?} replacing {} fragment(s)",
                replacement_for.len()
            ),
        }
    }
}
