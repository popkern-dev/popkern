//! Typed repeated-apex partition plan and application sequence.

use super::*;
use popkern_core::topology::TopologyRole;

#[path = "repeated_apex_partition_execution_options.rs"]
mod execution_options;
pub(crate) use execution_options::RepeatedApexPartitionExecutionOptions;

pub(super) struct RepeatedApexPartitionPlan {
    cap_face: FaceHandle,
    cap_loop: LoopHandle,
    cap_boundary: Vec<PartitionBoundaryEdge>,
    retired_spans: Vec<RetiredApexSpan>,
}

impl RepeatedApexPartitionPlan {
    pub(super) fn build(
        brep: &BRep,
        cap_face: FaceHandle,
        loop_plan: &RepeatedApexLoopPlan,
        apex_vertex_handles: &[VertexHandle],
        face_edge_handles: &[EdgeHandle],
        n_fillets: usize,
    ) -> Result<Self, ApexBuildError> {
        if face_edge_handles.len() != n_fillets {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition needs {n_fillets} FaceEdges, got {}",
                    face_edge_handles.len()
                ),
            });
        }

        let mut cap_boundary = Vec::with_capacity(2 * n_fillets);
        for k in 0..n_fillets {
            let t_a = apex_vertex_handles[2 * k];
            let t_b = apex_vertex_handles[2 * k + 1];
            cap_boundary.push(PartitionBoundaryEdge::FilletEnd {
                fillet_index: k,
                t_a,
                t_b,
            });

            let next_k = (k + 1) % n_fillets;
            let expected_start = apex_vertex_handles[2 * k + 1];
            let expected_end = apex_vertex_handles[2 * next_k];
            let face_edge = face_edge_handles[k];
            let (cap_half_edge, neighbor_half_edge) = oriented_edge_half_edges_on_face(
                brep,
                cap_face,
                face_edge,
                expected_start,
                expected_end,
            )
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition: FaceEdge[{k}]={face_edge:?} is not a \
                     cap-side half-edge {expected_start:?}->{expected_end:?} on {cap_face:?}"
                ),
            })?;
            cap_boundary.push(PartitionBoundaryEdge::FaceEdge {
                fillet_index: k,
                edge: face_edge,
                cap_half_edge,
                neighbor_half_edge,
            });
        }

        let retired_spans = loop_plan
            .arcs
            .iter()
            .map(|arc| RetiredApexSpan {
                arc_index: arc.index,
                kind: arc.kind,
                half_edges: arc_half_edges_between(loop_plan, arc.start_pos, arc.end_pos),
            })
            .collect::<Vec<_>>();

        let partition = Self {
            cap_face,
            cap_loop: loop_plan.loop_handle,
            cap_boundary,
            retired_spans,
        };
        partition.validate(n_fillets)?;
        Ok(partition)
    }

    fn validate(&self, n_fillets: usize) -> Result<(), ApexBuildError> {
        if self.cap_boundary.len() != 2 * n_fillets {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition produced {} cap boundary entries; expected {}",
                    self.cap_boundary.len(),
                    2 * n_fillets
                ),
            });
        }
        for k in 0..n_fillets {
            match self.cap_boundary[2 * k] {
                PartitionBoundaryEdge::FilletEnd { fillet_index, .. } if fillet_index == k => {}
                _ => {
                    return Err(ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition cap boundary slot {} is not FilletEnd[{k}]",
                            2 * k
                        ),
                    });
                }
            }
            match self.cap_boundary[2 * k + 1] {
                PartitionBoundaryEdge::FaceEdge { fillet_index, .. } if fillet_index == k => {}
                _ => {
                    return Err(ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition cap boundary slot {} is not FaceEdge[{k}]",
                            2 * k + 1
                        ),
                    });
                }
            }
        }
        Ok(())
    }

    pub(super) fn diagnostic_summary(&self) -> String {
        let boundary = self
            .cap_boundary
            .iter()
            .map(|edge| match edge {
                PartitionBoundaryEdge::FilletEnd {
                    fillet_index,
                    t_a,
                    t_b,
                } => format!("FilletEnd[{fillet_index}]({t_a:?}->{t_b:?})"),
                PartitionBoundaryEdge::FaceEdge {
                    fillet_index,
                    edge,
                    cap_half_edge,
                    neighbor_half_edge,
                } => format!(
                    "FaceEdge[{fillet_index}]({edge:?}, cap={cap_half_edge:?}, neighbor={neighbor_half_edge:?})"
                ),
            })
            .collect::<Vec<_>>()
            .join(" -> ");
        let retired = self
            .retired_spans
            .iter()
            .map(|span| {
                format!(
                    "arc[{}] {:?}:{}he",
                    span.arc_index,
                    span.kind,
                    span.half_edges.len()
                )
            })
            .collect::<Vec<_>>()
            .join("; ");
        format!(
            "cap={:?}/{:?}; boundary=[{}]; retired=[{}]",
            self.cap_face, self.cap_loop, boundary, retired
        )
    }

    pub(super) fn apply(
        &self,
        brep: &mut BRep,
        plan: &ApexPlan,
        apex_vertex_handles: &[VertexHandle],
        fillet_surfaces: &[FilletSurfaceResult],
        n_fillets: usize,
        execution_options: RepeatedApexPartitionExecutionOptions,
    ) -> Result<
        (
            FaceHandle,
            Vec<EdgeHandle>,
            Vec<crate::types::ApexConstructionDiagnostic>,
        ),
        ApexBuildError,
    > {
        let mut diagnostics = Vec::new();
        let mut face_edge_arcs: Vec<Option<FaceEdgeOwnership>> = vec![None; n_fillets];
        let mut base_span: Option<&RetiredApexSpan> = None;
        for span in &self.retired_spans {
            match span.kind {
                RepeatedApexArcKind::FaceEdge { fillet_index } => {
                    let ownership = FaceEdgeOwnership::from_retired_span(
                        brep,
                        span,
                        self.cap_face,
                        fillet_index,
                        n_fillets,
                    )?;
                    face_edge_arcs[fillet_index] = Some(ownership);
                }
                RepeatedApexArcKind::Base => {
                    base_span = Some(span);
                }
                RepeatedApexArcKind::Unresolved => {
                    return Err(ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition cannot apply unresolved arc[{}]",
                            span.arc_index
                        ),
                    });
                }
            }
        }

        let face_edge_arcs = face_edge_arcs
            .into_iter()
            .enumerate()
            .map(|(k, arc)| {
                arc.ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition missing FaceEdge arc for fillet {k}"
                    ),
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let base_residual_plan = base_span
            .map(|span| plan_retired_base_residual_loop(brep, span))
            .transpose()?;

        let mut fillet_end_edges = Vec::with_capacity(n_fillets);
        let mut cap_cycle = Vec::with_capacity(2 * n_fillets);
        let mut fillet_end_ownership = Vec::with_capacity(n_fillets);

        for k in 0..n_fillets {
            let t_a = apex_vertex_handles[2 * k];
            let t_b = apex_vertex_handles[2 * k + 1];
            let (edge, cap_half_edge, cylinder_candidate_half_edge) =
                insert_unwired_edge_pair(brep, t_a, t_b);
            let fillet_index = planned_fillet_index_for_sector(plan, k)?;
            let curve = &plan.target_topology.apex_edges[2 * k].curve;
            attach_fillet_end_curve_to_edge(brep, edge, curve, fillet_index, fillet_surfaces)?;
            let target_cylinder_face =
                target_cylinder_face_for_fillet(brep, fillet_index, fillet_surfaces)?;

            let previous = (k + n_fillets - 1) % n_fillets;
            fillet_end_ownership.push(FilletEndOwnership {
                fillet_index,
                edge,
                cap_half_edge,
                cylinder_candidate_half_edge,
                target_cylinder_face,
                retired_apex_cell: RetiredApexCell {
                    previous_face_edge_index: previous,
                    next_face_edge_index: k,
                    t_a_to_apex: face_edge_arcs[previous].t_a_next_to_apex,
                    apex_to_t_b: face_edge_arcs[k].apex_to_t_b,
                },
            });
            cap_cycle.push(cap_half_edge);
            cap_cycle.push(face_edge_arcs[k].face_edge_cap_half_edge);
            fillet_end_edges.push(edge);
        }

        validate_partition_ownership_records(brep, &face_edge_arcs, &fillet_end_ownership)?;
        let mut boundary_outcome =
            PartitionBoundaryOutcome::from_current_cap_cycle(cap_cycle.clone());
        if !execution_options.apply_cylinder_trim_rewrite {
            if execution_options.require_cylinder_trim_rewrite_ready {
                validate_cylinder_trim_rewrite_readiness(brep, &fillet_end_ownership)?;
            }
            if execution_options.require_cylinder_trim_rebuild_plan {
                validate_cylinder_trim_rebuild_plans(brep, &fillet_end_ownership)?;
            }
        }
        if execution_options.apply_cylinder_trim_rewrite {
            let rewrite_result = apply_cylinder_trim_rewrites_transactionally(
                brep,
                self.cap_face,
                self.cap_loop,
                &face_edge_arcs,
                &fillet_end_ownership,
                execution_options.require_v_free_repeated_apex,
            )?;
            boundary_outcome = rewrite_result.boundary_outcome;
            diagnostics.extend(rewrite_result.diagnostics);
        }
        if !execution_options.apply_cylinder_trim_rewrite {
            if execution_options.adopt_cylinder_fillet_end_endpoints {
                adopt_fillet_end_endpoints_into_target_cylinders(brep, &fillet_end_ownership)?;
                refresh_fillet_end_curves_after_endpoint_adoption(
                    brep,
                    plan,
                    &fillet_end_ownership,
                    fillet_surfaces,
                )?;
            }
            if execution_options.require_shared_cylinder_fillet_end {
                validate_target_cylinder_loop_ownership_ready(brep, &fillet_end_ownership)?;
            }
        }
        wire_existing_loop(
            brep,
            self.cap_loop,
            self.cap_face,
            &boundary_outcome.cap_cycle,
        )?;

        let existing_target_cylinders_own_all_strips =
            fillet_end_ownership.iter().all(|ownership| {
                boundary_outcome.uses_existing_target_cylinder_strip(ownership.fillet_index)
            });
        for ownership in &fillet_end_ownership {
            if boundary_outcome.uses_existing_target_cylinder_strip(ownership.fillet_index) {
                validate_existing_target_cylinder_strip(
                    brep,
                    ownership,
                    plan.apex_vertex,
                    fillet_surfaces,
                )?;
                assign_fillet_strip_surface(
                    brep,
                    ownership.target_cylinder_face,
                    ownership.fillet_index,
                    fillet_surfaces,
                )?;
                let live_seam_fallback = boundary_outcome
                    .fillet_end_cap_side(ownership.fillet_index)
                    .and_then(|cap_side| {
                        brep.half_edges
                            .get(cap_side)
                            .map(|half_edge| half_edge.twin)
                    });
                let residual_loop = residual_strip_loop_with_duplicated_seam(
                    brep,
                    ownership,
                    boundary_outcome.strip_side_override(ownership.fillet_index),
                    live_seam_fallback,
                )?;
                let residual_face =
                    insert_face_with_outer_loop(brep, self.cap_face, &residual_loop)?;
                brep.faces[residual_face].role = TopologyRole::ConstructionResidual;
                brep.faces[residual_face].surface = None;
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace §14.9/partition] retired FilletEnd[{}] apex cell into \
                         unsurfaced residual face {residual_face:?}; production strip is target \
                         cylinder {:?}",
                        ownership.fillet_index, ownership.target_cylinder_face
                    );
                }
            } else {
                let strip_loop = boundary_outcome.strip_loop(ownership.fillet_index, ownership);
                protect_source_loops_before_partition_loop_claim(brep, &strip_loop)?;
                let strip_face = insert_face_with_outer_loop(brep, self.cap_face, &strip_loop)?;
                assign_fillet_strip_surface(
                    brep,
                    strip_face,
                    ownership.fillet_index,
                    fillet_surfaces,
                )?;
            }
        }

        if let Some(base_span) = base_span {
            let base_loop_half_edges = if base_span
                .half_edges
                .iter()
                .all(|half_edge| brep.half_edges.contains_key(*half_edge))
            {
                let apex = base_span
                    .half_edges
                    .first()
                    .and_then(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
                    .ok_or_else(|| ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition base span arc[{}] has no live \
                             start",
                            base_span.arc_index
                        ),
                    })?;
                let (seam_edge, seam_a, seam_b) = insert_unwired_edge_pair(brep, apex, apex);
                brep.edges[seam_edge].role = TopologyRole::ConstructionResidual;
                let mut live_loop_half_edges = base_span.half_edges.clone();
                live_loop_half_edges.push(seam_a);
                live_loop_half_edges.push(seam_b);
                live_loop_half_edges
            } else {
                base_residual_plan
                    .as_ref()
                    .ok_or_else(|| ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition base span arc[{}] was consumed \
                                 before durable residual authoring data was captured",
                            base_span.arc_index
                        ),
                    })?
                    .author_loop(brep)?
            };
            let base_face =
                insert_face_with_outer_loop(brep, self.cap_face, &base_loop_half_edges)?;
            brep.faces[base_face].role = TopologyRole::ConstructionResidual;
            brep.faces[base_face].surface = None;
            if base_span
                .half_edges
                .iter()
                .any(|half_edge| !brep.half_edges.contains_key(*half_edge))
            {
                mark_loop_edges_construction_residual(brep, &base_loop_half_edges);
            }
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/partition] retired base apex span into face {base_face:?}"
                );
            }
        }

        if existing_target_cylinders_own_all_strips {
            demote_retired_apex_residual_faces(
                brep,
                self.cap_face,
                plan.apex_vertex,
                &self.retired_spans,
                &fillet_end_ownership,
            );
        }

        // Carrier absorption (Stage-D seam-ownership inventory): planar
        // production faces take back the coplanar scaffold fragments
        // that separate them from their true seam partners.
        absorb_coplanar_residual_scaffold(brep);
        // Retired-subcomplex role coherence (see helper docs): must
        // run after every residual face demotion above.
        demote_residual_interior_edges(brep);
        // Handover completion: park genuinely orphaned loop claimers in
        // closed residual homes (must be the LAST loop-structure pass —
        // mid-transaction displaced claims are deliberate hand-offs, and
        // the default cube path's PK-074 pocket completion consumes
        // displaced claims even later, at dispatch: gate to the strict
        // V-free contract, whose scaffold this pass makes
        // loop-partitionable so the carrier-absorption welds can apply).
        if execution_options.require_v_free_repeated_apex {
            rehome_all_displaced_claimers(brep)?;
        }
        brep.normalize_edge_representatives();
        repair_stale_shell_anchors(brep)?;
        boundary_outcome.refresh_from_current_cap_loop(brep, self.cap_face)?;
        if boundary_outcome.matches_current_cap_cycle(&cap_cycle) {
            assert_normalized_cap_boundary(brep, self.cap_face, apex_vertex_handles, n_fillets)?;
        } else {
            assert_cap_boundary_matches_selected_cycle(
                brep,
                self.cap_face,
                &boundary_outcome.cap_cycle,
            )?;
            for (slot, ownership) in fillet_end_ownership.iter().enumerate() {
                let Some(cap_side) = boundary_outcome.fillet_end_cap_side(ownership.fillet_index)
                else {
                    continue;
                };
                let edge = brep
                    .half_edges
                    .get(cap_side)
                    .map(|half_edge| half_edge.edge);
                if let Some(edge) = edge {
                    fillet_end_edges[slot] = edge;
                }
            }
        }
        materialize_unowned_half_edge_residual_faces(brep, self.cap_face)?;
        validate_repeated_apex_normalization_structure(brep)?;
        validate_partition_cap_to_cylinder_seams(brep, self.cap_face, n_fillets)?;
        Ok((self.cap_face, fillet_end_edges, diagnostics))
    }
}

enum PartitionBoundaryEdge {
    FilletEnd {
        fillet_index: usize,
        t_a: VertexHandle,
        t_b: VertexHandle,
    },
    FaceEdge {
        fillet_index: usize,
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        neighbor_half_edge: HalfEdgeHandle,
    },
}

pub(super) struct RetiredApexSpan {
    pub(super) arc_index: usize,
    pub(super) kind: RepeatedApexArcKind,
    pub(super) half_edges: Vec<HalfEdgeHandle>,
}
