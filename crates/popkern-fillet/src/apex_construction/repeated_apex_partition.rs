//! Planned repeated-apex cap partitioning for constructive Gregory apexes.
//!
//! This module owns the non-serial graph rewrite used when a repeated-apex cap
//! loop cannot be normalized by any order of FilletEnd `mef`s. It builds a
//! typed partition from live FaceEdge half-edges, rewires the cap loop to the
//! planned FilletEnd/FaceEdge cycle, and validates topology plus semantic
//! cap-to-cylinder seam ownership before accepting.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveOnEdge, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, Vertex, VertexHandle,
};

use crate::fillet_surface::FilletSurfaceResult;

use super::constructive_geometry::{
    assign_fillet_strip_surface, attach_fillet_end_curve_to_edge, planned_fillet_index_for_sector,
};
use super::constructive_topology::{
    assert_normalized_cap_boundary, validate_repeated_apex_normalization_structure,
};
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::repeated_apex_loop_plan::{
    RepeatedApexArcKind, RepeatedApexLoopPlan, plan_repeated_apex_loop_normalization,
};
use super::repeated_apex_partition_topology::{
    insert_boundary_half_edge, insert_face_with_outer_loop, insert_unwired_edge_pair,
    oriented_edge_half_edges_on_face, validate_partition_cap_to_cylinder_seams, wire_existing_loop,
};
use super::repeated_apex_side_support::{
    SideSupportRegionDiagnosticSnapshot, apply_selected_side_support_continuation,
};
use super::repeated_apex_trim_subcomplex::{
    TopologyBalanceRequirement, TopologyCountDelta, TrimRewriteSubcomplexBuilder,
    TrimRewriteSubcomplexPlan, topology_count_summary,
};

#[path = "repeated_apex_partition_anchor_diagnostics.rs"]
mod partition_anchor_diagnostics;
#[path = "repeated_apex_partition_boundary_outcome.rs"]
mod partition_boundary_outcome;
#[path = "repeated_apex_partition_cap_boundary_cycles.rs"]
mod partition_cap_boundary_cycles;
#[path = "repeated_apex_partition_cap_boundary_rewrite.rs"]
mod partition_cap_boundary_rewrite;
#[path = "repeated_apex_partition_cylinder_anchors.rs"]
mod partition_cylinder_anchors;
#[path = "repeated_apex_partition_cylinder_geometry.rs"]
mod partition_cylinder_geometry;
#[path = "repeated_apex_partition_cylinder_loop_adoption.rs"]
mod partition_cylinder_loop_adoption;
#[path = "repeated_apex_partition_cylinder_trim_fragments.rs"]
mod partition_cylinder_trim_fragments;
#[path = "repeated_apex_partition_cylinder_trim_readiness.rs"]
mod partition_cylinder_trim_readiness;
#[path = "repeated_apex_partition_cylinder_trim_rebuild.rs"]
mod partition_cylinder_trim_rebuild;
#[path = "repeated_apex_partition_cylinder_trim_topology.rs"]
mod partition_cylinder_trim_topology;
#[path = "repeated_apex_partition_exact_trim_curves.rs"]
mod partition_exact_trim_curves;
#[path = "repeated_apex_partition_far_side_trim.rs"]
mod partition_far_side_trim;
#[path = "repeated_apex_partition_generator_anchors.rs"]
mod partition_generator_anchors;
#[path = "repeated_apex_partition_ownership.rs"]
mod partition_ownership;
#[path = "repeated_apex_partition_rebuild_diagnostics.rs"]
mod partition_rebuild_diagnostics;
#[path = "repeated_apex_partition_shared_fillet_end_trim.rs"]
mod partition_shared_fillet_end_trim;
#[path = "repeated_apex_partition_support_cycle_primitives.rs"]
mod partition_support_cycle_primitives;
#[path = "repeated_apex_partition_support_intersections.rs"]
mod partition_support_intersections;
#[path = "repeated_apex_partition_support_loop_rewrite.rs"]
mod partition_support_loop_rewrite;
#[path = "repeated_apex_partition_support_owners.rs"]
mod partition_support_owners;
#[path = "repeated_apex_partition_support_split_cycles.rs"]
mod partition_support_split_cycles;
#[path = "repeated_apex_partition_support_trim_materialization.rs"]
mod partition_support_trim_materialization;
#[path = "repeated_apex_partition_support_trim_sides.rs"]
mod partition_support_trim_sides;
#[path = "repeated_apex_partition_support_trim_validation.rs"]
mod partition_support_trim_validation;

use partition_anchor_diagnostics::generator_anchor_blocker;
use partition_boundary_outcome::PartitionBoundaryOutcome;
use partition_cap_boundary_cycles::{
    assert_cap_boundary_matches_selected_cycle, collect_loop_half_edges,
    intended_cap_boundary_cycle, reversed_twin_cycle, validate_half_edge_cycle,
};
use partition_cap_boundary_rewrite::plan_cap_boundary_rewrite_for_cylinder_trim;
use partition_cylinder_anchors::{
    coincident_vertex_on_face_outer_loop, containing_outer_loop_edge_by_curve,
    generator_anchor_on_cylinder_loop, materialize_support_boundary_generator_anchor,
    outer_loop_path_between_vertices, path_off_cylinder_vertices,
    support_topology_generator_anchor_pair, vertex_is_on_face_outer_loop,
    vertices_are_materialized_generator_anchors,
};
use partition_cylinder_geometry::{
    CylinderTrimSurface, arc_half_edges_between, author_missing_on_cylinder_line_trim_for_endpoint,
    cylinder_radial_vector, half_edge_end_vertex, point_lies_on_cylinder_face,
    target_cylinder_face_for_fillet,
};
use partition_cylinder_loop_adoption::adopt_fillet_end_endpoints_into_target_cylinders;
use partition_cylinder_trim_fragments::{
    FragmentAuthoringStatus, ReplacementRunReownershipPlan, ReplacementRunSpanOwnership,
    ReplacementRunSupportCandidate, RetiredCylinderTrimFragment, RetiredCylinderTrimRun,
};
use partition_cylinder_trim_rebuild::{
    SemanticCylinderTrimRebuildPlan, SemanticCylinderTrimSegment,
};
use partition_cylinder_trim_topology::{
    CylinderLoopOrientationCandidate, CylinderLoopOrientationKind, CylinderTrimRewriteTopologyPlan,
    FarSideTrimStrategy, FilletEndSeamRewrite, OrientedCylinderTrimSegment,
    PendingCylinderLoopReplacement,
};
use partition_exact_trim_curves::{
    exact_cylinder_arc_trim_curve_between, exact_cylinder_trim_curve_between,
    exact_support_cylinder_intersection_curve_between,
};
use partition_far_side_trim::materialized_far_side_trim_sides;
use partition_generator_anchors::materialize_missing_generator_anchors;
use partition_ownership::{
    FaceEdgeOwnership, FilletEndOwnership, RetiredApexCell, validate_partition_ownership_records,
    validate_target_cylinder_loop_ownership_ready,
};
use partition_rebuild_diagnostics::{
    adjacent_face_across_half_edge, face_surface_kind, fragment_replacement_authoring_status,
    replacement_run_curve_blocker, replacement_run_exact_curve_kind,
    replacement_run_reownership_blocker, replacement_run_reownership_plan,
    retired_trim_runs_by_adjacent_face,
};
use partition_shared_fillet_end_trim::materialized_shared_fillet_end_cylinder_half_edge;
use partition_support_cycle_primitives::{
    closed_half_edge_cycles_from_unwired_set, half_edge_is_degenerate,
    support_cycle_open_chain_diagnostic, support_split_half_edge_end_vertex,
    vertices_are_coincident,
};
use partition_support_intersections::{
    plane_plane_cylinder_joint_points, support_plane_faces_containing_point,
};
use partition_support_loop_rewrite::materialize_closed_support_trim_cycles;
use partition_support_owners::{
    GeneratorTrimSupportOwnerQuery, support_face_for_generator_trim_side,
    target_cylinder_support_faces, target_cylinder_support_faces_with_semantic,
    unique_far_side_support_face_for_trim_side,
};
use partition_support_split_cycles::materialize_support_split_cycles;
use partition_support_trim_materialization::materialized_generator_trim_side;
use partition_support_trim_sides::{
    PreMaterializedTrimSide, SupportOwnedTrimSide, SupportOwnedTrimSource,
    pre_materialize_far_side_trim_sides, pre_materialized_trim_sides,
};
use partition_support_trim_validation::validate_replacement_loop_support_side_ownership;

pub(super) fn try_build_planned_repeated_apex_partition(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    n_fillets: usize,
) -> Result<Option<(FaceHandle, Vec<EdgeHandle>)>, ApexBuildError> {
    let Some(loop_plan) = plan_repeated_apex_loop_normalization(
        brep,
        cap_face,
        plan.apex_vertex,
        apex_vertex_handles,
        n_fillets,
    ) else {
        return Ok(None);
    };
    if !loop_plan.is_commercial_grade_candidate(n_fillets) {
        return Ok(None);
    }

    let partition_plan = RepeatedApexPartitionPlan::build(
        brep,
        cap_face,
        &loop_plan,
        apex_vertex_handles,
        face_edge_handles,
        n_fillets,
    )?;
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace §14.9/partition] planned graph partition: {}",
            partition_plan.diagnostic_summary()
        );
    }

    partition_plan
        .apply(brep, plan, apex_vertex_handles, fillet_surfaces, n_fillets)
        .map(Some)
}

struct RepeatedApexPartitionPlan {
    cap_face: FaceHandle,
    cap_loop: popkern_core::topology::LoopHandle,
    cap_boundary: Vec<PartitionBoundaryEdge>,
    retired_spans: Vec<RetiredApexSpan>,
}

impl RepeatedApexPartitionPlan {
    fn build(
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

    fn diagnostic_summary(&self) -> String {
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

    fn apply(
        &self,
        brep: &mut BRep,
        plan: &ApexPlan,
        apex_vertex_handles: &[VertexHandle],
        fillet_surfaces: &[FilletSurfaceResult],
        n_fillets: usize,
    ) -> Result<(FaceHandle, Vec<EdgeHandle>), ApexBuildError> {
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
        if std::env::var_os("POPKERN_A2_APPLY_CYLINDER_TRIM_REWRITE").is_some() {
            boundary_outcome = apply_cylinder_trim_rewrites_transactionally(
                brep,
                self.cap_face,
                self.cap_loop,
                &face_edge_arcs,
                &fillet_end_ownership,
            )?;
        }
        if std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REWRITE_READY").is_some() {
            validate_cylinder_trim_rewrite_readiness(brep, &fillet_end_ownership)?;
        }
        if std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REBUILD_PLAN").is_some() {
            validate_cylinder_trim_rebuild_plans(brep, &fillet_end_ownership)?;
        }
        if std::env::var_os("POPKERN_A2_ADOPT_CYLINDER_FILLET_END_ENDPOINTS").is_some() {
            adopt_fillet_end_endpoints_into_target_cylinders(brep, &fillet_end_ownership)?;
        }
        if std::env::var_os("POPKERN_A2_REQUIRE_SHARED_CYLINDER_FILLET_END").is_some() {
            validate_target_cylinder_loop_ownership_ready(brep, &fillet_end_ownership)?;
        }
        wire_existing_loop(
            brep,
            self.cap_loop,
            self.cap_face,
            &boundary_outcome.cap_cycle,
        )?;

        for ownership in &fillet_end_ownership {
            let strip_face = insert_face_with_outer_loop(
                brep,
                self.cap_face,
                &boundary_outcome.strip_loop(ownership.fillet_index, ownership),
            )?;
            assign_fillet_strip_surface(brep, strip_face, ownership.fillet_index, fillet_surfaces)?;
        }

        if let Some(base_span) = base_span {
            let apex = base_span
                .half_edges
                .first()
                .and_then(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition base span arc[{}] has no live start",
                        base_span.arc_index
                    ),
                })?;
            let (_seam_edge, seam_a, seam_b) = insert_unwired_edge_pair(brep, apex, apex);
            let mut base_loop_half_edges = base_span.half_edges.clone();
            base_loop_half_edges.push(seam_a);
            base_loop_half_edges.push(seam_b);
            let base_face =
                insert_face_with_outer_loop(brep, self.cap_face, &base_loop_half_edges)?;
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace §14.9/partition] retired base apex span into face {base_face:?}"
                );
            }
        }

        brep.normalize_edge_representatives();
        if boundary_outcome.matches_current_cap_cycle(&cap_cycle) {
            assert_normalized_cap_boundary(brep, self.cap_face, apex_vertex_handles, n_fillets)?;
        } else {
            assert_cap_boundary_matches_selected_cycle(
                brep,
                self.cap_face,
                &boundary_outcome.cap_cycle,
            )?;
        }
        validate_repeated_apex_normalization_structure(brep)?;
        validate_partition_cap_to_cylinder_seams(brep, self.cap_face, n_fillets)?;
        Ok((self.cap_face, fillet_end_edges))
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

struct RetiredApexSpan {
    arc_index: usize,
    kind: RepeatedApexArcKind,
    half_edges: Vec<HalfEdgeHandle>,
}

fn apply_materialized_cylinder_loop_replacements(
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
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
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

fn validate_cylinder_trim_rebuild_plans(
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

fn validate_cylinder_trim_rewrite_readiness(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    let rebuild_plans = plan_semantic_cylinder_trim_rebuilds(brep, fillet_ends)?;
    for rebuild_plan in rebuild_plans {
        rebuild_plan.validate_rewrite_readiness()?;
    }
    Ok(())
}

fn apply_cylinder_trim_rewrites_transactionally(
    brep: &mut BRep,
    cap_face: FaceHandle,
    cap_loop: popkern_core::topology::LoopHandle,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Result<PartitionBoundaryOutcome, ApexBuildError> {
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
        let subcomplex_plan = plan_trim_subcomplex_for_cylinder_rewrite(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
        )?;
        validate_cylinder_trim_rewrite_topology_preflight(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            &subcomplex_plan,
        )?;
        let pre_materialized = pre_materialize_far_side_trim_sides(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
        )?;
        let side_support_continuation =
            apply_selected_side_support_continuation(brep, &subcomplex_plan)?;
        apply_materialized_cylinder_loop_replacements(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
            &subcomplex_plan,
            &pre_materialized,
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
        boundary_outcome.plan_strip_loop_overrides(brep, fillet_ends)?;
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/cylinder-trim-rewrite] topology_plans=[{}]; side_support_continuation={}; {}; {}",
                topology_plans
                    .iter()
                    .map(CylinderTrimRewriteTopologyPlan::diagnostic_label)
                    .collect::<Vec<_>>()
                    .join("; "),
                side_support_continuation.diagnostic_label(),
                cap_boundary_plan.diagnostic_label(),
                boundary_outcome.diagnostic_label()
            );
        }
        cap_boundary_plan.apply_adjacent_loop_rewrites(brep)?;
        Ok(boundary_outcome)
    })();
    if result.is_err() {
        *brep = rollback_snapshot;
    }
    result
}

fn validate_cylinder_trim_rewrite_topology_preflight(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    subcomplex_plan: &TrimRewriteSubcomplexPlan,
) -> Result<(), ApexBuildError> {
    let side_support_diagnostics =
        SideSupportRegionDiagnosticSnapshot::build_read_only(brep, subcomplex_plan);
    debug_assert!(side_support_diagnostics.validate_fresh(brep));
    let before = topology_count_summary(brep);
    let mut candidate = brep.clone();
    let pre_materialized = pre_materialize_far_side_trim_sides(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
    )?;
    let side_support_continuation =
        apply_selected_side_support_continuation(&mut candidate, subcomplex_plan)?;
    match apply_materialized_cylinder_loop_replacements(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
        subcomplex_plan,
        &pre_materialized,
    ) {
        Ok(()) => {
            let after = topology_count_summary(&candidate);
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace §14.9/topology-preflight] accepted before {}; after {}; delta {}; side_support_continuation={}; {}; {}",
                    before.diagnostic_label(),
                    after.diagnostic_label(),
                    TopologyCountDelta::between(before, after).diagnostic_label(),
                    side_support_continuation.diagnostic_label(),
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
        if !fragment.materialized_anchor_replacement {
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

fn plan_semantic_cylinder_trim_rebuilds(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<Vec<SemanticCylinderTrimRebuildPlan>, ApexBuildError> {
    fillet_ends
        .iter()
        .map(|fillet_end| SemanticCylinderTrimRebuildPlan::build(brep, fillet_end))
        .collect()
}

impl SemanticCylinderTrimSegment {
    fn diagnostic_label(&self) -> String {
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
