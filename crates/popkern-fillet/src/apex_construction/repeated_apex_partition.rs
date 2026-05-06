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
#[path = "repeated_apex_partition_cylinder_anchors.rs"]
mod partition_cylinder_anchors;
#[path = "repeated_apex_partition_cylinder_geometry.rs"]
mod partition_cylinder_geometry;
#[path = "repeated_apex_partition_exact_trim_curves.rs"]
mod partition_exact_trim_curves;
#[path = "repeated_apex_partition_ownership.rs"]
mod partition_ownership;
#[path = "repeated_apex_partition_rebuild_diagnostics.rs"]
mod partition_rebuild_diagnostics;
#[path = "repeated_apex_partition_support_intersections.rs"]
mod partition_support_intersections;
#[path = "repeated_apex_partition_support_owners.rs"]
mod partition_support_owners;

use partition_anchor_diagnostics::generator_anchor_blocker;
use partition_boundary_outcome::PartitionBoundaryOutcome;
use partition_cap_boundary_cycles::{
    assert_cap_boundary_matches_selected_cycle, collect_loop_half_edges,
    intended_cap_boundary_cycle, reversed_twin_cycle, validate_half_edge_cycle,
};
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
use partition_exact_trim_curves::{
    exact_cylinder_arc_trim_curve_between, exact_cylinder_trim_curve_between,
    exact_support_cylinder_intersection_curve_between,
};
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
use partition_support_intersections::{
    plane_plane_cylinder_joint_points, support_plane_faces_containing_point,
};
use partition_support_owners::{
    GeneratorTrimSupportOwnerQuery, support_face_for_generator_trim_side,
    target_cylinder_support_faces, target_cylinder_support_faces_with_semantic,
    unique_far_side_support_face_for_trim_side,
};

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

fn adopt_fillet_end_endpoints_into_target_cylinders(
    brep: &mut BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
        let start = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        let end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has \
                     stale twin",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            }
        })?;
        adopt_vertex_into_target_cylinder_loop(
            brep,
            fillet_end.target_cylinder_face,
            start,
            fillet_end.fillet_index,
        )?;
        adopt_vertex_into_target_cylinder_loop(
            brep,
            fillet_end.target_cylinder_face,
            end,
            fillet_end.fillet_index,
        )?;
    }
    Ok(())
}

fn materialize_missing_generator_anchors(
    brep: &mut BRep,
    cap_face: FaceHandle,
    fillet_ends: &[FilletEndOwnership],
    face_edges: &[FaceEdgeOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
        let allowed_support_faces = target_cylinder_support_faces(
            brep,
            fillet_end.target_cylinder_face,
            cap_face,
            face_edges,
        );
        let start = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        let end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has \
                     stale twin",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            }
        })?;
        let start_point = brep
            .vertices
            .get(start)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] start endpoint {start:?} \
                         is stale while materializing generator anchors",
                    fillet_end.fillet_index
                ),
            })?;
        let end_point = brep
            .vertices
            .get(end)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] end endpoint {end:?} is \
                         stale while materializing generator anchors",
                    fillet_end.fillet_index
                ),
            })?;
        let start_has_anchor =
            generator_anchor_on_cylinder_loop(brep, fillet_end.target_cylinder_face, start_point)
                .is_some();
        let end_has_anchor =
            generator_anchor_on_cylinder_loop(brep, fillet_end.target_cylinder_face, end_point)
                .is_some();
        if start_has_anchor && end_has_anchor {
            continue;
        }
        let (start_anchor, end_anchor) = support_topology_generator_anchor_pair(
            brep,
            fillet_end.target_cylinder_face,
            start_point,
            end_point,
            &allowed_support_faces,
        )
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot materialize support-aware \
                 generator anchors on target cylinder face {:?}: start {}; end {}",
                fillet_end.fillet_index,
                fillet_end.target_cylinder_face,
                generator_anchor_blocker(brep, fillet_end.target_cylinder_face, start_point),
                generator_anchor_blocker(brep, fillet_end.target_cylinder_face, end_point)
            ),
        })?;
        for (endpoint, has_anchor, anchor_point) in [
            (start, start_has_anchor, start_anchor),
            (end, end_has_anchor, end_anchor),
        ] {
            if has_anchor {
                continue;
            }
            let anchor = materialize_support_boundary_generator_anchor(
                brep,
                fillet_end.target_cylinder_face,
                endpoint,
                &anchor_point,
                &allowed_support_faces,
                fillet_end.fillet_index,
            )?;
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                let support_boundary_hits = support_plane_faces_containing_point(
                    brep,
                    fillet_end.target_cylinder_face,
                    anchor_point.point,
                )
                .into_iter()
                .map(|support_face| {
                    let boundary_edge =
                        containing_outer_loop_edge_by_curve(brep, support_face, anchor_point.point);
                    format!("{support_face:?}: boundary_edge={boundary_edge:?}")
                })
                .collect::<Vec<_>>()
                .join(", ");
                eprintln!(
                    "[a2-trace §14.9/cylinder-anchor] materialized support-aware FilletEnd[{}] \
                     anchor {anchor:?} at ({:.6},{:.6},{:.6}) for endpoint {endpoint:?} on {:?}; \
                     support hits [{}]",
                    fillet_end.fillet_index,
                    anchor_point.point.x,
                    anchor_point.point.y,
                    anchor_point.point.z,
                    fillet_end.target_cylinder_face,
                    support_boundary_hits
                );
            }
        }
    }
    Ok(())
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

#[derive(Clone)]
struct SupportOwnedTrimSide {
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    fillet_index: usize,
    half_edge: HalfEdgeHandle,
    support_side_owned: bool,
    support_source: Option<SupportOwnedTrimSource>,
    requires_support_ownership: bool,
}

#[derive(Clone, Copy)]
struct SupportOwnedTrimSource {
    loop_: LoopHandle,
    prev: HalfEdgeHandle,
    next: HalfEdgeHandle,
}

#[derive(Clone)]
struct PreMaterializedTrimSide {
    cylinder_face: FaceHandle,
    fillet_index: usize,
    start: VertexHandle,
    end: VertexHandle,
    trim_sides: Vec<SupportOwnedTrimSide>,
}

fn pre_materialized_trim_sides(
    pre_materialized: &[PreMaterializedTrimSide],
    cylinder_face: FaceHandle,
    fillet_index: usize,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<SupportOwnedTrimSide>> {
    pre_materialized
        .iter()
        .find(|entry| {
            entry.cylinder_face == cylinder_face
                && entry.fillet_index == fillet_index
                && entry.start == start
                && entry.end == end
        })
        .map(|entry| entry.trim_sides.clone())
}

fn pre_materialize_far_side_trim_sides(
    brep: &mut BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> Result<Vec<PreMaterializedTrimSide>, ApexBuildError> {
    let mut materialized = Vec::new();
    for (rebuild_plan, topology_plan) in rebuild_plans.iter().zip(topology_plans) {
        let Some(fragment) =
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)
        else {
            continue;
        };
        if !fragment.materialized_anchor_replacement {
            continue;
        }
        let allowed_support_faces =
            target_cylinder_support_faces(brep, rebuild_plan.cylinder_face, cap_face, face_edges);
        for segment in &topology_plan.loop_candidate.segments {
            if !segment.label.contains("far-side") {
                continue;
            }
            let trim_sides = materialized_far_side_trim_sides(
                brep,
                rebuild_plan.cylinder_face,
                segment.start,
                segment.end,
                segment.label,
                rebuild_plan.fillet_index,
                &allowed_support_faces,
            )?;
            materialized.push(PreMaterializedTrimSide {
                cylinder_face: rebuild_plan.cylinder_face,
                fillet_index: rebuild_plan.fillet_index,
                start: segment.start,
                end: segment.end,
                trim_sides,
            });
        }
    }
    Ok(materialized)
}

fn materialize_closed_support_trim_cycles(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    let rollback = brep.clone();
    let result = materialize_closed_support_trim_cycles_candidate(brep, face_edges, trim_sides);
    if let Err(err) = result {
        *brep = rollback;
        return Err(err);
    }
    validate_live_half_edge_cycles_after_support_rewrite(brep)?;
    Ok(())
}

fn validate_live_half_edge_cycles_after_support_rewrite(brep: &BRep) -> Result<(), ApexBuildError> {
    for (loop_handle, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = std::collections::BTreeSet::new();
        for _ in 0..=brep.half_edges.len() {
            let half_edge =
                brep.half_edges
                    .get(cursor)
                    .ok_or_else(|| ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: missing half-edge {cursor:?}"
                        ),
                    })?;
            if half_edge.loop_ != loop_handle {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: half-edge {cursor:?} belongs to \
                         {:?}",
                        half_edge.loop_
                    ),
                });
            }
            if !visited.insert(cursor) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: repeated half-edge {cursor:?}"
                    ),
                });
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        if cursor != start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support ownership rewrite produced invalid \
                     live loop {loop_handle:?}: cycle did not close from {start:?}"
                ),
            });
        }
    }
    Ok(())
}

fn repair_support_source_loops_after_cylinder_rewire(
    brep: &mut BRep,
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        let Some(source) = trim_side.support_source else {
            continue;
        };
        if !brep.loops.contains_key(source.loop_) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot repair \
                     support loop after cylinder trim ownership transfer: stale source loop \
                     {loop_handle:?}",
                    fillet_index = trim_side.fillet_index,
                    loop_handle = source.loop_
                ),
            });
        }
        if !brep.half_edges.contains_key(source.prev)
            || !brep.half_edges.contains_key(source.next)
            || !brep.half_edges.contains_key(trim_side.half_edge)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot repair \
                     support loop after cylinder trim ownership transfer: stale source links \
                     prev={prev:?} moved={moved:?} next={next:?}",
                    fillet_index = trim_side.fillet_index,
                    prev = source.prev,
                    moved = trim_side.half_edge,
                    next = source.next
                ),
            });
        }
        let prev_still_points_to_moved = brep.half_edges[source.prev].next == trim_side.half_edge;
        let next_still_points_to_moved = brep.half_edges[source.next].prev == trim_side.half_edge;
        if !prev_still_points_to_moved || !next_still_points_to_moved {
            continue;
        }
        let replacement = author_support_source_loop_replacement_span(brep, trim_side, source)?;
        brep.half_edges[source.prev].next = replacement;
        brep.half_edges[source.next].prev = replacement;
        brep.half_edges[replacement].prev = source.prev;
        brep.half_edges[replacement].next = source.next;
        brep.half_edges[replacement].loop_ = source.loop_;
        if brep.loops[source.loop_].half_edge == trim_side.half_edge {
            brep.loops[source.loop_].half_edge = replacement;
        }
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/support-loop-repair] FilletEnd[{}] replaced stolen cylinder \
                 side {:?} in source support loop {:?} with boundary span {:?} between {:?}->{:?}",
                trim_side.fillet_index,
                trim_side.half_edge,
                source.loop_,
                replacement,
                source.prev,
                source.next
            );
        }
    }
    Ok(())
}

fn author_support_source_loop_replacement_span(
    brep: &mut BRep,
    trim_side: &SupportOwnedTrimSide,
    source: SupportOwnedTrimSource,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let start = brep
        .half_edges
        .get(trim_side.half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale moved half-edge",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let end = brep
        .half_edges
        .get(source.next)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale next half-edge {:?}",
                trim_side.fillet_index, trim_side.half_edge, source.next
            ),
        })?;
    let curve = brep
        .half_edges
        .get(trim_side.half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: moved edge has no curve",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let (edge, replacement) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    if support_split_half_edge_end_vertex(brep, replacement) != Some(end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] authored support replacement \
                 span with wrong endpoint: expected {start:?}->{end:?}, got {replacement:?}",
                trim_side.fillet_index
            ),
        });
    }
    Ok(replacement)
}

fn materialize_closed_support_trim_cycles_candidate(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    trace_topology_counts(brep, "support-loop-entry");
    repair_support_source_loops_after_cylinder_rewire(brep, trim_sides)?;
    let mut by_support_face: std::collections::BTreeMap<FaceHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    let mut generated_trim_count_by_support_face =
        std::collections::BTreeMap::<FaceHandle, usize>::new();
    let mut support_faces_with_existing_boundary = std::collections::BTreeSet::new();
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        if brep
            .half_edges
            .get(trim_side.half_edge)
            .and_then(|half_edge| brep.loops.get(half_edge.loop_))
            .is_none_or(|loop_data| loop_data.face != trim_side.cylinder_face)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for {cylinder_half_edge:?}: it is not owned by target \
                     cylinder face {cylinder_face:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge,
                    cylinder_face = trim_side.cylinder_face
                ),
            });
        }
        let twin = brep
            .half_edges
            .get(trim_side.half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for stale half-edge {cylinder_half_edge:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge
                ),
            })?;
        if trim_side.support_side_owned {
            continue;
        }
        by_support_face
            .entry(trim_side.support_face)
            .or_default()
            .push(twin);
        *generated_trim_count_by_support_face
            .entry(trim_side.support_face)
            .or_default() += 1;
    }
    for face_edge in face_edges {
        if let Some(half_edges) = by_support_face.get_mut(&face_edge.support_face)
            && !half_edge_is_degenerate(brep, face_edge.support_half_edge)
        {
            half_edges.push(face_edge.support_half_edge);
            support_faces_with_existing_boundary.insert(face_edge.support_face);
        }
    }

    let mut support_groups = by_support_face.into_iter().collect::<Vec<_>>();
    support_groups.sort_by_key(|(support_face, _)| {
        let generated_count = generated_trim_count_by_support_face
            .get(support_face)
            .copied()
            .unwrap_or_default();
        (
            std::cmp::Reverse(generated_count),
            support_faces_with_existing_boundary.contains(support_face),
            *support_face,
        )
    });

    for (support_face, mut half_edges) in support_groups {
        trace_topology_counts(brep, "support-face-before");
        half_edges.retain(|half_edge| !half_edge_is_degenerate(brep, *half_edge));
        let support_split_plan = materialize_support_split_cycles(brep, support_face, &half_edges)?;
        if !support_split_plan.is_complete() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: generated trim sides {:?} do not form closed support \
                     cycle(s) yet. Interior support subdivision must add the missing support \
                     boundary spans before these cylinder trims can be exposed. {}; {}",
                    half_edges,
                    support_cycle_open_chain_diagnostic(brep, &half_edges),
                    support_split_plan.diagnostic()
                ),
            });
        }
        if support_split_plan.cycles.len() != 1 {
            let counts = topology_count_summary(brep);
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: got {} closed support cycles. Multi-cycle support \
                     subdivision requires an Euler-balanced residual face/ring operator; direct \
                     inner-loop wiring or duplicate-face insertion would violate shell topology. \
                     Current topology deficit before this support face is {}; one face-loop per \
                     closed support cycle would add only {}.",
                    support_split_plan.cycles.len(),
                    counts.deficit_to_genus_zero(),
                    support_split_plan.cycles.len()
                ),
            });
        }
        let cycle = &support_split_plan.cycles[0];
        validate_half_edge_cycle(brep, cycle, "materialized support trim cycle", support_face)?;
        let loop_handle = brep
            .faces
            .get(support_face)
            .map(|face| face.outer_loop)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot replace stale support face \
                     {support_face:?}"
                ),
            })?;
        wire_existing_loop(brep, loop_handle, support_face, cycle)?;
        trace_topology_counts(brep, "support-face-after");
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/support-loop-replace] support {support_face:?} -> \
                 rewired outer loop {:?} using {:?}",
                loop_handle, cycle
            );
        }
    }
    Ok(())
}

fn trace_topology_counts(brep: &BRep, label: &'static str) {
    if std::env::var_os("POPKERN_A2_TRACE").is_none() {
        return;
    }
    eprintln!("{}", topology_count_summary(brep).trace_label(label));
}

struct SupportSplitCyclePlan {
    cycles: Vec<Vec<HalfEdgeHandle>>,
    missing_spans: Vec<(VertexHandle, VertexHandle)>,
}

impl SupportSplitCyclePlan {
    fn is_complete(&self) -> bool {
        self.missing_spans.is_empty()
    }

    fn diagnostic(&self) -> String {
        if self.missing_spans.is_empty() {
            return "support split cycles are closed".to_string();
        }
        let spans = self
            .missing_spans
            .iter()
            .map(|(start, end)| format!("{start:?}->{end:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "support split planner requires authored two-sided span(s) [{}] before it can create \
             support-owned loops",
            spans
        )
    }
}

fn materialize_support_split_cycles(
    brep: &mut BRep,
    support_face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Result<SupportSplitCyclePlan, ApexBuildError> {
    if let Some(cycle) =
        support_replacement_outer_loop_from_boundary_order(brep, support_face, half_edges)
        && half_edge_cycle_uses_unique_edges(&cycle)
        && validate_half_edge_cycle(
            brep,
            &cycle,
            "boundary-ordered support replacement loop",
            support_face,
        )
        .is_ok()
    {
        return Ok(SupportSplitCyclePlan {
            cycles: vec![cycle],
            missing_spans: Vec::new(),
        });
    }
    let cycles = closed_half_edge_cycles_from_unwired_set(brep, half_edges);
    let covered = cycles.iter().map(Vec::len).sum::<usize>();
    if covered == half_edges.len() {
        return Ok(SupportSplitCyclePlan {
            cycles,
            missing_spans: Vec::new(),
        });
    }
    let missing_spans = missing_support_split_spans(brep, half_edges)
        .into_iter()
        .filter(|(start, end)| !vertices_are_coincident(brep, *start, *end))
        .collect::<Vec<_>>();
    if missing_spans.is_empty() {
        return Ok(SupportSplitCyclePlan {
            cycles: Vec::new(),
            missing_spans: Vec::new(),
        });
    }
    let mut supplemented = half_edges.to_vec();
    for (start, end) in &missing_spans {
        let half_edge = author_support_split_closing_span(brep, support_face, *start, *end)?;
        supplemented.push(half_edge);
    }
    let cycles = closed_half_edge_cycles_from_unwired_set(brep, &supplemented);
    let covered = cycles.iter().map(Vec::len).sum::<usize>();
    if covered == supplemented.len() {
        return Ok(SupportSplitCyclePlan {
            cycles,
            missing_spans: Vec::new(),
        });
    }
    trace_support_split_incomplete_cycles(
        brep,
        support_face,
        &missing_spans,
        &supplemented,
        &cycles,
    );
    Ok(SupportSplitCyclePlan {
        cycles: Vec::new(),
        missing_spans,
    })
}

fn author_support_split_closing_span(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    validate_support_split_span_geometry(brep, support_face, start, end)?;
    let rollback = brep.clone();
    let result = (|| {
        let start_point = brep
            .vertices
            .get(start)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author support closing span on \
                         {support_face:?}: stale start vertex {start:?}"
                ),
            })?;
        let end_point = brep
            .vertices
            .get(end)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author support closing span on \
                     {support_face:?}: stale end vertex {end:?}"
                ),
            })?;
        let curve = CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author degenerate support closing \
                     span {start:?}->{end:?} on {support_face:?}"
                ),
            })?;
        let before_edges = brep.edges.len();
        let before_half_edges = brep.half_edges.len();
        let before_faces = brep.faces.len();
        let before_loops = brep.loops.len();
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(curve);
        if brep.edges.len() != before_edges + 1
            || brep.half_edges.len() != before_half_edges + 1
            || brep.faces.len() != before_faces
            || brep.loops.len() != before_loops
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span {start:?}->{end:?} on \
                     {support_face:?} violated topology delta dV=0 dE=1 dHE=1 dF=0 dL=0"
                ),
            });
        }
        if !brep.boundary_edges.contains(&edge) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span {start:?}->{end:?} on \
                     {support_face:?} did not mark its unowned side as boundary"
                ),
            });
        }
        let actual_start = brep.half_edges[half_edge].origin;
        let actual_end = support_split_half_edge_end_vertex(brep, half_edge);
        if actual_start != start || actual_end != Some(end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support closing span on {support_face:?} \
                     has wrong orientation: expected {start:?}->{end:?}, got \
                     {actual_start:?}->{actual_end:?}"
                ),
            });
        }
        Ok(half_edge)
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn trace_support_split_incomplete_cycles(
    brep: &BRep,
    support_face: FaceHandle,
    missing_spans: &[(VertexHandle, VertexHandle)],
    supplemented: &[HalfEdgeHandle],
    cycles: &[Vec<HalfEdgeHandle>],
) {
    if std::env::var_os("POPKERN_A2_TRACE").is_none() {
        return;
    }
    let spans = missing_spans
        .iter()
        .map(|(start, end)| format!("{start:?}->{end:?}"))
        .collect::<Vec<_>>()
        .join(", ");
    let supplemented_edges = support_cycle_open_chain_diagnostic(brep, supplemented);
    let cycle_sizes = cycles.iter().map(Vec::len).collect::<Vec<_>>();
    eprintln!(
        "[a2-trace §14.9/support-split] support {support_face:?} supplemented spans [{spans}] \
         did not close all half-edges; supplemented=[{supplemented_edges}] cycle_sizes={cycle_sizes:?}"
    );
}

fn support_replacement_outer_loop_from_boundary_order(
    brep: &BRep,
    support_face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
) -> Option<Vec<HalfEdgeHandle>> {
    if half_edges.is_empty() {
        return None;
    }
    let boundary_vertices =
        super::constructive_topology::collect_outer_loop_vertices(brep, support_face)?;
    let boundary_index = boundary_vertices
        .iter()
        .enumerate()
        .map(|(index, vertex)| (*vertex, index))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mut ordered = half_edges
        .iter()
        .filter_map(|half_edge| {
            let start = brep.half_edges.get(*half_edge)?.origin;
            let end = half_edge_end_vertex(brep, *half_edge)?;
            let start_index = boundary_index.get(&start).copied()?;
            boundary_index.get(&end)?;
            Some((start_index, *half_edge))
        })
        .collect::<Vec<_>>();
    if ordered.len() != half_edges.len() {
        return None;
    }
    ordered.sort_by_key(|(start_index, half_edge)| (*start_index, *half_edge));
    let ordered_half_edges = ordered
        .iter()
        .map(|(_, half_edge)| *half_edge)
        .collect::<Vec<_>>();
    let mut cycle = Vec::new();
    for (index, half_edge) in ordered_half_edges.iter().enumerate() {
        cycle.push(*half_edge);
        let end = half_edge_end_vertex(brep, *half_edge)?;
        let next = ordered_half_edges[(index + 1) % ordered_half_edges.len()];
        let next_start = brep.half_edges.get(next)?.origin;
        if end == next_start {
            continue;
        }
        let path =
            support_outer_loop_forward_path_between_vertices(brep, support_face, end, next_start)?;
        cycle.extend(path);
    }
    Some(cycle)
}

fn half_edge_cycle_uses_unique_edges(cycle: &[HalfEdgeHandle]) -> bool {
    let unique = cycle
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    unique.len() == cycle.len()
}

fn support_outer_loop_forward_path_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    support_outer_loop_paths_between_vertices(brep, support_face, start, end)
        .into_iter()
        .next()
}

fn support_outer_loop_paths_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Vec<Vec<HalfEdgeHandle>> {
    let Some(face) = brep.faces.get(support_face) else {
        return Vec::new();
    };
    let Ok(loop_edges) = collect_loop_half_edges(brep, face.outer_loop) else {
        return Vec::new();
    };
    let n = loop_edges.len();
    if n == 0 {
        return Vec::new();
    }
    let mut paths = Vec::new();
    for start_index in 0..n {
        if brep
            .half_edges
            .get(loop_edges[start_index])
            .map(|he| he.origin)
            != Some(start)
        {
            continue;
        }
        let mut path = Vec::new();
        for offset in 0..n {
            let index = (start_index + offset) % n;
            path.push(loop_edges[index]);
            let next_index = (index + 1) % n;
            if brep
                .half_edges
                .get(loop_edges[next_index])
                .map(|he| he.origin)
                == Some(end)
            {
                paths.push(path);
                break;
            }
        }
    }
    paths
}

fn missing_support_split_spans(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<(VertexHandle, VertexHandle)> {
    let mut starts = std::collections::BTreeMap::<VertexHandle, usize>::new();
    let mut ends = std::collections::BTreeMap::<VertexHandle, usize>::new();
    for half_edge in half_edges {
        let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) else {
            continue;
        };
        let Some(end) = half_edge_end_vertex(brep, *half_edge) else {
            continue;
        };
        *starts.entry(start).or_default() += 1;
        *ends.entry(end).or_default() += 1;
    }
    let mut open_starts = starts
        .iter()
        .filter_map(|(vertex, count)| {
            let end_count = ends.get(vertex).copied().unwrap_or_default();
            (count > &end_count).then_some((*vertex, count - end_count))
        })
        .collect::<Vec<_>>();
    let mut open_ends = ends
        .iter()
        .filter_map(|(vertex, count)| {
            let start_count = starts.get(vertex).copied().unwrap_or_default();
            (count > &start_count).then_some((*vertex, count - start_count))
        })
        .collect::<Vec<_>>();
    let mut spans = Vec::new();
    for (end_vertex, end_count) in open_ends.drain(..) {
        for _ in 0..end_count {
            let Some((start_vertex, start_count)) = open_starts.first_mut() else {
                return spans;
            };
            spans.push((end_vertex, *start_vertex));
            *start_count -= 1;
            if *start_count == 0 {
                open_starts.remove(0);
            }
        }
    }
    spans
}

fn validate_support_split_span_geometry(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<(), ApexBuildError> {
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let Some(SurfaceGeometry::Plane { origin, normal, .. }) = brep
        .faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 support face {support_face:?} has no plane surface"
            ),
        });
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    let start_distance = (start_point - *origin).dot(normal).abs();
    let end_distance = (end_point - *origin).dot(normal).abs();
    if start_distance > tol || end_distance > tol {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 endpoints are off support face {support_face:?} by {:.3e}/{:.3e}",
                start_distance, end_distance
            ),
        });
    }
    if start_point.distance_to(&end_point) <= tol {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author support split span {start:?}->{end:?}: \
                 degenerate endpoints"
            ),
        });
    }
    Ok(())
}

fn half_edge_is_degenerate(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(start) = brep.half_edges.get(half_edge).map(|data| data.origin) else {
        return true;
    };
    let Some(end) = half_edge_end_vertex(brep, half_edge) else {
        return true;
    };
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return true;
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return true;
    };
    start_point.distance_to(&end_point) <= brep.tolerance.linear().max(1.0e-9)
}

fn vertices_are_coincident(brep: &BRep, start: VertexHandle, end: VertexHandle) -> bool {
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return false;
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return false;
    };
    start_point.distance_to(&end_point) <= brep.tolerance.linear().max(1.0e-9)
}

fn support_cycle_open_chain_diagnostic(brep: &BRep, half_edges: &[HalfEdgeHandle]) -> String {
    half_edges
        .iter()
        .map(|half_edge| {
            let start = brep
                .half_edges
                .get(*half_edge)
                .map(|data| data.origin)
                .unwrap_or_default();
            let end = support_split_half_edge_end_vertex(brep, *half_edge).unwrap_or_default();
            format!("{half_edge:?}:{start:?}->{end:?}")
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn support_split_half_edge_end_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<VertexHandle> {
    half_edge_end_vertex(brep, half_edge)
        .or_else(|| boundary_curve_end_vertex(brep, half_edge))
        .or_else(|| {
            let next = brep.half_edges.get(half_edge)?.next;
            brep.half_edges
                .get(next)
                .map(|next_half_edge| next_half_edge.origin)
        })
}

fn boundary_curve_end_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    if !brep.boundary_edges.contains(&half_edge_data.edge) {
        return None;
    }
    let curve = brep.edges.get(half_edge_data.edge)?.curve.as_ref()?;
    let start_point = brep.vertices.get(half_edge_data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let end_point = if start_point.distance_to(&p0) <= start_point.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&end_point);
            (distance <= tol).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
}

fn closed_half_edge_cycles_from_unwired_set(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<Vec<HalfEdgeHandle>> {
    let mut by_start: std::collections::BTreeMap<VertexHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    for half_edge in half_edges {
        if let Some(start) = brep.half_edges.get(*half_edge).map(|data| data.origin) {
            by_start.entry(start).or_default().push(*half_edge);
        }
    }
    if by_start.values().any(|edges| edges.len() != 1) {
        return Vec::new();
    }
    let input = half_edges
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut visited = std::collections::BTreeSet::new();
    let mut cycles = Vec::new();
    for start_half_edge in half_edges {
        if visited.contains(start_half_edge) {
            continue;
        }
        let mut cycle = Vec::new();
        let mut cursor = *start_half_edge;
        for _ in 0..=half_edges.len() {
            if !input.contains(&cursor) || !visited.insert(cursor) {
                return Vec::new();
            }
            cycle.push(cursor);
            let Some(end) = support_split_half_edge_end_vertex(brep, cursor) else {
                return Vec::new();
            };
            let Some(next) = support_split_next_half_edge_from_end(brep, &by_start, end) else {
                return Vec::new();
            };
            if next == *start_half_edge {
                cycles.push(cycle);
                break;
            }
            cursor = next;
        }
    }
    if visited.len() == half_edges.len() {
        cycles
    } else {
        Vec::new()
    }
}

fn support_split_next_half_edge_from_end(
    brep: &BRep,
    by_start: &std::collections::BTreeMap<VertexHandle, Vec<HalfEdgeHandle>>,
    end: VertexHandle,
) -> Option<HalfEdgeHandle> {
    if let Some(next) = by_start.get(&end).and_then(|edges| edges.first()).copied() {
        return Some(next);
    }
    let coincident = by_start
        .iter()
        .filter_map(|(start, edges)| {
            vertices_are_coincident(brep, end, *start).then_some(edges.first().copied())
        })
        .flatten()
        .collect::<Vec<_>>();
    if coincident.len() == 1 {
        Some(coincident[0])
    } else {
        None
    }
}

fn validate_replacement_loop_support_side_ownership(
    brep: &BRep,
    cylinder_face: FaceHandle,
    fillet_index: usize,
    cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let mut missing = Vec::new();
    for half_edge in cycle {
        let Some(twin) = brep.half_edges.get(*half_edge).map(|data| data.twin) else {
            missing.push(format!("{half_edge:?}: stale cylinder half-edge"));
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(twin) else {
            missing.push(format!("{half_edge:?}: stale twin {twin:?}"));
            continue;
        };
        let Some(twin_loop) = brep.loops.get(twin_data.loop_) else {
            missing.push(format!(
                "{half_edge:?}/twin {twin:?}: no support-side loop owns the new trim side; {}",
                replacement_trim_support_split_diagnostic(brep, cylinder_face, *half_edge)
            ));
            continue;
        };
        if twin_loop.face == cylinder_face {
            missing.push(format!(
                "{half_edge:?}/twin {twin:?}: twin is still owned by target cylinder face \
                 {cylinder_face:?}"
            ));
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot accept \
             materialized target-cylinder loop yet: support-side ownership is incomplete for \
             replacement trim side(s): {}. The next construction slice must split or create the \
             adjacent support loops/faces that own these exact cylinder trims before exposing the \
             rewritten cylinder face.",
            missing.join("; ")
        ),
    })
}

fn replacement_trim_support_split_diagnostic(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> String {
    let Some(start) = brep.half_edges.get(half_edge).map(|data| data.origin) else {
        return "stale trim start prevents support split planning".to_string();
    };
    let Some(end) = half_edge_end_vertex(brep, half_edge) else {
        return "stale trim end prevents support split planning".to_string();
    };
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return format!("stale trim start vertex {start:?} prevents support split planning");
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return format!("stale trim end vertex {end:?} prevents support split planning");
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut support_faces = Vec::new();
    for (face, face_data) in brep.faces.iter() {
        if face == cylinder_face {
            continue;
        }
        let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
            face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let start_distance = (start_point - *origin).dot(normal).abs();
        let end_distance = (end_point - *origin).dot(normal).abs();
        if start_distance > tol || end_distance > tol {
            continue;
        }
        let start_on_boundary = vertex_is_on_face_outer_loop(brep, face, start);
        let end_on_boundary = vertex_is_on_face_outer_loop(brep, face, end);
        let start_edge = containing_outer_loop_edge_by_curve(brep, face, start_point);
        let end_edge = containing_outer_loop_edge_by_curve(brep, face, end_point);
        support_faces.push(format!(
            "{face:?}: start_boundary={start_on_boundary}, end_boundary={end_on_boundary}, \
             start_edge={start_edge:?}, end_edge={end_edge:?}"
        ));
    }
    if support_faces.is_empty() {
        return format!(
            "no real support plane contains trim endpoints {start:?}->{end:?}; construction must \
             refuse rather than inventing a support owner"
        );
    }
    format!(
        "candidate support split face(s) [{}]; anchors that are neither existing boundary vertices \
         nor splittable boundary-edge points require an interior support-face subdivision before \
         the cylinder loop can be exposed",
        support_faces.join("; ")
    )
}

fn existing_shared_fillet_end_half_edge_for_segment(
    brep: &BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let (_edge, cap_half_edge, cylinder_half_edge) = rebuild_plan.shared_fillet_end_half_edges()?;
    [cap_half_edge, cylinder_half_edge]
        .into_iter()
        .find(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .is_some_and(|data| data.origin == start)
                && half_edge_end_vertex(brep, *half_edge)
                    .is_some_and(|actual_end| actual_end == end)
        })
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot find existing shared \
                 FilletEnd half-edge {:?}->{:?} while replacing target-cylinder loop",
                rebuild_plan.fillet_index, start, end
            ),
        })
}

fn materialized_shared_fillet_end_cylinder_half_edge(
    brep: &mut BRep,
    rebuild_plan: &SemanticCylinderTrimRebuildPlan,
    topology_plan: &CylinderTrimRewriteTopologyPlan,
    start: VertexHandle,
    end: VertexHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let existing =
        existing_shared_fillet_end_half_edge_for_segment(brep, rebuild_plan, start, end)?;
    if !matches!(
        topology_plan.seam_rewrite,
        FilletEndSeamRewrite::ReuseCurrentTwinPair { .. }
    ) {
        return Ok(existing);
    }
    let original_edge = brep
        .half_edges
        .get(existing)
        .map(|half_edge| half_edge.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot duplicate stale shared \
                 cylinder seam {existing:?}",
                rebuild_plan.fillet_index
            ),
        })?;
    let curve = brep
        .edges
        .get(original_edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot duplicate shared cylinder \
                 seam {existing:?}: original edge {original_edge:?} has no curve",
                rebuild_plan.fillet_index
            ),
        })?;
    let (edge, half_edge) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    if support_split_half_edge_end_vertex(brep, half_edge) != Some(end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] duplicated shared cylinder seam \
                 has wrong orientation: expected {start:?}->{end:?}, got {half_edge:?}",
                rebuild_plan.fillet_index
            ),
        });
    }
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace §14.9/cylinder-loop-replace] FilletEnd[{}] duplicated shared cylinder \
             seam {:?}->{:?}: original {:?}, cylinder boundary {:?}",
            rebuild_plan.fillet_index, start, end, existing, half_edge
        );
    }
    Ok(half_edge)
}

fn materialized_generator_trim_side(
    brep: &mut BRep,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<
    (
        HalfEdgeHandle,
        FaceHandle,
        bool,
        Option<SupportOwnedTrimSource>,
    ),
    ApexBuildError,
> {
    let curve =
        exact_cylinder_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index)?;
    match materialize_support_face_mef_trim(
        brep,
        support_face,
        start,
        end,
        curve.clone(),
        label,
        fillet_index,
    ) {
        Ok((half_edge, actual_support_face, support_source)) => {
            return Ok((half_edge, actual_support_face, true, Some(support_source)));
        }
        Err(err) => {
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} through support-face mef on {support_face:?}: \
                     {err}"
                );
            }
        }
    }
    let (edge, half_edge, _twin) = insert_unwired_edge_pair(brep, start, end);
    brep.edges[edge].curve = Some(curve);
    Ok((half_edge, support_face, false, None))
}

fn materialize_support_face_mef_trim(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<(HalfEdgeHandle, FaceHandle, SupportOwnedTrimSource), ApexBuildError> {
    let rollback = brep.clone();
    let result = (|| {
        let (actual_support_face, support_start, support_end) =
            support_residual_face_split_half_edges(brep, support_face, start, end).ok_or_else(
                || ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                         {label} {start:?}->{end:?} by support-face mef: no residual support \
                         face has both endpoints on one outer loop. {}",
                        support_family_endpoint_diagnostic(brep, support_face, start, end)
                    ),
                },
            )?;
        let result =
            brep.mef(support_start, support_end)
                .map_err(|err| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} by support-face mef on {actual_support_face:?}: \
                     {err:?}"
                    ),
                })?;
        brep.edges[result.edge].curve = Some(curve);
        let support_source = SupportOwnedTrimSource {
            loop_: brep.half_edges[result.half_edge].loop_,
            prev: brep.half_edges[result.half_edge].prev,
            next: brep.half_edges[result.half_edge].next,
        };
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} by support-face mef on {actual_support_face:?}; cylinder \
                 side {:?}, support side {:?}",
                result.half_edge, result.half_edge_twin
            );
        }
        Ok((result.half_edge, actual_support_face, support_source))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn support_family_endpoint_diagnostic(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> String {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    let reports = candidate_faces
        .into_iter()
        .map(|face| {
            let start_on_loop = vertex_is_on_face_outer_loop(brep, face, start);
            let end_on_loop = vertex_is_on_face_outer_loop(brep, face, end);
            let start_edge = brep.vertices.get(start).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            let end_edge = brep.vertices.get(end).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            format!(
                "{face:?}: start_on_loop={start_on_loop} end_on_loop={end_on_loop} \
                 start_edge={start_edge:?} end_edge={end_edge:?}"
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("support family endpoint status [{reports}]")
}

fn support_residual_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    candidate_faces.into_iter().find_map(|candidate_face| {
        support_face_split_half_edges(brep, candidate_face, start, end)
            .map(|(support_start, support_end)| (candidate_face, support_start, support_end))
    })
}

fn support_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face = brep.faces.get(support_face)?;
    let loop_handle = face.outer_loop;
    let loop_data = brep.loops.get(loop_handle)?;
    let mut start_edges = Vec::new();
    let mut end_edges = Vec::new();
    let mut cursor = loop_data.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        if half_edge.origin == start {
            start_edges.push(cursor);
        }
        if half_edge.origin == end {
            end_edges.push(cursor);
        }
        cursor = half_edge.next;
        if cursor == loop_data.half_edge {
            break;
        }
    }
    match (start_edges.as_slice(), end_edges.as_slice()) {
        ([support_start], [support_end]) if support_start != support_end => {
            Some((*support_start, *support_end))
        }
        _ => None,
    }
}

fn coplanar_residual_support_faces(brep: &BRep, support_face: FaceHandle) -> Vec<FaceHandle> {
    let Some(source_face) = brep.faces.get(support_face) else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: source_origin,
        normal: source_normal,
        ..
    }) = source_face
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == support_face || face_data.shell != source_face.shell {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let normals_parallel = source_normal.cross(normal).length() <= tol;
            let same_plane = (*origin - *source_origin).dot(source_normal).abs() <= tol;
            (normals_parallel && same_plane).then_some(face)
        })
        .collect()
}

fn unique_support_face_for_trim_endpoints(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<FaceHandle, ApexBuildError> {
    let start_point = brep
        .vertices
        .get(start)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose support \
                 owner for {label}: stale start vertex {start:?}"
            ),
        })?;
    let end_point = brep
        .vertices
        .get(end)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose support \
                 owner for {label}: stale end vertex {end:?}"
            ),
        })?;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face || !allowed_support_faces.contains(&face) {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let start_distance = (start_point - *origin).dot(normal).abs();
            let end_distance = (end_point - *origin).dot(normal).abs();
            (start_distance <= tol && end_distance <= tol).then_some(face)
        })
        .collect::<Vec<_>>();
    candidates.sort();
    candidates.dedup();
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose unique \
             support owner for {label} {start:?}->{end:?}: candidates {:?} from allowed faces \
             {:?}",
            candidates, allowed_support_faces
        ),
    })
}

fn materialized_far_side_trim_sides(
    brep: &mut BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<Vec<SupportOwnedTrimSide>, ApexBuildError> {
    if let Ok(support_face) = unique_support_face_for_trim_endpoints(
        brep,
        cylinder_face,
        start,
        end,
        label,
        fillet_index,
        allowed_support_faces,
    ) {
        let (half_edge, actual_support_face, support_side_owned, support_source) =
            materialized_generator_trim_side(
                brep,
                cylinder_face,
                support_face,
                start,
                end,
                label,
                fillet_index,
            )?;
        return Ok(vec![SupportOwnedTrimSide {
            cylinder_face,
            support_face: actual_support_face,
            fillet_index,
            half_edge,
            support_side_owned,
            support_source,
            requires_support_ownership: true,
        }]);
    }

    let (edge, half_edge, twin) = insert_unwired_edge_pair(brep, start, end);
    match exact_cylinder_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index) {
        Ok(curve) => {
            match unique_far_side_support_face_for_trim_side(
                brep,
                cylinder_face,
                half_edge,
                allowed_support_faces,
            ) {
                Ok(support_face) => {
                    brep.edges[edge].curve = Some(curve);
                    return Ok(vec![SupportOwnedTrimSide {
                        cylinder_face,
                        support_face,
                        fillet_index,
                        half_edge,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    }]);
                }
                Err(_) => {
                    brep.half_edges.remove(half_edge);
                    if brep.half_edges.contains_key(twin) {
                        brep.half_edges.remove(twin);
                    }
                    brep.edges.remove(edge);
                }
            }
        }
        Err(_) => {
            brep.half_edges.remove(half_edge);
            if brep.half_edges.contains_key(twin) {
                brep.half_edges.remove(twin);
            }
            brep.edges.remove(edge);
        }
    }

    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let start_supports = support_plane_faces_containing_point(brep, cylinder_face, start_point)
        .into_iter()
        .filter(|support_face| allowed_support_faces.contains(support_face))
        .collect::<Vec<_>>();
    let end_supports = support_plane_faces_containing_point(brep, cylinder_face, end_point)
        .into_iter()
        .filter(|support_face| allowed_support_faces.contains(support_face))
        .collect::<Vec<_>>();
    let mut candidates = Vec::new();
    for start_support in &start_supports {
        for end_support in &end_supports {
            for joint_point in
                plane_plane_cylinder_joint_points(brep, cylinder_face, *start_support, *end_support)
            {
                let score =
                    joint_point.distance_to(&start_point) + joint_point.distance_to(&end_point);
                candidates.push((*start_support, *end_support, joint_point, score));
            }
        }
    }
    candidates.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    if let Some((start_support, end_support, joint_point, _)) = candidates.into_iter().next() {
        let joint = brep.vertices.insert(Vertex {
            position: joint_point,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: Some(brep.tolerance.linear()),
        });
        let (edge_a, half_edge_a, _twin_a) = insert_unwired_edge_pair(brep, start, joint);
        let (edge_b, half_edge_b, _twin_b) = insert_unwired_edge_pair(brep, joint, end);
        let curve_a = exact_support_cylinder_intersection_curve_between(
            brep,
            cylinder_face,
            start,
            joint,
            start_point,
            joint_point,
            label,
            fillet_index,
        );
        let curve_b = exact_support_cylinder_intersection_curve_between(
            brep,
            cylinder_face,
            joint,
            end,
            joint_point,
            end_point,
            label,
            fillet_index,
        );
        match (curve_a, curve_b) {
            (Ok(curve_a), Ok(curve_b)) => {
                brep.edges[edge_a].curve = Some(curve_a);
                brep.edges[edge_b].curve = Some(curve_b);
                if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                    eprintln!(
                        "[a2-trace §14.9/cylinder-loop-replace] split far-side FilletEnd[{}] \
                         through joint {joint:?} using supports {:?}->{:?}",
                        fillet_index, start_support, end_support
                    );
                }
                return Ok(vec![
                    SupportOwnedTrimSide {
                        cylinder_face,
                        support_face: start_support,
                        fillet_index,
                        half_edge: half_edge_a,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    },
                    SupportOwnedTrimSide {
                        cylinder_face,
                        support_face: end_support,
                        fillet_index,
                        half_edge: half_edge_b,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    },
                ]);
            }
            (left, right) => {
                let left_status = match left {
                    Ok(_) => "ok".to_string(),
                    Err(err) => format!("{err:?}"),
                };
                let right_status = match right {
                    Ok(_) => "ok".to_string(),
                    Err(err) => format!("{err:?}"),
                };
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] found \
                         plane-plane-cylinder joint {joint:?} for {label}, but could not author \
                         both exact support/cylinder sub-spans: left {start:?}->{joint:?} \
                         via {start_support:?}: {left_status}; right {joint:?}->{end:?} via \
                         {end_support:?}: {right_status}"
                    ),
                });
            }
        }
    }
    if let Ok(curve) =
        exact_cylinder_arc_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index)
    {
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(curve);
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/cylinder-loop-replace] FilletEnd[{fillet_index}] authored \
                 cylinder-only far-side boundary arc {start:?}->{end:?}; support transition \
                 has no plane-plane-cylinder joint"
            );
        }
        return Ok(vec![SupportOwnedTrimSide {
            cylinder_face,
            support_face: start_supports
                .first()
                .copied()
                .or_else(|| end_supports.first().copied())
                .unwrap_or_default(),
            fillet_index,
            half_edge,
            support_side_owned: false,
            support_source: None,
            requires_support_ownership: false,
        }]);
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot split exact \
             far-side cylinder trim {start:?}->{end:?}: no plane-plane-cylinder transition \
             point from start supports {:?} to end supports {:?}",
            start_supports, end_supports
        ),
    })
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

fn plan_cap_boundary_rewrite_for_cylinder_trim(
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
        let neighbor_face = neighbor_loop
            .and_then(|loop_handle| brep.loops.get(loop_handle).map(|loop_data| loop_data.face));
        if neighbor_face == Some(cap_face) {
            continue;
        }
        rewrites.push(AdjacentLoopRewriteRequirement {
            stolen_cap_side: *stolen_cap_side,
            replacement_neighbor_side,
            neighbor_loop,
            neighbor_face,
            current_neighbor_cycle: if let Some(loop_handle) = neighbor_loop {
                collect_loop_half_edges(brep, loop_handle)?
            } else {
                Vec::new()
            },
            reauthored_neighbor_cycle: if let Some(loop_handle) = neighbor_loop {
                let current = collect_loop_half_edges(brep, loop_handle)?;
                let reauthored = reversed_twin_cycle(brep, &current, cap_face)?;
                validate_half_edge_cycle(
                    brep,
                    &reauthored,
                    "reauthored adjacent boundary",
                    neighbor_face.unwrap_or(cap_face),
                )?;
                reauthored
            } else {
                Vec::new()
            },
            reason: if neighbor_loop.is_some() {
                "live adjacent loop must be reauthored because replacing a stolen half-edge with \
                 its twin reverses local loop direction"
            } else {
                "new FilletEnd twin is not yet assigned to a live adjacent loop; caller must \
                 attach its twin to the later fillet strip loop"
            },
        });
    }
    Ok(rewrites)
}

#[derive(Clone)]
struct SemanticCylinderTrimRebuildPlan {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    start_anchor: VertexHandle,
    cap_start: VertexHandle,
    cap_end: VertexHandle,
    end_anchor: VertexHandle,
    segments: Vec<SemanticCylinderTrimSegment>,
    retired_fragments: Vec<RetiredCylinderTrimFragment>,
}

#[derive(Clone)]
enum SemanticCylinderTrimSegment {
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

#[derive(Clone)]
struct RetiredCylinderTrimFragment {
    path: Vec<HalfEdgeHandle>,
    vertices: Vec<VertexHandle>,
    adjacent_face_runs: Vec<RetiredCylinderTrimRun>,
    adjacent_faces: Vec<FaceHandle>,
    off_surface_vertices: Vec<VertexHandle>,
    materialized_anchor_replacement: bool,
}

#[derive(Clone)]
struct RetiredCylinderTrimRun {
    adjacent_face: FaceHandle,
    path: Vec<HalfEdgeHandle>,
    start: VertexHandle,
    end: VertexHandle,
    vertices: Vec<VertexHandle>,
    off_surface_vertices: Vec<VertexHandle>,
}

struct ReplacementRunReownershipPlan {
    spans: Vec<ReplacementRunSpanOwnership>,
}

struct ReplacementRunSpanOwnership {
    start: VertexHandle,
    end: VertexHandle,
    candidate_faces: Vec<ReplacementRunSupportCandidate>,
}

struct ReplacementRunSupportCandidate {
    face: FaceHandle,
    curve_kind: &'static str,
}

struct FragmentAuthoringStatus {
    authorable: bool,
    diagnostic: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FarSideTrimStrategy {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CylinderLoopOrientationKind {
    CurrentSharedFilletEnd,
    ReauthoredSharedFilletEnd,
}

#[derive(Clone)]
struct CylinderLoopOrientationCandidate {
    kind: CylinderLoopOrientationKind,
    far_side_strategy: FarSideTrimStrategy,
    segments: Vec<OrientedCylinderTrimSegment>,
}

struct CylinderTrimRewriteTopologyPlan {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    loop_candidate: CylinderLoopOrientationCandidate,
    seam_rewrite: FilletEndSeamRewrite,
}

struct PendingCylinderLoopReplacement {
    fillet_index: usize,
    cylinder_face: FaceHandle,
    kind: CylinderLoopOrientationKind,
    cycle: Vec<HalfEdgeHandle>,
}

enum FilletEndSeamRewrite {
    ReuseCurrentTwinPair {
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        current_cylinder_half_edge: HalfEdgeHandle,
    },
    ReauthorCapBoundaryAndCylinderSide {
        edge: EdgeHandle,
        cap_half_edge: HalfEdgeHandle,
        current_cylinder_half_edge: HalfEdgeHandle,
    },
}

struct CapBoundaryRewritePlan {
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
    reason: &'static str,
}

#[derive(Clone)]
struct OrientedCylinderTrimSegment {
    label: &'static str,
    start: VertexHandle,
    end: VertexHandle,
}

impl SemanticCylinderTrimRebuildPlan {
    fn build(brep: &BRep, fillet_end: &FilletEndOwnership) -> Result<Self, ApexBuildError> {
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
                 rebuild: no generator anchor on target cylinder face {:?} for end endpoint \
                 {:?}: {}",
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

    fn oriented_cylinder_cycle(
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

    fn current_cylinder_order_is_constructive(&self) -> bool {
        self.retired_fragments.first().is_some_and(|fragment| {
            !fragment.materialized_anchor_replacement && fragment.off_surface_vertices.is_empty()
        })
    }

    fn validate_rewrite_readiness(&self) -> Result<(), ApexBuildError> {
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

    fn validate_replacement_curve_authoring_ready(
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

    fn fragment_for_strategy(
        &self,
        strategy: FarSideTrimStrategy,
    ) -> Option<&RetiredCylinderTrimFragment> {
        match strategy {
            FarSideTrimStrategy::Forward => self.retired_fragments.first(),
            FarSideTrimStrategy::Reverse => self.retired_fragments.get(1),
        }
    }

    fn plan_cylinder_trim_rewrite_topology(
        &self,
        brep: &BRep,
    ) -> Result<CylinderTrimRewriteTopologyPlan, ApexBuildError> {
        let loop_candidate = self
            .selected_cylinder_loop_orientation_candidate(brep)?
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

    fn shared_fillet_end_half_edges(
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
        let Some(strategy) = self.selected_far_side_strategy(brep) else {
            return Ok(None);
        };
        let candidates = self.cylinder_loop_orientation_candidates(brep, strategy)?;
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
    ) -> Result<Vec<CylinderLoopOrientationCandidate>, ApexBuildError> {
        let mut candidates = Vec::new();
        candidates.push(CylinderLoopOrientationCandidate {
            kind: CylinderLoopOrientationKind::CurrentSharedFilletEnd,
            far_side_strategy: strategy,
            segments: self.oriented_cylinder_cycle_for_strategy(brep, strategy)?,
        });
        if strategy == FarSideTrimStrategy::Reverse {
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

    fn forward_off_surface_vertices(&self) -> &[VertexHandle] {
        self.retired_fragments
            .first()
            .map(|fragment| fragment.off_surface_vertices.as_slice())
            .unwrap_or(&[])
    }

    fn reverse_path_status(&self) -> String {
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

    fn segment_summary(&self) -> String {
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

impl SemanticCylinderTrimSegment {
    fn endpoints(&self) -> (VertexHandle, VertexHandle) {
        match self {
            SemanticCylinderTrimSegment::Generator { start, end }
            | SemanticCylinderTrimSegment::SharedFilletEnd { start, end, .. }
            | SemanticCylinderTrimSegment::FarSideCylinderTrim { start, end, .. } => (*start, *end),
        }
    }

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

impl RetiredCylinderTrimFragment {
    fn from_path(brep: &BRep, face: FaceHandle, path: Vec<HalfEdgeHandle>) -> Self {
        let vertices = path
            .iter()
            .filter_map(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .collect();
        let adjacent_face_runs = retired_trim_runs_by_adjacent_face(brep, face, &path);
        let adjacent_faces = path
            .iter()
            .filter_map(|half_edge| {
                let twin = brep.half_edges.get(*half_edge)?.twin;
                let twin_loop = brep.half_edges.get(twin)?.loop_;
                let twin_face = brep.loops.get(twin_loop)?.face;
                (twin_face != face).then_some(twin_face)
            })
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        let off_surface_vertices = path_off_cylinder_vertices(brep, face, &path);
        Self {
            path,
            vertices,
            adjacent_face_runs,
            adjacent_faces,
            off_surface_vertices,
            materialized_anchor_replacement: false,
        }
    }

    fn materialized_anchor_replacement(start: VertexHandle, end: VertexHandle) -> Self {
        Self {
            path: Vec::new(),
            vertices: vec![start, end],
            adjacent_face_runs: Vec::new(),
            adjacent_faces: Vec::new(),
            off_surface_vertices: Vec::new(),
            materialized_anchor_replacement: true,
        }
    }

    fn diagnostic_label(&self) -> String {
        if self.materialized_anchor_replacement {
            return format!(
                "materialized-anchor replacement {:?}->{:?}",
                self.vertices.first().copied().unwrap_or_default(),
                self.vertices.get(1).copied().unwrap_or_default()
            );
        }
        format!(
            "{} half-edge(s), vertices {:?}, adjacent faces {:?}, runs [{}], \
             off-cylinder vertices {:?}",
            self.path.len(),
            self.vertices,
            self.adjacent_faces,
            self.adjacent_face_runs
                .iter()
                .map(RetiredCylinderTrimRun::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.off_surface_vertices
        )
    }
}

impl RetiredCylinderTrimRun {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}: {} half-edge(s), {:?}->{:?}, vertices {:?}, off-cylinder {:?}, action {}",
            self.adjacent_face,
            self.path.len(),
            self.start,
            self.end,
            self.vertices,
            self.off_surface_vertices,
            if self.requires_replacement() {
                "replace"
            } else {
                "reuse"
            }
        )
    }

    fn requires_replacement(&self) -> bool {
        !self.off_surface_vertices.is_empty()
    }
}

impl ReplacementRunReownershipPlan {
    fn diagnostic_label(&self) -> String {
        self.spans
            .iter()
            .map(ReplacementRunSpanOwnership::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

impl ReplacementRunSpanOwnership {
    fn diagnostic_label(&self) -> String {
        let candidates = self
            .candidate_faces
            .iter()
            .map(|candidate| format!("{:?}/{}", candidate.face, candidate.curve_kind))
            .collect::<Vec<_>>()
            .join("|");
        format!(
            "{:?}->{:?} candidates [{}]",
            self.start, self.end, candidates
        )
    }
}

impl CylinderTrimRewriteTopologyPlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] cylinder {:?}: {}, {}",
            self.fillet_index,
            self.cylinder_face,
            self.loop_candidate.diagnostic_label(),
            self.seam_rewrite.diagnostic_label()
        )
    }
}

impl FilletEndSeamRewrite {
    fn diagnostic_label(&self) -> String {
        match self {
            FilletEndSeamRewrite::ReuseCurrentTwinPair {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => format!(
                "reuse current FilletEnd cap/cylinder twin pair edge={edge:?} \
                 cap_he={cap_half_edge:?} current_cylinder_he={current_cylinder_half_edge:?}"
            ),
            FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                edge,
                cap_half_edge,
                current_cylinder_half_edge,
            } => format!(
                "reauthor cap boundary and cylinder side for FilletEnd seam edge={edge:?} \
                 cap_he={cap_half_edge:?} current_cylinder_he={current_cylinder_half_edge:?}"
            ),
        }
    }
}

impl CapBoundaryRewritePlan {
    fn apply_adjacent_loop_rewrites(&self, brep: &mut BRep) -> Result<(), ApexBuildError> {
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
            if rewrite.reauthored_neighbor_cycle.is_empty() {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition adjacent loop rewrite has no \
                         reauthored cycle: {}",
                        rewrite.diagnostic_label()
                    ),
                });
            }
            wire_existing_loop(brep, loop_handle, face, &rewrite.reauthored_neighbor_cycle)?;
        }
        Ok(())
    }

    fn boundary_outcome(
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

    fn validate_adjacent_loop_rewrite_readiness(&self) -> Result<(), ApexBuildError> {
        if !self.required {
            return Ok(());
        }
        let missing_adjacent_cycles = self
            .adjacent_loop_rewrites
            .iter()
            .filter(|rewrite| {
                rewrite.neighbor_loop.is_some() && rewrite.reauthored_neighbor_cycle.is_empty()
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

    fn validate_rewrite_postconditions(
        &self,
        brep: &BRep,
        face_edges: &[FaceEdgeOwnership],
        topology_plans: &[CylinderTrimRewriteTopologyPlan],
    ) -> Result<(), ApexBuildError> {
        if !self.required {
            return Ok(());
        }
        for plan in topology_plans {
            let FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide {
                cap_half_edge,
                current_cylinder_half_edge,
                ..
            } = plan.seam_rewrite
            else {
                continue;
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
        }
        for face_edge in face_edges {
            if !self
                .reauthored_cap_cycle
                .contains(&face_edge.support_half_edge)
            {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cap-boundary rewrite for FaceEdge[{}] \
                         does not include support half-edge {:?} as the reauthored cap side",
                        face_edge.fillet_index, face_edge.support_half_edge
                    ),
                });
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
        }
        Ok(())
    }

    fn diagnostic_label(&self) -> String {
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
             neighbor_face={:?}, current_neighbor_cycle={:?}, reauthored_neighbor_cycle={:?}: {}",
            self.stolen_cap_side,
            self.replacement_neighbor_side,
            self.neighbor_loop,
            self.neighbor_face,
            self.current_neighbor_cycle,
            self.reauthored_neighbor_cycle,
            self.reason
        )
    }
}

impl CylinderLoopOrientationCandidate {
    fn requires_fillet_end_reorientation(&self) -> bool {
        self.kind == CylinderLoopOrientationKind::ReauthoredSharedFilletEnd
    }

    fn validate_contiguous(&self, fillet_index: usize) -> Result<(), String> {
        for i in 0..self.segments.len() {
            let current = &self.segments[i];
            let next = &self.segments[(i + 1) % self.segments.len()];
            if current.end != next.start {
                return Err(format!(
                    "FilletEnd[{fillet_index}] segment {} ends at {:?}, but segment {} starts \
                     at {:?}",
                    current.label, current.end, next.label, next.start
                ));
            }
        }
        Ok(())
    }

    fn diagnostic_label(&self) -> String {
        let segments = self
            .segments
            .iter()
            .map(OrientedCylinderTrimSegment::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "{:?}/{:?} reorient_fillet_end={} segments=[{}]",
            self.kind,
            self.far_side_strategy,
            self.requires_fillet_end_reorientation(),
            segments
        )
    }
}

impl OrientedCylinderTrimSegment {
    fn diagnostic_label(&self) -> String {
        format!("{} {:?}->{:?}", self.label, self.start, self.end)
    }
}

fn adopt_vertex_into_target_cylinder_loop(
    brep: &mut BRep,
    target_cylinder_face: FaceHandle,
    vertex: VertexHandle,
    fillet_index: usize,
) -> Result<(), ApexBuildError> {
    if vertex_is_on_face_outer_loop(brep, target_cylinder_face, vertex) {
        return Ok(());
    }
    let point = brep
        .vertices
        .get(vertex)
        .map(|vertex_data| vertex_data.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] endpoint vertex \
                 {vertex:?} is stale"
            ),
        })?;
    if let Some(coincident) = coincident_vertex_on_face_outer_loop(
        brep,
        target_cylinder_face,
        point,
        brep.tolerance.linear(),
    ) {
        if coincident != vertex {
            brep.weld_vertex_pair(vertex, coincident).map_err(|err| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] failed to \
                         adopt coincident cylinder-loop vertex {coincident:?} into {vertex:?}: {err:?}"
                    ),
                }
            })?;
        }
        return Ok(());
    }

    let edge = if let Some(edge) =
        containing_outer_loop_edge_by_curve(brep, target_cylinder_face, point)
    {
        edge
    } else if let Some(edge) =
        author_missing_on_cylinder_line_trim_for_endpoint(brep, target_cylinder_face, point)?
    {
        edge
    } else {
        let surface_status = if point_lies_on_cylinder_face(brep, target_cylinder_face, point) {
            "point lies on the target cylinder surface but no current trim edge owns it"
        } else {
            "point is not on the target cylinder surface"
        };
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] endpoint {vertex:?} \
                 at ({:.6},{:.6},{:.6}) does not lie on any curve-owned outer-loop edge of \
                 target cylinder face {target_cylinder_face:?}; {surface_status}",
                point.x, point.y, point.z
            ),
        });
    };
    let tolerance = brep.tolerance;
    let split = brep
        .split_edge_at_point(edge, point, &tolerance)
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed to split \
                 target cylinder edge {edge:?} at endpoint {vertex:?}: {err:?}"
            ),
        })?;
    brep.weld_vertex_pair(vertex, split.vertex)
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed to adopt \
                 split vertex {:?} into endpoint {vertex:?}: {err:?}",
                split.vertex
            ),
        })?;
    Ok(())
}
