//! Planned repeated-apex cap partitioning for constructive Gregory apexes.
//!
//! This module owns the non-serial graph rewrite used when a repeated-apex cap
//! loop cannot be normalized by any order of FilletEnd `mef`s. It builds a
//! typed partition from live FaceEdge half-edges, rewires the cap loop to the
//! planned FilletEnd/FaceEdge cycle, and validates topology plus semantic
//! cap-to-cylinder seam ownership before accepting.

use popkern_core::geometry::curve_inversion::invert_curve_at_point;
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, Vertex, VertexHandle};

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
    insert_face_with_outer_loop, insert_unwired_edge_pair, oriented_edge_half_edges_on_face,
    validate_partition_cap_to_cylinder_seams, wire_existing_loop,
};
use super::repeated_apex_trim_subcomplex::{
    TopologyBalanceRequirement, TopologyCountDelta, TrimEdgeConstruction, TrimOwnerPair,
    TrimRewriteEdgePlan, TrimRewriteSubcomplexPlan, topology_count_summary,
    trim_rewrite_support_groups,
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

#[derive(Clone)]
struct FaceEdgeOwnership {
    fillet_index: usize,
    edge: EdgeHandle,
    apex_to_t_b: HalfEdgeHandle,
    face_edge_cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    t_a_next_to_apex: HalfEdgeHandle,
}

impl FaceEdgeOwnership {
    fn from_retired_span(
        brep: &BRep,
        span: &RetiredApexSpan,
        cap_face: FaceHandle,
        fillet_index: usize,
        n_fillets: usize,
    ) -> Result<Self, ApexBuildError> {
        if span.half_edges.len() != 3 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition expected FaceEdge arc[{fillet_index}] to \
                     have 3 half-edges, got {}",
                    span.half_edges.len()
                ),
            });
        }
        let face_edge_cap_half_edge = span.half_edges[1];
        let face = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .and_then(|half_edge| {
                brep.loops
                    .get(half_edge.loop_)
                    .map(|loop_data| loop_data.face)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     middle half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        if face != cap_face {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] middle \
                     half-edge is on {face:?}, expected cap face {cap_face:?}"
                ),
            });
        }
        let edge = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     cap half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        let support_half_edge = brep
            .half_edges
            .get(face_edge_cap_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     cap half-edge {face_edge_cap_half_edge:?}"
                ),
            })?;
        let support_face = brep
            .half_edges
            .get(support_half_edge)
            .and_then(|half_edge| {
                brep.loops
                    .get(half_edge.loop_)
                    .map(|loop_data| loop_data.face)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] has stale \
                     support half-edge {support_half_edge:?}"
                ),
            })?;
        if support_face == cap_face {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc[{fillet_index}] support side \
                     is still on cap face {cap_face:?}"
                ),
            });
        }
        if fillet_index >= n_fillets {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge arc fillet index {fillet_index} \
                     out of range for {n_fillets} fillets"
                ),
            });
        }
        Ok(Self {
            fillet_index,
            edge,
            apex_to_t_b: span.half_edges[0],
            face_edge_cap_half_edge,
            support_half_edge,
            support_face,
            t_a_next_to_apex: span.half_edges[2],
        })
    }
}

struct FilletEndOwnership {
    fillet_index: usize,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    cylinder_candidate_half_edge: HalfEdgeHandle,
    target_cylinder_face: FaceHandle,
    retired_apex_cell: RetiredApexCell,
}

struct RetiredApexCell {
    previous_face_edge_index: usize,
    next_face_edge_index: usize,
    t_a_to_apex: HalfEdgeHandle,
    apex_to_t_b: HalfEdgeHandle,
}

struct PartitionBoundaryOutcome {
    cap_cycle: Vec<HalfEdgeHandle>,
    strip_side_overrides: std::collections::BTreeMap<usize, HalfEdgeHandle>,
    strip_loop_overrides: std::collections::BTreeMap<usize, Vec<HalfEdgeHandle>>,
}

impl PartitionBoundaryOutcome {
    fn from_current_cap_cycle(cap_cycle: Vec<HalfEdgeHandle>) -> Self {
        Self {
            cap_cycle,
            strip_side_overrides: std::collections::BTreeMap::new(),
            strip_loop_overrides: std::collections::BTreeMap::new(),
        }
    }

    fn strip_loop(
        &self,
        fillet_index: usize,
        ownership: &FilletEndOwnership,
    ) -> Vec<HalfEdgeHandle> {
        if let Some(loop_half_edges) = self.strip_loop_overrides.get(&fillet_index) {
            return loop_half_edges.clone();
        }
        let shared = self
            .strip_side_overrides
            .get(&fillet_index)
            .copied()
            .unwrap_or(ownership.cylinder_candidate_half_edge);
        vec![
            shared,
            ownership.retired_apex_cell.t_a_to_apex,
            ownership.retired_apex_cell.apex_to_t_b,
        ]
    }

    fn diagnostic_label(&self) -> String {
        let overrides = self
            .strip_side_overrides
            .iter()
            .map(|(fillet_index, half_edge)| format!("FilletEnd[{fillet_index}]={half_edge:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        let strip_loops = self
            .strip_loop_overrides
            .iter()
            .map(|(fillet_index, half_edges)| format!("FilletEnd[{fillet_index}]={half_edges:?}"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "selected cap_cycle={:?}; strip-side overrides=[{}]; strip-loop overrides=[{}]",
            self.cap_cycle, overrides, strip_loops
        )
    }

    fn matches_current_cap_cycle(&self, current_cap_cycle: &[HalfEdgeHandle]) -> bool {
        self.cap_cycle == current_cap_cycle
    }

    fn plan_strip_loop_overrides(
        &mut self,
        brep: &BRep,
        fillet_ends: &[FilletEndOwnership],
    ) -> Result<(), ApexBuildError> {
        for ownership in fillet_ends {
            let candidates = self.strip_loop_candidates(brep, ownership)?;
            let mut rejected = Vec::new();
            let mut selected = None;
            for candidate in candidates {
                match validate_half_edge_cycle(
                    brep,
                    &candidate,
                    "selected FilletEnd strip boundary",
                    ownership.target_cylinder_face,
                ) {
                    Ok(()) => {
                        selected = Some(candidate);
                        break;
                    }
                    Err(err) => rejected.push(format!("{candidate:?}: {err:?}")),
                }
            }
            let Some(selected) = selected else {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition target-cylinder trim rewrite cannot yet \
                         mutate: no contiguous strip loop candidate for FilletEnd[{}]; rejected {}",
                        ownership.fillet_index,
                        rejected.join("; ")
                    ),
                });
            };
            self.strip_loop_overrides
                .insert(ownership.fillet_index, selected);
        }
        Ok(())
    }

    fn strip_loop_candidates(
        &self,
        brep: &BRep,
        ownership: &FilletEndOwnership,
    ) -> Result<Vec<Vec<HalfEdgeHandle>>, ApexBuildError> {
        let shared = self
            .strip_side_overrides
            .get(&ownership.fillet_index)
            .copied()
            .unwrap_or(ownership.cylinder_candidate_half_edge);
        let ta_to_apex = ownership.retired_apex_cell.t_a_to_apex;
        let apex_to_tb = ownership.retired_apex_cell.apex_to_t_b;
        let apex_to_ta = brep
            .half_edges
            .get(ta_to_apex)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] strip candidate references \
                     stale T_a->apex half-edge {ta_to_apex:?}",
                    ownership.fillet_index
                ),
            })?;
        let tb_to_apex = brep
            .half_edges
            .get(apex_to_tb)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] strip candidate references \
                     stale apex->T_b half-edge {apex_to_tb:?}",
                    ownership.fillet_index
                ),
            })?;
        Ok(vec![
            vec![shared, ta_to_apex, apex_to_tb],
            vec![shared, tb_to_apex, apex_to_ta],
        ])
    }
}

fn validate_partition_ownership_records(
    brep: &BRep,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    for face_edge in face_edges {
        let cap_twin = brep
            .half_edges
            .get(face_edge.face_edge_cap_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] cap half-edge {:?} is stale",
                    face_edge.fillet_index, face_edge.face_edge_cap_half_edge
                ),
            })?;
        if cap_twin != face_edge.support_half_edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] support half-edge mismatch: \
                     recorded {:?}, actual {:?}",
                    face_edge.fillet_index, face_edge.support_half_edge, cap_twin
                ),
            });
        }
        let support_face = brep
            .half_edges
            .get(face_edge.support_half_edge)
            .and_then(|half_edge| {
                brep.loops
                    .get(half_edge.loop_)
                    .map(|loop_data| loop_data.face)
            })
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] support half-edge {:?} is stale",
                    face_edge.fillet_index, face_edge.support_half_edge
                ),
            })?;
        if support_face != face_edge.support_face {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] support face changed from \
                     {:?} to {:?}",
                    face_edge.fillet_index, face_edge.support_face, support_face
                ),
            });
        }
        let cap_edge = brep
            .half_edges
            .get(face_edge.face_edge_cap_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] cap half-edge {:?} is stale",
                    face_edge.fillet_index, face_edge.face_edge_cap_half_edge
                ),
            })?;
        if cap_edge != face_edge.edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FaceEdge[{}] edge changed from {:?} to {:?}",
                    face_edge.fillet_index, face_edge.edge, cap_edge
                ),
            });
        }
    }

    for fillet_end in fillet_ends {
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
        let cap_twin = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        if cap_twin != fillet_end.cylinder_candidate_half_edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cylinder candidate mismatch: \
                     recorded {:?}, actual {:?}",
                    fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge, cap_twin
                ),
            });
        }
        let cap_edge = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        if cap_edge != fillet_end.edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] edge changed from {:?} to {:?}",
                    fillet_end.fillet_index, fillet_end.edge, cap_edge
                ),
            });
        }
        if !brep
            .faces
            .get(fillet_end.target_cylinder_face)
            .and_then(|face| face.surface.as_ref())
            .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} \
                     is not a live cylinder face",
                    fillet_end.fillet_index, fillet_end.target_cylinder_face
                ),
            });
        }
        let cylinder_start = brep
            .half_edges
            .get(fillet_end.cylinder_candidate_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cylinder candidate {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge
                ),
            })?;
        let cylinder_end = half_edge_end_vertex(brep, fillet_end.cylinder_candidate_half_edge)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cylinder candidate {:?} \
                     has stale twin",
                    fillet_end.fillet_index, fillet_end.cylinder_candidate_half_edge
                ),
            })?;
        if cylinder_start != cap_end || cylinder_end != cap_start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cylinder candidate is not \
                     the reverse of the cap half-edge",
                    fillet_end.fillet_index
                ),
            });
        }
        let retired_a_start = brep
            .half_edges
            .get(fillet_end.retired_apex_cell.t_a_to_apex)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] retired T_a->apex half-edge \
                     {:?} is stale",
                    fillet_end.fillet_index, fillet_end.retired_apex_cell.t_a_to_apex
                ),
            })?;
        let retired_apex = half_edge_end_vertex(brep, fillet_end.retired_apex_cell.t_a_to_apex)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] retired T_a->apex half-edge \
                     {:?} has stale twin",
                    fillet_end.fillet_index, fillet_end.retired_apex_cell.t_a_to_apex
                ),
            })?;
        let retired_b_start = brep
            .half_edges
            .get(fillet_end.retired_apex_cell.apex_to_t_b)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] retired apex->T_b half-edge \
                     {:?} is stale",
                    fillet_end.fillet_index, fillet_end.retired_apex_cell.apex_to_t_b
                ),
            })?;
        let retired_b_end = half_edge_end_vertex(brep, fillet_end.retired_apex_cell.apex_to_t_b)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] retired apex->T_b half-edge \
                     {:?} has stale twin",
                    fillet_end.fillet_index, fillet_end.retired_apex_cell.apex_to_t_b
                ),
            })?;
        if retired_a_start != cap_start
            || retired_apex != retired_b_start
            || retired_b_end != cap_end
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] retired apex cell is not \
                     T_a -> apex -> T_b for cap endpoints {:?}->{:?}",
                    fillet_end.fillet_index, cap_start, cap_end
                ),
            });
        }
        if fillet_end.retired_apex_cell.previous_face_edge_index >= face_edges.len()
            || fillet_end.retired_apex_cell.next_face_edge_index >= face_edges.len()
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] has invalid retired apex \
                     cell indices {} and {} for {} FaceEdges",
                    fillet_end.fillet_index,
                    fillet_end.retired_apex_cell.previous_face_edge_index,
                    fillet_end.retired_apex_cell.next_face_edge_index,
                    face_edges.len()
                ),
            });
        }
    }
    Ok(())
}

fn validate_target_cylinder_loop_ownership_ready(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
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
        let Some(vertices) = super::constructive_topology::collect_outer_loop_vertices(
            brep,
            fillet_end.target_cylinder_face,
        ) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] target cylinder face {:?} \
                     has no walkable outer loop",
                    fillet_end.fillet_index, fillet_end.target_cylinder_face
                ),
            });
        };
        if !vertices.contains(&cap_start) || !vertices.contains(&cap_end) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot yet author shared \
                     cap/cylinder ownership: target cylinder face {:?} outer loop does not \
                     contain both cap endpoints {:?}->{:?}",
                    fillet_end.fillet_index, fillet_end.target_cylinder_face, cap_start, cap_end
                ),
            });
        }
    }
    Ok(())
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
) -> Result<(), ApexBuildError> {
    let mut all_support_owned_trim_sides = Vec::new();
    let mut pending_cylinder_loop_replacements = Vec::new();
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
        let mut cycle = Vec::with_capacity(topology_plan.loop_candidate.segments.len());
        let mut support_owned_trim_sides = Vec::new();
        for segment in &topology_plan.loop_candidate.segments {
            if segment.label.contains("shared FilletEnd") {
                let half_edge = existing_shared_fillet_end_half_edge_for_segment(
                    brep,
                    rebuild_plan,
                    segment.start,
                    segment.end,
                )?;
                cycle.push(half_edge);
            } else if segment.label.contains("far-side") {
                let trim_sides = materialized_far_side_trim_sides(
                    brep,
                    rebuild_plan.cylinder_face,
                    segment.start,
                    segment.end,
                    segment.label,
                    rebuild_plan.fillet_index,
                    &allowed_support_faces,
                )?;
                support_owned_trim_sides.extend(trim_sides.iter().cloned());
                cycle.extend(trim_sides.iter().map(|trim_side| trim_side.half_edge));
            } else {
                let support_face = support_face_for_generator_trim_side(
                    brep,
                    rebuild_plan.cylinder_face,
                    cap_face,
                    segment.start,
                    segment.end,
                    segment.label,
                    rebuild_plan.fillet_index,
                    &allowed_support_faces,
                )?;
                let (half_edge, actual_support_face, support_side_owned) =
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
}

fn target_cylinder_support_faces(
    brep: &BRep,
    cylinder_face: FaceHandle,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
) -> std::collections::BTreeSet<FaceHandle> {
    let semantic_support_faces = face_edges
        .iter()
        .map(|face_edge| face_edge.support_face)
        .collect::<std::collections::BTreeSet<_>>();
    let mut support_faces = std::collections::BTreeSet::new();
    let Some(face_data) = brep.faces.get(cylinder_face) else {
        return support_faces;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return support_faces;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return support_faces;
    }
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        if let Some(adjacent_face) = adjacent_face_across_half_edge(brep, cylinder_face, cursor)
            && adjacent_face != cap_face
            && face_is_planar_support(brep, adjacent_face)
        {
            support_faces.insert(adjacent_face);
        }
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    if support_faces.is_empty() {
        return semantic_support_faces;
    }
    support_faces
}

fn face_is_planar_support(brep: &BRep, face: FaceHandle) -> bool {
    matches!(
        brep.faces
            .get(face)
            .and_then(|face_data| face_data.surface.as_ref())
            .map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Plane { .. })
    )
}

fn materialize_closed_support_trim_cycles(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<(), ApexBuildError> {
    trace_topology_counts(brep, "support-loop-entry");
    let mut by_support_face: std::collections::BTreeMap<FaceHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    for trim_side in trim_sides {
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
    }
    for face_edge in face_edges {
        if let Some(half_edges) = by_support_face.get_mut(&face_edge.support_face) {
            if !half_edge_is_degenerate(brep, face_edge.support_half_edge) {
                half_edges.push(face_edge.support_half_edge);
            }
        }
    }

    for (support_face, mut half_edges) in by_support_face {
        trace_topology_counts(brep, "support-face-before");
        coalesce_degenerate_support_open_vertices(brep, &half_edges)?;
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

fn coalesce_degenerate_support_open_vertices(
    brep: &mut BRep,
    half_edges: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    for (start, end) in missing_support_split_spans(brep, half_edges) {
        if start == end || !vertices_are_coincident(brep, start, end) {
            continue;
        }
        brep.weld_vertex_pair(start, end)
            .map_err(|err| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot coalesce degenerate support open-chain \
                 vertices {start:?}/{end:?}: {err:?}"
                ),
            })?;
    }
    Ok(())
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
    brep: &BRep,
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
        validate_support_split_span_geometry(brep, support_face, *start, *end)?;
        let Some(path) = support_outer_loop_path_between_vertices(brep, support_face, *start, *end)
        else {
            return Ok(SupportSplitCyclePlan {
                cycles: Vec::new(),
                missing_spans,
            });
        };
        supplemented.extend(path);
    }
    let cycles = closed_half_edge_cycles_from_unwired_set(brep, &supplemented);
    let covered = cycles.iter().map(Vec::len).sum::<usize>();
    if covered == supplemented.len() {
        return Ok(SupportSplitCyclePlan {
            cycles,
            missing_spans: Vec::new(),
        });
    }
    Ok(SupportSplitCyclePlan {
        cycles: Vec::new(),
        missing_spans,
    })
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

fn support_outer_loop_path_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let face = brep.faces.get(support_face)?;
    let loop_data = brep.loops.get(face.outer_loop)?;
    let loop_start = loop_data.half_edge;
    let mut candidate_starts = Vec::new();
    let mut cursor = loop_start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.origin == start {
            candidate_starts.push(cursor);
        }
        cursor = half_edge.next;
        if cursor == loop_start {
            break;
        }
    }
    candidate_starts
        .into_iter()
        .filter_map(|path_start| {
            let mut path = Vec::new();
            let mut cursor = path_start;
            for _ in 0..=brep.half_edges.len() {
                let half_edge = brep.half_edges.get(cursor)?;
                if brep.loops.get(half_edge.loop_)?.face != support_face {
                    return None;
                }
                path.push(cursor);
                if half_edge_end_vertex(brep, cursor)? == end {
                    return Some(path);
                }
                cursor = half_edge.next;
                if cursor == path_start {
                    break;
                }
            }
            None
        })
        .min_by_key(Vec::len)
}

fn support_outer_loop_forward_path_between_vertices(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let face = brep.faces.get(support_face)?;
    let loop_data = brep.loops.get(face.outer_loop)?;
    let loop_start = loop_data.half_edge;
    let mut cursor = loop_start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.origin == start {
            let mut path = Vec::new();
            let mut path_cursor = cursor;
            for _ in 0..=brep.half_edges.len() {
                let path_half_edge = brep.half_edges.get(path_cursor)?;
                if brep.loops.get(path_half_edge.loop_)?.face != support_face {
                    return None;
                }
                path.push(path_cursor);
                if half_edge_end_vertex(brep, path_cursor)? == end {
                    return Some(path);
                }
                path_cursor = path_half_edge.next;
                if path_cursor == cursor {
                    break;
                }
            }
            return None;
        }
        cursor = half_edge.next;
        if cursor == loop_start {
            break;
        }
    }
    None
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
            let end = half_edge_end_vertex(brep, *half_edge).unwrap_or_default();
            format!("{half_edge:?}:{start:?}->{end:?}")
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn support_face_for_generator_trim_side(
    brep: &BRep,
    cylinder_face: FaceHandle,
    cap_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<FaceHandle, ApexBuildError> {
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face || face == cap_face {
                return None;
            }
            if !allowed_support_faces.contains(&face) {
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
            if start_distance > tol || end_distance > tol {
                return None;
            }
            let start_on_boundary = vertex_is_on_face_outer_loop(brep, face, start);
            let end_on_boundary = vertex_is_on_face_outer_loop(brep, face, end);
            (start_on_boundary || end_on_boundary).then_some(face)
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose a unique \
             support owner for generator trim {label} {start:?}->{end:?}: boundary support \
             candidates {:?} from allowed faces {:?}. Generator trim support ownership must come \
             from the real support face adjacent to the target cylinder, not from a coplanar \
             duplicate.",
            candidates, allowed_support_faces
        ),
    })
}

fn unique_far_side_support_face_for_trim_side(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<FaceHandle, ApexBuildError> {
    let start = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot find far-side support owner for stale \
                 trim half-edge {half_edge:?}"
            ),
        })?;
    let end =
        half_edge_end_vertex(brep, half_edge).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot find far-side support owner for trim \
             half-edge {half_edge:?}: stale twin"
            ),
        })?;
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            if !allowed_support_faces.contains(&face) {
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
            "planned repeated-apex partition cannot choose a unique far-side support face for \
             generated trim {half_edge:?} {start:?}->{end:?}: adjacent-support candidates {:?} \
             from allowed faces {:?}. Far-side support ownership must stay within the two \
             FaceEdge support owners adjacent to the FilletEnd.",
            candidates, allowed_support_faces
        ),
    })
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
            let Some(end) = half_edge_end_vertex(brep, cursor) else {
                return Vec::new();
            };
            let Some(next) = by_start.get(&end).and_then(|edges| edges.first()).copied() else {
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

fn materialized_generator_trim_side(
    brep: &mut BRep,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<(HalfEdgeHandle, FaceHandle, bool), ApexBuildError> {
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
        Ok((half_edge, actual_support_face)) => {
            return Ok((half_edge, actual_support_face, true));
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
    Ok((half_edge, support_face, false))
}

fn materialize_support_face_mef_trim(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<(HalfEdgeHandle, FaceHandle), ApexBuildError> {
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
        if std::env::var_os("POPKERN_A2_TRACE").is_some() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} by support-face mef on {actual_support_face:?}; cylinder \
                 side {:?}, support side {:?}",
                result.half_edge, result.half_edge_twin
            );
        }
        Ok((result.half_edge, actual_support_face))
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

fn exact_cylinder_trim_curve_between(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let start_point = brep
        .vertices
        .get(start)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: stale start vertex {start:?}"
            ),
        })?;
    let end_point = brep
        .vertices
        .get(end)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: stale end vertex {end:?}"
            ),
        })?;
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, fillet_index)?;
    if cylinder
        .validate_generator_segment(brep, fillet_index, start, end)
        .is_ok()
    {
        return CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                     generator trim {label}: degenerate endpoints {start:?}->{end:?}"
                ),
            });
    }
    exact_support_cylinder_intersection_curve_between(
        brep,
        cylinder_face,
        start,
        end,
        start_point,
        end_point,
        label,
        fillet_index,
    )
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
        let (half_edge, actual_support_face, support_side_owned) =
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
        }]);
    }

    let (edge, half_edge, _twin) = insert_unwired_edge_pair(brep, start, end);
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
                    }]);
                }
                Err(_) => {
                    brep.half_edges.remove(half_edge);
                    if let Some(twin) = brep.edges.get(edge).map(|edge_data| {
                        brep.half_edges
                            .get(edge_data.half_edge)
                            .map(|half_edge_data| half_edge_data.twin)
                            .unwrap_or_default()
                    }) && brep.half_edges.contains_key(twin)
                    {
                        brep.half_edges.remove(twin);
                    }
                    brep.edges.remove(edge);
                }
            }
        }
        Err(_) => {
            brep.half_edges.remove(half_edge);
            if let Some(twin) = brep.edges.get(edge).map(|edge_data| {
                brep.half_edges
                    .get(edge_data.half_edge)
                    .map(|half_edge_data| half_edge_data.twin)
                    .unwrap_or_default()
            }) && brep.half_edges.contains_key(twin)
            {
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
    for (start_support, end_support, joint_point, _) in candidates {
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
                    },
                    SupportOwnedTrimSide {
                        cylinder_face,
                        support_face: end_support,
                        fillet_index,
                        half_edge: half_edge_b,
                        support_side_owned: false,
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
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot split exact \
             far-side cylinder trim {start:?}->{end:?}: no plane-plane-cylinder transition \
             point from start supports {:?} to end supports {:?}",
            start_supports, end_supports
        ),
    })
}

#[allow(clippy::too_many_arguments)]
fn exact_support_cylinder_intersection_curve_between(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    start_point: Point3,
    end_point: Point3,
    label: &'static str,
    fillet_index: usize,
) -> Result<CurveOnEdge, ApexBuildError> {
    let cylinder_surface = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                 cylinder trim {label}: target cylinder face {cylinder_face:?} has no surface"
            ),
        })?;
    let mut blockers = Vec::new();
    for (support_face, face_data) in brep.faces.iter() {
        if support_face == cylinder_face {
            continue;
        }
        let Some(support_surface) = face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = support_surface else {
            continue;
        };
        let tol = brep.tolerance.linear().max(1.0e-9);
        let start_plane_distance = (start_point - *origin).dot(normal).abs();
        let end_plane_distance = (end_point - *origin).dot(normal).abs();
        if start_plane_distance > tol || end_plane_distance > tol {
            blockers.push(format!(
                "{support_face:?}: plane distances {:.3e}/{:.3e}",
                start_plane_distance, end_plane_distance
            ));
            continue;
        }
        if line_segment_lies_on_cylinder_and_plane(
            brep,
            cylinder_face,
            start_point,
            end_point,
            *origin,
            *normal,
        ) {
            return CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                         line support/cylinder trim for {label}: degenerate endpoints \
                         {start:?}->{end:?}"
                    ),
                });
        }
        let anchors = EdgeAnchors::Open {
            start: start_point,
            end: end_point,
        };
        if let Ok(triple) = compute_edge_geometry(
            cylinder_surface,
            support_surface,
            anchors,
            None,
            &brep.tolerance,
            &EdgeTripleOptions::default(),
        ) {
            let trim_start = triple.curve_3d.point_at(triple.t_start);
            let trim_end = triple.curve_3d.point_at(triple.t_end);
            if trim_start
                .distance_to(&start_point)
                .max(trim_end.distance_to(&end_point))
                <= brep.tolerance.linear().max(1.0e-9)
            {
                return Ok(CurveOnEdge {
                    geometry: triple.curve_3d,
                    t_start: triple.t_start,
                    t_end: triple.t_end,
                });
            }
            blockers.push(format!(
                "{support_face:?}: edge-triple endpoint mismatch {:.3e}/{:.3e}",
                trim_start.distance_to(&start_point),
                trim_end.distance_to(&end_point)
            ));
            continue;
        }
        let result = intersect_surfaces(
            cylinder_surface,
            (0.0, std::f64::consts::TAU),
            (-10.0, 10.0),
            support_surface,
            (-10.0, 10.0),
            (-10.0, 10.0),
            &brep.tolerance,
        )
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed analytic \
                 support/cylinder intersection for {label} on {support_face:?}: {err:?}"
            ),
        })?;
        for intersection in result.intersections {
            let SurfaceSurfaceIntersection::Curve(curve_on_surface) = intersection else {
                continue;
            };
            let range = curve_on_surface.curve_3d.parameter_range();
            let candidate = CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve_on_surface.curve_3d.clone()),
                t_start: range.0,
                t_end: range.1,
            };
            let Ok(t_start) = invert_curve_at_point(&candidate, start_point, &brep.tolerance)
            else {
                continue;
            };
            let Ok(t_end) = invert_curve_at_point(&candidate, end_point, &brep.tolerance) else {
                continue;
            };
            if (t_start - t_end).abs() <= brep.tolerance.linear().max(1.0e-9) {
                continue;
            }
            return Ok(CurveOnEdge {
                geometry: candidate.geometry,
                t_start,
                t_end,
            });
        }
        blockers.push(format!(
            "{support_face:?}: analytic intersection did not invert both endpoints {:?}->{:?}",
            start, end
        ));
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author exact \
             support/cylinder curve for {label} {start:?}->{end:?}; tried support planes [{}]",
            blockers.join("; ")
        ),
    })
}

fn support_plane_faces_containing_point(
    brep: &BRep,
    cylinder_face: FaceHandle,
    point: Point3,
) -> Vec<FaceHandle> {
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            ((point - *origin).dot(normal).abs() <= tol).then_some(face)
        })
        .collect()
}

fn line_segment_lies_on_cylinder_and_plane(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: Point3,
    end: Point3,
    plane_origin: Point3,
    plane_normal: Vec3,
) -> bool {
    let direction = end - start;
    if direction.length() <= brep.tolerance.linear().max(1.0e-12) {
        return false;
    }
    let tol = brep.tolerance.linear().max(1.0e-9);
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let point = start + direction * t;
        if (point - plane_origin).dot(&plane_normal).abs() > tol {
            return false;
        }
        if !point_lies_on_cylinder_face(brep, cylinder_face, point) {
            return false;
        }
    }
    true
}

fn plane_plane_cylinder_joint_points(
    brep: &BRep,
    cylinder_face: FaceHandle,
    plane_a: FaceHandle,
    plane_b: FaceHandle,
) -> Vec<Point3> {
    let Some(SurfaceGeometry::Cylinder {
        origin: cylinder_origin,
        axis,
        radius,
    }) = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_a_origin,
        normal: plane_a_normal,
        ..
    }) = brep
        .faces
        .get(plane_a)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_b_origin,
        normal: plane_b_normal,
        ..
    }) = brep
        .faces
        .get(plane_b)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let n1 = *plane_a_normal;
    let n2 = *plane_b_normal;
    let direction = n1.cross(&n2);
    let denom = direction.dot(&direction);
    if denom <= 1.0e-24 {
        return Vec::new();
    }
    let c1 = n1.dot(&(*plane_a_origin - Point3::new(0.0, 0.0, 0.0)));
    let c2 = n2.dot(&(*plane_b_origin - Point3::new(0.0, 0.0, 0.0)));
    let point_vec = (n2.cross(&direction) * c1 + direction.cross(&n1) * c2) * (1.0 / denom);
    let line_point = Point3::new(point_vec.x, point_vec.y, point_vec.z);
    line_cylinder_intersections(line_point, direction, *cylinder_origin, *axis, *radius)
        .into_iter()
        .filter(|point| point_lies_on_cylinder_face(brep, cylinder_face, *point))
        .filter(|point| {
            let tol = brep.tolerance.linear().max(1.0e-9);
            ((*point - *plane_a_origin).dot(plane_a_normal)).abs() <= tol
                && ((*point - *plane_b_origin).dot(plane_b_normal)).abs() <= tol
        })
        .collect()
}

fn line_cylinder_intersections(
    line_point: Point3,
    line_direction: Vec3,
    cylinder_origin: Point3,
    cylinder_axis: Vec3,
    radius: f64,
) -> Vec<Point3> {
    let Some(axis_unit) = cylinder_axis.normalized() else {
        return Vec::new();
    };
    let rel = line_point - cylinder_origin;
    let direction_perp = line_direction - axis_unit * line_direction.dot(&axis_unit);
    let rel_perp = rel - axis_unit * rel.dot(&axis_unit);
    let a = direction_perp.dot(&direction_perp);
    let b = 2.0 * rel_perp.dot(&direction_perp);
    let c = rel_perp.dot(&rel_perp) - radius * radius;
    if a.abs() <= 1.0e-24 {
        return Vec::new();
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -1.0e-12 {
        return Vec::new();
    }
    if discriminant.abs() <= 1.0e-12 {
        let t = -b / (2.0 * a);
        return vec![line_point + line_direction * t];
    }
    let root = discriminant.sqrt();
    let t0 = (-b - root) / (2.0 * a);
    let t1 = (-b + root) / (2.0 * a);
    vec![
        line_point + line_direction * t0,
        line_point + line_direction * t1,
    ]
}

fn vertices_are_materialized_generator_anchors(
    brep: &BRep,
    target_cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> bool {
    [start, end].into_iter().all(|vertex| {
        brep.vertices.get(vertex).is_some_and(|data| {
            !vertex_is_on_face_outer_loop(brep, target_cylinder_face, vertex)
                && point_lies_on_cylinder_face(brep, target_cylinder_face, data.position)
        })
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
        validate_cylinder_trim_rewrite_topology_preflight(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
        )?;
        apply_materialized_cylinder_loop_replacements(
            brep,
            cap_face,
            face_edges,
            &rebuild_plans,
            &topology_plans,
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
                "[a2-trace §14.9/cylinder-trim-rewrite] topology_plans=[{}]; {}; {}",
                topology_plans
                    .iter()
                    .map(CylinderTrimRewriteTopologyPlan::diagnostic_label)
                    .collect::<Vec<_>>()
                    .join("; "),
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
) -> Result<(), ApexBuildError> {
    let subcomplex_plan = plan_trim_subcomplex_for_cylinder_rewrite(
        brep,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
    )?;
    let before = topology_count_summary(brep);
    let mut candidate = brep.clone();
    match apply_materialized_cylinder_loop_replacements(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
    ) {
        Ok(()) => {
            let after = topology_count_summary(&candidate);
            if after.deficit_to_genus_zero() > before.deficit_to_genus_zero() {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition target-cylinder trim rewrite preflight \
                         would increase the unresolved closed-shell topology deficit before the \
                         cap rewrite can complete the subcomplex: before {}; after {}; delta {}. \
                         Required balance: {}. Planned trim subcomplex: {}. The rewrite must author an Euler-balanced cylinder/support/cap \
                         subcomplex, not expose a loop set that validation would need to repair.",
                        before.diagnostic_label(),
                        after.diagnostic_label(),
                        TopologyCountDelta::between(before, after).diagnostic_label(),
                        TopologyBalanceRequirement::from_deficit_increase(
                            after.deficit_to_genus_zero() - before.deficit_to_genus_zero()
                        )
                        .diagnostic_label(),
                        subcomplex_plan.diagnostic_label()
                    ),
                });
            }
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!(
                    "[a2-trace §14.9/topology-preflight] accepted before {}; after {}; delta {}; {}",
                    before.diagnostic_label(),
                    after.diagnostic_label(),
                    TopologyCountDelta::between(before, after).diagnostic_label(),
                    subcomplex_plan.diagnostic_label()
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
                     {}; delta {}. Required balance: {}. Planned trim subcomplex: {}. This is a construction blocker, not a validation blocker: the \
                     implementation must plan the Euler-balanced replacement subcomplex before \
                     changing the live BRep.",
                    before.diagnostic_label(),
                    after.diagnostic_label(),
                    TopologyCountDelta::between(before, after).diagnostic_label(),
                    TopologyBalanceRequirement::from_deficit_increase(
                        after.deficit_to_genus_zero() - before.deficit_to_genus_zero()
                    )
                    .diagnostic_label(),
                    subcomplex_plan.diagnostic_label()
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
    let mut edge_plans = Vec::new();
    let mut expected_edge_only_delta = TopologyCountDelta::zero();
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
                edge_plans.push(TrimRewriteEdgePlan {
                    fillet_index: topology_plan.fillet_index,
                    owner_pair: TrimOwnerPair::CylinderCap {
                        cylinder_face: topology_plan.cylinder_face,
                        cap_face,
                    },
                    construction: TrimEdgeConstruction::ExistingSeam,
                    start: segment.start,
                    end: segment.end,
                });
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
                expected_edge_only_delta.add_assign(TopologyCountDelta::between(before, after));
                for trim_side in trim_sides {
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
                    edge_plans.push(TrimRewriteEdgePlan {
                        fillet_index: topology_plan.fillet_index,
                        owner_pair: TrimOwnerPair::CylinderSupport {
                            cylinder_face: trim_side.cylinder_face,
                            support_face: trim_side.support_face,
                        },
                        construction: TrimEdgeConstruction::NewExactSupportCylinderTrim,
                        start,
                        end,
                    });
                }
            } else {
                let support_face = support_face_for_generator_trim_side(
                    brep,
                    rebuild_plan.cylinder_face,
                    cap_face,
                    segment.start,
                    segment.end,
                    segment.label,
                    rebuild_plan.fillet_index,
                    &allowed_support_faces,
                )?;
                edge_plans.push(TrimRewriteEdgePlan {
                    fillet_index: topology_plan.fillet_index,
                    owner_pair: TrimOwnerPair::CylinderSupport {
                        cylinder_face: topology_plan.cylinder_face,
                        support_face,
                    },
                    construction: TrimEdgeConstruction::NewExactGeneratorTrim,
                    start: segment.start,
                    end: segment.end,
                });
                expected_edge_only_delta
                    .add_assign(TopologyCountDelta::from_components(0, 1, 0, 0, 0));
            }
        }
    }
    Ok(TrimRewriteSubcomplexPlan {
        support_groups: trim_rewrite_support_groups(&edge_plans),
        edge_plans,
        expected_edge_only_delta,
    })
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

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    let loop_data = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot collect stale loop {loop_handle:?}"
            ),
        })?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition loop {loop_handle:?} starts at stale half-edge \
                 {start:?}"
            ),
        });
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        let next = brep
            .half_edges
            .get(cursor)
            .map(|half_edge| half_edge.next)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition loop {loop_handle:?} references stale \
                     half-edge {cursor:?}"
                ),
            })?;
        cursor = next;
        if cursor == start {
            return Ok(cycle);
        }
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition loop {loop_handle:?} did not close within {} \
             half-edge steps",
            brep.half_edges.len()
        ),
    })
}

fn intended_cap_boundary_cycle(
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
) -> Vec<HalfEdgeHandle> {
    let mut cycle = Vec::with_capacity(2 * fillet_ends.len());
    for (fillet_end, face_edge) in fillet_ends.iter().zip(face_edges) {
        cycle.push(fillet_end.cap_half_edge);
        cycle.push(face_edge.face_edge_cap_half_edge);
    }
    cycle
}

fn reversed_twin_cycle(
    brep: &BRep,
    current_cap_cycle: &[HalfEdgeHandle],
    cap_face: FaceHandle,
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    current_cap_cycle
        .iter()
        .rev()
        .map(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .map(|half_edge_data| half_edge_data.twin)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cap-boundary rewrite references stale \
                         half-edge {half_edge:?} on cap face {cap_face:?}"
                    ),
                })
        })
        .collect()
}

fn validate_half_edge_cycle(
    brep: &BRep,
    cycle: &[HalfEdgeHandle],
    label: &'static str,
    cap_face: FaceHandle,
) -> Result<(), ApexBuildError> {
    if cycle.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition {label} for cap face {cap_face:?} is empty"
            ),
        });
    }
    for i in 0..cycle.len() {
        let current = cycle[i];
        let next = cycle[(i + 1) % cycle.len()];
        let current_end =
            half_edge_end_vertex(brep, current).ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition {label} for cap face {cap_face:?} has stale \
                     current half-edge {current:?}"
                ),
            })?;
        let next_start = brep
            .half_edges
            .get(next)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition {label} for cap face {cap_face:?} has stale \
                     next half-edge {next:?}"
                ),
            })?;
        if current_end != next_start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition {label} for cap face {cap_face:?} is not \
                     contiguous: {current:?} ends at {current_end:?}, but {next:?} starts at \
                     {next_start:?}"
                ),
            });
        }
    }
    Ok(())
}

fn assert_cap_boundary_matches_selected_cycle(
    brep: &BRep,
    cap_face: FaceHandle,
    selected_cap_cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let cap_loop = brep
        .faces
        .get(cap_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition selected cap boundary references stale cap face \
                 {cap_face:?}"
            ),
        })?
        .outer_loop;
    let actual = collect_loop_half_edges(brep, cap_loop)?;
    if actual == selected_cap_cycle {
        return Ok(());
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition selected cap boundary mismatch: got {:?}, expected \
             selected cycle {:?}",
            actual, selected_cap_cycle
        ),
    })
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
    ReuseCurrentTwinPair,
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
                    if !brep
                        .edges
                        .get(*edge)
                        .is_some_and(|edge_data| edge_data.curve.is_some())
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
            FilletEndSeamRewrite::ReuseCurrentTwinPair
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
            FilletEndSeamRewrite::ReuseCurrentTwinPair => {
                "reuse current FilletEnd cap/cylinder twin pair".to_string()
            }
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
            let FilletEndSeamRewrite::ReauthorCapBoundaryAndCylinderSide { cap_half_edge, .. } =
                plan.seam_rewrite
            else {
                continue;
            };
            strip_side_overrides.insert(plan.fillet_index, cap_half_edge);
        }
        PartitionBoundaryOutcome {
            cap_cycle: self.reauthored_cap_cycle.clone(),
            strip_side_overrides,
            strip_loop_overrides: std::collections::BTreeMap::new(),
        }
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

struct CylinderTrimSurface {
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
}

impl CylinderTrimSurface {
    fn from_face(
        brep: &BRep,
        face: FaceHandle,
        fillet_index: usize,
    ) -> Result<Self, ApexBuildError> {
        let Some(surface) = brep
            .faces
            .get(face)
            .and_then(|face_data| face_data.surface.as_ref())
        else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] target cylinder \
                     face {face:?} has no surface"
                ),
            });
        };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: _,
        } = surface.geometry
        else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] target face \
                     {face:?} is not a cylinder"
                ),
            });
        };
        Ok(Self { origin, axis })
    }

    fn validate_generator_segment(
        &self,
        brep: &BRep,
        fillet_index: usize,
        start: VertexHandle,
        end: VertexHandle,
    ) -> Result<(), ApexBuildError> {
        let start_point = brep.vertices[start].position;
        let end_point = brep.vertices[end].position;
        let start_radial =
            cylinder_radial_vector(start_point, self.origin, self.axis).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot classify \
                     generator trim {start:?}->{end:?}: degenerate cylinder axis"
                    ),
                }
            })?;
        let end_radial =
            cylinder_radial_vector(end_point, self.origin, self.axis).ok_or_else(|| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot classify \
                     generator trim {start:?}->{end:?}: degenerate cylinder axis"
                    ),
                }
            })?;
        let Some(start_unit) = start_radial.normalized() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} has degenerate start radial vector"
                ),
            });
        };
        let Some(end_unit) = end_radial.normalized() else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} has degenerate end radial vector"
                ),
            });
        };
        if start_unit.dot(&end_unit) < 1.0 - 1.0e-5 {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] generator trim \
                     {start:?}->{end:?} does not keep a constant cylinder radial direction"
                ),
            });
        }
        Ok(())
    }
}

fn retired_trim_runs_by_adjacent_face(
    brep: &BRep,
    cylinder_face: FaceHandle,
    path: &[HalfEdgeHandle],
) -> Vec<RetiredCylinderTrimRun> {
    let mut runs = Vec::new();
    let mut current_face = None;
    let mut current_path = Vec::new();
    for half_edge in path {
        let Some(adjacent_face) = adjacent_face_across_half_edge(brep, cylinder_face, *half_edge)
        else {
            continue;
        };
        if current_face.is_some_and(|face| face != adjacent_face) && !current_path.is_empty() {
            runs.push(RetiredCylinderTrimRun::from_path(
                brep,
                cylinder_face,
                current_face.expect("current face exists"),
                std::mem::take(&mut current_path),
            ));
        }
        current_face = Some(adjacent_face);
        current_path.push(*half_edge);
    }
    if let Some(adjacent_face) = current_face
        && !current_path.is_empty()
    {
        runs.push(RetiredCylinderTrimRun::from_path(
            brep,
            cylinder_face,
            adjacent_face,
            current_path,
        ));
    }
    runs
}

impl RetiredCylinderTrimRun {
    fn from_path(
        brep: &BRep,
        cylinder_face: FaceHandle,
        adjacent_face: FaceHandle,
        path: Vec<HalfEdgeHandle>,
    ) -> Self {
        let start = path
            .first()
            .and_then(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .unwrap_or_default();
        let end = path
            .last()
            .and_then(|half_edge| half_edge_end_vertex(brep, *half_edge))
            .unwrap_or_default();
        let vertices = path
            .iter()
            .filter_map(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .collect();
        let off_surface_vertices = path_off_cylinder_vertices(brep, cylinder_face, &path);
        Self {
            adjacent_face,
            path,
            start,
            end,
            vertices,
            off_surface_vertices,
        }
    }
}

fn adjacent_face_across_half_edge(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    let twin_loop = brep.half_edges.get(twin)?.loop_;
    let twin_face = brep.loops.get(twin_loop)?.face;
    (twin_face != cylinder_face).then_some(twin_face)
}

fn face_surface_kind(brep: &BRep, face: FaceHandle) -> &'static str {
    match brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "Plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "Cone",
        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
        Some(SurfaceGeometry::Torus { .. }) => "Torus",
        Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "Gregory",
        None => "None",
    }
}

fn replacement_run_exact_curve_kind(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> Option<&'static str> {
    let support_surface = brep
        .faces
        .get(run.adjacent_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?;
    let SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    } = support_surface
    else {
        return None;
    };
    let start_point = brep.vertices.get(run.start)?.position;
    let end_point = brep.vertices.get(run.end)?.position;
    let start_plane_distance = (start_point - *plane_origin).dot(plane_normal).abs();
    let end_plane_distance = (end_point - *plane_origin).dot(plane_normal).abs();
    let tol = brep.tolerance.linear().max(1.0e-9);
    if start_plane_distance > tol || end_plane_distance > tol {
        return None;
    }
    if !point_lies_on_cylinder_face(brep, cylinder_face, start_point)
        || !point_lies_on_cylinder_face(brep, cylinder_face, end_point)
    {
        return None;
    }
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, 0).ok()?;
    if cylinder
        .validate_generator_segment(brep, 0, run.start, run.end)
        .is_ok()
    {
        return Some("Line");
    }
    let SurfaceGeometry::Cylinder { axis, .. } = brep
        .faces
        .get(cylinder_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?
    else {
        return None;
    };
    if plane_normal.dot(axis).abs() > 1.0e-9 {
        return Some("Ellipse");
    }
    None
}

fn fragment_replacement_authoring_status(
    brep: &BRep,
    cylinder_face: FaceHandle,
    fragment: &RetiredCylinderTrimFragment,
) -> FragmentAuthoringStatus {
    let mut labels = Vec::new();
    let mut authorable = true;
    for run in &fragment.adjacent_face_runs {
        if !run.requires_replacement() {
            labels.push(format!(
                "{:?} {:?}->{:?}: reuse",
                run.adjacent_face, run.start, run.end
            ));
            continue;
        }
        if let Some(kind) = replacement_run_exact_curve_kind(brep, cylinder_face, run) {
            labels.push(format!(
                "{:?} {:?}->{:?}: exact {kind}",
                run.adjacent_face, run.start, run.end
            ));
            continue;
        }
        if let Some(plan) = replacement_run_reownership_plan(brep, cylinder_face, run) {
            labels.push(format!(
                "{:?} {:?}->{:?}: split/re-own {}",
                run.adjacent_face,
                run.start,
                run.end,
                plan.diagnostic_label()
            ));
            continue;
        }
        authorable = false;
        labels.push(format!(
            "{:?} {:?}->{:?}: blocked ({})",
            run.adjacent_face,
            run.start,
            run.end,
            replacement_run_curve_blocker(brep, cylinder_face, run)
        ));
    }
    FragmentAuthoringStatus {
        authorable,
        diagnostic: labels.join("; "),
    }
}

fn replacement_run_reownership_plan(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> Option<ReplacementRunReownershipPlan> {
    if run.path.len() <= 1 {
        return None;
    }
    let mut spans = Vec::new();
    for half_edge in &run.path {
        let start = brep.half_edges.get(*half_edge)?.origin;
        let end = half_edge_end_vertex(brep, *half_edge)?;
        let candidate_faces = replacement_span_support_candidates(brep, cylinder_face, start, end);
        if candidate_faces.is_empty() {
            return None;
        }
        spans.push(ReplacementRunSpanOwnership {
            start,
            end,
            candidate_faces,
        });
    }
    Some(ReplacementRunReownershipPlan { spans })
}

fn replacement_run_reownership_blocker(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> String {
    if run.path.len() <= 1 {
        return "run has one span, so split/re-own decomposition does not apply".to_string();
    }
    let span_reports = run
        .path
        .iter()
        .filter_map(|half_edge| {
            let start = brep.half_edges.get(*half_edge)?.origin;
            let end = half_edge_end_vertex(brep, *half_edge)?;
            let candidates = replacement_span_support_candidates(brep, cylinder_face, start, end);
            let candidate_label = if candidates.is_empty() {
                "no exact support-plane candidates".to_string()
            } else {
                candidates
                    .iter()
                    .map(|candidate| format!("{:?}/{}", candidate.face, candidate.curve_kind))
                    .collect::<Vec<_>>()
                    .join("|")
            };
            Some(format!("{start:?}->{end:?}: {candidate_label}"))
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!("reownership spans [{span_reports}]")
}

fn replacement_span_support_candidates(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Vec<ReplacementRunSupportCandidate> {
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            let SurfaceGeometry::Plane { .. } = face_data.surface.as_ref()?.geometry else {
                return None;
            };
            let curve_kind =
                replacement_span_exact_curve_kind(brep, cylinder_face, face, start, end)?;
            Some(ReplacementRunSupportCandidate { face, curve_kind })
        })
        .collect()
}

fn replacement_span_exact_curve_kind(
    brep: &BRep,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<&'static str> {
    let support_surface = brep
        .faces
        .get(support_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?;
    let SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    } = support_surface
    else {
        return None;
    };
    let start_point = brep.vertices.get(start)?.position;
    let end_point = brep.vertices.get(end)?.position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    if (start_point - *plane_origin).dot(plane_normal).abs() > tol
        || (end_point - *plane_origin).dot(plane_normal).abs() > tol
    {
        return None;
    }
    if !point_lies_on_cylinder_face(brep, cylinder_face, start_point)
        || !point_lies_on_cylinder_face(brep, cylinder_face, end_point)
    {
        return None;
    }
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, 0).ok()?;
    if cylinder
        .validate_generator_segment(brep, 0, start, end)
        .is_ok()
    {
        return Some("Line");
    }
    let SurfaceGeometry::Cylinder { axis, .. } = brep
        .faces
        .get(cylinder_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?
    else {
        return None;
    };
    if plane_normal.dot(axis).abs() > 1.0e-9 {
        return Some("Ellipse");
    }
    None
}

fn replacement_run_curve_blocker(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> String {
    let Some(start_point) = brep.vertices.get(run.start).map(|vertex| vertex.position) else {
        return format!("start vertex {:?} is stale", run.start);
    };
    let Some(end_point) = brep.vertices.get(run.end).map(|vertex| vertex.position) else {
        return format!("end vertex {:?} is stale", run.end);
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    }) = brep
        .faces
        .get(run.adjacent_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return "support surface is not a plane".to_string();
    };
    let Some(SurfaceGeometry::Cylinder {
        origin: cylinder_origin,
        axis,
        ..
    }) = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return "target face is not a cylinder".to_string();
    };
    let plane_start = (start_point - *plane_origin).dot(plane_normal).abs();
    let plane_end = (end_point - *plane_origin).dot(plane_normal).abs();
    let plane_axis_dot = plane_normal.dot(axis).abs();
    let radial_alignment = match (
        cylinder_radial_vector(start_point, *cylinder_origin, *axis)
            .and_then(|radial| radial.normalized()),
        cylinder_radial_vector(end_point, *cylinder_origin, *axis)
            .and_then(|radial| radial.normalized()),
    ) {
        (Some(start_radial), Some(end_radial)) => start_radial.dot(&end_radial),
        _ => f64::NAN,
    };
    format!(
        "plane distances start={plane_start:.3e}, end={plane_end:.3e}, \
         |plane_normal.cylinder_axis|={plane_axis_dot:.3e}, radial_alignment={radial_alignment:.6}"
    )
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

fn generator_anchor_on_cylinder_loop(
    brep: &BRep,
    face: FaceHandle,
    endpoint: Point3,
) -> Option<VertexHandle> {
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return None;
    };
    let endpoint_radial = cylinder_radial_vector(endpoint, origin, axis)?;
    if (endpoint_radial.length() - radius).abs()
        > brep.tolerance.linear().max(radius.abs() * 1.0e-6)
    {
        return None;
    }
    let endpoint_unit = endpoint_radial.normalized()?;
    let mut vertices = super::constructive_topology::collect_outer_loop_vertices(brep, face)?;
    let extra_vertices = brep
        .vertices
        .iter()
        .filter_map(|(vertex, data)| {
            (!vertices.contains(&vertex)
                && (data.outgoing == HalfEdgeHandle::default()
                    || point_lies_on_cylinder_face(brep, face, data.position)))
            .then_some(vertex)
        })
        .collect::<Vec<_>>();
    vertices.extend(extra_vertices);
    vertices
        .into_iter()
        .filter_map(|vertex| {
            let point = brep.vertices.get(vertex)?.position;
            if point.distance_to(&endpoint) <= brep.tolerance.linear() {
                return None;
            }
            if !point_lies_on_cylinder_face(brep, face, point) {
                return None;
            }
            let radial = cylinder_radial_vector(point, origin, axis)?;
            let radial_unit = radial.normalized()?;
            if radial_unit.dot(&endpoint_unit) < 1.0 - 1.0e-5 {
                return None;
            }
            Some((vertex, point.distance_to(&endpoint)))
        })
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(vertex, _)| vertex)
}

fn materialize_support_boundary_generator_anchor(
    brep: &mut BRep,
    target_cylinder_face: FaceHandle,
    endpoint: VertexHandle,
    anchor: &SupportTopologyAnchorCandidate,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
    fillet_index: usize,
) -> Result<VertexHandle, ApexBuildError> {
    if let Some(existing) = coincident_vertex_on_face_outer_loop(
        brep,
        anchor.support_face,
        anchor.point,
        brep.tolerance.linear(),
    ) {
        if !point_lies_on_cylinder_face(
            brep,
            target_cylinder_face,
            brep.vertices[existing].position,
        ) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot reuse \
                     support boundary anchor {existing:?} for endpoint {endpoint:?}: it is not \
                     on target cylinder face {target_cylinder_face:?}"
                ),
            });
        }
        return Ok(existing);
    }
    let mut support_edges = Vec::new();
    if let Some(edge) = containing_outer_loop_edge_by_curve(brep, anchor.support_face, anchor.point)
    {
        support_edges.push((anchor.support_face, edge));
    }
    for support_face in allowed_support_faces {
        if *support_face == anchor.support_face {
            continue;
        }
        let Some(face_data) = brep.faces.get(*support_face) else {
            continue;
        };
        let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
            face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        if (anchor.point - *origin).dot(normal).abs() > brep.tolerance.linear().max(1.0e-9) {
            continue;
        }
        if let Some(edge) = containing_outer_loop_edge_by_curve(brep, *support_face, anchor.point) {
            support_edges.push((*support_face, edge));
        }
    }
    support_edges.sort_by_key(|(support_face, edge)| (*support_face, *edge));
    support_edges.dedup_by_key(|(_, edge)| *edge);
    if support_edges.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot materialize \
                 support-aware generator anchor for endpoint {endpoint:?}: computed point \
                 ({:.6},{:.6},{:.6}) lies on support face {:?} but not on a splittable \
                 curve-owned support boundary edge",
                anchor.point.x, anchor.point.y, anchor.point.z, anchor.support_face
            ),
        });
    }
    if std::env::var_os("POPKERN_A2_TRACE").is_some() {
        eprintln!(
            "[a2-trace §14.9/cylinder-anchor] support boundary promotion for endpoint \
             {endpoint:?} on {:?}: {:?}",
            target_cylinder_face, support_edges
        );
    }
    let tolerance = brep.tolerance;
    let mut primary_vertex = None;
    for (support_face, edge) in support_edges {
        let split = brep
            .split_edge_at_point(edge, anchor.point, &tolerance)
            .map_err(|err| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] failed to split \
                     support boundary edge {edge:?} on {support_face:?} for endpoint \
                     {endpoint:?}: {err:?}"
                ),
            })?;
        if let Some(primary) = primary_vertex {
            brep.weld_vertex_pair(primary, split.vertex)
                .map_err(|err| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] failed to \
                         weld support-promoted anchor {:?} into {:?}: {err:?}",
                        split.vertex, primary
                    ),
                })?;
        } else {
            primary_vertex = Some(split.vertex);
        }
    }
    let split_vertex = primary_vertex.ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] did not produce a support \
             boundary anchor for endpoint {endpoint:?}"
        ),
    })?;
    if !point_lies_on_cylinder_face(
        brep,
        target_cylinder_face,
        brep.vertices[split_vertex].position,
    ) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] split support \
                 boundary at {split_vertex:?}, but the resulting anchor is not on \
                 target cylinder face {target_cylinder_face:?}",
                split_vertex = split_vertex
            ),
        });
    }
    Ok(split_vertex)
}

#[derive(Clone)]
struct SupportTopologyAnchorCandidate {
    point: Point3,
    support_face: FaceHandle,
    endpoint_distance: f64,
}

fn support_topology_generator_anchor_pair(
    brep: &BRep,
    face: FaceHandle,
    start_endpoint: Point3,
    end_endpoint: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Option<(
    SupportTopologyAnchorCandidate,
    SupportTopologyAnchorCandidate,
)> {
    let start_candidates =
        support_topology_anchor_candidates(brep, face, start_endpoint, allowed_support_faces);
    let end_candidates =
        support_topology_anchor_candidates(brep, face, end_endpoint, allowed_support_faces);
    let mut pairs = Vec::new();
    for start in &start_candidates {
        for end in &end_candidates {
            let direct_support = start.support_face == end.support_face;
            let transition_support = !plane_plane_cylinder_joint_points(
                brep,
                face,
                start.support_face,
                end.support_face,
            )
            .is_empty();
            if !direct_support && !transition_support {
                continue;
            }
            let score = start.endpoint_distance
                + end.endpoint_distance
                + if direct_support { 0.0 } else { 1.0e-6 };
            pairs.push((score, start.clone(), end.clone()));
        }
    }
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    pairs.into_iter().next().map(|(_, start, end)| (start, end))
}

fn support_topology_anchor_candidates(
    brep: &BRep,
    face: FaceHandle,
    endpoint: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<SupportTopologyAnchorCandidate> {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return Vec::new();
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return Vec::new();
    };
    let Some(radial) = cylinder_radial_vector(endpoint, origin, axis) else {
        return Vec::new();
    };
    if (radial.length() - radius).abs() > brep.tolerance.linear().max(radius.abs() * 1.0e-6) {
        return Vec::new();
    }
    let Some(radial_unit) = radial.normalized() else {
        return Vec::new();
    };
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return Vec::new();
    }
    let tol = brep.tolerance.linear().max(1.0e-9);
    let line_origin = origin + radial_unit * radius;
    let mut candidates = Vec::new();
    for (support_face, face_data) in brep.faces.iter() {
        if support_face == face {
            continue;
        }
        if !allowed_support_faces.contains(&support_face) {
            continue;
        }
        let Some(SurfaceGeometry::Plane {
            origin: plane_origin,
            normal,
            ..
        }) = face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let denom = axis.dot(normal);
        if denom.abs() <= 1.0e-12 {
            continue;
        }
        let lambda = (*plane_origin - line_origin).dot(normal) / denom;
        let point = line_origin + axis * lambda;
        if point.distance_to(&endpoint) <= tol {
            continue;
        }
        if !point_lies_on_cylinder_face(brep, face, point) {
            continue;
        }
        if (point - *plane_origin).dot(normal).abs() > tol {
            continue;
        }
        if candidates
            .iter()
            .any(|candidate: &SupportTopologyAnchorCandidate| {
                candidate.point.distance_to(&point) <= tol && candidate.support_face == support_face
            })
        {
            continue;
        }
        candidates.push(SupportTopologyAnchorCandidate {
            point,
            support_face,
            endpoint_distance: point.distance_to(&endpoint),
        });
    }
    candidates.sort_by(|a, b| {
        a.endpoint_distance
            .partial_cmp(&b.endpoint_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    candidates
}

fn generator_anchor_blocker(brep: &BRep, face: FaceHandle, endpoint: Point3) -> String {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return "target face has no surface".to_string();
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return "target face is not a cylinder".to_string();
    };
    let Some(endpoint_radial) = cylinder_radial_vector(endpoint, origin, axis) else {
        return "endpoint has degenerate radial vector to cylinder axis".to_string();
    };
    let endpoint_surface_error = (endpoint_radial.length() - radius).abs();
    let Some(endpoint_unit) = endpoint_radial.normalized() else {
        return "endpoint radial vector cannot be normalized".to_string();
    };
    let Some(vertices) = super::constructive_topology::collect_outer_loop_vertices(brep, face)
    else {
        return "target cylinder face has no walkable outer loop".to_string();
    };
    let mut reports = vertices
        .into_iter()
        .filter_map(|vertex| {
            let point = brep.vertices.get(vertex)?.position;
            let radial = cylinder_radial_vector(point, origin, axis)?;
            let radial_distance = radial.length();
            let radial_unit = radial.normalized()?;
            let surface_error = (radial_distance - radius).abs();
            let radial_alignment = radial_unit.dot(&endpoint_unit);
            Some((vertex, surface_error, radial_alignment, point))
        })
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| {
        let a_score = a.1 + (1.0 - a.2).abs() * radius.abs().max(1.0e-9);
        let b_score = b.1 + (1.0 - b.2).abs() * radius.abs().max(1.0e-9);
        a_score
            .partial_cmp(&b_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let nearest = reports
        .iter()
        .take(4)
        .map(|(vertex, surface_error, radial_alignment, point)| {
            format!(
                "{vertex:?}@({:.6},{:.6},{:.6}) surface_error={surface_error:.3e}, \
                 radial_alignment={radial_alignment:.6}",
                point.x, point.y, point.z
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!(
        "endpoint cylinder surface_error={endpoint_surface_error:.3e}; nearest radial loop \
         candidates [{}]. A commercial rebuild must materialize anchor vertices on the stored \
         cylinder surface before wiring generator trims.",
        nearest
    )
}

fn outer_loop_path_between_vertices(
    brep: &BRep,
    face: FaceHandle,
    start_vertex: VertexHandle,
    end_vertex: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let mut cursor = loop_data.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.origin == start_vertex {
            let mut path = Vec::new();
            let mut path_cursor = cursor;
            for _ in 0..=brep.half_edges.len() {
                let path_half_edge = brep.half_edges.get(path_cursor)?;
                path.push(path_cursor);
                let end = half_edge_end_vertex(brep, path_cursor)?;
                if end == end_vertex {
                    return Some(path);
                }
                path_cursor = path_half_edge.next;
                if path_cursor == cursor {
                    return None;
                }
            }
            return None;
        }
        cursor = half_edge.next;
        if cursor == loop_data.half_edge {
            break;
        }
    }
    None
}

fn path_off_cylinder_vertices(
    brep: &BRep,
    face: FaceHandle,
    path: &[HalfEdgeHandle],
) -> Vec<VertexHandle> {
    path.iter()
        .filter_map(|half_edge| {
            let vertex = brep.half_edges.get(*half_edge)?.origin;
            let point = brep.vertices.get(vertex)?.position;
            (!point_lies_on_cylinder_face(brep, face, point)).then_some(vertex)
        })
        .collect()
}

fn vertex_is_on_face_outer_loop(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    super::constructive_topology::collect_outer_loop_vertices(brep, face)
        .is_some_and(|vertices| vertices.contains(&vertex))
}

fn coincident_vertex_on_face_outer_loop(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: f64,
) -> Option<VertexHandle> {
    let vertices = super::constructive_topology::collect_outer_loop_vertices(brep, face)?;
    let tol2 = tolerance * tolerance;
    vertices.into_iter().find(|vertex| {
        brep.vertices.get(*vertex).is_some_and(|vertex_data| {
            let delta = vertex_data.position - point;
            delta.dot(&delta) <= tol2
        })
    })
}

fn containing_outer_loop_edge_by_curve(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
) -> Option<EdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        let edge = brep.edges.get(half_edge.edge)?;
        if let Some(curve) = edge.curve.as_ref()
            && invert_curve_at_point(curve, point, &brep.tolerance).is_ok()
            && point_is_not_edge_endpoint(brep, half_edge.edge, point)
        {
            return Some(half_edge.edge);
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    None
}

fn author_missing_on_cylinder_line_trim_for_endpoint(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
) -> Result<Option<EdgeHandle>, ApexBuildError> {
    let Some(half_edge) = containing_on_cylinder_line_chord_half_edge(brep, face, point) else {
        return Ok(None);
    };
    let edge = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge_data| half_edge_data.edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} is stale"
            ),
        })?;
    let start = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge_data| brep.vertices.get(half_edge_data.origin))
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has stale start vertex"
            ),
        })?;
    let end_vertex =
        half_edge_end_vertex(brep, half_edge).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has stale twin"
            ),
        })?;
    let end = brep
        .vertices
        .get(end_vertex)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 end vertex {end_vertex:?} is stale"
            ),
        })?;
    let curve =
        CurveOnEdge::line_segment(start, end, brep.tolerance.linear()).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot author missing cylinder trim curve: \
                 half-edge {half_edge:?} has degenerate endpoints"
                ),
            }
        })?;
    brep.edges[edge].curve = Some(curve);
    Ok(Some(edge))
}

fn containing_on_cylinder_line_chord_half_edge(
    brep: &BRep,
    face: FaceHandle,
    point: Point3,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        let start_point = brep.vertices.get(half_edge.origin)?.position;
        let end_vertex = half_edge_end_vertex(brep, cursor)?;
        let end_point = brep.vertices.get(end_vertex)?.position;
        if point_lies_on_cylinder_face(brep, face, start_point)
            && point_lies_on_cylinder_face(brep, face, end_point)
            && point_lies_on_segment(point, start_point, end_point, brep.tolerance.linear())
        {
            return Some(cursor);
        }
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    None
}

fn point_lies_on_segment(point: Point3, start: Point3, end: Point3, tolerance: f64) -> bool {
    let edge = end - start;
    let len2 = edge.dot(&edge);
    if len2 <= tolerance.max(1.0e-12).powi(2) {
        return false;
    }
    let t = (point - start).dot(&edge) / len2;
    if t <= tolerance || t >= 1.0 - tolerance {
        return false;
    }
    let projected = start + edge * t;
    point.distance_to(&projected) <= tolerance.max(edge.length() * 1.0e-8)
}

fn point_is_not_edge_endpoint(brep: &BRep, edge: EdgeHandle, point: Point3) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    let Some(twin) = brep.half_edges.get(half_edge.twin) else {
        return false;
    };
    let start = brep.vertices[half_edge.origin].position;
    let end = brep.vertices[twin.origin].position;
    let tol = brep.tolerance.linear();
    point.distance_to(&start) > tol && point.distance_to(&end) > tol
}

fn point_lies_on_cylinder_face(brep: &BRep, face: FaceHandle, point: Point3) -> bool {
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
    else {
        return false;
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return false;
    };
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return false;
    }
    let rel = point - origin;
    let axial = rel.dot(&axis) / axis_len2;
    let center = origin + axis * axial;
    let radial = point.distance_to(&center);
    (radial - radius).abs() <= brep.tolerance.linear().max(radius.abs() * 1.0e-6)
}

fn cylinder_radial_vector(
    point: Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
) -> Option<popkern_core::geometry::Vec3> {
    let axis_len2 = axis.dot(&axis);
    if axis_len2 <= 1.0e-24 {
        return None;
    }
    let rel = point - origin;
    let axial = rel.dot(&axis) / axis_len2;
    let center = origin + axis * axial;
    Some(point - center)
}

fn target_cylinder_face_for_fillet(
    brep: &BRep,
    fillet_index: usize,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<FaceHandle, ApexBuildError> {
    let target = fillet_surfaces
        .get(fillet_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("missing fillet surface for FilletEnd[{fillet_index}]"),
        })?;
    brep.faces
        .iter()
        .filter_map(|(face_handle, face)| {
            let surface = face.surface.as_ref()?;
            cylinder_surfaces_match(&surface.geometry, &target.surface).then_some(face_handle)
        })
        .next()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] has no live \
                 cylinder face matching its fillet surface"
            ),
        })
}

fn cylinder_surfaces_match(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    let (
        SurfaceGeometry::Cylinder {
            origin: a_origin,
            axis: a_axis,
            radius: a_radius,
        },
        SurfaceGeometry::Cylinder {
            origin: b_origin,
            axis: b_axis,
            radius: b_radius,
        },
    ) = (a, b)
    else {
        return false;
    };
    if (a_radius - b_radius).abs() > a_radius.abs().max(b_radius.abs()) * 1.0e-5 + 1.0e-9 {
        return false;
    }
    let Some(a_axis_unit) = a_axis.normalized() else {
        return false;
    };
    let Some(b_axis_unit) = b_axis.normalized() else {
        return false;
    };
    if a_axis_unit.dot(&b_axis_unit).abs() < 1.0 - 1.0e-5 {
        return false;
    }
    let delta = *a_origin - *b_origin;
    let axial = delta.dot(&b_axis_unit);
    let perpendicular = delta - b_axis_unit * axial;
    perpendicular.length() <= 1.0e-5
}

fn half_edge_end_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin_data| twin_data.origin)
}

fn arc_half_edges_between(
    loop_plan: &RepeatedApexLoopPlan,
    start_pos: usize,
    end_pos: usize,
) -> Vec<HalfEdgeHandle> {
    let mut out = Vec::new();
    let n = loop_plan.half_edges.len();
    if n == 0 {
        return out;
    }
    let mut pos = start_pos % n;
    loop {
        out.push(loop_plan.half_edges[pos]);
        pos = (pos + 1) % n;
        if pos == end_pos {
            break;
        }
        if out.len() > n {
            break;
        }
    }
    out
}
