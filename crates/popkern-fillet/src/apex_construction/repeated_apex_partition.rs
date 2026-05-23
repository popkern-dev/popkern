//! Planned repeated-apex cap partitioning for constructive Gregory apexes.
//!
//! This module owns the non-serial graph rewrite used when a repeated-apex cap
//! loop cannot be normalized by any order of FilletEnd `mef`s. It builds a
//! typed partition from live FaceEdge half-edges, rewires the cap loop to the
//! planned FilletEnd/FaceEdge cycle, and validates topology plus semantic
//! cap-to-cylinder seam ownership before accepting.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveOnEdge, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdge, HalfEdgeHandle, LoopHandle, SharedTrimEdgeInput,
    TopologyRole, Vertex, VertexHandle, commit_shared_trim_edge_transactionally,
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
    SideSupportRegionContinuationOutcome, SideSupportRegionDiagnosticSnapshot,
    apply_selected_side_support_continuation,
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
#[path = "repeated_apex_partition_cap_boundary_rewrite_curve_refresh.rs"]
mod partition_cap_boundary_rewrite_curve_refresh;
#[path = "repeated_apex_partition_cap_boundary_rewrite_diagnostics.rs"]
mod partition_cap_boundary_rewrite_diagnostics;
#[path = "repeated_apex_partition_cylinder_anchors.rs"]
mod partition_cylinder_anchors;
#[path = "repeated_apex_partition_cylinder_geometry.rs"]
mod partition_cylinder_geometry;
#[path = "repeated_apex_partition_cylinder_loop_adoption.rs"]
mod partition_cylinder_loop_adoption;
#[path = "repeated_apex_partition_cylinder_loop_replacement.rs"]
mod partition_cylinder_loop_replacement;
#[path = "repeated_apex_partition_cylinder_trim_fragments.rs"]
mod partition_cylinder_trim_fragments;
#[path = "repeated_apex_partition_cylinder_trim_readiness.rs"]
mod partition_cylinder_trim_readiness;
#[path = "repeated_apex_partition_cylinder_trim_rebuild.rs"]
mod partition_cylinder_trim_rebuild;
#[path = "repeated_apex_partition_cylinder_trim_topology.rs"]
mod partition_cylinder_trim_topology;
#[path = "repeated_apex_partition_cylinder_trim_transaction.rs"]
mod partition_cylinder_trim_transaction;
#[path = "repeated_apex_partition_exact_trim_curves.rs"]
mod partition_exact_trim_curves;
#[path = "repeated_apex_partition_far_side_trim.rs"]
mod partition_far_side_trim;
#[path = "repeated_apex_partition_generator_anchors.rs"]
mod partition_generator_anchors;
#[path = "repeated_apex_partition_ownership.rs"]
mod partition_ownership;
#[path = "repeated_apex_partition_plan.rs"]
mod partition_plan;
#[path = "repeated_apex_partition_plan_curve_refresh.rs"]
mod partition_plan_curve_refresh;
#[path = "repeated_apex_partition_planner.rs"]
mod partition_planner;
#[path = "repeated_apex_partition_rebuild_diagnostics.rs"]
mod partition_rebuild_diagnostics;
#[path = "repeated_apex_partition_residuals.rs"]
mod partition_residuals;
#[path = "repeated_apex_partition_shared_fillet_end_trim.rs"]
mod partition_shared_fillet_end_trim;
#[path = "repeated_apex_partition_support_bridge_trim_claim.rs"]
mod partition_support_bridge_trim_claim;
#[path = "repeated_apex_partition_support_cycle_primitives.rs"]
mod partition_support_cycle_primitives;
#[path = "repeated_apex_partition_support_family.rs"]
mod partition_support_family;
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
use partition_cap_boundary_rewrite::protect_source_loops_before_partition_loop_claim;
use partition_cap_boundary_rewrite_curve_refresh::author_materialized_shared_fillet_end_trims_for_plan;
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
use partition_cylinder_loop_replacement::{
    MaterializedCylinderLoopReplacementRequest, apply_materialized_cylinder_loop_replacements,
};
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
use partition_cylinder_trim_transaction::{
    apply_cylinder_trim_rewrites_transactionally, validate_cylinder_trim_rebuild_plans,
    validate_cylinder_trim_rewrite_readiness,
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
pub(super) use partition_plan::RepeatedApexPartitionExecutionOptions;
use partition_plan::RetiredApexSpan;
use partition_plan_curve_refresh::refresh_fillet_end_curves_after_endpoint_adoption;
use partition_planner::try_plan_repeated_apex_partition;
use partition_rebuild_diagnostics::{
    adjacent_face_across_half_edge, face_surface_kind, fragment_replacement_authoring_status,
    replacement_run_curve_blocker, replacement_run_exact_curve_kind,
    replacement_run_reownership_blocker, replacement_run_reownership_plan,
    retired_trim_runs_by_adjacent_face,
};
use partition_residuals::{
    demote_retired_apex_residual_faces, mark_loop_edges_construction_residual,
    plan_retired_base_residual_loop, residual_strip_loop_with_duplicated_seam,
    validate_existing_target_cylinder_strip,
};
use partition_shared_fillet_end_trim::materialized_shared_fillet_end_cylinder_half_edge;
use partition_support_bridge_trim_claim::{
    SideSupportBridgeTrimClaimRequest, claim_side_support_bridge_boundary_for_cylinder_trim,
};
use partition_support_cycle_primitives::{
    closed_half_edge_cycles_from_unwired_set, half_edge_is_degenerate,
    support_cycle_open_chain_diagnostic, support_split_half_edge_end_vertex,
    vertices_are_coincident,
};
use partition_support_family::{coplanar_residual_support_faces, support_family_faces};
use partition_support_intersections::{
    plane_plane_cylinder_joint_points, support_plane_faces_containing_point,
};
use partition_support_loop_rewrite::{
    materialize_authored_support_boundary_residual_cycles, materialize_closed_support_trim_cycles,
    repair_support_source_loops_after_cylinder_rewire,
};
use partition_support_owners::{
    GeneratorTrimSupportOwnerQuery, support_face_for_generator_trim_side,
    target_cylinder_support_faces, target_cylinder_support_faces_with_semantic,
    unique_far_side_support_face_for_trim_side,
};
use partition_support_split_cycles::materialize_support_split_cycles;
use partition_support_trim_materialization::{
    GeneratorTrimMaterializationRequest, materialized_generator_trim_side,
};
use partition_support_trim_sides::{
    PreMaterializedTrimSide, SupportOwnedTrimSide, SupportOwnedTrimSource,
    pre_materialize_far_side_trim_sides, pre_materialized_trim_sides,
};
use partition_support_trim_validation::validate_replacement_loop_support_side_ownership;

pub(super) type RepeatedApexNormalization = (
    FaceHandle,
    Vec<EdgeHandle>,
    Vec<crate::types::ApexConstructionDiagnostic>,
);
pub(super) type OptionalRepeatedApexNormalization = Option<RepeatedApexNormalization>;

pub(super) fn try_build_planned_repeated_apex_partition(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    execution_options: RepeatedApexPartitionExecutionOptions,
) -> Result<OptionalRepeatedApexNormalization, ApexBuildError> {
    let n_fillets = plan.cap_kind.n_fillets();
    let Some(partition_plan) = try_plan_repeated_apex_partition(
        brep,
        cap_face,
        plan.apex_vertex,
        apex_vertex_handles,
        face_edge_handles,
        n_fillets,
    )?
    else {
        return Ok(None);
    };
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/partition] planned graph partition: {}",
            partition_plan.diagnostic_summary()
        );
    }

    partition_plan
        .apply(
            brep,
            plan,
            apex_vertex_handles,
            fillet_surfaces,
            n_fillets,
            execution_options,
        )
        .map(Some)
}
