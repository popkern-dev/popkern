use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use serde::{Deserialize, Serialize};

#[path = "types_repair_artifact.rs"]
mod types_repair_artifact;

pub use types_repair_artifact::build_blend_repair_action_artifact;

use super::types_exact_boundary::{CapContinuityReport, ExactBoundaryCapReport};
use super::types_repair_exact_boundary::ExactBoundaryCapDiagnostic;

/// Structured diagnostics for repair and hardening actions performed during a
/// blend operation.
#[derive(Debug, Clone, Default)]
pub struct BlendRepairReport {
    /// PK-141: trims resolved through the curved-support region path
    /// (gates the stale-vertex settle — it must not touch bodies the
    /// region machinery never shaped).
    pub region_resolved_trims: usize,
    /// Whole-blend support-face seam relocation diagnostics.
    pub support_seam_relocation: SupportSeamRelocationReport,
    /// Straight edge curves reauthored to terminate at their owning topology
    /// vertices after late blend topology rewrites.
    pub edge_curve_endpoint_coherence_repairs: usize,
    /// Per-apex construction route decisions for vertex-blend/direct-apex route tracking.
    pub apex_routes: Vec<ApexConstructionRouteDecision>,
    /// Structured apex construction diagnostics and declined repair candidates.
    pub apex_diagnostics: Vec<ApexConstructionDiagnostic>,
    /// PK-148 terminal-arc-trim settlement at open fillet chain ends.
    pub terminal_arc_trim: super::TerminalArcTrimReport,
    /// PK-152 blend×blend corner-seam settlement at adjacent-fillet
    /// pair corners.
    pub blend_seam_corner: super::BlendSeamCornerReport,
    /// PK-153 retirement of unglued construction-residual scaffolding
    /// (apex cells, displaced fragments, holding faces with no
    /// production-twinned side) before the result is returned.
    pub residual_scaffolding: crate::residual_scaffolding::ResidualScaffoldingRetirementReport,
    /// PK-244 serving-boundary orientation-flag reconciliation: faces
    /// whose stored flag contradicted their measured loop winding
    /// (or, for pole-degenerate charts, the tangency witness) and
    /// were flipped to the honest value.
    pub orientation_flag_reconciliations: usize,
    /// PK-250 serving-boundary self-intersection retries: direct-apex
    /// corners whose SERVED region measured a face-face
    /// self-intersection — or, since PK-260, whose direct-apex
    /// construction recorded a serving-shape decline (see
    /// [`SERVING_SHAPE_DECLINE_PREFIX`]) — forcing a one-shot
    /// pipeline re-run with the vertex excluded from direct-apex
    /// promotion so the vertex-blend dispatch serves the corner instead.
    pub serving_self_intersection_retries: Vec<ServingSelfIntersectionRetry>,
    /// PK-251 serving-boundary footprint repairs: direct-apex corners
    /// whose served cap crossed a strip support and were repaired IN
    /// PLACE (interior control points pushed out of the strip
    /// cylinder at the validator's own overlap samples) so the
    /// direct-apex serve survives. An uncleared repair falls through
    /// to the PK-250 exclusion retry.
    pub serving_footprint_repairs: Vec<ServingFootprintRepair>,
}

/// Marker prefix of the serving-shape gate's decline reason (PK-260).
/// The direct-apex dispatch formats its
/// [`ApexConstructionRoute::DirectApexFallbackToVertexBlend`] reason with
/// this prefix when a body defect (glued surfaceless residuals or
/// dead boundary sides) survives every serving-shape rung; the
/// serving retry matches it to re-run the pipeline with the corner
/// excluded, because the mid-pipeline fallback runs the vertex-blend
/// builder on direct-apex-convention anchors and serves a degraded
/// body (measured 8.5e-2 edge/surface deviations on the oblique
/// base-corner class) while the full re-run re-places perp-walk anchors
/// and serves clean.
pub const SERVING_SHAPE_DECLINE_PREFIX: &str = "serving-shape defect persists";

/// Construction route used for one smooth multi-fillet apex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApexConstructionRoute {
    /// direct-apex was not requested; the `insert_vertex_blend`
    /// corner-blend rung ladder handled the apex.
    VertexBlend,
    /// direct-apex planner and builder completed; vertex-blend skipped.
    DirectApex,
    /// direct-apex was requested but declined; vertex-blend handled it.
    DirectApexFallbackToVertexBlend,
}

/// Per-apex vertex-blend/direct-apex route diagnostic.
#[derive(Debug, Clone)]
pub struct ApexConstructionRouteDecision {
    /// Original apex vertex being resolved.
    pub vertex: VertexHandle,
    /// Original fillet edges incident to this apex.
    pub incident_edges: Vec<EdgeHandle>,
    /// Route that produced the final topology for this apex.
    pub route: ApexConstructionRoute,
    /// Faces authored by direct-apex, empty for vertex-blend and fallback routes.
    pub built_faces: Vec<FaceHandle>,
    /// Number of corner patch reports produced by this route.
    pub corner_patch_count: usize,
    /// Number of exact-boundary sector ownership records extracted for this route.
    pub exact_boundary_sector_count: usize,
    /// Number of exact-boundary curves prepared without approximation.
    pub exact_boundary_prepared_curve_count: usize,
    /// Number of generated inner spokes prepared without topology mutation.
    pub exact_boundary_inner_spoke_count: usize,
    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub exact_boundary_ribbon_count: usize,
    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub exact_boundary_solved_surface_count: usize,
    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub exact_boundary_pcurve_pair_count: usize,
    /// Candidate-only performance counters and budget evidence for exact-boundary caps.
    pub exact_boundary_performance: ExactBoundaryCandidatePerformanceReport,
    /// Versioned exact-boundary report for GUI and programmatic consumers.
    pub exact_boundary_report: ExactBoundaryCapReport,
    /// direct-apex decline reason when [`ApexConstructionRoute::DirectApexFallbackToVertexBlend`] applies.
    pub fallback_reason: Option<String>,
}

/// Durable exact-boundary candidate counters for route reports, diagnostics,
/// future runtime gating, and UI/API explanations.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactBoundaryCandidatePerformanceReport {
    /// Number of extracted sector ownership records.
    pub sector_count: usize,
    /// Number of exact curves prepared without approximation.
    pub prepared_curve_count: usize,
    /// Number of generated inner spokes prepared without topology mutation.
    pub inner_spoke_count: usize,
    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub support_ribbon_count: usize,
    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub solved_surface_count: usize,
    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub pcurve_pair_count: usize,
    /// Number of support-ribbon samples consumed by candidate preparation.
    pub support_sample_count: usize,
    /// Number of plane support ribbons evaluated.
    pub plane_support_ribbon_count: usize,
    /// Number of cylindrical support ribbons evaluated.
    pub cylinder_support_ribbon_count: usize,
    /// Number of conical support ribbons evaluated.
    pub cone_support_ribbon_count: usize,
    /// Number of spherical support ribbons evaluated.
    pub sphere_support_ribbon_count: usize,
    /// Number of trusted NURBS support ribbons evaluated.
    pub nurbs_support_ribbon_count: usize,
    /// Total NURBS control points allocated for solved exact-boundary sector surfaces.
    pub surface_control_point_count: usize,
    /// Number of materialized sector faces accepted by validation.
    pub validation_sector_count: usize,
    /// Number of pcurve pairs accepted by validation.
    pub validation_pcurve_pair_count: usize,
    /// Number of generated spoke edges accepted by validation.
    pub validation_spoke_edge_count: usize,
    /// Elapsed candidate-preparation and validation time in microseconds.
    pub elapsed_micros: Option<u128>,
    /// First deterministic decline code suitable for a future runtime skip policy.
    pub deterministic_skip_code: Option<String>,
}

/// Structured apex-construction diagnostic emitted by construction planning and
/// repair preflights.
#[derive(Debug, Clone, PartialEq)]
pub enum ApexConstructionDiagnostic {
    /// Exact-boundary cap candidate preparation declined before it could be
    /// committed as the active cap primitive.
    ExactBoundaryCap(ExactBoundaryCapDiagnostic),
    /// Target-cylinder trim replacement could not materialize an exact,
    /// ownership-balanced far-side trim.
    CylinderTrimMaterialization(CylinderTrimMaterializationDiagnostic),
    /// Native Gregory cap/support seam continuity did not meet the commercial
    /// G1 gate after construction.
    CapSupportContinuity(CapSupportContinuityDiagnostic),
    /// Side-support continuation topology was applied before repeated-apex
    /// target-cylinder trim rewrites.
    SideSupportContinuation(SideSupportContinuationDiagnostic),
    /// KP20 support-cell macro piece was split into sub-cells along U and V;
    /// records whether cap-vertex-aware parameters or the uniform fallback
    /// drove each axis, plus the fallback reason when applicable.
    SupportCellSplit(SupportCellSplitDiagnostic),
    /// Coherence-invariant validation outcome at an operator boundary.
    /// Carries the [`popkern_core::topology::CoherenceReport`] populated
    /// by `BRep::validate_coherence`, plus a stage label that names the
    /// operator whose return triggered the validation. An empty report
    /// means every tracked invariant held at the requested scope.
    CoherenceValidation(CoherenceValidationDiagnostic),
}

/// Coherence-validation diagnostic carried inside
/// [`ApexConstructionDiagnostic::CoherenceValidation`].
#[derive(Debug, Clone, PartialEq)]
pub struct CoherenceValidationDiagnostic {
    /// Operator boundary that emitted this diagnostic (e.g.
    /// `"corner_patch_support_completion"`).
    pub stage: String,
    /// Snapshot of the coherence report at the operator boundary.
    pub report: popkern_core::topology::CoherenceReport,
}

/// Per-axis status of a cap-vertex-aware support-cell split.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportCellSplitAxisStatus {
    /// Caller did not supply a cap-vertex projection for this axis;
    /// uniform `i / cells` split was used.
    UniformNoProjection,
    /// Cap-vertex projection produced a sequence that matched
    /// `i / cells` to within tolerance; uniform formula used directly.
    UniformEquivalent,
    /// Cap-vertex-aware sequence applied to the split.
    VertexAware,
    /// Cap-vertex projection produced a sequence that failed admissibility
    /// (`fallback_reason` carries the failure mode); uniform formula used.
    UniformFallback,
}

/// Why a cap-vertex-aware split fell back to the uniform formula.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportCellSplitFallbackReason {
    /// One or more sub-cells would have arc-length below `MIN_SUBCELL_ARC`.
    MinArc,
    /// Newton iteration did not converge for at least one projection.
    DidNotConverge,
    /// More projected parameters than `cells_axis - 1`; current grid
    /// candidate cannot align every cap-vertex.
    Overcount,
    /// Caller supplied an empty projection set (no cap-vertices on this side).
    NoProjection,
}

/// One axis of a support-cell split's diagnostic facts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportCellSplitAxisDiagnostic {
    /// Status enum for this axis.
    pub status: SupportCellSplitAxisStatus,
    /// Reason this axis fell back to uniform (`None` when not a fallback).
    pub fallback_reason: Option<SupportCellSplitFallbackReason>,
    /// Cap-vertices the caller supplied for this axis.
    pub cap_vertex_count: usize,
    /// Cap-vertices that produced a usable interior projection
    /// (Converged in the open interval, not corner-coincident, not duplicate).
    pub projected_count: usize,
    /// Interior splits filled by uniform subdivision because the projection
    /// undercount left gaps (`cells_axis - 1 - projected_count` when this
    /// path is taken; `0` otherwise).
    pub uniform_fill_count: usize,
    /// Minimum sub-cell arc-length observed under the chosen sequence.
    pub min_subcell_arc: f64,
    /// Maximum `|projected_t - i / cells_axis|` across interior splits — the
    /// raw "distance from uniform" before the equivalence shortcut is applied.
    pub uniform_equivalence_max_delta: f64,
}

/// Cap-vertex-aware split outcome for one macro piece (Stage A3 of the
/// cap-vertex-aware split design).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportCellSplitDiagnostic {
    /// Sector index identifying the macro piece's parent sector.
    pub sector_index: usize,
    /// Quadrant within the sector (canonical `Kp20PieceId` index 0..3).
    pub piece_index: usize,
    /// Sub-cell count along U.
    pub cells_u: usize,
    /// Sub-cell count along V.
    pub cells_v: usize,
    /// Cap-vertices that did not lie on any of the four boundary sides.
    pub cap_vertex_interior_count: usize,
    /// U-axis split facts.
    pub u_axis: SupportCellSplitAxisDiagnostic,
    /// V-axis split facts.
    pub v_axis: SupportCellSplitAxisDiagnostic,
}

/// Applied side-support continuation diagnostic for repeated-apex topology
/// ownership.
#[derive(Debug, Clone, PartialEq)]
pub struct SideSupportContinuationDiagnostic {
    /// Candidate stage that emitted the diagnostic.
    pub stage: String,
    /// Stable diagnostic code.
    pub code: String,
    /// Cap face whose repeated-apex trim rewrite required side-support
    /// ownership continuation.
    pub cap_face: FaceHandle,
    /// Semantic topology primitive that was applied.
    pub primitive: String,
    /// Whether a support-ownership bridge payload was materialized.
    pub applied: bool,
    /// Number of bridge faces authored by the primitive.
    pub bridge_face_count: usize,
    /// Number of bridge loops authored by the primitive.
    pub bridge_loop_count: usize,
}

/// Native cap/support continuity diagnostic for a constructed apex cap.
#[derive(Debug, Clone, PartialEq)]
pub struct CapSupportContinuityDiagnostic {
    /// Semantic target, measured result, and decline code for this continuity gate.
    pub continuity: CapContinuityReport,
    /// Cap face whose support seam was measured.
    pub cap_face: FaceHandle,
    /// Adjacent support or fillet face across the worst seam.
    pub support_face: FaceHandle,
    /// Topological edge shared by the cap and support face.
    pub edge: EdgeHandle,
    /// Scope that performed the measurement, for example candidate or live topology.
    pub measurement_scope: String,
    /// Normalized seam parameter where the selected failing residual was measured.
    pub sample_parameter: Option<f64>,
    /// Cap-side pcurve UV at the selected failing sample, when available.
    pub cap_uv: Option<(f64, f64)>,
    /// Support-side pcurve UV at the selected failing sample, when available.
    pub support_uv: Option<(f64, f64)>,
    /// Distance from the cap pcurve image to the owning edge curve at the sample.
    pub cap_curve_gap: Option<f64>,
    /// Distance from the support pcurve image to the owning edge curve at the sample.
    pub support_curve_gap: Option<f64>,
    /// Worst sampled unit-normal parallelism deficit.
    pub worst_g1_coplanarity: f64,
    /// Commercial gate used for the measurement.
    pub gate: f64,
    /// Number of cap/support edges measured for the apex.
    pub edges_measured: usize,
}

/// Materialized target-cylinder trim replacement diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub struct CylinderTrimMaterializationDiagnostic {
    /// Fillet index whose target cylinder trim was inspected.
    pub fillet_index: usize,
    /// Target cylinder face being rewritten.
    pub cylinder_face: FaceHandle,
    /// Start vertex of the declined far-side trim.
    pub start: VertexHandle,
    /// End vertex of the declined far-side trim.
    pub end: VertexHandle,
    /// Selected far-side trim strategy.
    pub strategy: CylinderTrimFarSideStrategy,
    /// Support-plane candidates containing the start endpoint.
    pub start_support_faces: Vec<FaceHandle>,
    /// Support-plane candidates containing the end endpoint.
    pub end_support_faces: Vec<FaceHandle>,
    /// Stable decline reason.
    pub reason: CylinderTrimMaterializationDeclineReason,
    /// Whether exact curve authoring was available for the selected span.
    pub exact_authorable: bool,
}

/// Far-side trim direction selected for target-cylinder loop replacement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderTrimFarSideStrategy {
    /// Forward loop path between target-cylinder anchors.
    Forward,
    /// Reverse loop path between target-cylinder anchors.
    Reverse,
}

/// Stable reason a materialized cylinder trim replacement declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderTrimMaterializationDeclineReason {
    /// Endpoints belong to different support planes and no
    /// plane-plane-cylinder transition point exists.
    NoPlanePlaneCylinderTransition,
    /// The preflight declined for a different currently unclassified reason.
    Other,
}

/// Aggregate support-face seam relocation diagnostics.
#[derive(Debug, Clone, Default)]
pub struct SupportSeamRelocationReport {
    /// Number of support seams authored and committed.
    pub authored: usize,
    /// Number of candidate attempts declined before commit.
    pub declined: usize,
    /// Number of candidate attempts inspected.
    pub candidates: usize,
    /// Maximum endpoint movement applied or attempted.
    pub max_endpoint_motion: f64,
    /// Maximum sampled curve/cap/support disagreement observed.
    pub max_sampled_disagreement: f64,
    /// Transaction rollback reason, when the cap-neighborhood job rolled back.
    pub rollback_reason: Option<SupportSeamRollbackReason>,
    /// Per-candidate decisions retained for API and GUI diagnostics.
    pub decisions: Vec<SupportSeamDecision>,
}

impl SupportSeamRelocationReport {
    pub(crate) fn merge(&mut self, other: SupportSeamRelocationReport) {
        self.authored += other.authored;
        self.declined += other.declined;
        self.candidates += other.candidates;
        self.max_endpoint_motion = self.max_endpoint_motion.max(other.max_endpoint_motion);
        self.max_sampled_disagreement = self
            .max_sampled_disagreement
            .max(other.max_sampled_disagreement);
        self.rollback_reason = self.rollback_reason.or(other.rollback_reason);
        self.decisions.extend(other.decisions);
    }

    pub(crate) fn push(&mut self, decision: SupportSeamDecision) {
        self.candidates += 1;
        if decision.accepted {
            self.authored += 1;
        } else {
            self.declined += 1;
        }
        self.max_endpoint_motion = self.max_endpoint_motion.max(decision.max_endpoint_motion);
        self.max_sampled_disagreement = self
            .max_sampled_disagreement
            .max(decision.max_sampled_disagreement);
        self.decisions.push(decision);
    }
}

/// One support-face seam relocation decision.
#[derive(Debug, Clone)]
pub struct SupportSeamDecision {
    /// Candidate half-edge on the cap face.
    pub half_edge: HalfEdgeHandle,
    /// Gregory cap face being repaired.
    pub cap_face: FaceHandle,
    /// Adjacent support face, when the candidate reached support classification.
    pub support_face: Option<FaceHandle>,
    /// Canonical cap boundary side attempted.
    pub boundary_side: Option<BoundarySide>,
    /// Whether this candidate committed.
    pub accepted: bool,
    /// Decline reason for rejected candidates.
    pub reason: Option<SupportSeamDeclineReason>,
    /// Maximum endpoint movement applied or attempted for this candidate.
    pub max_endpoint_motion: f64,
    /// Maximum sampled curve/cap/support disagreement for this candidate.
    pub max_sampled_disagreement: f64,
}

/// Reason a support-face seam relocation candidate declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportSeamDeclineReason {
    /// Candidate was not adjacent to a planar support face.
    NotPlanarSupportCandidate,
    /// Candidate could not be matched to a canonical Gregory boundary.
    NoCanonicalBoundary,
    /// Canonical boundary endpoints were degenerate.
    DegenerateBoundary,
    /// Endpoint relocation exceeded the local motion gate.
    ExcessiveEndpointMotion,
    /// Authored cap boundary did not match relocated edge vertices.
    CapBoundaryMismatch,
    /// Adjacent support face could not be lifted to contain the seam.
    SupportFaceLiftFailed,
    /// Sampled curve/cap/support agreement failed.
    SampledValidationFailed,
    /// Moved vertices were incoherent with incident edge curves.
    MovedVertexCurveIncoherent,
    /// Candidate would remove the local cap-cylinder neighbor evidence.
    CapCylinderNeighborLost,
    /// Candidate would reduce global Gregory cap-cylinder neighbor evidence.
    GlobalCapCylinderNeighborLost,
}

/// Stable serialized label for a support-seam relocation decline.
pub fn support_seam_decline_reason_label(reason: SupportSeamDeclineReason) -> &'static str {
    match reason {
        SupportSeamDeclineReason::NotPlanarSupportCandidate => "not_planar_support_candidate",
        SupportSeamDeclineReason::NoCanonicalBoundary => "no_canonical_boundary",
        SupportSeamDeclineReason::DegenerateBoundary => "degenerate_boundary",
        SupportSeamDeclineReason::ExcessiveEndpointMotion => "excessive_endpoint_motion",
        SupportSeamDeclineReason::CapBoundaryMismatch => "cap_boundary_mismatch",
        SupportSeamDeclineReason::SupportFaceLiftFailed => "support_face_lift_failed",
        SupportSeamDeclineReason::SampledValidationFailed => "sampled_validation_failed",
        SupportSeamDeclineReason::MovedVertexCurveIncoherent => "moved_vertex_curve_incoherent",
        SupportSeamDeclineReason::CapCylinderNeighborLost => "cap_cylinder_neighbor_lost",
        SupportSeamDeclineReason::GlobalCapCylinderNeighborLost => {
            "global_cap_cylinder_neighbor_lost"
        }
    }
}

/// Reason a support-face seam relocation job rolled back after candidate work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportSeamRollbackReason {
    /// Resulting topology failed loop closure validation.
    LoopClosureFailed,
    /// Existing Gregory cap-cylinder boundary hints would have been lost.
    CapCylinderHintsLost,
    /// Existing Gregory cap-cylinder neighbor evidence would have been lost.
    CapCylinderNeighborsLost,
}

/// Stable serialized label for a support-seam relocation rollback.
pub fn support_seam_rollback_reason_label(reason: SupportSeamRollbackReason) -> &'static str {
    match reason {
        SupportSeamRollbackReason::LoopClosureFailed => "loop_closure_failed",
        SupportSeamRollbackReason::CapCylinderHintsLost => "cap_cylinder_hints_lost",
        SupportSeamRollbackReason::CapCylinderNeighborsLost => "cap_cylinder_neighbors_lost",
    }
}

#[path = "types_repair_serving.rs"]
mod types_repair_serving;
pub use types_repair_serving::{ServingFootprintRepair, ServingSelfIntersectionRetry};
