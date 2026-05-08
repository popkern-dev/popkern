//! Phase 6c-A2.4 Gregory cap surface assembly.
//!
//! Bridges from A2's [`ApexPlan`] + post-stitch BRep state to the
//! legacy [`apply_three_quad_decomposition_to_blend`] entry point
//! (the K-P promotion ladder).  The bridge:
//!
//! 1. Walks `plan.target_topology.outer_loop_edge_indices` and emits
//!    one [`BoundarySegment`] per apex_edge — `FilletEnd` for
//!    [`TargetEdgeRole::CapVsCyl`], `FaceEdge` for
//!    [`TargetEdgeRole::CapVsPlane`].
//! 2. Resolves the `face` field of each `FaceEdge` by walking the
//!    BRep edge's twin half-edge.
//! 3. Calls the legacy K-P promotion ladder via
//!    `apply_three_quad_decomposition_to_blend`.
//!
//! Plus topological helpers used by `construct_apex_inner` to
//! identify the cap face after apex teardown and merge any unmerged
//! corner sub-faces left by the V-spoke kef sequence on Gregory
//! caps with mixed-host mef sites.
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md`
//! §A2.4.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::error::ApexBuildError;
use super::gregory_cap_boundary::{
    author_support_compatible_face_edge_tangents, find_face_across_edge,
    inject_physical_hex_boundary_cps, missing_adjacent_face_error, oriented_cap_edge_curve,
    validate_physical_sector_handles, validate_planned_gregory_loop_roles,
};
pub(super) use super::gregory_cap_topology::{
    find_face_bounded_by_all_apex_edges, merge_cap_side_accumulators,
};
use super::plan::ApexPlan;
use super::target_topology::TargetCapFace;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::construct::three_quad_apply::{
    ThreeQuadFailure, apply_three_quad_decomposition_to_blend,
};
use crate::types::CornerPatchInfo;
use crate::vertex_classifier::VertexBlendClassification;

#[derive(Debug, Clone, Copy)]
pub(super) struct GregoryCapSector {
    pub sector_index: usize,
    pub fillet_index: usize,
    pub fillet_end_edge: EdgeHandle,
    pub face_edge: EdgeHandle,
}

/// **A2.4 Gregory cap surface assembly** — bridge from A2's plan +
/// post-stitch BRep state to legacy `three_quad_apply::
/// apply_three_quad_decomposition_to_blend`.
///
/// Builds a `BlendBoundary` from `plan.target_topology.apex_edges`
/// (each apex_edge becomes a `BoundarySegment::FilletEnd` or
/// `FaceEdge` based on its role + curve), then invokes the
/// existing K-P promotion ladder.  The ladder produces N or 4N
/// sub-quad Gregory patches whose surfaces are assigned to the
/// poked sub-faces.
#[allow(clippy::too_many_arguments)]
pub(super) fn assemble_gregory_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    apex_position: Point3,
    _apex_vertex_handles: &[VertexHandle],
    apex_edge_handles: &[Option<EdgeHandle>],
    sectors: &[GregoryCapSector],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<Vec<CornerPatchInfo>, ApexBuildError> {
    if fillet_surfaces.is_empty() {
        return Err(ApexBuildError::PlanMismatch {
            reason: "Gregory cap surface assembly needs fillet_surfaces; caller \
                     used construct_apex (legacy entry) instead of \
                     construct_apex_with_surfaces"
                .into(),
        });
    }

    let outer_indices = match &plan.target_topology.cap_face {
        TargetCapFace::Gregory {
            outer_loop_edge_indices,
        } => outer_loop_edge_indices,
        other => {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "assemble_gregory_cap_surface called with non-Gregory cap_face: {other:?}"
                ),
            });
        }
    };
    if outer_indices.len() < 6 || !outer_indices.len().is_multiple_of(2) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory cap outer loop must have 2N >= 6 edges, got {}",
                outer_indices.len()
            ),
        });
    }
    let n = outer_indices.len() / 2;
    validate_physical_sector_handles(sectors, n, apex_edge_handles)?;

    validate_planned_gregory_loop_roles(plan, outer_indices, n)?;

    // Build boundary segments from physical sector-owned BRep edges, not from
    // planner curve copies. The cap half-edge determines orientation.
    let mut segments: Vec<BoundarySegment> = Vec::with_capacity(2 * n);
    for sector_index in 0..n {
        let sector = sectors
            .iter()
            .find(|sector| sector.sector_index == sector_index)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("Gregory sectors missing sector index {sector_index}"),
            })?;
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: sector.fillet_index,
            curve: oriented_cap_edge_curve(brep, cap_face, sector.fillet_end_edge, "FilletEnd")?,
        });
        let face = find_face_across_edge(brep, sector.face_edge, cap_face).ok_or_else(|| {
            missing_adjacent_face_error(brep, cap_face, sector.face_edge, "FaceEdge", sector_index)
        })?;
        segments.push(BoundarySegment::FaceEdge {
            face,
            curve: oriented_cap_edge_curve(brep, cap_face, sector.face_edge, "FaceEdge")?,
        });
    }

    let mut boundary = BlendBoundary {
        vertex: plan.apex_vertex,
        vertex_position: apex_position,
        segments,
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    match crate::blend::sphere_fit_boundary::apply_sphere_fit_substitution(
        &mut boundary,
        fillet_surfaces,
        brep,
        tolerance,
    ) {
        Ok(()) => {
            boundary.sphere_fit_substituted = true;
        }
        Err(reason) => {
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                eprintln!("[a2-trace] Gregory sphere-fit boundary substitution skipped: {reason}");
            }
        }
    }
    inject_physical_hex_boundary_cps(&mut boundary)?;
    if boundary.sphere_fit_substituted {
        author_support_compatible_face_edge_tangents(&mut boundary, fillet_surfaces, brep);
    }

    // Classification — A2 fires only at SameSideConvex apices today
    // (cube/pyramid mixed-radius corners).  A2.5 will lift this to
    // honor the planner's own classification once OppositeSide /
    // saddle handling lands.
    let classification = VertexBlendClassification::SameSideConvex;

    // Drive the legacy K-P promotion ladder with all rungs enabled.
    // The ladder picks the highest-quality rung that converges; for
    // cube_mixed N=3 the typical landing is biquintic-full or two-
    // stripe with K-P binonic promotion above.  These flags mirror
    // the legacy default-on dispatch when the ladder is fully
    // exercised; subsequent passes can plumb FilletOptions through
    // for user-driven knobs.
    let result = apply_three_quad_decomposition_to_blend(
        brep,
        cap_face,
        &boundary,
        fillet_surfaces,
        plan.apex_vertex,
        n,
        classification,
        true, // use_deg5
        true, // use_biquintic_full
        true, // use_biseptic
        true, // use_two_stripe
        true, // use_kp_binonic
        true, // use_kp_polar
    );

    match result {
        Ok(patches) => Ok(patches),
        Err(ThreeQuadFailure::Overlap) => Err(ApexBuildError::PlanMismatch {
            reason: "Gregory three-quad decomposition self-intersects (overlap remediation)".into(),
        }),
        Err(ThreeQuadFailure::Other(s)) => Err(ApexBuildError::PlanMismatch {
            reason: format!("Gregory three-quad decomposition failed: {s}"),
        }),
    }
}
