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
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md`
//! §A2.4.

use popkern_core::geometry::{CurveEval, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::cap_kind::CapKind;
use super::constructive_geometry::cylinder_arc_through_edge_endpoints;
use super::error::ApexBuildError;
use super::gregory_cap_boundary::{
    author_support_compatible_face_edge_tangents, find_face_across_edge,
    half_edge_for_face_on_edge, inject_physical_hex_boundary_cps, missing_adjacent_face_error,
    oriented_cap_edge_curve, validate_physical_sector_handles, validate_planned_gregory_loop_roles,
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

    // The K-P cap builder consumes sectors in planner order: FilletEnd[k]
    // followed by the support FaceEdge[k]. Physical BRep loop traversal can
    // rotate after repeated-apex topology normalization, so boundary assembly
    // stays tied to the sector map and validates/repairs individual curve
    // orientation instead of reordering sectors.
    let segments = gregory_boundary_segments_from_sector_order(
        brep,
        cap_face,
        plan,
        sectors,
        fillet_surfaces,
        n,
    )?;
    trace_gregory_boundary_segments(&segments);

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
    if !apply_planned_gregory_boundary(&mut boundary, plan, n)? {
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
                if crate::a2_trace::enabled() {
                    eprintln!(
                        "[a2-trace] Gregory sphere-fit boundary substitution skipped: {reason}"
                    );
                }
            }
        }
    }
    inject_physical_hex_boundary_cps(&mut boundary)?;
    if boundary.sphere_fit_substituted {
        author_support_compatible_face_edge_tangents(&mut boundary, fillet_surfaces, brep);
    }

    let classification = VertexBlendClassification::SameSideConvex;

    if let Some(patches) = super::pns_cap_candidate::try_boundary_preserving_pns_cap(
        brep,
        cap_face,
        &boundary,
        fillet_surfaces,
        plan.apex_vertex,
        n,
        classification,
        sectors,
    ) {
        return Ok(patches);
    }

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

fn apply_planned_gregory_boundary(
    boundary: &mut BlendBoundary,
    plan: &ApexPlan,
    n: usize,
) -> Result<bool, ApexBuildError> {
    let CapKind::Gregory {
        per_fillet_terminations,
        per_fillet_caps,
    } = &plan.cap_kind
    else {
        return Ok(false);
    };
    if per_fillet_terminations.len() != n || per_fillet_caps.len() != n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory cap boundary plan expected {n} terminations and cap spheres, got {}/{}",
                per_fillet_terminations.len(),
                per_fillet_caps.len()
            ),
        });
    }
    if boundary.segments.len() != 2 * n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory cap boundary has {} segments for {n} fillets",
                boundary.segments.len()
            ),
        });
    }

    for k in 0..n {
        let [start, end] = per_fillet_terminations[k];
        let fillet_segment =
            boundary
                .segments
                .get_mut(2 * k)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!("Gregory boundary missing FilletEnd segment {k}"),
                })?;
        let BoundarySegment::FilletEnd { curve, .. } = fillet_segment else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!("Gregory boundary segment {} is not a FilletEnd", 2 * k),
            });
        };
        *curve = align_planned_curve_to_endpoints(curve, start, end, 1.0e-6)
            .unwrap_or_else(|| line_curve(start, end));

        let next = (k + 1) % n;
        let bridge_start = end;
        let bridge_end = per_fillet_terminations[next][0];
        let face_segment =
            boundary
                .segments
                .get_mut(2 * k + 1)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!("Gregory boundary missing FaceEdge segment {k}"),
                })?;
        let BoundarySegment::FaceEdge { curve, .. } = face_segment else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!("Gregory boundary segment {} is not a FaceEdge", 2 * k + 1),
            });
        };
        *curve = line_curve(bridge_start, bridge_end);
    }
    boundary.per_fillet_caps = per_fillet_caps.clone();
    boundary.sphere_fit_substituted = true;
    Ok(true)
}

fn gregory_boundary_segments_from_sector_order(
    brep: &BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    n: usize,
) -> Result<Vec<BoundarySegment>, ApexBuildError> {
    let mut segments = Vec::with_capacity(2 * n);
    for sector_index in 0..n {
        let sector = sectors
            .iter()
            .find(|sector| sector.sector_index == sector_index)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("Gregory sectors missing sector index {sector_index}"),
            })?;
        segments.push(gregory_fillet_segment(
            brep,
            cap_face,
            plan,
            sector,
            fillet_surfaces,
        )?);
        segments.push(gregory_face_segment(brep, cap_face, sector)?);
    }
    Ok(segments)
}

fn gregory_fillet_segment(
    brep: &BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sector: &GregoryCapSector,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Result<BoundarySegment, ApexBuildError> {
    let planned_fillet_curve = plan
        .target_topology
        .apex_edges
        .get(2 * sector.sector_index)
        .map(|edge| &edge.curve)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory sector {} missing planned FilletEnd curve",
                sector.sector_index
            ),
        })?;
    Ok(BoundarySegment::FilletEnd {
        fillet_index: sector.fillet_index,
        curve: oriented_fillet_end_boundary_curve(
            brep,
            cap_face,
            sector.fillet_end_edge,
            sector.fillet_index,
            planned_fillet_curve,
            fillet_surfaces,
        )?,
    })
}

fn gregory_face_segment(
    brep: &BRep,
    cap_face: FaceHandle,
    sector: &GregoryCapSector,
) -> Result<BoundarySegment, ApexBuildError> {
    let face = find_face_across_edge(brep, sector.face_edge, cap_face).ok_or_else(|| {
        missing_adjacent_face_error(
            brep,
            cap_face,
            sector.face_edge,
            "FaceEdge",
            sector.sector_index,
        )
    })?;
    Ok(BoundarySegment::FaceEdge {
        face,
        curve: oriented_cap_edge_curve(brep, cap_face, sector.face_edge, "FaceEdge")?,
    })
}

fn oriented_fillet_end_boundary_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    fillet_index: usize,
    planned_curve: &popkern_core::geometry::NurbsCurve,
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
) -> Result<popkern_core::geometry::NurbsCurve, ApexBuildError> {
    let cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory FilletEnd edge {edge:?} is not bounded by cap face {cap_face:?}"
            ),
        }
    })?;
    let start = brep
        .half_edges
        .get(cap_half_edge)
        .and_then(|half_edge| brep.vertices.get(half_edge.origin))
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory FilletEnd[{fillet_index}] edge {edge:?} has stale cap-side start"
            ),
        })?;
    let end = brep
        .half_edges
        .get(cap_half_edge)
        .and_then(|half_edge| brep.half_edges.get(half_edge.twin))
        .and_then(|twin| brep.vertices.get(twin.origin))
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory FilletEnd[{fillet_index}] edge {edge:?} has stale cap-side end"
            ),
        })?;
    let surface = fillet_surfaces
        .get(fillet_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("Gregory FilletEnd[{fillet_index}] has no fillet surface"),
        })?
        .surface
        .clone();
    if let Some((geometry, t_start, t_end)) =
        cylinder_arc_through_edge_endpoints(&surface, start, end)
    {
        let curve = geometry.to_nurbs_in_range((t_start, t_end));
        if curve_matches_endpoints(&curve, start, end, brep.tolerance.linear().max(1.0e-6)) {
            return Ok(curve);
        }
        if crate::a2_trace::enabled() {
            let (u0, u1) = curve.parameter_range();
            let p0 = curve.point_at(u0);
            let p1 = curve.point_at(u1);
            eprintln!(
                "[a2-trace §A2.4/gregory-boundary] FilletEnd[{fillet_index}] analytic cylinder \
                 arc on {edge:?} does not match cap endpoints: curve distances {:.3e}/{:.3e}; \
                 using planner curve aligned to cap half-edge",
                p0.distance_to(&start),
                p1.distance_to(&end)
            );
        }
    }
    if let Some(curve) = align_planned_curve_to_endpoints(
        planned_curve,
        start,
        end,
        brep.tolerance.linear().max(1.0e-6),
    ) {
        return Ok(curve);
    }
    oriented_cap_edge_curve(brep, cap_face, edge, "FilletEnd")
}

fn align_planned_curve_to_endpoints(
    curve: &popkern_core::geometry::NurbsCurve,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> Option<popkern_core::geometry::NurbsCurve> {
    if curve_matches_endpoints(curve, start, end, tolerance) {
        return Some(curve.clone());
    }
    let reversed = reverse_nurbs_curve(curve);
    curve_matches_endpoints(&reversed, start, end, tolerance).then_some(reversed)
}

fn curve_matches_endpoints(
    curve: &popkern_core::geometry::NurbsCurve,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> bool {
    let (t0, t1) = curve.parameter_range();
    curve.point_at(t0).distance_to(&start) <= tolerance
        && curve.point_at(t1).distance_to(&end) <= tolerance
}

fn reverse_nurbs_curve(
    curve: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::NurbsCurve {
    let first = *curve.knots.first().unwrap_or(&0.0);
    let last = *curve.knots.last().unwrap_or(&1.0);
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|k| first + last - *k)
        .collect();
    knots.sort_by(f64::total_cmp);
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    popkern_core::geometry::NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}

fn line_curve(start: Point3, end: Point3) -> popkern_core::geometry::NurbsCurve {
    popkern_core::geometry::NurbsCurve {
        degree: 1,
        control_points: vec![start, end],
        knots: vec![0.0, 0.0, 1.0, 1.0],
        weights: vec![1.0, 1.0],
    }
}

fn trace_gregory_boundary_segments(segments: &[BoundarySegment]) {
    if !crate::a2_trace::enabled() {
        return;
    }
    for (index, segment) in segments.iter().enumerate() {
        let (label, curve) = match segment {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => (format!("FilletEnd[{fillet_index}]"), curve),
            BoundarySegment::FaceEdge { face, curve } => (format!("FaceEdge[{face:?}]"), curve),
        };
        let (t0, t1) = curve.parameter_range();
        let start = curve.point_at(t0);
        let end = curve.point_at(t1);
        eprintln!(
            "[a2-trace §A2.4/gregory-boundary] segment[{index}] {label}: \
             ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6})",
            start.x, start.y, start.z, end.x, end.y, end.z
        );
    }
}
