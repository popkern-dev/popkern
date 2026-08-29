//! Cap surface assignment for constructive apex building.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{
    CornerPatchInfo, ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapReport,
};

use super::constructive_geometry::planned_fillet_index_for_sector;
use super::error::ApexBuildError;
use super::exact_boundary_cap::facade::{
    ExactBoundaryCapCandidateStatus, run_exact_boundary_cap_candidate,
};
use super::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;
use super::gregory_cap::GregoryCapSector;
use super::plan::ApexPlan;

/// Gregory cap assignment output, including exact-boundary preflight metrics.
pub(super) struct GregoryCapAssignment {
    pub corner_patches: Vec<CornerPatchInfo>,
    pub exact_boundary_sectors: Vec<ExactBoundaryCapSectorInput>,
    pub exact_boundary_committed_faces: Vec<FaceHandle>,
    pub exact_boundary_prepared_curve_count: usize,
    pub exact_boundary_inner_spoke_count: usize,
    pub exact_boundary_ribbon_count: usize,
    pub exact_boundary_solved_surface_count: usize,
    pub exact_boundary_pcurve_pair_count: usize,
    pub exact_boundary_performance: ExactBoundaryCandidatePerformanceReport,
    pub exact_boundary_report: ExactBoundaryCapReport,
    pub exact_boundary_diagnostics: Vec<crate::types::ApexConstructionDiagnostic>,
}

/// Assign `SurfaceGeometry::Sphere` to the cap face.
pub(super) fn assign_sphere_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    center: Point3,
    radius: f64,
) -> Result<(), ApexBuildError> {
    let face = brep
        .faces
        .get_mut(cap_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("assign_sphere_cap_surface: cap_face {cap_face:?} no longer exists"),
        })?;
    face.surface = Some(SurfaceOnFace {
        geometry: {
            if popkern_core::trace::trace_enabled("sphere-fit") {
                eprintln!("[corner-method] apex-sphere r={radius:.3}");
            }
            SurfaceGeometry::Sphere { center, radius }
        },
        orientation: true,
    });
    Ok(())
}

/// Delegate Gregory cap surface construction to `assemble_gregory_cap_surface`.
#[allow(clippy::too_many_arguments)]
pub(super) fn assign_gregory_cap_surface(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    apex_position: Point3,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[Option<EdgeHandle>],
    fillet_end_handles: &[EdgeHandle],
    sectors: &[GregoryCapSector],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<GregoryCapAssignment, ApexBuildError> {
    let n = plan.cap_kind.n_fillets();
    if face_edge_handles.len() != n || fillet_end_handles.len() != n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "assign_gregory_cap_surface: expected N={n} face edges + N fillet ends, \
                 got {} face_edge_handles, {} fillet_end_handles",
                face_edge_handles.len(),
                fillet_end_handles.len()
            ),
        });
    }
    validate_gregory_sectors(n, sectors)?;

    let mut apex_edge_handles: Vec<Option<EdgeHandle>> = Vec::with_capacity(2 * n);
    for k in 0..n {
        apex_edge_handles.push(Some(fillet_end_handles[k]));
        // None at welded shared corners — no FaceEdge exists there.
        apex_edge_handles.push(face_edge_handles[k]);
    }
    let face_edge_support_surfaces =
        face_edge_support_surfaces(brep, cap_face, face_edge_handles, tolerance.linear());

    let exact_boundary_candidate = run_exact_boundary_cap_candidate(brep, cap_face, plan, sectors)?;
    if exact_boundary_candidate.status == ExactBoundaryCapCandidateStatus::Accepted {
        return Ok(GregoryCapAssignment {
            corner_patches: exact_boundary_corner_patches(plan, &exact_boundary_candidate),
            exact_boundary_prepared_curve_count: exact_boundary_candidate.prepared_curve_count(),
            exact_boundary_inner_spoke_count: exact_boundary_candidate.inner_spoke_count(),
            exact_boundary_ribbon_count: exact_boundary_candidate.ribbon_count(),
            exact_boundary_solved_surface_count: exact_boundary_candidate.solved_surface_count(),
            exact_boundary_pcurve_pair_count: exact_boundary_candidate.pcurve_pair_count(),
            exact_boundary_sectors: exact_boundary_candidate.sectors,
            exact_boundary_committed_faces: exact_boundary_candidate.committed_faces,
            exact_boundary_performance: exact_boundary_candidate.performance,
            exact_boundary_report: exact_boundary_candidate.report,
            exact_boundary_diagnostics: exact_boundary_candidate.diagnostics,
        });
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace cap-assign] exact-boundary declined; falling through to \
             assemble_gregory_cap_surface for cap {cap_face:?}"
        );
    }
    let surface_assembly = crate::apex_construction::gregory_cap::assemble_gregory_cap_surface(
        brep,
        cap_face,
        plan,
        apex_position,
        apex_vertex_handles,
        &apex_edge_handles,
        sectors,
        fillet_surfaces,
        tolerance,
    )?;
    let corner_patches = surface_assembly.corner_patches;
    restore_face_edge_support_surfaces(brep, &face_edge_support_surfaces, tolerance.linear());
    author_plane_side_face_edge_pcurves(brep, face_edge_handles);
    let exact_boundary_prepared_curve_count = exact_boundary_candidate.prepared_curve_count();
    let exact_boundary_inner_spoke_count = exact_boundary_candidate.inner_spoke_count();
    let exact_boundary_ribbon_count = exact_boundary_candidate.ribbon_count();
    let exact_boundary_solved_surface_count = exact_boundary_candidate.solved_surface_count();
    let exact_boundary_pcurve_pair_count = exact_boundary_candidate.pcurve_pair_count();
    let mut exact_boundary_diagnostics = exact_boundary_candidate.diagnostics;
    exact_boundary_diagnostics.extend(surface_assembly.diagnostics);
    Ok(GregoryCapAssignment {
        corner_patches,
        exact_boundary_prepared_curve_count,
        exact_boundary_inner_spoke_count,
        exact_boundary_ribbon_count,
        exact_boundary_solved_surface_count,
        exact_boundary_pcurve_pair_count,
        exact_boundary_sectors: exact_boundary_candidate.sectors,
        exact_boundary_committed_faces: exact_boundary_candidate.committed_faces,
        exact_boundary_performance: exact_boundary_candidate.performance,
        exact_boundary_report: exact_boundary_candidate.report,
        exact_boundary_diagnostics,
    })
}

/// Semantic ownership records for an ACCEPTED exact-boundary cap
/// candidate (PK-119).
///
/// The exact-boundary route (one-patch Coons wedge, or the multipatch
/// partition fallback) commits per-sector NURBS faces whose ownership
/// previously lived only on the apex route report — corner-patch
/// consumers (continuity bookkeeping, fixtures, GUI selectors) could
/// not find the cap without route-report archaeology. Emit one
/// `CornerPatchInfo` per committed face:
///
/// - `method` is `Coons` — the construction is a Coons wedge solve
///   with ribbon derivative rows, not a Gregory assembly; the
///   post-pipeline method reclassifier preserves it (like PnS).
/// - `achieved_continuity` is `G1`: candidate acceptance gates the
///   cap/support unit-normal coplanarity deficit (≤5e-2) along every
///   accepted seam in addition to the G0 residual gates. Known
///   evidence-density limitation on the pyramid4 class is tracked as
///   PK-125. A measured continuity report, when present, wins.
/// - `exact_boundary_cap` quarantines these faces from vertex-blend
///   pcurve-completion passes — acceptance already gated every
///   cap/support pcurve pair transactionally.
fn exact_boundary_corner_patches(
    plan: &ApexPlan,
    candidate: &super::exact_boundary_cap::facade::ExactBoundaryCapCandidateResult,
) -> Vec<CornerPatchInfo> {
    let achieved = candidate
        .report
        .continuity
        .achieved
        .unwrap_or(crate::types::ContinuityLevel::G1);
    candidate
        .committed_faces
        .iter()
        .map(|&face| CornerPatchInfo {
            vertex: plan.apex_vertex,
            n_sides: plan.cap_kind.n_fillets(),
            method: crate::types::CornerMethod::Coons,
            achieved_continuity: achieved,
            classification: crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
            face,
            bridge_fallbacks: Vec::new(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: None,
            cap_sphere: None,
            setback_cap: false,
            exact_boundary_cap: true,
            sub_patch_provenance: None,
        })
        .collect()
}

fn face_edge_support_surfaces(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edge_handles: &[Option<EdgeHandle>],
    tolerance: f64,
) -> Vec<(FaceHandle, SurfaceOnFace)> {
    face_edge_handles
        .iter()
        .flatten()
        .filter_map(|edge| {
            let support_face = face_across_edge(brep, *edge, cap_face)?;
            let surface = brep.faces.get(support_face)?.surface.as_ref()?.clone();
            edge_curve_fits_surface(brep, *edge, &surface.geometry, tolerance)
                .then_some((support_face, surface))
        })
        .collect()
}

fn restore_face_edge_support_surfaces(
    brep: &mut BRep,
    surfaces: &[(FaceHandle, SurfaceOnFace)],
    tolerance: f64,
) {
    for (face, surface) in surfaces {
        if !edge_owner_surface_still_needed(brep, *face, &surface.geometry, tolerance) {
            continue;
        }
        let Some(face_data) = brep.faces.get_mut(*face) else {
            continue;
        };
        face_data.surface = Some(surface.clone());
    }
}

fn edge_owner_surface_still_needed(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    brep.faces.get(face).is_some_and(|_| {
        face_outer_edges(brep, face)
            .into_iter()
            .any(|edge| edge_curve_fits_surface(brep, edge, surface, tolerance))
    })
}

fn face_outer_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut out = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            break;
        };
        out.push(half_edge.edge);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    out
}

fn face_across_edge(brep: &BRep, edge: EdgeHandle, cap_face: FaceHandle) -> Option<FaceHandle> {
    let edge_data = brep.edges.get(edge)?;
    let first = edge_data.half_edge;
    let second = brep.half_edges.get(first)?.twin;
    [first, second].into_iter().find_map(|half_edge| {
        let face = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?.face;
        (face != cap_face).then_some(face)
    })
}

fn edge_curve_fits_surface(
    brep: &BRep,
    edge: EdgeHandle,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    let Some(curve) = brep.edges.get(edge).and_then(|edge| edge.curve.as_ref()) else {
        return false;
    };
    let (t0, t1) = curve.geometry.parameter_range();
    (0..12).all(|index| {
        let s = index as f64 / 11.0;
        let point = curve.geometry.point_at(t0 + s * (t1 - t0));
        let Some((u, v)) = project_point_to_uv(surface, &point) else {
            return false;
        };
        surface.point_at(u, v).distance_to(&point) <= tolerance.max(1.0e-6)
    })
}

pub(super) fn author_plane_side_face_edge_pcurves(
    brep: &mut BRep,
    face_edge_handles: &[Option<EdgeHandle>],
) {
    for &edge in face_edge_handles.iter().flatten() {
        let Some(edge_data) = brep.edges.get(edge).cloned() else {
            continue;
        };
        let Some(curve) = edge_data.curve.as_ref() else {
            continue;
        };
        let Some(primary) = brep.half_edges.get(edge_data.half_edge).cloned() else {
            continue;
        };
        for (half_edge, forward) in [(edge_data.half_edge, true), (primary.twin, false)] {
            let Some(face) = brep
                .half_edges
                .get(half_edge)
                .and_then(|half_edge_data| brep.loops.get(half_edge_data.loop_))
                .map(|loop_data| loop_data.face)
            else {
                continue;
            };
            let Some(surface) = brep
                .faces
                .get(face)
                .and_then(|face_data| face_data.surface.as_ref())
            else {
                continue;
            };
            if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
                continue;
            }
            let Some(pcurve) = projected_curve_pcurve(&curve.geometry, &surface.geometry, forward)
            else {
                continue;
            };
            brep.half_edges[half_edge].pcurve = Some(pcurve);
        }
    }
}

fn projected_curve_pcurve(
    curve: &impl CurveEval,
    surface: &SurfaceGeometry,
    forward: bool,
) -> Option<Pcurve> {
    let (t_start, t_end) = curve.parameter_range();
    let samples = 33usize;
    let mut uv_points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for i in 0..samples {
        let s = i as f64 / (samples - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let point = curve.point_at(t_start + curve_s * (t_end - t_start));
        let (u, v) = project_point_to_uv(surface, &point)?;
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &params, 5)
        .ok()
        .map(Into::into)
}

pub(super) fn gregory_sectors_from_handles(
    plan: &ApexPlan,
    fillet_end_handles: &[EdgeHandle],
    face_edge_handles: &[Option<EdgeHandle>],
) -> Result<Vec<GregoryCapSector>, ApexBuildError> {
    if fillet_end_handles.len() != face_edge_handles.len() {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "gregory sectors need equal FilletEnd/FaceEdge counts, got {} and {}",
                fillet_end_handles.len(),
                face_edge_handles.len()
            ),
        });
    }
    fillet_end_handles
        .iter()
        .zip(face_edge_handles)
        .enumerate()
        .map(|(sector_index, (&fillet_end_edge, &face_edge))| {
            if face_edge.is_none() && !plan.cap_kind.corner_is_welded(sector_index) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory sector {sector_index} has no FaceEdge but its corner is not \
                         declared welded"
                    ),
                });
            }
            Ok(GregoryCapSector {
                sector_index,
                fillet_index: planned_fillet_index_for_sector(plan, sector_index)?,
                fillet_end_edge,
                face_edge,
            })
        })
        .collect()
}

fn validate_gregory_sectors(
    n_fillets: usize,
    sectors: &[GregoryCapSector],
) -> Result<(), ApexBuildError> {
    if sectors.len() != n_fillets {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "expected {n_fillets} Gregory cap sectors, got {}",
                sectors.len()
            ),
        });
    }
    let mut sector_indices_seen = vec![false; n_fillets];
    let mut fillets_seen = vec![false; n_fillets];
    let mut fillet_edges_seen = std::collections::BTreeSet::new();
    let mut face_edges_seen = std::collections::BTreeSet::new();
    for sector in sectors {
        let Some(seen) = sector_indices_seen.get_mut(sector.sector_index) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory cap sector has out-of-range sector index {} for N={n_fillets}",
                    sector.sector_index
                ),
            });
        };
        if *seen {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory cap sectors repeat sector index {}",
                    sector.sector_index
                ),
            });
        }
        *seen = true;
        let Some(seen) = fillets_seen.get_mut(sector.fillet_index) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory cap sector has out-of-range fillet index {} for N={n_fillets}",
                    sector.fillet_index
                ),
            });
        };
        if *seen {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory cap sector repeats fillet index {}",
                    sector.fillet_index
                ),
            });
        }
        *seen = true;
        if !fillet_edges_seen.insert(sector.fillet_end_edge) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory cap sectors repeat FilletEnd edge {:?}",
                    sector.fillet_end_edge
                ),
            });
        }
        if let Some(face_edge) = sector.face_edge
            && !face_edges_seen.insert(face_edge)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!("Gregory cap sectors repeat FaceEdge {face_edge:?}"),
            });
        }
    }
    Ok(())
}
