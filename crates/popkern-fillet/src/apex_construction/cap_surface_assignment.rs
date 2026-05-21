//! Cap surface assignment for constructive apex building.

use std::time::Instant;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{CornerPatchInfo, ExactBoundaryCandidatePerformanceReport};

use super::constructive_geometry::planned_fillet_index_for_sector;
use super::error::ApexBuildError;
use super::exact_boundary_cap::boundary::{
    ExactBoundaryCurveBudget, prepare_exact_boundary_curves,
};
use super::exact_boundary_cap::commit::commit_exact_boundary_sector_surfaces;
use super::exact_boundary_cap::corner::{
    ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes,
};
use super::exact_boundary_cap::facts::attach_shared_repair_fact_ids;
use super::exact_boundary_cap::pcurves::{
    ExactBoundaryPcurveBudget, prepare_exact_boundary_pcurves,
};
use super::exact_boundary_cap::performance::{
    finalize_performance_report, record_pcurve_report, record_validation_report,
    support_sample_count, surface_control_point_count,
};
use super::exact_boundary_cap::ribbon::{
    ExactBoundaryRibbonBudget, prepare_exact_boundary_ribbons,
};
use super::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, extract_exact_boundary_cap_sectors,
};
use super::exact_boundary_cap::solve::{ExactBoundarySolveBudget, solve_exact_boundary_surfaces};
use super::exact_boundary_cap::topology::materialize_exact_boundary_sector_topology;
use super::exact_boundary_cap::validation::{
    ExactBoundaryValidationBudget, validate_exact_boundary_candidate,
};
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
        geometry: SurfaceGeometry::Sphere { center, radius },
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
    face_edge_handles: &[EdgeHandle],
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
        apex_edge_handles.push(Some(face_edge_handles[k]));
    }

    let exact_boundary_started = Instant::now();
    let exact_boundary_transaction_snapshot = brep.clone();
    let mut exact_boundary_sectors =
        extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors)
            .map_err(|diagnostic| diagnostic.into_build_error())?;
    let mut exact_boundary_performance = ExactBoundaryCandidatePerformanceReport {
        sector_count: exact_boundary_sectors.len(),
        ..Default::default()
    };
    let prepared_boundary_curves = prepare_exact_boundary_curves(
        brep,
        &exact_boundary_sectors,
        ExactBoundaryCurveBudget::default(),
    )
    .map_err(|diagnostic| diagnostic.into_build_error())?;
    let exact_boundary_prepared_curve_count = prepared_boundary_curves.curves.len();
    exact_boundary_performance.prepared_curve_count = exact_boundary_prepared_curve_count;
    let mut exact_boundary_diagnostics = Vec::new();
    let prepared_inner_spokes = match prepare_exact_boundary_inner_spokes(
        &exact_boundary_sectors,
        &prepared_boundary_curves,
        ExactBoundaryInnerSpokeBudget::default(),
    ) {
        Ok(prepared_spokes) => Some(prepared_spokes),
        Err(diagnostic) => {
            exact_boundary_diagnostics.push(diagnostic.into_apex_diagnostic());
            None
        }
    };
    let exact_boundary_inner_spoke_count = prepared_inner_spokes
        .as_ref()
        .map(|prepared_spokes| prepared_spokes.sectors.len() * 2)
        .unwrap_or(0);
    exact_boundary_performance.inner_spoke_count = exact_boundary_inner_spoke_count;
    let prepared_ribbons = match prepare_exact_boundary_ribbons(
        brep,
        &exact_boundary_sectors,
        &prepared_boundary_curves,
        ExactBoundaryRibbonBudget::default(),
    ) {
        Ok(prepared_ribbons) => Some(prepared_ribbons),
        Err(diagnostic) => {
            exact_boundary_diagnostics.push(diagnostic.into_apex_diagnostic());
            None
        }
    };
    let exact_boundary_ribbon_count = prepared_ribbons
        .as_ref()
        .map(|prepared_ribbons| prepared_ribbons.ribbons.len())
        .unwrap_or(0);
    exact_boundary_performance.support_ribbon_count = exact_boundary_ribbon_count;
    exact_boundary_performance.support_sample_count = prepared_ribbons
        .as_ref()
        .map(support_sample_count)
        .unwrap_or(0);
    let prepared_solved_surfaces = match (prepared_inner_spokes.as_ref(), prepared_ribbons.as_ref())
    {
        (Some(prepared_inner_spokes), Some(prepared_ribbons)) => {
            match solve_exact_boundary_surfaces(
                prepared_inner_spokes,
                prepared_ribbons,
                ExactBoundarySolveBudget::default(),
            ) {
                Ok(solved) => Some(solved),
                Err(diagnostic) => {
                    exact_boundary_diagnostics.push(diagnostic.into_apex_diagnostic());
                    None
                }
            }
        }
        _ => None,
    };
    let exact_boundary_solved_surface_count = prepared_solved_surfaces
        .as_ref()
        .map(|solved| solved.validated_surface_count())
        .unwrap_or(0);
    exact_boundary_performance.solved_surface_count = exact_boundary_solved_surface_count;
    exact_boundary_performance.surface_control_point_count = prepared_solved_surfaces
        .as_ref()
        .map(surface_control_point_count)
        .unwrap_or(0);
    let mut exact_boundary_committed = false;
    let mut exact_boundary_committed_faces = Vec::new();
    let exact_boundary_pcurve_pair_count = match (
        prepared_inner_spokes.as_ref(),
        prepared_ribbons.as_ref(),
        prepared_solved_surfaces.as_ref(),
    ) {
        (Some(prepared_inner_spokes), Some(prepared_ribbons), Some(prepared_solved_surfaces)) => {
            match prepare_exact_boundary_pcurves(
                prepared_inner_spokes,
                prepared_ribbons,
                prepared_solved_surfaces,
                ExactBoundaryPcurveBudget::default(),
            ) {
                Ok(prepared_pcurves) => {
                    let pair_count = prepared_pcurves.validated_pair_count();
                    record_pcurve_report(&mut exact_boundary_performance, &prepared_pcurves);
                    match materialize_exact_boundary_sector_topology(
                        brep,
                        cap_face,
                        &exact_boundary_sectors,
                        prepared_inner_spokes,
                    ) {
                        Ok(topology) => {
                            let _topology_evidence =
                                (topology.pole_vertex, topology.spoke_edges.len());
                            match validate_exact_boundary_candidate(
                                brep,
                                &topology,
                                &topology.sectors,
                                prepared_inner_spokes,
                                prepared_solved_surfaces,
                                &prepared_pcurves,
                                ExactBoundaryValidationBudget::default(),
                            ) {
                                Ok(validation_report) => {
                                    record_validation_report(
                                        &mut exact_boundary_performance,
                                        validation_report,
                                    );
                                    exact_boundary_sectors = topology.sectors;
                                    exact_boundary_committed_faces = exact_boundary_sectors
                                        .iter()
                                        .map(|sector| sector.cap_face)
                                        .collect();
                                    if let Err(diagnostic) = commit_exact_boundary_sector_surfaces(
                                        brep,
                                        cap_face,
                                        &exact_boundary_sectors,
                                        prepared_solved_surfaces,
                                        &prepared_pcurves,
                                    ) {
                                        exact_boundary_diagnostics
                                            .push(diagnostic.into_apex_diagnostic());
                                    } else {
                                        exact_boundary_committed = true;
                                    }
                                }
                                Err(diagnostic) => {
                                    exact_boundary_diagnostics
                                        .push(diagnostic.into_apex_diagnostic());
                                }
                            }
                        }
                        Err(diagnostic) => {
                            exact_boundary_diagnostics.push(diagnostic.into_apex_diagnostic());
                        }
                    }
                    pair_count
                }
                Err(diagnostic) => {
                    exact_boundary_diagnostics.push(diagnostic.into_apex_diagnostic());
                    0
                }
            }
        }
        _ => 0,
    };
    if !exact_boundary_diagnostics.is_empty() || exact_boundary_pcurve_pair_count == 0 {
        *brep = exact_boundary_transaction_snapshot;
        exact_boundary_committed = false;
        exact_boundary_committed_faces.clear();
    }
    attach_shared_repair_fact_ids(brep, &mut exact_boundary_diagnostics);
    finalize_performance_report(
        &mut exact_boundary_performance,
        exact_boundary_started,
        &exact_boundary_diagnostics,
    );
    if exact_boundary_committed {
        return Ok(GregoryCapAssignment {
            corner_patches: Vec::new(),
            exact_boundary_sectors,
            exact_boundary_committed_faces,
            exact_boundary_prepared_curve_count,
            exact_boundary_inner_spoke_count,
            exact_boundary_ribbon_count,
            exact_boundary_solved_surface_count,
            exact_boundary_pcurve_pair_count,
            exact_boundary_performance,
            exact_boundary_diagnostics,
        });
    }
    let corner_patches = crate::apex_construction::gregory_cap::assemble_gregory_cap_surface(
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
    author_plane_side_face_edge_pcurves(brep, face_edge_handles);
    Ok(GregoryCapAssignment {
        corner_patches,
        exact_boundary_sectors,
        exact_boundary_committed_faces,
        exact_boundary_prepared_curve_count,
        exact_boundary_inner_spoke_count,
        exact_boundary_ribbon_count,
        exact_boundary_solved_surface_count,
        exact_boundary_pcurve_pair_count,
        exact_boundary_performance,
        exact_boundary_diagnostics,
    })
}

pub(super) fn author_plane_side_face_edge_pcurves(
    brep: &mut BRep,
    face_edge_handles: &[EdgeHandle],
) {
    for &edge in face_edge_handles {
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
    face_edge_handles: &[EdgeHandle],
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
        if !face_edges_seen.insert(sector.face_edge) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!("Gregory cap sectors repeat FaceEdge {:?}", sector.face_edge),
            });
        }
    }
    Ok(())
}
