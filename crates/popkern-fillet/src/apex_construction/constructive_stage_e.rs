//! Stage E surface assignment for constructive apex building.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::CornerPatchInfo;

use super::constructive_geometry::planned_fillet_index_for_sector;
use super::error::ApexBuildError;
use super::gregory_cap::GregoryCapSector;
use super::plan::ApexPlan;

/// Stage E (Sphere): assign `SurfaceGeometry::Sphere` to the cap face.
pub(super) fn stage_e_assign_sphere(
    brep: &mut BRep,
    cap_face: FaceHandle,
    center: Point3,
    radius: f64,
) -> Result<(), ApexBuildError> {
    let face = brep
        .faces
        .get_mut(cap_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("stage_e_assign_sphere: cap_face {cap_face:?} no longer exists"),
        })?;
    face.surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere { center, radius },
        orientation: true,
    });
    Ok(())
}

/// Stage E (Gregory): delegate to `assemble_gregory_cap_surface`.
#[allow(clippy::too_many_arguments)]
pub(super) fn stage_e_assign_gregory(
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
) -> Result<Vec<CornerPatchInfo>, ApexBuildError> {
    let n = plan.cap_kind.n_fillets();
    if face_edge_handles.len() != n || fillet_end_handles.len() != n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "stage_e_assign_gregory: expected N={n} face edges + N fillet ends, \
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

    crate::apex_construction::gregory_cap::assemble_gregory_cap_surface(
        brep,
        cap_face,
        plan,
        apex_position,
        apex_vertex_handles,
        &apex_edge_handles,
        sectors,
        fillet_surfaces,
        tolerance,
    )
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
