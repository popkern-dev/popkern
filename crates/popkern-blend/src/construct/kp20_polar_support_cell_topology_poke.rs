//! Topology poke contract for KP20 support-cell primitive emission.
#![allow(dead_code)]

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::topology::{BRep, FaceHandle};

use super::Kp20ExpandedSupportPokeShape;
use crate::construct::kp20_polar_expanded_support_binding::{
    Kp20ExpandedSupportSurfaceAssignment, Kp20ExpandedSupportSurfaceSlot,
    attach_kp20_expanded_support_surfaces,
    plan_kp20_expanded_support_surface_bindings_with_sources,
};
use crate::construct::kp20_polar_support_cell_inventory::Kp20SupportCellInventory;

use crate::construct::kp20_polar_support_cell_boundary_reconcile as support_cell_boundary_reconcile;
use crate::construct::kp20_polar_support_cell_pcurve_author as support_cell_pcurve_author;
use crate::construct::kp20_polar_support_cell_topology_pattern as support_cell_topology_pattern;
#[path = "kp20_polar_support_cell_unit_domain.rs"]
mod support_cell_unit_domain;

use crate::construct::kp20_polar_support_cell_boundary_trim::face_outer_half_edges;
use support_cell_topology_pattern::Kp20SupportCellTopologyPattern;
use support_cell_unit_domain::unit_domain_loop_candidates;

#[derive(Debug, Clone)]
pub(crate) struct Kp20SupportCellPokeResult {
    pub(crate) support_cell_faces: Vec<FaceHandle>,
    pub(crate) shape: Kp20ExpandedSupportPokeShape,
    pub(crate) surface_slots: Vec<Kp20ExpandedSupportSurfaceSlot>,
}

#[derive(Debug, Clone)]
pub(crate) struct Kp20SupportCellCommitResult {
    pub(crate) support_cell_faces: Vec<FaceHandle>,
    pub(crate) support_surface_indices: Vec<usize>,
}

pub(crate) fn commit_kp20_support_cell_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    if sector_count == 0 || !support_surfaces.len().is_multiple_of(sector_count) {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell commit expected a support surface inventory divisible by sector count {sector_count}, got {}",
            support_surfaces.len(),
        )));
    }
    let support_split_sectors = (0..sector_count).collect::<Vec<_>>();
    let cells_per_target_sector = support_surfaces.len() / sector_count;
    commit_kp20_support_cell_patchwork_checked(
        brep,
        face,
        central_point,
        sector_count,
        support_surfaces,
        &support_split_sectors,
        cells_per_target_sector,
    )
}

pub(crate) fn commit_kp20_support_cell_patchwork_for_inventory(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
    inventory: &Kp20SupportCellInventory,
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    let expected_surfaces = inventory.target_sectors.len() * inventory.cells_per_target_sector;
    if inventory.sector_count != sector_count || support_surfaces.len() != expected_surfaces {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell inventory mismatch: ring={sector_count}, inventory={}, expected_surfaces={expected_surfaces}, got={}",
            inventory.sector_count,
            support_surfaces.len(),
        )));
    }
    commit_kp20_support_cell_patchwork_checked(
        brep,
        face,
        central_point,
        sector_count,
        support_surfaces,
        &inventory.target_sectors,
        inventory.cells_per_target_sector,
    )
}

fn commit_kp20_support_cell_patchwork_checked(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
    support_split_sectors: &[usize],
    cells_per_target_sector: usize,
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    let before = brep.clone();
    match commit_kp20_support_cell_patchwork_inner(
        brep,
        face,
        central_point,
        sector_count,
        support_surfaces,
        support_split_sectors,
        cells_per_target_sector,
    ) {
        Ok(result) => Ok(result),
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn commit_kp20_support_cell_patchwork_inner(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
    support_split_sectors: &[usize],
    cells_per_target_sector: usize,
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    if sector_count == 0 || support_split_sectors.is_empty() {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell commit expected target sectors within ring sector count {sector_count}, got {} targets for {} surfaces",
            support_split_sectors.len(),
            support_surfaces.len(),
        )));
    }
    let Some(pattern) =
        Kp20SupportCellTopologyPattern::for_cells_per_sector(cells_per_target_sector)
    else {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell commit does not support {cells_per_target_sector} cells per target sector",
        )));
    };
    pattern.validate_support_surface_inventory(support_surfaces.len())?;
    let poke = poke_kp20_support_cell_patchwork_with_pattern(
        brep,
        face,
        central_point,
        sector_count,
        support_split_sectors,
        pattern,
    )?;
    let source_map =
        match_support_cell_surfaces_to_slots(brep, pattern, &poke.surface_slots, support_surfaces)?;
    let binding_plan = plan_kp20_expanded_support_surface_bindings_with_sources(
        &poke.surface_slots,
        support_surfaces.len(),
        &source_map,
    )?;
    attach_kp20_expanded_support_surfaces(brep, &binding_plan, support_surfaces)?;
    crate::construct::post_pipeline::determinism_trace("post_attach_surfaces", brep);
    support_cell_pcurve_author::author_optimized_support_cell_pcurves(
        brep,
        &poke.surface_slots,
        support_surfaces,
        &source_map,
    )?;
    crate::construct::post_pipeline::determinism_trace("post_pcurve_author", brep);
    support_cell_boundary_reconcile::reconcile_support_cell_full_boundaries(
        brep,
        &poke.surface_slots,
    )?;
    crate::construct::post_pipeline::determinism_trace("post_commit_reconcile", brep);
    let support_surface_indices = binding_plan
        .bindings
        .iter()
        .map(|binding| binding.support_surface_index)
        .collect();
    Ok(Kp20SupportCellCommitResult {
        support_cell_faces: poke.support_cell_faces,
        support_surface_indices,
    })
}

fn match_support_cell_surfaces_to_slots(
    brep: &BRep,
    pattern: Kp20SupportCellTopologyPattern,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surfaces: &[GregoryPatch],
) -> Result<Vec<usize>, CoreError> {
    let mut used = vec![false; support_surfaces.len()];
    let mut source_map = Vec::with_capacity(slots.len());
    for (slot_index, slot) in slots.iter().enumerate() {
        let mut best = None;
        let source_range = pattern.source_range_for_slot(slot_index, support_surfaces.len());
        for (source_index, surface) in support_surfaces.iter().enumerate() {
            if used[source_index] || !source_range.contains(&source_index) {
                continue;
            }
            let residual = support_surface_slot_residual(brep, slot.face, surface)?;
            if best.is_none_or(|(_, best_residual)| residual < best_residual) {
                best = Some((source_index, residual));
            }
        }
        let Some((source_index, _)) = best else {
            return Err(CoreError::PreconditionFailed(
                "kp20 support-cell binding ran out of support surfaces".into(),
            ));
        };
        used[source_index] = true;
        source_map.push(source_index);
    }
    Ok(source_map)
}

fn support_surface_slot_residual(
    brep: &BRep,
    face: FaceHandle,
    surface: &GregoryPatch,
) -> Result<f64, CoreError> {
    let half_edges = face_outer_half_edges(brep, face)?;
    let candidates = unit_domain_loop_candidates(half_edges.len());
    if candidates.is_empty() {
        return Ok(f64::INFINITY);
    }
    Ok(candidates
        .iter()
        .map(|candidate| {
            half_edges
                .iter()
                .zip(candidate.iter())
                .map(|(&half_edge, &(u, v))| {
                    let position = brep.vertices[brep.half_edges[half_edge].origin].position;
                    position.distance_to(&surface.point_at(u, v))
                })
                .sum::<f64>()
        })
        .fold(f64::INFINITY, f64::min))
}

pub(crate) fn poke_kp20_support_cell_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_split_sectors: &[usize],
) -> Result<Kp20SupportCellPokeResult, CoreError> {
    let pattern = Kp20SupportCellTopologyPattern::current_four_cell_support_grid();
    poke_kp20_support_cell_patchwork_with_pattern(
        brep,
        face,
        central_point,
        support_split_sectors.len(),
        support_split_sectors,
        pattern,
    )
}

fn poke_kp20_support_cell_patchwork_with_pattern(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_split_sectors: &[usize],
    pattern: Kp20SupportCellTopologyPattern,
) -> Result<Kp20SupportCellPokeResult, CoreError> {
    let before = brep.clone();
    match poke_kp20_support_cell_patchwork_inner(
        brep,
        face,
        central_point,
        support_split_sectors,
        sector_count,
        pattern,
    ) {
        Ok(result) => {
            if let Err(error) = brep.validate() {
                *brep = before;
                Err(CoreError::PreconditionFailed(format!(
                    "kp20 support-cell poke produced invalid BRep: {error}"
                )))
            } else {
                Ok(result)
            }
        }
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn poke_kp20_support_cell_patchwork_inner(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_split_sectors: &[usize],
    sector_count: usize,
    pattern: Kp20SupportCellTopologyPattern,
) -> Result<Kp20SupportCellPokeResult, CoreError> {
    let macro_poke = crate::construct::poke_4n::poke_n_gon_face_into_expected_4n_quads(
        brep,
        face,
        central_point,
        sector_count,
    )?;
    let n = macro_poke.sub_faces.len() / 4;
    let mut support_cell_faces = Vec::with_capacity(4 * support_split_sectors.len());
    for &sector in support_split_sectors {
        if sector >= n {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 support-cell poke sector {sector} out of range for n={n}"
            )));
        }
        let support_face = macro_poke.sub_faces[4 * sector + 3];
        if pattern == Kp20SupportCellTopologyPattern::FourCellSupportGrid {
            support_cell_faces.extend(pattern.split_sector_support_face(brep, support_face)?);
        } else {
            support_cell_faces.push(support_face);
        }
    }
    if pattern != Kp20SupportCellTopologyPattern::FourCellSupportGrid {
        support_cell_faces = pattern.split_sector_support_faces(brep, &support_cell_faces)?;
    }
    let shape = support_cell_poke_shape(brep, &support_cell_faces)?;
    let surface_slots = support_cell_surface_slots(brep, &support_cell_faces)?;
    validate_support_cell_surface_slots(&shape, &surface_slots)?;
    Ok(Kp20SupportCellPokeResult {
        support_cell_faces,
        shape,
        surface_slots,
    })
}

fn support_cell_poke_shape(
    brep: &BRep,
    support_cell_faces: &[FaceHandle],
) -> Result<Kp20ExpandedSupportPokeShape, CoreError> {
    let mut quad_count = 0;
    let mut pentagon_count = 0;
    for &face in support_cell_faces {
        match face_outer_half_edges(brep, face)?.len() {
            4 => quad_count += 1,
            5.. => pentagon_count += 1,
            len => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 support-cell poke emitted unsupported {len}-edge support cell"
                )));
            }
        }
    }
    Ok(Kp20ExpandedSupportPokeShape {
        face_count: support_cell_faces.len(),
        quad_count,
        pentagon_count,
    })
}

fn support_cell_surface_slots(
    brep: &BRep,
    support_cell_faces: &[FaceHandle],
) -> Result<Vec<Kp20ExpandedSupportSurfaceSlot>, CoreError> {
    let mut slots = Vec::with_capacity(support_cell_faces.len());
    for &face in support_cell_faces {
        let boundary_edge_count = face_outer_half_edges(brep, face)?.len();
        let assignment = match boundary_edge_count {
            4 => Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch,
            5.. => Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace,
            _ => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 support-cell slot expected at least 4 support-cell edges, got {boundary_edge_count}"
                )));
            }
        };
        slots.push(Kp20ExpandedSupportSurfaceSlot {
            face,
            boundary_edge_count,
            assignment,
        });
    }
    Ok(slots)
}

fn validate_support_cell_surface_slots(
    shape: &Kp20ExpandedSupportPokeShape,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<(), CoreError> {
    if slots.len() != shape.face_count
        || slots
            .iter()
            .filter(|slot| {
                slot.assignment == Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch
            })
            .count()
            != shape.quad_count
        || slots
            .iter()
            .filter(|slot| {
                slot.assignment == Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
            })
            .count()
            != shape.pentagon_count
    {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell slots do not match rectangular topology: faces={} slots={} pentagons={}",
            shape.face_count,
            slots.len(),
            shape.pentagon_count,
        )));
    }
    Ok(())
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_topology_poke_tests.rs"]
mod tests;
