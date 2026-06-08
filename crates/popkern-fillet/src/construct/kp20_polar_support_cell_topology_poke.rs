//! Topology poke contract for KP20 support-cell primitive emission.
#![allow(dead_code)]

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::Kp20ExpandedSupportPokeShape;
use crate::construct::kp20_polar_expanded_support_binding::{
    Kp20ExpandedSupportSurfaceAssignment, Kp20ExpandedSupportSurfaceSlot,
    attach_kp20_expanded_support_surfaces,
    plan_kp20_expanded_support_surface_bindings_with_sources,
};

#[path = "kp20_polar_support_cell_boundary_reconcile.rs"]
mod support_cell_boundary_reconcile;
#[path = "kp20_polar_support_cell_pcurve_author.rs"]
mod support_cell_pcurve_author;

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
    support_surfaces: &[popkern_core::geometry::gregory::GregoryPatch],
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    let before = brep.clone();
    match commit_kp20_support_cell_patchwork_inner(brep, face, central_point, support_surfaces) {
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
    support_surfaces: &[popkern_core::geometry::gregory::GregoryPatch],
) -> Result<Kp20SupportCellCommitResult, CoreError> {
    if support_surfaces.is_empty() || support_surfaces.len() % 4 != 0 {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell commit expected 4N support surfaces, got {}",
            support_surfaces.len(),
        )));
    }
    let sector_count = support_surfaces.len() / 4;
    let support_split_sectors = (0..sector_count).collect::<Vec<_>>();
    let poke = poke_kp20_support_cell_patchwork(brep, face, central_point, &support_split_sectors)?;
    let source_map =
        match_support_cell_surfaces_to_slots(brep, &poke.surface_slots, support_surfaces)?;
    let binding_plan = plan_kp20_expanded_support_surface_bindings_with_sources(
        &poke.surface_slots,
        support_surfaces.len(),
        &source_map,
    )?;
    attach_kp20_expanded_support_surfaces(brep, &binding_plan, support_surfaces)?;
    support_cell_pcurve_author::author_optimized_support_cell_pcurves(
        brep,
        &poke.surface_slots,
        support_surfaces,
        &source_map,
    )?;
    support_cell_boundary_reconcile::reconcile_support_cell_full_boundaries(
        brep,
        &poke.surface_slots,
    )?;
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
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surfaces: &[popkern_core::geometry::gregory::GregoryPatch],
) -> Result<Vec<usize>, CoreError> {
    let mut used = vec![false; support_surfaces.len()];
    let mut source_map = Vec::with_capacity(slots.len());
    for (slot_index, slot) in slots.iter().enumerate() {
        let mut best = None;
        let sector_start = 4 * (slot_index / 4);
        let sector_end = (sector_start + 4).min(support_surfaces.len());
        for (source_index, surface) in support_surfaces.iter().enumerate() {
            if used[source_index] || source_index < sector_start || source_index >= sector_end {
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
    surface: &popkern_core::geometry::gregory::GregoryPatch,
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

fn unit_domain_loop_candidates(edge_count: usize) -> Vec<Vec<(f64, f64)>> {
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let bases = if edge_count == 4 {
        vec![corners.to_vec()]
    } else if edge_count == 5 {
        (0..4)
            .map(|split| {
                let mut out = Vec::with_capacity(5);
                for index in 0..4 {
                    out.push(corners[index]);
                    if index == split {
                        let next = corners[(index + 1) % 4];
                        out.push((
                            (corners[index].0 + next.0) * 0.5,
                            (corners[index].1 + next.1) * 0.5,
                        ));
                    }
                }
                out
            })
            .collect()
    } else {
        Vec::new()
    };
    let mut candidates = Vec::new();
    for base in bases {
        for shift in 0..base.len() {
            let rotated = (0..base.len())
                .map(|index| base[(index + shift) % base.len()])
                .collect::<Vec<_>>();
            candidates.push(rotated.clone());
            candidates.push(rotated.into_iter().rev().collect());
        }
    }
    candidates
}

pub(crate) fn poke_kp20_support_cell_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_split_sectors: &[usize],
) -> Result<Kp20SupportCellPokeResult, CoreError> {
    let before = brep.clone();
    match poke_kp20_support_cell_patchwork_inner(brep, face, central_point, support_split_sectors) {
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
) -> Result<Kp20SupportCellPokeResult, CoreError> {
    let macro_poke = crate::construct::poke_4n::poke_n_gon_face_into_expected_4n_quads(
        brep,
        face,
        central_point,
        support_split_sectors.len(),
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
        support_cell_faces.extend(split_support_face_into_four_cells(brep, support_face)?);
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

fn split_support_face_into_four_cells(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<Vec<FaceHandle>, CoreError> {
    let lanes = split_quad_like_face_into_two_lanes(brep, face)?;
    let mut cells = Vec::with_capacity(4);
    for lane in lanes {
        cells.extend(split_quad_like_face_into_two_lanes(brep, lane)?);
    }
    Ok(cells)
}

fn split_quad_like_face_into_two_lanes(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<[FaceHandle; 2], CoreError> {
    let hes = face_outer_half_edges(brep, face)?;
    match hes.len() {
        4 => {
            let first_edge = brep.half_edges[hes[1]].edge;
            let opposite_edge = brep.half_edges[hes[3]].edge;
            brep.split_edge(first_edge, 0.5)?;
            brep.split_edge(opposite_edge, 0.5)?;
        }
        5 => {
            let opposite_edge = brep.half_edges[hes[3]].edge;
            brep.split_edge(opposite_edge, 0.5)?;
        }
        6..=8 => return split_presplit_face_into_quad_like_pair(brep, face),
        len => {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 support-cell poke expected quad or one-pre-split face, got {len} half-edges",
            )));
        }
    }
    let hes = face_outer_half_edges(brep, face)?;
    if hes.len() != 6 {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell poke expected 6 half-edges after splits, got {}",
            hes.len()
        )));
    }
    let mef = brep.mef(hes[2], hes[5])?;
    Ok([mef.face, face])
}

fn split_presplit_face_into_quad_like_pair(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<[FaceHandle; 2], CoreError> {
    let hes = face_outer_half_edges(brep, face)?;
    for first in 0..hes.len() {
        for second in 0..hes.len() {
            if first == second {
                continue;
            }
            let mut trial = brep.clone();
            let Ok(mef) = trial.mef(hes[first], hes[second]) else {
                continue;
            };
            let Ok(left) = face_outer_half_edges(&trial, face) else {
                continue;
            };
            let Ok(right) = face_outer_half_edges(&trial, mef.face) else {
                continue;
            };
            if (4..=5).contains(&left.len())
                && (4..=5).contains(&right.len())
                && trial.validate().is_ok()
            {
                *brep = trial;
                return Ok([mef.face, face]);
            }
        }
    }
    Err(CoreError::PreconditionFailed(format!(
        "kp20 support-cell poke could not split pre-split {}-edge support face into quad-like lanes",
        hes.len(),
    )))
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
            5 => pentagon_count += 1,
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
            5 => Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace,
            _ => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 support-cell slot expected 4- or 5-edge support cell, got {boundary_edge_count}"
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

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 support-cell poke face invalid".into()))?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 support-cell poke loop invalid".into()))?
        .half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 support-cell poke loop did not close".into(),
    ))
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_topology_poke_tests.rs"]
mod tests;
