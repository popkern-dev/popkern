use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::kp20_polar_expanded_support_binding::{
    Kp20ExpandedSupportSurfaceAssignment, Kp20ExpandedSupportSurfaceSlot,
    attach_kp20_expanded_support_surfaces, author_kp20_expanded_support_pcurves,
    plan_kp20_expanded_support_shared_trim_bindings,
};
#[path = "kp20_polar_support_cell_topology_poke.rs"]
mod support_cell_topology_poke;
pub(crate) use support_cell_topology_poke::commit_kp20_support_cell_patchwork_for_inventory;

#[derive(Debug, Clone)]
pub(crate) struct Kp20ExpandedSupportPokeResult {
    pub(crate) sub_faces: Vec<FaceHandle>,
    pub(crate) shape: Kp20ExpandedSupportPokeShape,
    pub(crate) surface_slots: Vec<Kp20ExpandedSupportSurfaceSlot>,
}

#[derive(Debug, Clone)]
pub(crate) struct Kp20ExpandedSupportCommitResult {
    pub(crate) sub_faces: Vec<FaceHandle>,
    pub(crate) support_surface_indices: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportPokeShape {
    pub(crate) face_count: usize,
    pub(crate) quad_count: usize,
    pub(crate) pentagon_count: usize,
}

impl Kp20ExpandedSupportPokeShape {
    pub(crate) fn requires_trimmed_surface_assignment(&self) -> bool {
        self.pentagon_count > 0
    }
}

pub(crate) fn poke_kp20_expanded_support_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_split_sectors: &[usize],
) -> Result<Kp20ExpandedSupportPokeResult, CoreError> {
    let before = brep.clone();
    match poke_kp20_expanded_support_patchwork_inner(
        brep,
        face,
        central_point,
        support_split_sectors,
    ) {
        Ok(result) => {
            if let Err(error) = brep.validate() {
                *brep = before;
                Err(CoreError::PreconditionFailed(format!(
                    "kp20 expanded support poke produced invalid BRep: {error}"
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

pub(crate) fn commit_kp20_expanded_support_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_surfaces: &[GregoryPatch],
) -> Result<Kp20ExpandedSupportCommitResult, CoreError> {
    let before = brep.clone();
    match commit_kp20_expanded_support_patchwork_inner(brep, face, central_point, support_surfaces)
    {
        Ok(result) => Ok(result),
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn commit_kp20_expanded_support_patchwork_inner(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_surfaces: &[GregoryPatch],
) -> Result<Kp20ExpandedSupportCommitResult, CoreError> {
    if support_surfaces.is_empty() || !support_surfaces.len().is_multiple_of(7) {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support commit expected 7N support surfaces, got {}",
            support_surfaces.len(),
        )));
    }
    let sector_count = support_surfaces.len() / 7;
    let support_split_sectors = (0..sector_count).collect::<Vec<_>>();
    let poke =
        poke_kp20_expanded_support_patchwork(brep, face, central_point, &support_split_sectors)?;
    if !poke.shape.requires_trimmed_surface_assignment() {
        return Err(CoreError::PreconditionFailed(
            "kp20 expanded support commit expected trimmed macro topology slots".into(),
        ));
    }
    let binding_plan = plan_kp20_expanded_support_shared_trim_bindings(
        &poke.surface_slots,
        support_surfaces.len(),
    )?;
    attach_kp20_expanded_support_surfaces(brep, &binding_plan, support_surfaces)?;
    author_kp20_expanded_support_pcurves(brep, &binding_plan)?;
    let support_surface_indices = binding_plan
        .bindings
        .iter()
        .map(|binding| binding.support_surface_index)
        .collect();
    Ok(Kp20ExpandedSupportCommitResult {
        sub_faces: poke.sub_faces,
        support_surface_indices,
    })
}

fn poke_kp20_expanded_support_patchwork_inner(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    support_split_sectors: &[usize],
) -> Result<Kp20ExpandedSupportPokeResult, CoreError> {
    let macro_poke = super::poke_4n::poke_n_gon_face_into_4n_quads(brep, face, central_point)?;
    let n = macro_poke.sub_faces.len() / 4;
    let mut support_lanes = vec![None; n];
    for &sector in support_split_sectors {
        if sector >= n {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 expanded support poke sector {sector} out of range for n={n}"
            )));
        }
        let macro_face = macro_poke.sub_faces[4 * sector + 3];
        support_lanes[sector] = Some(split_support_face_into_lanes(brep, macro_face)?);
    }
    let mut sub_faces =
        Vec::with_capacity(macro_poke.sub_faces.len() + 2 * support_split_sectors.len());
    for (sector, lane) in support_lanes.iter().enumerate().take(n) {
        for piece in 0..4 {
            let macro_face = macro_poke.sub_faces[4 * sector + piece];
            if let (3, Some([boundary, interior])) = (piece, *lane) {
                sub_faces.push(boundary);
                sub_faces.push(interior);
            } else {
                sub_faces.extend(split_fragmented_macro_face_into_quads(brep, macro_face)?);
            }
        }
    }
    let sub_faces = quadrangulate_fragmented_faces(brep, sub_faces)?;
    let shape = expanded_support_poke_shape(brep, &sub_faces)?;
    let surface_slots = expanded_support_surface_slots(brep, &sub_faces)?;
    validate_expanded_support_surface_slots(&shape, &surface_slots)?;
    Ok(Kp20ExpandedSupportPokeResult {
        sub_faces,
        shape,
        surface_slots,
    })
}

fn expanded_support_poke_shape(
    brep: &BRep,
    sub_faces: &[FaceHandle],
) -> Result<Kp20ExpandedSupportPokeShape, CoreError> {
    let mut quad_count = 0;
    let mut pentagon_count = 0;
    for &sub_face in sub_faces {
        match face_outer_half_edges(brep, sub_face)?.len() {
            4 => quad_count += 1,
            5 => pentagon_count += 1,
            len => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 expanded support poke emitted unsupported {len}-edge face"
                )));
            }
        }
    }
    Ok(Kp20ExpandedSupportPokeShape {
        face_count: sub_faces.len(),
        quad_count,
        pentagon_count,
    })
}

fn expanded_support_surface_slots(
    brep: &BRep,
    sub_faces: &[FaceHandle],
) -> Result<Vec<Kp20ExpandedSupportSurfaceSlot>, CoreError> {
    let mut slots = Vec::with_capacity(sub_faces.len());
    for &face in sub_faces {
        let boundary_edge_count = face_outer_half_edges(brep, face)?.len();
        let assignment = match boundary_edge_count {
            4 => Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch,
            5 => Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace,
            len => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 expanded support surface slot cannot classify {len}-edge face"
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

fn validate_expanded_support_surface_slots(
    shape: &Kp20ExpandedSupportPokeShape,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<(), CoreError> {
    let rectangular_slots = slots
        .iter()
        .filter(|slot| {
            slot.assignment == Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch
        })
        .count();
    let trimmed_slots = slots
        .iter()
        .filter(|slot| slot.assignment == Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace)
        .count();
    if slots.len() != shape.face_count
        || rectangular_slots != shape.quad_count
        || trimmed_slots != shape.pentagon_count
    {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support surface slots do not match topology shape: faces={} slots={} quads={} rectangular_slots={} pentagons={} trimmed_slots={}",
            shape.face_count,
            slots.len(),
            shape.quad_count,
            rectangular_slots,
            shape.pentagon_count,
            trimmed_slots,
        )));
    }
    Ok(())
}

fn split_support_face_into_lanes(
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
        len => {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 expanded support poke expected quad or one-pre-split support face, got {len} half-edges",
            )));
        }
    }
    let hes = face_outer_half_edges(brep, face)?;
    if hes.len() != 6 {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support poke expected 6 half-edges after support splits, got {}",
            hes.len()
        )));
    }
    let mef = brep.mef(hes[2], hes[5])?;
    Ok([mef.face, face])
}

fn split_fragmented_macro_face_into_quads(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<Vec<FaceHandle>, CoreError> {
    match face_outer_half_edges(brep, face)?.len() {
        4 => Ok(vec![face]),
        5 => split_pentagon_face_into_quads(brep, face).map(|faces| faces.to_vec()),
        len => Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support poke expected quad or fragmented pentagon macro face, got {len} half-edges",
        ))),
    }
}

fn split_pentagon_face_into_quads(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<[FaceHandle; 2], CoreError> {
    let hes5 = face_outer_half_edges(brep, face)?;
    for &half_edge in &hes5 {
        let mut candidate = brep.clone();
        let split_edge = candidate.half_edges[half_edge].edge;
        if candidate.split_edge(split_edge, 0.5).is_err() {
            continue;
        }
        let Ok(hes6) = face_outer_half_edges(&candidate, face) else {
            continue;
        };
        for first in 0..hes6.len() {
            for second in 0..hes6.len() {
                if first == second {
                    continue;
                }
                let mut trial = candidate.clone();
                let Ok(mef) = trial.mef(hes6[first], hes6[second]) else {
                    continue;
                };
                if face_outer_half_edges(&trial, face).is_ok_and(|hes| hes.len() == 4)
                    && face_outer_half_edges(&trial, mef.face).is_ok_and(|hes| hes.len() == 4)
                    && trial.validate().is_ok()
                {
                    *brep = trial;
                    return Ok([mef.face, face]);
                }
            }
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 expanded support poke could not split fragmented macro pentagon into quads".into(),
    ))
}

fn quadrangulate_fragmented_faces(
    brep: &mut BRep,
    faces: Vec<FaceHandle>,
) -> Result<Vec<FaceHandle>, CoreError> {
    let mut out = Vec::with_capacity(faces.len());
    for face in faces {
        match face_outer_half_edges(brep, face)?.len() {
            4 => out.push(face),
            5 => out.extend(split_pentagon_face_into_quads(brep, face)?),
            len => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 expanded support poke final quadrangulation found {len}-edge face"
                )));
            }
        }
    }
    Ok(out)
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 expanded support poke face invalid".into()))?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 expanded support poke loop invalid".into()))?
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
        "kp20 expanded support poke loop did not close".into(),
    ))
}

#[cfg(test)]
#[path = "kp20_polar_expanded_support_topology_poke_tests.rs"]
mod tests;
