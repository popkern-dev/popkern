use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::face_outer_half_edges;

const MAX_PRESPLIT_SUPPORT_FACE_BOUNDARY_EDGES: usize = 32;

#[allow(clippy::enum_variant_names)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kp20SupportCellTopologyPattern {
    FourCellSupportGrid,
    EightCellSupportGrid,
    SixteenCellSupportGrid,
}

impl Kp20SupportCellTopologyPattern {
    pub(crate) fn current_four_cell_support_grid() -> Self {
        Self::FourCellSupportGrid
    }

    pub(crate) fn eight_cell_support_grid() -> Self {
        Self::EightCellSupportGrid
    }

    pub(crate) fn for_cells_per_sector(cells_per_sector: usize) -> Option<Self> {
        match cells_per_sector {
            4 => Some(Self::FourCellSupportGrid),
            8 => Some(Self::EightCellSupportGrid),
            16 => Some(Self::SixteenCellSupportGrid),
            _ => None,
        }
    }

    pub(crate) fn cells_per_sector(self) -> usize {
        match self {
            Self::FourCellSupportGrid => 4,
            Self::EightCellSupportGrid => 8,
            Self::SixteenCellSupportGrid => 16,
        }
    }

    pub(crate) fn validate_support_surface_inventory(
        self,
        support_surface_count: usize,
    ) -> Result<usize, CoreError> {
        let cells_per_sector = self.cells_per_sector();
        if support_surface_count == 0 || !support_surface_count.is_multiple_of(cells_per_sector) {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 support-cell commit expected {cells_per_sector}N support surfaces, got {support_surface_count}",
            )));
        }
        Ok(support_surface_count / cells_per_sector)
    }

    pub(crate) fn source_range_for_slot(
        self,
        slot_index: usize,
        support_surface_count: usize,
    ) -> std::ops::Range<usize> {
        let cells_per_sector = self.cells_per_sector();
        let sector_start = cells_per_sector * (slot_index / cells_per_sector);
        let sector_end = (sector_start + cells_per_sector).min(support_surface_count);
        sector_start..sector_end
    }

    pub(crate) fn split_sector_support_face(
        self,
        brep: &mut BRep,
        face: FaceHandle,
    ) -> Result<Vec<FaceHandle>, CoreError> {
        self.split_sector_support_faces(brep, &[face])
    }

    pub(crate) fn split_sector_support_faces(
        self,
        brep: &mut BRep,
        faces: &[FaceHandle],
    ) -> Result<Vec<FaceHandle>, CoreError> {
        match self {
            Self::FourCellSupportGrid => split_support_faces_into_binary_grid(brep, faces, 2),
            Self::EightCellSupportGrid => split_support_faces_into_binary_grid(brep, faces, 3),
            Self::SixteenCellSupportGrid => split_support_faces_into_binary_grid(brep, faces, 4),
        }
    }
}

fn split_support_faces_into_binary_grid(
    brep: &mut BRep,
    faces: &[FaceHandle],
    split_levels: usize,
) -> Result<Vec<FaceHandle>, CoreError> {
    let mut cells = faces.to_vec();
    for _ in 0..split_levels {
        let mut next = Vec::with_capacity(cells.len() * 2);
        for cell in cells {
            next.extend(split_quad_like_face_into_two_lanes(brep, cell)?);
        }
        cells = next;
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
        6..=MAX_PRESPLIT_SUPPORT_FACE_BOUNDARY_EDGES => {
            return split_presplit_face_into_quad_like_pair(brep, face);
        }
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
    let mut best = None;
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
            if !presplit_face_split_is_quad_like(hes.len(), left.len(), right.len())
                || trial.validate().is_err()
            {
                continue;
            }
            let score = PresplitFaceSplitScore::from_loop_lengths(left.len(), right.len());
            if best
                .as_ref()
                .is_none_or(|(best_score, _, _)| score < *best_score)
            {
                best = Some((score, trial, mef.face));
            }
        }
    }
    if let Some((_, trial, new_face)) = best {
        *brep = trial;
        return Ok([new_face, face]);
    }
    Err(CoreError::PreconditionFailed(format!(
        "kp20 support-cell poke could not split pre-split {}-edge support face into quad-like lanes",
        hes.len(),
    )))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PresplitFaceSplitScore {
    max_boundary_edges: usize,
    imbalance: usize,
}

impl PresplitFaceSplitScore {
    fn from_loop_lengths(left: usize, right: usize) -> Self {
        Self {
            max_boundary_edges: left.max(right),
            imbalance: left.abs_diff(right),
        }
    }
}

fn presplit_face_split_is_quad_like(parent: usize, left: usize, right: usize) -> bool {
    left >= 4 && right >= 4 && left < parent && right < parent && left + right == parent + 2
}
