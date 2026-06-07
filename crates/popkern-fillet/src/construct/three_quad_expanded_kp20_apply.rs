//! Expanded KP20 support patchwork application.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::promotion_ladder::PromotionResult;
use super::three_quad_apply::ThreeQuadApplyResult;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerPatchInfo,
    CornerSubPatchProvenance, CornerSubPatchRole,
};
use crate::vertex_classifier::VertexBlendClassification;

pub(super) struct ExpandedKp20ApplyInput<'a> {
    pub brep: &'a mut BRep,
    pub blend_face: FaceHandle,
    pub central_point: Point3,
    pub support_surfaces: &'a [GregoryPatch],
    pub promotion_result: PromotionResult,
    pub vertex: VertexHandle,
    pub n_sides: usize,
    pub classification: VertexBlendClassification,
    pub bridge_fallbacks: &'a [&'static str],
    pub relaxed_coplanarity: Option<f64>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

pub(super) fn expanded_kp20_support_commit_is_accepted(
    promotion_result: PromotionResult,
    support_surfaces: &[GregoryPatch],
    n_sides: usize,
) -> bool {
    matches!(promotion_result, PromotionResult::ExpandedKp20Support)
        && n_sides >= 3
        && support_surfaces.len() == 7 * n_sides
}

pub(super) fn kp20_support_cell_commit_is_accepted(
    promotion_result: PromotionResult,
    support_surfaces: &[GregoryPatch],
    n_sides: usize,
) -> bool {
    matches!(promotion_result, PromotionResult::Kp20SupportCell)
        && n_sides >= 3
        && support_surfaces.len() == 4 * n_sides
}

pub(super) fn try_apply_expanded_kp20_support_patchwork(
    input: ExpandedKp20ApplyInput<'_>,
) -> Result<Option<ThreeQuadApplyResult>, String> {
    if !matches!(input.promotion_result, PromotionResult::ExpandedKp20Support)
        || input.n_sides < 3
        || input.support_surfaces.len() != 7 * input.n_sides
    {
        return Ok(None);
    }

    let commit =
        super::kp20_polar_expanded_support_topology_poke::commit_kp20_expanded_support_patchwork(
            input.brep,
            input.blend_face,
            input.central_point,
            input.support_surfaces,
        )
        .map_err(|error| format!("expanded KP20 support patchwork commit failed: {error}"))?;

    if commit.sub_faces.len() != commit.support_surface_indices.len() {
        return Err(format!(
            "expanded KP20 support commit emitted {} faces for {} surface bindings",
            commit.sub_faces.len(),
            commit.support_surface_indices.len()
        ));
    }

    let patches = commit
        .sub_faces
        .iter()
        .zip(commit.support_surface_indices.iter())
        .map(|(&face, &source_index)| CornerPatchInfo {
            vertex: input.vertex,
            n_sides: input.n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification: input.classification,
            face,
            bridge_fallbacks: input.bridge_fallbacks.to_vec(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: input.relaxed_coplanarity,
            cap_sphere: None,
            sub_patch_provenance: Some(expanded_kp20_provenance(source_index, input.n_sides)),
        })
        .collect();

    Ok(Some(ThreeQuadApplyResult {
        patches,
        diagnostics: input.diagnostics,
    }))
}

pub(super) fn try_apply_kp20_specialized_patchwork(
    input: ExpandedKp20ApplyInput<'_>,
) -> Result<Option<ThreeQuadApplyResult>, String> {
    match input.promotion_result {
        PromotionResult::ExpandedKp20Support
            if expanded_kp20_support_commit_is_accepted(
                input.promotion_result,
                input.support_surfaces,
                input.n_sides,
            ) =>
        {
            try_apply_expanded_kp20_support_patchwork(input)
        }
        PromotionResult::Kp20SupportCell
            if kp20_support_cell_commit_is_accepted(
                input.promotion_result,
                input.support_surfaces,
                input.n_sides,
            ) =>
        {
            try_apply_kp20_support_cell_patchwork(input)
        }
        _ => Ok(None),
    }
}

pub(super) fn try_apply_kp20_support_cell_patchwork(
    input: ExpandedKp20ApplyInput<'_>,
) -> Result<Option<ThreeQuadApplyResult>, String> {
    if !kp20_support_cell_commit_is_accepted(
        input.promotion_result,
        input.support_surfaces,
        input.n_sides,
    ) {
        return Ok(None);
    }

    let commit =
        super::kp20_polar_expanded_support_topology_poke::commit_kp20_support_cell_patchwork(
            input.brep,
            input.blend_face,
            input.central_point,
            input.support_surfaces,
        )
        .map_err(|error| format!("KP20 support-cell commit failed: {error}"))?;

    if commit.support_cell_faces.len() != commit.support_surface_indices.len() {
        return Err(format!(
            "KP20 support-cell commit emitted {} faces for {} surface bindings",
            commit.support_cell_faces.len(),
            commit.support_surface_indices.len()
        ));
    }

    let patches = commit
        .support_cell_faces
        .iter()
        .zip(commit.support_surface_indices.iter())
        .map(|(&face, &source_index)| CornerPatchInfo {
            vertex: input.vertex,
            n_sides: input.n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification: input.classification,
            face,
            bridge_fallbacks: input.bridge_fallbacks.to_vec(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: input.relaxed_coplanarity,
            cap_sphere: None,
            sub_patch_provenance: Some(CornerSubPatchProvenance::kp20_support_cell(
                source_index / 4,
                input.n_sides,
                source_index,
                input.support_surfaces.len(),
            )),
        })
        .collect();

    Ok(Some(ThreeQuadApplyResult {
        patches,
        diagnostics: input.diagnostics,
    }))
}

fn expanded_kp20_provenance(source_index: usize, sector_count: usize) -> CornerSubPatchProvenance {
    let sector_index = source_index / 7;
    let local_index = source_index % 7;
    let patch_count = 7 * sector_count;
    match local_index {
        0..=2 => CornerSubPatchProvenance::kp20_expanded_macro_fragment(
            sector_index,
            sector_count,
            source_index,
            patch_count,
            CornerSubPatchRole::NearOriginLeft,
        ),
        3 => CornerSubPatchProvenance::kp20_expanded_macro_fragment(
            sector_index,
            sector_count,
            source_index,
            patch_count,
            CornerSubPatchRole::NearOriginRight,
        ),
        4 => CornerSubPatchProvenance::kp20_expanded_macro_fragment(
            sector_index,
            sector_count,
            source_index,
            patch_count,
            CornerSubPatchRole::NearMidpointRight,
        ),
        5 => CornerSubPatchProvenance::kp20_expanded_support_lane(
            sector_index,
            sector_count,
            source_index,
            patch_count,
            true,
        ),
        _ => CornerSubPatchProvenance::kp20_expanded_support_lane(
            sector_index,
            sector_count,
            source_index,
            patch_count,
            false,
        ),
    }
}
