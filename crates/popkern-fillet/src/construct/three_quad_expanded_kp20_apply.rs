//! Expanded KP20 support patchwork application.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::kp20_polar_support_cell_inventory::{
    Kp20SupportCellInventory, support_cell_inventory_from_provenance,
};
use super::promotion_ladder::PromotionResult;
use super::three_quad_apply::ThreeQuadApplyResult;
use crate::types::{
    ApexConstructionDiagnostic, ContinuityLevel, CornerMethod, CornerPatchInfo,
    CornerSubPatchProvenance, CornerSubPatchRole,
};
use crate::vertex_classifier::VertexBlendClassification;
use support_cell_transaction::{
    Kp20SupportCellApplyTransactionResult, Kp20SupportCellValidationFailure,
    run_kp20_support_cell_apply_transaction, support_cell_assembled_quality_fact,
};

#[path = "three_quad_kp20_support_cell_transaction.rs"]
mod support_cell_transaction;

pub(super) struct ExpandedKp20ApplyInput<'a> {
    pub brep: &'a mut BRep,
    pub blend_face: FaceHandle,
    pub central_point: Point3,
    pub support_surfaces: &'a [GregoryPatch],
    pub support_cell_floor_surfaces: Option<&'a [GregoryPatch]>,
    pub support_cell_provenance: Option<&'a [CornerSubPatchProvenance]>,
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
    provenance: Option<&[CornerSubPatchProvenance]>,
) -> bool {
    matches!(promotion_result, PromotionResult::Kp20SupportCell)
        && n_sides >= 3
        && kp20_support_cell_inventory(support_surfaces, n_sides, provenance).is_some()
}

fn kp20_support_cell_inventory(
    support_surfaces: &[GregoryPatch],
    n_sides: usize,
    provenance: Option<&[CornerSubPatchProvenance]>,
) -> Option<Kp20SupportCellInventory> {
    if let Some(provenance) = provenance {
        return support_cell_inventory_from_provenance(n_sides, support_surfaces.len(), provenance)
            .ok();
    }
    Kp20SupportCellInventory::full_ring(n_sides, support_surfaces.len())
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
                input.support_cell_provenance,
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
        input.support_cell_provenance,
    ) {
        return Ok(None);
    }

    let promoted_inventory = kp20_support_cell_inventory(
        input.support_surfaces,
        input.n_sides,
        input.support_cell_provenance,
    )
    .ok_or_else(|| {
        format!(
            "KP20 support-cell commit received unsupported {} surfaces for {} sectors",
            input.support_surfaces.len(),
            input.n_sides
        )
    })?;
    let promoted_cells_per_target_sector = promoted_inventory.cells_per_target_sector;
    let promoted_surfaces = input.support_surfaces;
    let promoted_diagnostics = input.diagnostics.clone();
    let fallback = input.support_cell_floor_surfaces.and_then(|surfaces| {
        kp20_support_cell_inventory(surfaces, input.n_sides, None)
            .map(|inventory| (surfaces, inventory))
    });
    let fallback_surface_count = fallback.as_ref().map(|(surfaces, _)| surfaces.len());
    let fallback_cells_per_target_sector = fallback
        .as_ref()
        .map(|(_, inventory)| inventory.cells_per_target_sector);
    let vertex = input.vertex;
    let blend_face = input.blend_face;
    let central_point = input.central_point;
    let n_sides = input.n_sides;
    let classification = input.classification;
    let bridge_fallbacks = input.bridge_fallbacks;
    let relaxed_coplanarity = input.relaxed_coplanarity;
    let transaction = run_kp20_support_cell_apply_transaction(
        input.brep,
        input.n_sides,
        promoted_surfaces.len(),
        promoted_cells_per_target_sector,
        fallback_surface_count,
        |brep| {
            apply_kp20_support_cell_surfaces(
                brep,
                Kp20SupportCellSurfaceApply {
                    blend_face,
                    central_point,
                    support_surfaces: promoted_surfaces,
                    support_cell_provenance: input.support_cell_provenance,
                    inventory: promoted_inventory.clone(),
                    vertex,
                    n_sides,
                    classification,
                    bridge_fallbacks,
                    relaxed_coplanarity,
                    diagnostics: promoted_diagnostics,
                },
            )
        },
        |brep| {
            brep.validate()
                .map_err(Kp20SupportCellValidationFailure::from)
        },
        |brep, diagnostic| {
            let Some((floor_surfaces, floor_inventory)) = fallback else {
                return Err("KP20 support-cell fallback surfaces missing".to_string());
            };
            let mut fallback_diagnostics = input.diagnostics.clone();
            fallback_diagnostics.push(diagnostic);
            apply_kp20_support_cell_surfaces(
                brep,
                Kp20SupportCellSurfaceApply {
                    blend_face,
                    central_point,
                    support_surfaces: floor_surfaces,
                    support_cell_provenance: None,
                    inventory: floor_inventory,
                    vertex,
                    n_sides,
                    classification,
                    bridge_fallbacks,
                    relaxed_coplanarity,
                    diagnostics: fallback_diagnostics,
                },
            )
        },
    )?;
    match transaction {
        Kp20SupportCellApplyTransactionResult::Promoted(mut result) => {
            result.diagnostics.push(support_cell_assembled_quality_fact(
                n_sides,
                "promoted",
                promoted_surfaces.len(),
                promoted_cells_per_target_sector,
                promoted_surfaces.len(),
                fallback_surface_count,
                "valid",
                None,
            ));
            Ok(Some(result))
        }
        Kp20SupportCellApplyTransactionResult::Fallback(mut result) => {
            let fallback_surface_count = fallback_surface_count.unwrap_or(result.patches.len());
            result.diagnostics.push(support_cell_assembled_quality_fact(
                n_sides,
                "fallback_four_cell_floor",
                fallback_surface_count,
                fallback_cells_per_target_sector.unwrap_or(fallback_surface_count / n_sides),
                promoted_surfaces.len(),
                Some(fallback_surface_count),
                "valid_after_promoted_decline",
                None,
            ));
            Ok(Some(result))
        }
    }
}

struct Kp20SupportCellSurfaceApply<'a> {
    blend_face: FaceHandle,
    central_point: Point3,
    support_surfaces: &'a [GregoryPatch],
    support_cell_provenance: Option<&'a [CornerSubPatchProvenance]>,
    inventory: Kp20SupportCellInventory,
    vertex: VertexHandle,
    n_sides: usize,
    classification: VertexBlendClassification,
    bridge_fallbacks: &'a [&'static str],
    relaxed_coplanarity: Option<f64>,
    diagnostics: Vec<ApexConstructionDiagnostic>,
}

fn apply_kp20_support_cell_surfaces(
    brep: &mut BRep,
    apply: Kp20SupportCellSurfaceApply<'_>,
) -> Result<ThreeQuadApplyResult, String> {
    let commit =
        super::kp20_polar_expanded_support_topology_poke::commit_kp20_support_cell_patchwork_for_inventory(
            brep,
            apply.blend_face,
            apply.central_point,
            apply.n_sides,
            apply.support_surfaces,
            &apply.inventory,
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
            vertex: apply.vertex,
            n_sides: apply.n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification: apply.classification,
            face,
            bridge_fallbacks: apply.bridge_fallbacks.to_vec(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: apply.relaxed_coplanarity,
            cap_sphere: None,
            sub_patch_provenance: Some(
                apply
                    .support_cell_provenance
                    .and_then(|items| items.get(source_index))
                    .copied()
                    .unwrap_or_else(|| {
                        CornerSubPatchProvenance::kp20_support_cell(
                            source_index / apply.inventory.cells_per_target_sector,
                            apply.n_sides,
                            source_index,
                            apply.support_surfaces.len(),
                        )
                    }),
            ),
        })
        .collect();
    Ok(ThreeQuadApplyResult {
        patches,
        diagnostics: apply.diagnostics,
    })
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
