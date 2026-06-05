//! KP20 polar comparison candidate assembly helpers.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    Kp20SectorPieceMaps, TensorBorderBi3, assemble_kp20_bi6_macro_cap, build_kp20_guide_ring,
    sample_kp20_macro_patch_from_piece_maps, sample_kp20_macro_patch_from_tau_maps,
};

use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_sampled,
};
use super::super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_macro_patches,
};

pub(super) fn mapped_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    maps: Kp20SectorPieceMaps,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let sampled_macro_patches = sampled_macro_patches_from_maps(star_point, tensor_border, maps)?;
    candidate_summary_from_sampled(tensor_border, sampled_macro_patches)
}

pub(super) fn lifted_mapped_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    maps: Kp20SectorPieceMaps,
) -> Result<Kp20SupportNormalLiftReport, String> {
    let sampled_macro_patches = sampled_macro_patches_from_maps(star_point, tensor_border, maps)?;
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    lift_candidate_from_macro_patches(star_point, tensor_border, cap.macro_patches)
}

pub(super) fn tau_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .enumerate()
        .map(|(sector, guide_sector)| {
            sample_kp20_macro_patch_from_tau_maps(guide_sector, guide.n, sector)
        })
        .collect::<Result<Vec<_>, _>>()?;
    candidate_summary_from_sampled(tensor_border, sampled_macro_patches)
}

fn sampled_macro_patches_from_maps(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    maps: Kp20SectorPieceMaps,
) -> Result<Vec<Kp20MacroPatch>, String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    Ok(guide
        .sectors
        .iter()
        .map(|sector| sample_kp20_macro_patch_from_piece_maps(sector, maps))
        .collect())
}
