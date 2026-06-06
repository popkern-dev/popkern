//! Primary KP20 polar macro candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6CapInput, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport,
    Kp20GuideBuildInput, Kp20PieceMapDomainReport, Kp20PieceMapStitchReport,
    RectangularGregoryBiNonic, TensorBorderBi3, assemble_kp20_bi6_cap, assemble_kp20_bi6_macro_cap,
    build_kp20_guide_ring, elevate_bi6_grid_to_bi9, kp20_piece_map_domain_report,
    kp20_piece_map_stitch_report, kp20_positive_nonoverlap_c2_sector_piece_maps,
    measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_piece_maps,
};

use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, gregory_patches_from_macro_patches,
    validation_is_clean,
};
use super::labels::{piece_counts_label, validation_label};

pub(super) struct Kp20PrimaryMappedCandidate {
    pub(super) patches: Vec<GregoryPatch>,
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    domain: Kp20PieceMapDomainReport,
    stitch: Kp20PieceMapStitchReport,
}

impl Kp20PrimaryMappedCandidate {
    pub(super) fn is_acceptable(&self) -> bool {
        validation_is_clean(self.validation) && self.summary.fold_count == 0
    }

    pub(super) fn detail(&self) -> String {
        format!(
            "map_family=positive_nonoverlap_c2;status=sampled;folds={};piece_folds={};validation:{};inside_triangle={};non_positive_orientation_maps={};c2_edge={:.3e};origin_strip_edge={:.3e};origin_strip_reversed={};min_center_sep={:.3e}",
            self.summary.fold_count,
            piece_counts_label(self.summary.piece_fold_counts),
            validation_label(self.validation),
            self.domain.stays_inside_triangle(),
            self.domain.non_positive_orientation_maps,
            self.stitch.max_internal_c2_edge_mismatch,
            self.stitch.max_origin_strip_edge_mismatch,
            self.stitch.origin_strip_reversed,
            self.stitch.min_piece_center_separation,
        )
    }
}

pub(super) fn try_positive_nonoverlap_c2_macro_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20PrimaryMappedCandidate, String> {
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    let domain = kp20_piece_map_domain_report(maps);
    if !domain.stays_inside_triangle() || domain.non_positive_orientation_maps > 0 {
        return Err(format!(
            "positive_nonoverlap_c2 map invalid: inside_triangle={} non_positive_orientation_maps={}",
            domain.stays_inside_triangle(),
            domain.non_positive_orientation_maps,
        ));
    }
    let stitch = kp20_piece_map_stitch_report(kp20_positive_nonoverlap_c2_sector_piece_maps());
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_positive_nonoverlap_c2_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let (summary, validation) =
        candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
    Ok(Kp20PrimaryMappedCandidate {
        patches: gregory_patches_from_macro_patches(&cap.macro_patches),
        summary,
        validation,
        domain,
        stitch,
    })
}

pub(crate) fn try_kp20_polar_macro_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Vec<GregoryPatch>, String> {
    let cap = assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let report = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    if !report.passes(1e-8) {
        return Err(format!("validation residuals exceed tolerance: {report:?}"));
    }

    Ok(cap
        .patches_bi6
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect())
}
