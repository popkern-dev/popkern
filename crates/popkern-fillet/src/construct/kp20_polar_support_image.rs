//! Support-image reparameterization candidates for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20PatchCorner,
    Kp20PieceCornerMaps, Kp20SectorPieceMaps, TensorBorderBi3, build_kp20_guide_ring,
    kp20_default_sector_piece_maps, sample_kp20_macro_patch_from_piece_maps,
};

use super::kp20_polar_comparison::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_sampled, is_better_summary,
    validation_is_clean,
};

pub(super) fn support_image_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let default_sampled = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(sector, kp20_default_sector_piece_maps())
        })
        .collect::<Vec<_>>();
    let baseline = candidate_summary_from_sampled(tensor_border, default_sampled)?;
    let mut best: Option<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport)> = None;
    for strength in [0.85, 0.92, 1.08, 1.15] {
        let sampled_macro_patches = guide
            .sectors
            .iter()
            .enumerate()
            .map(|(sector, guide_sector)| {
                let maps = support_image_piece_maps(
                    kp20_default_sector_piece_maps(),
                    sector,
                    strength,
                    &baseline.0.fold_sites,
                );
                sample_kp20_macro_patch_from_piece_maps(guide_sector, maps)
            })
            .collect::<Vec<_>>();
        let candidate = candidate_summary_from_sampled(tensor_border, sampled_macro_patches)?;
        if !validation_is_clean(candidate.1) {
            continue;
        }
        if is_better_summary(&candidate, Some(&baseline))
            && is_better_summary(&candidate, best.as_ref())
        {
            best = Some(candidate);
        }
    }
    best.ok_or_else(|| {
        "support-image reparameterization found no residual-clean improvement".to_string()
    })
}

fn support_image_piece_maps(
    mut maps: Kp20SectorPieceMaps,
    sector: usize,
    strength: f64,
    fold_sites: &[Kp20FoldSite],
) -> Kp20SectorPieceMaps {
    for (piece, corner_maps) in &mut maps {
        if fold_sites
            .iter()
            .any(|site| site.sector == sector && site.piece == *piece)
        {
            *corner_maps = reparameterized_corner_maps(*corner_maps, strength);
        }
    }
    maps
}

fn reparameterized_corner_maps(
    corner_maps: Kp20PieceCornerMaps,
    strength: f64,
) -> Kp20PieceCornerMaps {
    let base = corner_maps[0].1;
    let center = [
        base.origin[0] + 0.5 * base.du[0] + 0.5 * base.dv[0],
        base.origin[1] + 0.5 * base.du[1] + 0.5 * base.dv[1],
    ];
    let origin = [
        center[0] + strength * (base.origin[0] - center[0]),
        center[1] + strength * (base.origin[1] - center[1]),
    ];
    square_maps(
        origin,
        [base.du[0] * strength, base.du[1] * strength],
        [base.dv[0] * strength, base.dv[1] * strength],
    )
}

fn square_maps(origin: [f64; 2], du: [f64; 2], dv: [f64; 2]) -> Kp20PieceCornerMaps {
    [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new(origin, du, dv),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new(add2(origin, du), scale2(du, -1.0), dv),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new(add2(origin, dv), du, scale2(dv, -1.0)),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new(
                add2(add2(origin, du), dv),
                scale2(du, -1.0),
                scale2(dv, -1.0),
            ),
        ),
    ]
}

const fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

const fn scale2(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}
