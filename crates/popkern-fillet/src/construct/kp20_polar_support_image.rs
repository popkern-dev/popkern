//! Support-image reparameterization candidates for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput,
    Kp20MacroPatch, Kp20PatchCorner, Kp20PieceCornerMaps, Kp20SectorPieceMaps, TensorBorderBi3,
    assemble_kp20_bi6_macro_cap, build_kp20_guide_ring, kp20_default_sector_piece_maps,
    sample_kp20_macro_patch_from_piece_maps,
};

use super::kp20_polar_comparison::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_cap, is_better_summary,
    validation_is_clean,
};
use super::kp20_polar_support_lift::{
    Kp20SupportNormalLiftReport, lift_candidate_from_support_image_macro_patches,
};

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportImageCandidate {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) macro_patches: Vec<Kp20MacroPatch>,
    pub(super) params: Kp20SupportImageParams,
}

struct Kp20SupportImageSearchCandidate {
    image: Kp20SupportImageCandidate,
    lift: Kp20SupportNormalLiftReport,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SupportImageParams {
    uniform: f64,
    du_scale: f64,
    dv_scale: f64,
    shear: f64,
}

pub(super) fn support_image_and_lift_candidate_summaries(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> (
    Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
    Result<Kp20SupportNormalLiftReport, String>,
) {
    match support_image_composed_candidate(star_point, tensor_border) {
        Ok(candidate) => {
            let summary = Ok((candidate.image.summary, candidate.image.validation));
            let lift = Ok(candidate.lift);
            (summary, lift)
        }
        Err(err) => (Err(err.clone()), Err(err)),
    }
}

fn support_image_composed_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20SupportImageSearchCandidate, String> {
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
    let baseline = candidate_from_sampled(
        tensor_border,
        default_sampled,
        Kp20SupportImageParams::identity(),
    )?;
    let mut image_candidates = Vec::new();
    let mut evaluated = 0_usize;
    let mut rejected = 0_usize;
    for params in Kp20SupportImageParams::coordinate_search_family() {
        evaluated += 1;
        let sampled_macro_patches = guide
            .sectors
            .iter()
            .enumerate()
            .map(|(sector, guide_sector)| {
                let maps = support_image_piece_maps(
                    kp20_default_sector_piece_maps(),
                    sector,
                    params,
                    &baseline.summary.fold_sites,
                );
                sample_kp20_macro_patch_from_piece_maps(guide_sector, maps)
            })
            .collect::<Vec<_>>();
        let candidate = candidate_from_sampled(tensor_border, sampled_macro_patches, params)?;
        if !validation_is_clean(candidate.validation) {
            rejected += 1;
            continue;
        }
        let candidate_summary = (candidate.summary.clone(), candidate.validation);
        let baseline_summary = (baseline.summary.clone(), baseline.validation);
        if is_better_summary(&candidate_summary, Some(&baseline_summary)) {
            image_candidates.push(candidate);
        } else {
            rejected += 1;
        }
    }
    image_candidates.sort_by(|a, b| compare_summaries(&a.summary, &b.summary));
    let mut best: Option<Kp20SupportImageSearchCandidate> = None;
    let mut composed_evaluated = 0_usize;
    for candidate in image_candidates.into_iter().take(1) {
        composed_evaluated += 1;
        let lift = lift_candidate_from_support_image_macro_patches(
            star_point,
            tensor_border,
            candidate.macro_patches.clone(),
        )?;
        if !validation_is_clean(lift.validation) {
            rejected += 1;
            continue;
        }
        let candidate_lift_summary = (lift.summary.clone(), lift.validation);
        let best_lift_summary = best
            .as_ref()
            .map(|candidate| (candidate.lift.summary.clone(), candidate.lift.validation));
        if is_better_summary(&candidate_lift_summary, best_lift_summary.as_ref()) {
            best = Some(Kp20SupportImageSearchCandidate {
                image: candidate,
                lift,
            });
        } else {
            rejected += 1;
        }
    }
    best.map(|mut candidate| {
        candidate
            .image
            .summary
            .fold_details
            .push(format!("map_params={}", candidate.image.params.label()));
        candidate.image.summary.fold_details.push(format!(
            "map_search=evaluated:{evaluated},composed:{composed_evaluated},accepted:1,rejected:{rejected}"
        ));
        candidate.lift.summary.fold_details.push(format!(
            "map_params={};map_search=evaluated:{evaluated},composed:{composed_evaluated},accepted:1,rejected:{rejected}",
            candidate.image.params.label()
        ));
        candidate
    })
    .ok_or_else(|| {
        "support-image reparameterization found no residual-clean improvement".to_string()
    })
}

fn compare_summaries(a: &Kp20SamplerFoldSummary, b: &Kp20SamplerFoldSummary) -> std::cmp::Ordering {
    a.fold_count.cmp(&b.fold_count).then_with(|| {
        b.min_jacobian_normal_cos
            .partial_cmp(&a.min_jacobian_normal_cos)
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}

fn candidate_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<Kp20MacroPatch>,
    params: Kp20SupportImageParams,
) -> Result<Kp20SupportImageCandidate, String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let (summary, validation) =
        candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
    Ok(Kp20SupportImageCandidate {
        summary,
        validation,
        macro_patches: cap.macro_patches,
        params,
    })
}

fn support_image_piece_maps(
    mut maps: Kp20SectorPieceMaps,
    sector: usize,
    params: Kp20SupportImageParams,
    fold_sites: &[Kp20FoldSite],
) -> Kp20SectorPieceMaps {
    for (piece, corner_maps) in &mut maps {
        if fold_sites
            .iter()
            .any(|site| site.sector == sector && site.piece == *piece)
        {
            *corner_maps = reparameterized_corner_maps(*corner_maps, params);
        }
    }
    maps
}

fn reparameterized_corner_maps(
    corner_maps: Kp20PieceCornerMaps,
    params: Kp20SupportImageParams,
) -> Kp20PieceCornerMaps {
    let base = corner_maps[0].1;
    let center = [
        base.origin[0] + 0.5 * base.du[0] + 0.5 * base.dv[0],
        base.origin[1] + 0.5 * base.du[1] + 0.5 * base.dv[1],
    ];
    let du = [
        params.uniform * params.du_scale * base.du[0] + params.shear * base.dv[0],
        params.uniform * params.du_scale * base.du[1] + params.shear * base.dv[1],
    ];
    let dv = [
        params.shear * base.du[0] + params.uniform * params.dv_scale * base.dv[0],
        params.shear * base.du[1] + params.uniform * params.dv_scale * base.dv[1],
    ];
    let origin = [
        center[0] - 0.5 * du[0] - 0.5 * dv[0],
        center[1] - 0.5 * du[1] - 0.5 * dv[1],
    ];
    square_maps(origin, du, dv)
}

impl Kp20SupportImageParams {
    fn identity() -> Self {
        Self {
            uniform: 1.0,
            du_scale: 1.0,
            dv_scale: 1.0,
            shear: 0.0,
        }
    }

    fn coordinate_search_family() -> [Self; 9] {
        [
            Self {
                uniform: 0.92,
                du_scale: 1.0,
                dv_scale: 1.0,
                shear: 0.04,
            },
            Self {
                uniform: 0.89,
                du_scale: 1.0,
                dv_scale: 1.0,
                shear: 0.04,
            },
            Self {
                uniform: 0.95,
                du_scale: 1.0,
                dv_scale: 1.0,
                shear: 0.04,
            },
            Self {
                uniform: 0.92,
                du_scale: 0.96,
                dv_scale: 1.0,
                shear: 0.04,
            },
            Self {
                uniform: 0.92,
                du_scale: 1.04,
                dv_scale: 1.0,
                shear: 0.04,
            },
            Self {
                uniform: 0.92,
                du_scale: 1.0,
                dv_scale: 0.96,
                shear: 0.04,
            },
            Self {
                uniform: 0.92,
                du_scale: 1.0,
                dv_scale: 1.04,
                shear: 0.04,
            },
            Self {
                uniform: 0.92,
                du_scale: 1.0,
                dv_scale: 1.0,
                shear: 0.02,
            },
            Self {
                uniform: 0.92,
                du_scale: 1.0,
                dv_scale: 1.0,
                shear: 0.06,
            },
        ]
    }

    fn label(self) -> String {
        format!(
            "uniform:{:.3},du:{:.3},dv:{:.3},shear:{:.3}",
            self.uniform, self.du_scale, self.dv_scale, self.shear
        )
    }
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
