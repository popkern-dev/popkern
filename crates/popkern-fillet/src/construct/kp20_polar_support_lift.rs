//! Transactional support-normal strip adjustment for KP20 polar macro candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    Kp20PieceId, TensorBorderBi3, assemble_kp20_bi6_macro_cap, build_kp20_guide_ring,
    sample_kp20_macro_patch_from_default_maps,
};
use popkern_core::geometry::{Point3, Vec3};

use super::kp20_polar_comparison::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_cap,
};

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportNormalLiftReport {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) accepted_steps: Vec<Kp20SupportNormalStep>,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportNormalStep {
    pub(super) piece_mask: [bool; 4],
    pub(super) scale: f64,
    pub(super) fold_count: usize,
    pub(super) piece_fold_counts: [usize; 4],
}

pub(super) fn normal_lift_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20SupportNormalLiftReport, String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(sample_kp20_macro_patch_from_default_maps)
        .collect::<Vec<_>>();
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let mut accepted = cap.macro_patches.clone();
    let mut accepted_summary = candidate_summary_from_cap(tensor_border, accepted.clone())?;
    let mut accepted_steps = Vec::new();
    let target_masks = [
        [false, true, false, false],
        [false, false, true, false],
        [false, false, false, true],
        [false, false, true, true],
    ];
    for _ in 0..4 {
        let mut best_step: Option<(
            Vec<Kp20MacroPatch>,
            (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
            f64,
            [bool; 4],
        )> = None;
        for scale in [-0.25, -0.1, 0.1, 0.25] {
            for piece_mask in target_masks {
                let mut candidate = accepted.clone();
                apply_support_normal_lift(
                    &mut candidate,
                    star_point,
                    tensor_border,
                    scale,
                    piece_mask,
                    &accepted_summary.0.fold_sites,
                );
                let summary = candidate_summary_from_cap(tensor_border, candidate.clone())?;
                if !validation_is_clean(summary.1) {
                    continue;
                }
                if is_better_summary(&summary, Some(&accepted_summary))
                    && is_better_summary(&summary, best_step.as_ref().map(|(_, s, _, _)| s))
                {
                    best_step = Some((candidate, summary, scale, piece_mask));
                }
            }
        }
        let Some((candidate, summary, scale, piece_mask)) = best_step else {
            break;
        };
        accepted = candidate;
        accepted_summary = summary;
        accepted_steps.push(Kp20SupportNormalStep {
            piece_mask,
            scale,
            fold_count: accepted_summary.0.fold_count,
            piece_fold_counts: accepted_summary.0.piece_fold_counts,
        });
    }
    Ok(Kp20SupportNormalLiftReport {
        summary: accepted_summary.0,
        validation: accepted_summary.1,
        accepted_steps,
    })
}

fn is_better_summary(
    candidate: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    best: Option<&(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport)>,
) -> bool {
    let Some(best) = best else {
        return true;
    };
    candidate.0.fold_count < best.0.fold_count
        || (candidate.0.fold_count == best.0.fold_count
            && candidate.0.min_jacobian_normal_cos > best.0.min_jacobian_normal_cos)
}

fn validation_is_clean(validation: Kp20Bi6ValidationReport) -> bool {
    validation.max_tensor_border_residual <= 1.0e-8
        && validation.max_internal_c2_position_residual <= 1.0e-8
        && validation.max_internal_c2_first_residual <= 1.0e-8
        && validation.max_origin_g2_strip_residual <= 1.0e-8
}

fn apply_support_normal_lift(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    scale: f64,
    piece_mask: [bool; 4],
    fold_sites: &[Kp20FoldSite],
) {
    for (sector, macro_patch) in macro_patches.iter_mut().enumerate() {
        let row = tensor_border.sectors[sector][0];
        let spoke_a = row[0] - star_point;
        let spoke_b = row[3] - star_point;
        let normal = spoke_a
            .cross(&spoke_b)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let span = 0.5 * (spoke_a.length() + spoke_b.length());
        let displacement = normal * (scale * span);
        for id in Kp20PieceId::all() {
            if !piece_mask[id.index()] {
                continue;
            }
            let piece = macro_patch.piece_mut(id);
            for i in 0..7 {
                for j in 0..7 {
                    if !can_lift_coefficient(id, i, j) {
                        continue;
                    }
                    let weight =
                        lift_weight(i, j) * local_lift_weight(sector, id, i, j, fold_sites);
                    piece.cps[i][j] = piece.cps[i][j] + displacement * weight;
                }
            }
        }
    }
}

fn can_lift_coefficient(id: Kp20PieceId, i: usize, j: usize) -> bool {
    match id {
        Kp20PieceId::NearMidpointLeft | Kp20PieceId::NearMidpointRight => {
            i >= 3 && i != 4 && i != 6
        }
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => {
            i != 0 && i != 1 && i != 5 && !(i >= 4 && j <= 2)
        }
    }
}

fn lift_weight(i: usize, j: usize) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    (core::f64::consts::PI * u).sin().max(0.0) * (core::f64::consts::PI * v).sin().max(0.0)
}

fn local_lift_weight(
    sector: usize,
    piece: Kp20PieceId,
    i: usize,
    j: usize,
    fold_sites: &[Kp20FoldSite],
) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let site_weight = fold_sites
        .iter()
        .filter(|site| site.sector == sector && site.piece == piece)
        .map(|site| {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.3_f64 * 0.3_f64;
            (1.0 - r2 / radius2).max(0.0).powi(2)
        })
        .fold(0.0, f64::max);
    if site_weight > 0.0 {
        site_weight.max(0.35)
    } else if fold_sites
        .iter()
        .any(|site| site.sector == sector && site.piece == piece)
    {
        0.35
    } else {
        0.0
    }
}
