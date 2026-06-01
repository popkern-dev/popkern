//! Structured comparison diagnostics for KP20 polar macro sampling.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput,
    Kp20MacroPatch, Kp20PieceId, RectangularGregoryBiNonic, TensorBorderBi3,
    assemble_kp20_bi6_macro_cap, build_kp20_guide_ring, elevate_bi6_grid_to_bi9,
    measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_default_maps,
    sample_kp20_macro_patch_from_orientation_valid_maps, sample_kp20_macro_patch_from_tau_maps,
};
use popkern_core::geometry::{Point3, Vec3};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use crate::types::ApexConstructionDiagnostic;

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20SamplerFoldSummary {
    pub fold_count: usize,
    pub min_jacobian_normal_cos: f64,
    pub piece_fold_counts: [usize; 4],
}

pub(super) fn sampler_comparison_diagnostic(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    bounded_affine: Kp20SamplerFoldSummary,
) -> ApexConstructionDiagnostic {
    let orientation_valid = orientation_valid_candidate_summary(star_point, tensor_border);
    let normal_lift = normal_lift_candidate_summary(star_point, tensor_border);
    let sector_local_tau = tau_candidate_summary(star_point, tensor_border);
    let mut detail = format!(
        "bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={}",
        bounded_affine.fold_count,
        bounded_affine.min_jacobian_normal_cos,
        piece_counts_label(bounded_affine.piece_fold_counts)
    );
    append_candidate_detail(&mut detail, "orientation_valid", &orientation_valid);
    append_candidate_detail(&mut detail, "support_normal_lift", &normal_lift);
    append_candidate_detail(&mut detail, "sector_local_tau", &sector_local_tau);
    comparison_diagnostic(
        tensor_border.sectors.len(),
        detail,
        sector_local_tau
            .as_ref()
            .ok()
            .map(|(summary, _)| summary.min_jacobian_normal_cos),
    )
}

pub(super) fn summarize_intersections(
    intersections: &[(
        usize,
        popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    )],
) -> Kp20SamplerFoldSummary {
    let mut piece_fold_counts = [0_usize; 4];
    for (idx, _) in intersections {
        piece_fold_counts[idx % 4] += 1;
    }
    let min_jacobian_normal_cos = intersections
        .iter()
        .map(|(_, report)| report.min_jacobian_normal_cos)
        .fold(f64::INFINITY, f64::min);
    Kp20SamplerFoldSummary {
        fold_count: intersections.len(),
        min_jacobian_normal_cos,
        piece_fold_counts,
    }
}

fn tau_candidate_summary(
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

fn orientation_valid_candidate_summary(
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
        .map(sample_kp20_macro_patch_from_orientation_valid_maps)
        .collect::<Vec<_>>();
    candidate_summary_from_sampled(tensor_border, sampled_macro_patches)
}

fn candidate_summary_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<popkern_core::geometry::gregory::Kp20MacroPatch>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    candidate_summary_from_cap(tensor_border, cap.macro_patches)
}

fn normal_lift_candidate_summary(
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
        .map(sample_kp20_macro_patch_from_default_maps)
        .collect::<Vec<_>>();
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    let mut best = None;
    let piece_masks = [
        [true, true, true, true],
        [true, false, false, false],
        [false, true, false, false],
        [false, false, true, false],
        [false, false, false, true],
        [true, true, false, false],
        [false, false, true, true],
    ];
    for scale in [-0.25, -0.1, 0.1, 0.25] {
        for piece_mask in piece_masks {
            let mut candidate = cap.macro_patches.clone();
            apply_support_normal_lift(&mut candidate, star_point, tensor_border, scale, piece_mask);
            let summary = candidate_summary_from_cap(tensor_border, candidate)?;
            if is_better_summary(&summary, best.as_ref()) {
                best = Some(summary);
            }
        }
    }
    best.ok_or_else(|| "support-normal lift search had no candidates".to_string())
}

fn candidate_summary_from_cap(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = popkern_core::geometry::gregory::Kp20Bi6MacroCapOutput {
        n: tensor_border.sectors.len(),
        patches_bi6: flatten_macro_patches(&macro_patches),
        macro_patches,
        strip_projections: Vec::new(),
    };
    let validation = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    let candidate = cap
        .patches_bi6
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect::<Vec<_>>();
    let intersections = candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect::<Vec<_>>();
    Ok((summarize_intersections(&intersections), validation))
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

fn apply_support_normal_lift(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    scale: f64,
    piece_mask: [bool; 4],
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
                    let weight = lift_weight(i, j);
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

fn flatten_macro_patches(macro_patches: &[Kp20MacroPatch]) -> Vec<[[Point3; 7]; 7]> {
    let mut out = Vec::with_capacity(4 * macro_patches.len());
    for macro_patch in macro_patches {
        for id in Kp20PieceId::all() {
            out.push(macro_patch.piece(id).cps);
        }
    }
    out
}

fn append_candidate_detail(
    detail: &mut String,
    label: &str,
    summary: &Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String>,
) {
    match summary {
        Ok((folds, validation)) => {
            detail.push_str(&format!(
                ";{label}_status=sampled;{label}_folds={};{label}_min_jacobian_cos={:.3e};{label}_piece_folds={};{label}_tensor_border_residual={:.3e};{label}_internal_c2_position={:.3e};{label}_internal_c2_first={:.3e};{label}_origin_g2_strip={:.3e}",
                folds.fold_count,
                folds.min_jacobian_normal_cos,
                piece_counts_label(folds.piece_fold_counts),
                validation.max_tensor_border_residual,
                validation.max_internal_c2_position_residual,
                validation.max_internal_c2_first_residual,
                validation.max_origin_g2_strip_residual,
            ));
        }
        Err(err) => {
            detail.push_str(&format!(";{label}_status=declined;{label}_error={err}"));
        }
    }
}

fn comparison_diagnostic(
    n_sides: usize,
    detail: String,
    residual: Option<f64>,
) -> ApexConstructionDiagnostic {
    super::kp20_polar_diagnostics::dispatch_fact(
        "sampler_comparison",
        n_sides,
        Some(detail),
        residual,
        Some(SELF_INT_TOL),
    )
}

fn piece_counts_label(counts: [usize; 4]) -> String {
    let pieces = Kp20PieceId::all();
    format!(
        "{:?}:{}|{:?}:{}|{:?}:{}|{:?}:{}",
        pieces[0], counts[0], pieces[1], counts[1], pieces[2], counts[2], pieces[3], counts[3]
    )
}
