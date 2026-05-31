//! Structured comparison diagnostics for KP20 polar macro sampling.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20PieceId,
    RectangularGregoryBiNonic, TensorBorderBi3, assemble_kp20_bi6_macro_cap, build_kp20_guide_ring,
    elevate_bi6_grid_to_bi9, measure_kp20_bi6_macro_cap, sample_kp20_macro_patch_from_tau_maps,
};

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
    match tau_candidate_summary(star_point, tensor_border) {
        Ok((sector_local_tau, validation)) => comparison_diagnostic(
            tensor_border.sectors.len(),
            format!(
                "bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={};sector_local_tau_status=sampled;sector_local_tau_folds={};sector_local_tau_min_jacobian_cos={:.3e};sector_local_tau_piece_folds={};tau_tensor_border_residual={:.3e};tau_internal_c2_position={:.3e};tau_internal_c2_first={:.3e};tau_origin_g2_strip={:.3e}",
                bounded_affine.fold_count,
                bounded_affine.min_jacobian_normal_cos,
                piece_counts_label(bounded_affine.piece_fold_counts),
                sector_local_tau.fold_count,
                sector_local_tau.min_jacobian_normal_cos,
                piece_counts_label(sector_local_tau.piece_fold_counts),
                validation.max_tensor_border_residual,
                validation.max_internal_c2_position_residual,
                validation.max_internal_c2_first_residual,
                validation.max_origin_g2_strip_residual,
            ),
            Some(sector_local_tau.min_jacobian_normal_cos),
        ),
        Err(err) => comparison_diagnostic(
            tensor_border.sectors.len(),
            format!(
                "bounded_affine_folds={};bounded_affine_min_jacobian_cos={:.3e};bounded_affine_piece_folds={};sector_local_tau_status=declined;sector_local_tau_error={err}",
                bounded_affine.fold_count,
                bounded_affine.min_jacobian_normal_cos,
                piece_counts_label(bounded_affine.piece_fold_counts)
            ),
            None,
        ),
    }
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
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
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
