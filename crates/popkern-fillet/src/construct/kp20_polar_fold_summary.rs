//! Fold and validation summaries for KP20 polar macro candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId,
    RectangularGregoryBiNonic, TensorBorderBi3, assemble_kp20_bi6_macro_cap,
    elevate_bi6_grid_to_bi9, measure_kp20_bi6_macro_cap,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;

#[derive(Debug, Clone)]
pub(super) struct Kp20SamplerFoldSummary {
    pub fold_count: usize,
    pub min_jacobian_normal_cos: f64,
    pub piece_fold_counts: [usize; 4],
    pub fold_details: Vec<String>,
    pub(super) fold_sites: Vec<Kp20FoldSite>,
}

impl Kp20SamplerFoldSummary {
    pub(super) fn worst_fold_site(&self) -> Option<Kp20FoldSite> {
        self.fold_sites.iter().copied().min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    pub(super) fn fold_zone_counts(&self) -> Kp20FoldZoneCounts {
        let mut counts = Kp20FoldZoneCounts::default();
        for site in &self.fold_sites {
            if site.u >= 0.80 && site.v >= 0.80 {
                counts.origin_corner += 1;
            } else if site.u >= 0.80 || site.v >= 0.80 {
                counts.edge_shoulder += 1;
            } else {
                counts.interior += 1;
            }
        }
        counts
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20FoldSite {
    pub(super) sector: usize,
    pub(super) piece: Kp20PieceId,
    pub(super) u: f64,
    pub(super) v: f64,
    pub(super) min_jacobian_ratio: f64,
    pub(super) min_jacobian_normal_cos: f64,
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Kp20FoldZoneCounts {
    pub(super) origin_corner: usize,
    pub(super) edge_shoulder: usize,
    pub(super) interior: usize,
}

pub(super) fn summarize_intersections(
    intersections: &[(
        usize,
        popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    )],
) -> Kp20SamplerFoldSummary {
    let mut piece_fold_counts = [0_usize; 4];
    let mut fold_details = Vec::new();
    for (idx, report) in intersections {
        piece_fold_counts[idx % 4] += 1;
        let sector = idx / 4;
        let piece = Kp20PieceId::all()[idx % 4];
        fold_details.push(format!(
            "sector={sector},piece={piece:?},ratio={:.3e},cos={:.3e}",
            report.min_jacobian_ratio, report.min_jacobian_normal_cos
        ));
    }
    let min_jacobian_normal_cos = intersections
        .iter()
        .map(|(_, report)| report.min_jacobian_normal_cos)
        .fold(f64::INFINITY, f64::min);
    Kp20SamplerFoldSummary {
        fold_count: intersections.len(),
        min_jacobian_normal_cos,
        piece_fold_counts,
        fold_details,
        fold_sites: Vec::new(),
    }
}

pub(super) fn candidate_summary_from_sampled(
    tensor_border: &TensorBorderBi3,
    sampled_macro_patches: Vec<Kp20MacroPatch>,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches,
    })
    .map_err(|e| e.to_string())?;
    candidate_summary_from_cap(tensor_border, cap.macro_patches)
}

pub(super) fn is_better_summary(
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

pub(super) fn validation_is_clean(validation: Kp20Bi6ValidationReport) -> bool {
    validation.max_tensor_border_residual <= 1.0e-8
        && validation.max_internal_c2_position_residual <= 1.0e-8
        && validation.max_internal_c2_first_residual <= 1.0e-8
        && validation.max_origin_g2_strip_residual <= 1.0e-8
}

pub(super) fn gregory_patches_from_macro_patches(
    macro_patches: &[Kp20MacroPatch],
) -> Vec<GregoryPatch> {
    flatten_macro_patches(macro_patches)
        .into_iter()
        .map(elevate_bi6_grid_to_bi9)
        .map(|grid| GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid })))
        .collect()
}

pub(super) fn candidate_summary_from_cap(
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
    let candidate = gregory_patches_from_macro_patches(&cap.macro_patches);
    let intersections = candidate
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect::<Vec<_>>();
    Ok((
        summarize_candidate_intersections(&candidate, &intersections),
        validation,
    ))
}

pub(super) fn candidate_summary_from_changed_patches(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    unchanged_summary: &Kp20SamplerFoldSummary,
    changed_patches: &[[bool; 4]],
    max_candidate_folds: usize,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    candidate_summary_from_changed_patches_with_check(
        tensor_border,
        macro_patches,
        unchanged_summary,
        changed_patches,
        max_candidate_folds,
        FoldCheck::Exact,
    )
}

pub(super) fn fast_candidate_summary_from_changed_patches(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    unchanged_summary: &Kp20SamplerFoldSummary,
    changed_patches: &[[bool; 4]],
    max_candidate_folds: usize,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    candidate_summary_from_changed_patches_with_check(
        tensor_border,
        macro_patches,
        unchanged_summary,
        changed_patches,
        max_candidate_folds,
        FoldCheck::JacobianOnly,
    )
}

pub(super) fn fast_candidate_summary_from_changed_piece_sites(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    unchanged_summary: &Kp20SamplerFoldSummary,
    changed_piece_sites: &[(usize, Kp20PieceId)],
    max_candidate_folds: usize,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let mut changed_patches = vec![[false; 4]; tensor_border.sectors.len()];
    for &(sector, piece) in changed_piece_sites {
        changed_patches[sector][piece.index()] = true;
    }
    candidate_summary_from_changed_patches_with_check(
        tensor_border,
        macro_patches,
        unchanged_summary,
        &changed_patches,
        max_candidate_folds,
        FoldCheck::JacobianOnly,
    )
}

fn candidate_summary_from_changed_patches_with_check(
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    unchanged_summary: &Kp20SamplerFoldSummary,
    changed_patches: &[[bool; 4]],
    max_candidate_folds: usize,
    check: FoldCheck,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let cap = popkern_core::geometry::gregory::Kp20Bi6MacroCapOutput {
        n: tensor_border.sectors.len(),
        patches_bi6: flatten_macro_patches(&macro_patches),
        macro_patches,
        strip_projections: Vec::new(),
    };
    let validation = measure_kp20_bi6_macro_cap(&cap, tensor_border).map_err(|e| e.to_string())?;
    let mut fold_sites = unchanged_summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| !changed_patches[site.sector][site.piece.index()])
        .collect::<Vec<_>>();
    let patch_grids = cap.patches_bi6;
    for idx in changed_patch_indices(changed_patches, unchanged_summary) {
        let sector = idx / 4;
        let piece = Kp20PieceId::all()[idx % 4];
        let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
            grid: elevate_bi6_grid_to_bi9(patch_grids[idx]),
        }));
        let report = check.report(&patch);
        if report.intersects() {
            let sample = worst_jacobian_sample(&patch);
            fold_sites.push(Kp20FoldSite {
                sector,
                piece,
                u: sample.u,
                v: sample.v,
                min_jacobian_ratio: report.min_jacobian_ratio,
                min_jacobian_normal_cos: report.min_jacobian_normal_cos,
            });
            if fold_sites.len() > max_candidate_folds {
                break;
            }
        }
    }
    Ok((summary_from_fold_sites(fold_sites), validation))
}

#[derive(Debug, Clone, Copy)]
enum FoldCheck {
    Exact,
    JacobianOnly,
}

impl FoldCheck {
    fn report(
        self,
        patch: &GregoryPatch,
    ) -> popkern_core::geometry::gregory::GregorySelfIntersectionReport {
        match self {
            Self::Exact => patch.self_intersection_report(SELF_INT_TOL),
            Self::JacobianOnly => patch.jacobian_self_intersection_report(),
        }
    }
}

fn changed_patch_indices(
    changed_patches: &[[bool; 4]],
    unchanged_summary: &Kp20SamplerFoldSummary,
) -> Vec<usize> {
    let mut indices = Vec::new();
    for (sector, pieces) in changed_patches.iter().enumerate() {
        for piece in Kp20PieceId::all() {
            if pieces[piece.index()] {
                indices.push(sector * 4 + piece.index());
            }
        }
    }
    indices.sort_by(|a, b| {
        previous_fold_count(*b, unchanged_summary).cmp(&previous_fold_count(*a, unchanged_summary))
    });
    indices
}

fn previous_fold_count(idx: usize, summary: &Kp20SamplerFoldSummary) -> usize {
    let sector = idx / 4;
    let piece = Kp20PieceId::all()[idx % 4];
    summary
        .fold_sites
        .iter()
        .filter(|site| site.sector == sector && site.piece == piece)
        .count()
}

pub(super) fn summarize_candidate_intersections(
    candidate: &[GregoryPatch],
    intersections: &[(
        usize,
        popkern_core::geometry::gregory::GregorySelfIntersectionReport,
    )],
) -> Kp20SamplerFoldSummary {
    let mut summary = summarize_intersections(intersections);
    summary.fold_sites = intersections
        .iter()
        .map(|(idx, report)| {
            let sector = idx / 4;
            let piece = Kp20PieceId::all()[idx % 4];
            let sample = worst_jacobian_sample(&candidate[*idx]);
            Kp20FoldSite {
                sector,
                piece,
                u: sample.u,
                v: sample.v,
                min_jacobian_ratio: report.min_jacobian_ratio,
                min_jacobian_normal_cos: report.min_jacobian_normal_cos,
            }
        })
        .collect();
    summary.fold_details = fold_details_from_sites(&summary.fold_sites);
    summary
}

fn summary_from_fold_sites(fold_sites: Vec<Kp20FoldSite>) -> Kp20SamplerFoldSummary {
    let mut piece_fold_counts = [0_usize; 4];
    for site in &fold_sites {
        piece_fold_counts[site.piece.index()] += 1;
    }
    let min_jacobian_normal_cos = fold_sites
        .iter()
        .map(|site| site.min_jacobian_normal_cos)
        .fold(f64::INFINITY, f64::min);
    let fold_count = fold_sites.len();
    Kp20SamplerFoldSummary {
        fold_count,
        min_jacobian_normal_cos,
        piece_fold_counts,
        fold_details: fold_details_from_sites(&fold_sites),
        fold_sites,
    }
}

fn fold_details_from_sites(fold_sites: &[Kp20FoldSite]) -> Vec<String> {
    fold_sites
        .iter()
        .map(|site| {
            format!(
                "sector={},piece={:?},u={:.3},v={:.3},ratio={:.3e},cos={:.3e}",
                site.sector,
                site.piece,
                site.u,
                site.v,
                site.min_jacobian_ratio,
                site.min_jacobian_normal_cos
            )
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct WorstJacobianSample {
    u: f64,
    v: f64,
    min_jacobian_normal_cos: f64,
}

fn worst_jacobian_sample(patch: &GregoryPatch) -> WorstJacobianSample {
    let n = 8_usize;
    let mut samples = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let u = (i as f64 + 0.5) / n as f64;
            let v = (j as f64 + 0.5) / n as f64;
            let jacobian = patch.partial_u(u, v).cross(&patch.partial_v(u, v));
            let mag = jacobian.length();
            if mag >= 1.0e-14 {
                samples.push((u, v, jacobian, mag));
            }
        }
    }
    let Some((_, _, reference, max_mag)) = samples
        .iter()
        .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal))
        .copied()
    else {
        return WorstJacobianSample {
            u: 0.5,
            v: 0.5,
            min_jacobian_normal_cos: 1.0,
        };
    };
    let mag_threshold = 0.01 * max_mag;
    let mut worst = WorstJacobianSample {
        u: 0.5,
        v: 0.5,
        min_jacobian_normal_cos: 1.0,
    };
    for (u, v, jacobian, mag) in samples {
        if mag < mag_threshold {
            continue;
        }
        let cos = jacobian.dot(&reference) / (mag * max_mag);
        if cos < worst.min_jacobian_normal_cos {
            worst.u = u;
            worst.v = v;
            worst.min_jacobian_normal_cos = cos;
        }
    }
    worst
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fold_zone_counts_classify_sampler_edge_and_corner_sites() {
        let summary = Kp20SamplerFoldSummary {
            fold_count: 3,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [0; 4],
            fold_details: Vec::new(),
            fold_sites: vec![
                fold_site(0.8125, 0.9375),
                fold_site(0.8125, 0.1875),
                fold_site(0.5, 0.5),
            ],
        };

        let zones = summary.fold_zone_counts();

        assert_eq!(zones.origin_corner, 1);
        assert_eq!(zones.edge_shoulder, 1);
        assert_eq!(zones.interior, 1);
    }

    fn fold_site(u: f64, v: f64) -> Kp20FoldSite {
        Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearOriginLeft,
            u,
            v,
            min_jacobian_ratio: -0.1,
            min_jacobian_normal_cos: -0.1,
        }
    }
}
