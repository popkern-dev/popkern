//! Support-surface-owned KP20 polar candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6ValidationReport, Kp20MacroPatch, Kp20PieceId, TensorBorderBi3,
};
use popkern_core::geometry::{Point3, Vec3};

use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_changed_patches, is_better_summary, validation_is_clean,
};
use super::kp20_polar_support_image::Kp20SupportImageCandidate;

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceReport {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) stats: Kp20SupportSurfaceStats,
    pub(super) stage_trace: Vec<Kp20SupportSurfaceStageReport>,
    pub(super) decline_reason: &'static str,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceStats {
    pub(super) active_fold_pieces: usize,
    pub(super) orientation_rows: usize,
    pub(super) evaluated_candidates: usize,
    pub(super) exact_evaluated_candidates: usize,
    pub(super) accepted_candidates: usize,
    pub(super) rejected_candidates: usize,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20SupportSurfaceStageReport {
    pub(super) mode: &'static str,
    pub(super) evaluated_candidates: usize,
    pub(super) fast_finalist_folds: Option<usize>,
    pub(super) exact_evaluated: bool,
    pub(super) exact_folds: Option<usize>,
    pub(super) accepted: bool,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SupportSurfaceTrial {
    mode: Kp20SupportSurfaceMode,
    scale: f64,
    piece: Kp20PieceId,
}

#[derive(Debug, Clone, Copy)]
enum Kp20SupportSurfaceMode {
    SupportNormal,
    RadialBisector,
    SectorChord,
}

pub(super) fn support_surface_candidate_summary(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    support_image: &Kp20SupportImageCandidate,
) -> Result<Kp20SupportSurfaceReport, String> {
    let active_fold_pieces = active_fold_piece_count(&support_image.summary);
    let mut accepted_patches = support_image.macro_patches.clone();
    let mut accepted_summary = (support_image.summary.clone(), support_image.validation);
    let mut stats = Kp20SupportSurfaceStats {
        active_fold_pieces,
        orientation_rows: active_fold_pieces,
        evaluated_candidates: 0,
        exact_evaluated_candidates: 0,
        accepted_candidates: 0,
        rejected_candidates: 0,
    };
    let mut stage_trace = Vec::new();
    for mode in Kp20SupportSurfaceMode::stages() {
        let mut stage = Kp20SupportSurfaceStageReport {
            mode: mode.label(),
            evaluated_candidates: 0,
            fast_finalist_folds: None,
            exact_evaluated: false,
            exact_folds: None,
            accepted: false,
        };
        let mut best_fast: Option<(
            Vec<Kp20MacroPatch>,
            Vec<[bool; 4]>,
            (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        )> = None;
        for trial in support_surface_trials(mode, &accepted_summary.0) {
            stats.evaluated_candidates += 1;
            stage.evaluated_candidates += 1;
            let mut candidate = accepted_patches.clone();
            let changed_patches = apply_orientation_row_trial(
                &mut candidate,
                star_point,
                tensor_border,
                &accepted_summary.0.fold_sites,
                trial,
            );
            let fast_summary = fast_candidate_summary_from_changed_patches(
                tensor_border,
                candidate.clone(),
                &accepted_summary.0,
                &changed_patches,
                accepted_summary.0.fold_count,
            )?;
            if !validation_is_clean(fast_summary.1)
                || !is_better_summary(&fast_summary, Some(&accepted_summary))
            {
                stats.rejected_candidates += 1;
                continue;
            }
            if best_fast
                .as_ref()
                .is_some_and(|(_, _, best)| !is_better_summary(&fast_summary, Some(best)))
            {
                stats.rejected_candidates += 1;
                continue;
            }
            best_fast = Some((candidate, changed_patches, fast_summary));
        }
        let Some((candidate, changed_patches, fast_summary)) = best_fast else {
            stage_trace.push(stage);
            continue;
        };
        stage.fast_finalist_folds = Some(fast_summary.0.fold_count);
        stage.exact_evaluated = true;
        stats.exact_evaluated_candidates += 1;
        let summary = candidate_summary_from_changed_patches(
            tensor_border,
            candidate.clone(),
            &accepted_summary.0,
            &changed_patches,
            accepted_summary.0.fold_count,
        )?;
        stage.exact_folds = Some(summary.0.fold_count);
        if validation_is_clean(summary.1) && is_better_summary(&summary, Some(&accepted_summary)) {
            accepted_patches = candidate;
            accepted_summary = summary;
            stats.accepted_candidates += 1;
            stage.accepted = true;
        } else {
            stats.rejected_candidates += 1;
        }
        stage_trace.push(stage);
    }
    let decline_reason = if stats.accepted_candidates > 0 {
        "support_surface_orientation_rows_accepted"
    } else {
        "support_surface_no_improving_orientation_row"
    };
    Ok(Kp20SupportSurfaceReport {
        summary: accepted_summary.0,
        validation: accepted_summary.1,
        stats,
        stage_trace,
        decline_reason,
    })
}

fn active_fold_piece_count(summary: &Kp20SamplerFoldSummary) -> usize {
    active_fold_pieces(summary).len()
}

fn active_fold_pieces(summary: &Kp20SamplerFoldSummary) -> Vec<Kp20PieceId> {
    let mut active = [false; 4];
    for site in &summary.fold_sites {
        active[site.piece.index()] = true;
    }
    Kp20PieceId::all()
        .into_iter()
        .filter(|piece| active[piece.index()])
        .collect()
}

impl Kp20SupportSurfaceMode {
    fn stages() -> [Self; 3] {
        [Self::SupportNormal, Self::RadialBisector, Self::SectorChord]
    }

    fn label(self) -> &'static str {
        match self {
            Self::SupportNormal => "support_normal",
            Self::RadialBisector => "radial_bisector",
            Self::SectorChord => "sector_chord",
        }
    }
}

fn support_surface_trials(
    mode: Kp20SupportSurfaceMode,
    summary: &Kp20SamplerFoldSummary,
) -> Vec<Kp20SupportSurfaceTrial> {
    let mut trials = Vec::new();
    for scale in [0.06, -0.06] {
        for piece in active_fold_pieces(summary) {
            trials.push(Kp20SupportSurfaceTrial { mode, scale, piece });
        }
    }
    trials
}

fn apply_orientation_row_trial(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    fold_sites: &[Kp20FoldSite],
    trial: Kp20SupportSurfaceTrial,
) -> Vec<[bool; 4]> {
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    for site in fold_sites
        .iter()
        .filter(|site| site.piece == trial.piece)
        .copied()
    {
        let row = tensor_border.sectors[site.sector][0];
        let spoke_a = row[0] - star_point;
        let spoke_b = row[3] - star_point;
        let normal = spoke_a
            .cross(&spoke_b)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let direction = support_surface_direction(trial.mode, spoke_a, spoke_b, normal);
        let span = 0.5 * (spoke_a.length() + spoke_b.length());
        let piece = macro_patches[site.sector].piece_mut(site.piece);
        for i in 0..7 {
            for j in 0..7 {
                if !can_adjust_orientation_row(site.piece, i, j) {
                    continue;
                }
                let weight = orientation_row_weight(site, i, j);
                if weight == 0.0 {
                    continue;
                }
                piece.cps[i][j] = piece.cps[i][j] + direction * (trial.scale * span * weight);
                changed_patches[site.sector][site.piece.index()] = true;
            }
        }
    }
    changed_patches
}

fn support_surface_direction(
    mode: Kp20SupportSurfaceMode,
    spoke_a: Vec3,
    spoke_b: Vec3,
    normal: Vec3,
) -> Vec3 {
    match mode {
        Kp20SupportSurfaceMode::SupportNormal => normal,
        Kp20SupportSurfaceMode::RadialBisector => {
            (spoke_a + spoke_b).normalized().unwrap_or(normal)
        }
        Kp20SupportSurfaceMode::SectorChord => (spoke_b - spoke_a).normalized().unwrap_or(normal),
    }
}

fn can_adjust_orientation_row(piece: Kp20PieceId, i: usize, j: usize) -> bool {
    match piece {
        Kp20PieceId::NearMidpointLeft | Kp20PieceId::NearMidpointRight => {
            i == 3 && (2..=4).contains(&j)
        }
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => {
            j == 3 && (2..=4).contains(&i)
        }
    }
}

fn orientation_row_weight(site: Kp20FoldSite, i: usize, j: usize) -> f64 {
    let u = i as f64 / 6.0;
    let v = j as f64 / 6.0;
    let du = u - site.u;
    let dv = v - site.v;
    let r2 = du * du + dv * dv;
    let locality = (1.0 - r2 / (0.42_f64 * 0.42_f64)).max(0.0);
    let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
    locality * locality * severity.max(0.25)
}
