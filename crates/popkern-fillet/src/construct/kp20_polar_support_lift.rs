use super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_changed_patches,
    fast_candidate_summary_from_cap, fast_candidate_summary_from_changed_patches,
    is_better_summary, validation_is_clean,
};
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideBuildInput, Kp20MacroPatch,
    Kp20PieceId, TensorBorderBi3, assemble_kp20_bi6_macro_cap, build_kp20_guide_ring,
    sample_kp20_macro_patch_from_default_maps,
};
use popkern_core::geometry::{Point3, Vec3};
const MAX_SUPPORT_LIFT_CANDIDATES: usize = 64;
#[derive(Debug, Clone)]
pub(super) struct Kp20SupportNormalLiftReport {
    pub(super) summary: Kp20SamplerFoldSummary,
    pub(super) validation: Kp20Bi6ValidationReport,
    pub(super) macro_patches: Vec<Kp20MacroPatch>,
    pub(super) accepted_steps: Vec<Kp20SupportNormalStep>,
    pub(super) optimizer_stats: Kp20SupportOptimizerStats,
}
#[derive(Debug, Clone)]
pub(super) struct Kp20SupportNormalStep {
    pub(super) piece_mask: [bool; 4],
    pub(super) mode: Kp20SupportLiftMode,
    pub(super) scale: f64,
    pub(super) fold_count: usize,
    pub(super) piece_fold_counts: [usize; 4],
}
#[derive(Debug, Clone, Copy)]
pub(super) enum Kp20SupportLiftMode {
    SupportNormal,
    RadialBisector,
    SectorChord,
    FoldSiteLeastSquares,
}
#[derive(Debug, Default, Clone)]
pub(super) struct Kp20SupportOptimizerStats {
    pub(super) active_fold_sites: usize,
    pub(super) row_variable_count: usize,
    pub(super) evaluated_candidates: usize,
    pub(super) accepted_candidates: usize,
    pub(super) rejected_candidates: usize,
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
    lift_candidate_from_macro_patches(star_point, tensor_border, cap.macro_patches)
}

pub(super) fn lift_candidate_from_macro_patches(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Result<Kp20SupportNormalLiftReport, String> {
    lift_candidate_from_macro_patches_with_profile(
        star_point,
        tensor_border,
        macro_patches,
        Kp20LiftSearchProfile::General,
        MAX_SUPPORT_LIFT_CANDIDATES,
    )
}

pub(super) fn lift_candidate_from_support_image_macro_patches(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
) -> Result<Kp20SupportNormalLiftReport, String> {
    lift_candidate_from_macro_patches_with_profile(
        star_point,
        tensor_border,
        macro_patches,
        Kp20LiftSearchProfile::ComposedSupportImage,
        MAX_SUPPORT_LIFT_CANDIDATES,
    )
}

pub(super) fn lift_candidate_from_macro_patches_with_profile(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    macro_patches: Vec<Kp20MacroPatch>,
    profile: Kp20LiftSearchProfile,
    max_candidates: usize,
) -> Result<Kp20SupportNormalLiftReport, String> {
    let mut accepted = macro_patches;
    let mut accepted_summary = fast_candidate_summary_from_cap(tensor_border, accepted.clone())?;
    let mut accepted_steps = Vec::new();
    let mut optimizer_stats = Kp20SupportOptimizerStats::default();
    'lift_search: for _ in 0..4 {
        optimizer_stats.active_fold_sites = accepted_summary.0.fold_sites.len();
        optimizer_stats.row_variable_count =
            movable_coefficient_count(&accepted_summary.0.fold_sites);
        let mut best_step: Option<(
            Vec<Kp20MacroPatch>,
            (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
            Kp20SupportLiftMode,
            f64,
            [bool; 4],
        )> = None;
        for &mode in Kp20SupportLiftMode::search_order(profile) {
            let mut best_in_mode: Option<(
                Vec<Kp20MacroPatch>,
                (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
                Kp20SupportLiftMode,
                f64,
                [bool; 4],
            )> = None;
            'mode_search: for &scale in Kp20SupportLiftMode::scale_order(profile) {
                for piece_mask in
                    candidate_piece_masks(profile, mode, &accepted_summary.0.fold_sites)
                {
                    if optimizer_stats.evaluated_candidates >= max_candidates {
                        break 'lift_search;
                    }
                    let mut candidate = accepted.clone();
                    let changed_patches = apply_support_row_update(
                        &mut candidate,
                        star_point,
                        tensor_border,
                        mode,
                        scale,
                        piece_mask,
                        &accepted_summary.0.fold_sites,
                    );
                    optimizer_stats.evaluated_candidates += 1;
                    let fast_summary = fast_candidate_summary_from_changed_patches(
                        tensor_border,
                        candidate.clone(),
                        &accepted_summary.0,
                        &changed_patches,
                        accepted_summary.0.fold_count,
                    )?;
                    if !validation_is_clean(fast_summary.1) {
                        optimizer_stats.rejected_candidates += 1;
                        continue;
                    }
                    if !is_better_summary(&fast_summary, Some(&accepted_summary)) {
                        optimizer_stats.rejected_candidates += 1;
                        continue;
                    }
                    if best_in_mode
                        .as_ref()
                        .is_some_and(|(_, s, _, _, _)| !is_better_summary(&fast_summary, Some(s)))
                    {
                        optimizer_stats.rejected_candidates += 1;
                        continue;
                    }
                    let summary = candidate_summary_from_changed_patches(
                        tensor_border,
                        candidate.clone(),
                        &accepted_summary.0,
                        &changed_patches,
                        accepted_summary.0.fold_count,
                    )?;
                    if !validation_is_clean(summary.1)
                        || !is_better_summary(&summary, Some(&accepted_summary))
                    {
                        optimizer_stats.rejected_candidates += 1;
                        continue;
                    }
                    if best_in_mode
                        .as_ref()
                        .is_some_and(|(_, s, _, _, _)| !is_better_summary(&summary, Some(s)))
                    {
                        optimizer_stats.rejected_candidates += 1;
                        continue;
                    }
                    let drops_fold_count = summary.0.fold_count < accepted_summary.0.fold_count;
                    best_in_mode = Some((candidate, summary, mode, scale, piece_mask));
                    if profile == Kp20LiftSearchProfile::ComposedSupportImage && drops_fold_count {
                        break 'mode_search;
                    }
                }
            }
            if best_in_mode.as_ref().is_some_and(|(_, summary, _, _, _)| {
                best_step
                    .as_ref()
                    .is_none_or(|(_, best, _, _, _)| is_better_summary(summary, Some(best)))
            }) {
                best_step = best_in_mode;
            }
            if best_step.as_ref().is_some_and(|(_, summary, _, _, _)| {
                summary.0.fold_count < accepted_summary.0.fold_count
            }) {
                break;
            }
        }
        let Some((candidate, summary, mode, scale, piece_mask)) = best_step else {
            break;
        };
        accepted = candidate;
        accepted_summary = summary;
        optimizer_stats.accepted_candidates += 1;
        accepted_steps.push(Kp20SupportNormalStep {
            piece_mask,
            mode,
            scale,
            fold_count: accepted_summary.0.fold_count,
            piece_fold_counts: accepted_summary.0.piece_fold_counts,
        });
    }
    Ok(Kp20SupportNormalLiftReport {
        summary: accepted_summary.0,
        validation: accepted_summary.1,
        macro_patches: accepted,
        accepted_steps,
        optimizer_stats,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kp20LiftSearchProfile {
    General,
    ComposedSupportImage,
}

impl Kp20SupportLiftMode {
    fn search_order(profile: Kp20LiftSearchProfile) -> &'static [Self] {
        match profile {
            Kp20LiftSearchProfile::General => &[
                Self::SupportNormal,
                Self::FoldSiteLeastSquares,
                Self::RadialBisector,
                Self::SectorChord,
            ],
            Kp20LiftSearchProfile::ComposedSupportImage => &[
                Self::SupportNormal,
                Self::FoldSiteLeastSquares,
                Self::RadialBisector,
            ],
        }
    }

    fn scale_order(profile: Kp20LiftSearchProfile) -> &'static [f64] {
        match profile {
            Kp20LiftSearchProfile::General => &[-0.25, -0.1, 0.1, 0.25],
            Kp20LiftSearchProfile::ComposedSupportImage => &[0.25, 0.1, -0.1, -0.25],
        }
    }
}

fn candidate_piece_masks(
    profile: Kp20LiftSearchProfile,
    mode: Kp20SupportLiftMode,
    fold_sites: &[Kp20FoldSite],
) -> Vec<[bool; 4]> {
    let mut active = [false; 4];
    for site in fold_sites {
        active[site.piece.index()] = true;
    }
    if matches!(mode, Kp20SupportLiftMode::FoldSiteLeastSquares)
        || !matches!(mode, Kp20SupportLiftMode::SupportNormal)
    {
        return vec![active];
    }
    let mut masks = Vec::new();
    for piece in ordered_active_pieces(profile, fold_sites) {
        if active[piece.index()] {
            let mut mask = [false; 4];
            mask[piece.index()] = true;
            masks.push(mask);
        }
    }
    push_mask_if_active(&mut masks, active, [true, true, false, false]);
    push_mask_if_active(&mut masks, active, [false, false, true, true]);
    masks
}

fn ordered_active_pieces(
    profile: Kp20LiftSearchProfile,
    fold_sites: &[Kp20FoldSite],
) -> Vec<Kp20PieceId> {
    let mut pieces = Kp20PieceId::all().to_vec();
    if profile == Kp20LiftSearchProfile::ComposedSupportImage {
        pieces.sort_by(|a, b| {
            fold_sites
                .iter()
                .filter(|site| site.piece == *b)
                .count()
                .cmp(&fold_sites.iter().filter(|site| site.piece == *a).count())
                .then_with(|| a.index().cmp(&b.index()))
        });
    }
    pieces
}

fn push_mask_if_active(masks: &mut Vec<[bool; 4]>, active: [bool; 4], mask: [bool; 4]) {
    let useful = mask.iter().filter(|enabled| **enabled).count() >= 2
        && mask
            .into_iter()
            .enumerate()
            .all(|(idx, enabled)| !enabled || active[idx]);
    if useful && !masks.contains(&mask) {
        masks.push(mask);
    }
}

fn apply_support_row_update(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
    mode: Kp20SupportLiftMode,
    scale: f64,
    piece_mask: [bool; 4],
    fold_sites: &[Kp20FoldSite],
) -> Vec<[bool; 4]> {
    let mut changed_patches = vec![[false; 4]; macro_patches.len()];
    for (sector, macro_patch) in macro_patches.iter_mut().enumerate() {
        let row = tensor_border.sectors[sector][0];
        let spoke_a = row[0] - star_point;
        let spoke_b = row[3] - star_point;
        let normal = spoke_a
            .cross(&spoke_b)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let frame = SupportLiftFrame::new(mode, spoke_a, spoke_b, normal);
        let span = 0.5 * (spoke_a.length() + spoke_b.length());
        for id in Kp20PieceId::all() {
            if !piece_mask[id.index()] {
                continue;
            }
            if !has_fold_site(sector, id, fold_sites) {
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
                    if weight == 0.0 {
                        continue;
                    }
                    let displacement =
                        frame.displacement(sector, id, i, j, fold_sites) * (scale * span * weight);
                    piece.cps[i][j] = piece.cps[i][j] + displacement;
                    changed_patches[sector][id.index()] = true;
                }
            }
        }
    }
    changed_patches
}

fn has_fold_site(sector: usize, piece: Kp20PieceId, fold_sites: &[Kp20FoldSite]) -> bool {
    fold_sites
        .iter()
        .any(|site| site.sector == sector && site.piece == piece)
}

struct SupportLiftFrame {
    mode: Kp20SupportLiftMode,
    normal: Vec3,
    radial: Vec3,
    chord: Vec3,
}

impl SupportLiftFrame {
    fn new(mode: Kp20SupportLiftMode, spoke_a: Vec3, spoke_b: Vec3, normal: Vec3) -> Self {
        Self {
            mode,
            normal,
            radial: (spoke_a + spoke_b).normalized().unwrap_or(normal),
            chord: (spoke_b - spoke_a).normalized().unwrap_or(normal),
        }
    }

    fn displacement(
        &self,
        sector: usize,
        piece: Kp20PieceId,
        i: usize,
        j: usize,
        fold_sites: &[Kp20FoldSite],
    ) -> Vec3 {
        match self.mode {
            Kp20SupportLiftMode::SupportNormal => self.normal,
            Kp20SupportLiftMode::RadialBisector => self.radial,
            Kp20SupportLiftMode::SectorChord => self.chord,
            Kp20SupportLiftMode::FoldSiteLeastSquares => self
                .least_squares_displacement(sector, piece, i, j, fold_sites)
                .unwrap_or(self.normal),
        }
    }

    fn least_squares_displacement(
        &self,
        sector: usize,
        piece: Kp20PieceId,
        i: usize,
        j: usize,
        fold_sites: &[Kp20FoldSite],
    ) -> Option<Vec3> {
        let u = i as f64 / 6.0;
        let v = j as f64 / 6.0;
        let mut normal_weight = 0.0;
        let mut radial_weight = 0.0;
        let mut chord_weight = 0.0;
        for site in fold_sites
            .iter()
            .filter(|site| site.sector == sector && site.piece == piece)
        {
            let du = u - site.u;
            let dv = v - site.v;
            let r2 = du * du + dv * dv;
            let radius2 = 0.3_f64 * 0.3_f64;
            let locality = (1.0 - r2 / radius2).max(0.0);
            if locality == 0.0 {
                continue;
            }
            let severity = (1.0 - site.min_jacobian_normal_cos).max(0.0);
            let weight = locality * locality * severity;
            normal_weight += weight;
            radial_weight += weight * du;
            chord_weight += weight * dv;
        }
        let displacement =
            self.normal * normal_weight + self.radial * radial_weight + self.chord * chord_weight;
        displacement.normalized()
    }
}

fn movable_coefficient_count(fold_sites: &[Kp20FoldSite]) -> usize {
    let mut count = 0;
    for site in fold_sites {
        for i in 0..7 {
            for j in 0..7 {
                if can_lift_coefficient(site.piece, i, j) {
                    count += 1;
                }
            }
        }
    }
    count
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
    } else if has_fold_site(sector, piece, fold_sites) {
        0.35
    } else {
        0.0
    }
}
