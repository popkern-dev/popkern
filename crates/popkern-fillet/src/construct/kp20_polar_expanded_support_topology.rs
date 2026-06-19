//! Expanded support-patch topology measurement for KP20 limited lift.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6Patch, Kp20PieceId, RectangularGregoryBiNonic, TensorBorderBi3,
    elevate_bi6_grid_to_bi9,
};

use super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::kp20_polar_bi6_split::split_bi6_patch_u;
use super::kp20_polar_expanded_support_diagnostics::{
    Kp20ExpandedSupportFoldLocalization, localize_expanded_support_folds,
};
use super::kp20_polar_expanded_support_lane_shape::reauthored_support_lanes;
use super::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
use super::kp20_polar_support_lift::Kp20SupportNormalLiftReport;
use super::kp20_polar_support_surface_split::support_owner_piece;
use crate::types::{CornerSubPatchFamily, CornerSubPatchProvenance};

#[path = "kp20_polar_expanded_support_fragment_rebuild.rs"]
mod expanded_support_fragment_rebuild;
#[path = "kp20_polar_expanded_support_lane_trial.rs"]
mod expanded_support_lane_trial;
#[path = "kp20_polar_expanded_support_mode.rs"]
mod expanded_support_mode;
#[path = "kp20_polar_expanded_support_owner_candidate.rs"]
mod expanded_support_owner_candidate;
#[path = "kp20_polar_support_cell_primitive.rs"]
mod support_cell_primitive;

pub(super) use expanded_support_mode::Kp20ExpandedSupportMode;
use support_cell_primitive::Kp20SupportCellPrimitiveReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20ExpandedSupportTopologyReport {
    pub(super) evaluated: usize,
    pub(super) accepted: usize,
    pub(super) rejected: usize,
    pub(super) start_folds: usize,
    pub(super) measured_folds: usize,
    pub(super) retained_macro_folds: usize,
    pub(super) expanded_patch_folds: usize,
    pub(super) target_piece: Kp20PieceId,
    pub(super) target_sector: Option<usize>,
    pub(super) emitted_patch_count: usize,
    pub(super) measured_patch_count: usize,
    pub(super) topology_face_count: usize,
    pub(super) topology_pentagon_count: usize,
    pub(super) topology_alignment: &'static str,
    pub(super) split_parameter: f64,
    pub(super) fragment_rebuild: bool,
    pub(super) support_lane_reauthor: &'static str,
    pub(super) support_lane_reauthor_scale: Option<f64>,
    pub(super) support_cell_primitive: Kp20SupportCellPrimitiveReport,
    pub(super) emitted_patches: Vec<GregoryPatch>,
    pub(super) emitted_provenance: Vec<CornerSubPatchProvenance>,
    pub(super) fold_localization: Kp20ExpandedSupportFoldLocalization,
    pub(super) decline_reason: &'static str,
}

pub(super) fn measure_expanded_support_topology(
    _tensor_border: &TensorBorderBi3,
    baseline: Kp20SupportNormalLiftReport,
) -> Result<
    (
        Kp20SupportNormalLiftReport,
        Kp20ExpandedSupportTopologyReport,
    ),
    String,
> {
    let report = expanded_support_owner_candidate::measure_support_owner_candidates(&baseline);
    Ok((baseline, report))
}

#[derive(Debug, Clone)]
struct Kp20ExpandedSupportTopologyTrial {
    measured_folds: usize,
    split_parameter: f64,
    fragment_rebuild: bool,
    reauthor: Option<(Kp20ExpandedSupportMode, f64)>,
    fold_localization: Kp20ExpandedSupportFoldLocalization,
}

impl Kp20ExpandedSupportTopologyTrial {
    fn better_than(&self, other: &Self) -> bool {
        (
            self.measured_folds,
            self.fold_localization.support_lane_folds(),
        ) < (
            other.measured_folds,
            other.fold_localization.support_lane_folds(),
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn expanded_support_topology_trial(
    macro_patches: &[popkern_core::geometry::gregory::Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    support_piece: Kp20PieceId,
    split_parameter: f64,
    fragment_rebuild: bool,
    reauthor: Option<(Kp20ExpandedSupportMode, f64)>,
    retained_macro_folds: usize,
    macro_fragment_folds: usize,
) -> Kp20ExpandedSupportTopologyTrial {
    let support_lanes = expanded_support_lane_trial::expanded_support_lane_patchwork(
        macro_patches,
        targets,
        support_piece,
        split_parameter,
        reauthor,
    );
    let fold_localization =
        localize_expanded_support_folds(&support_lanes.patches, &support_lanes.provenance);
    let expanded_patch_folds = macro_fragment_folds + fold_localization.total_folds();
    Kp20ExpandedSupportTopologyTrial {
        measured_folds: retained_macro_folds + expanded_patch_folds,
        split_parameter,
        fragment_rebuild,
        reauthor,
        fold_localization,
    }
}

fn fragmented_macro_fold_count(
    macro_patches: &[popkern_core::geometry::gregory::Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    split_parameter: f64,
    fragment_rebuild: bool,
) -> usize {
    let mut folds = 0;
    for target in targets {
        let patch = *macro_patches[target.sector].piece(Kp20PieceId::NearOriginLeft);
        for fragment in macro_fragments(patch, split_parameter, fragment_rebuild) {
            if gregory_from_bi6(fragment)
                .self_intersection_report(SELF_INT_TOL)
                .intersects()
            {
                folds += 1;
            }
        }
    }
    folds
}

fn macro_fragments(
    patch: Kp20Bi6Patch,
    split_parameter: f64,
    fragment_rebuild: bool,
) -> [Kp20Bi6Patch; 3] {
    let (first, third) = split_bi6_patch_u(&patch, split_parameter);
    let (first, second) = split_bi6_patch_u(&first, 0.5);
    let fragments = [first, second, third];
    if fragment_rebuild {
        expanded_support_fragment_rebuild::rebuild_macro_fragments_from_fixed_boundary(fragments)
    } else {
        fragments
    }
}

pub(super) fn push_support_lane(
    patches: &mut Vec<GregoryPatch>,
    provenance: &mut Vec<CornerSubPatchProvenance>,
    patch: Kp20Bi6Patch,
    target: Kp20FoldSite,
    sector_count: usize,
    patch_count: usize,
    boundary_lane: bool,
) {
    let patch_index = patches.len();
    patches.push(gregory_from_bi6(patch));
    provenance.push(CornerSubPatchProvenance::kp20_expanded_support_lane(
        target.sector,
        sector_count,
        patch_index,
        patch_count,
        boundary_lane,
    ));
}

pub(super) struct Kp20ExpandedSupportPatchwork {
    patches: Vec<GregoryPatch>,
    provenance: Vec<CornerSubPatchProvenance>,
}

fn support_lane_reauthor_label(reauthor: Option<(Kp20ExpandedSupportMode, f64)>) -> &'static str {
    match reauthor.map(|(mode, _)| mode) {
        Some(mode) => mode.diagnostic_label(),
        None => "raw_support_split",
    }
}

fn emitted_candidate_patchwork(
    macro_patches: &[popkern_core::geometry::gregory::Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    support_piece: Kp20PieceId,
    split_parameter: f64,
    fragment_rebuild: bool,
    reauthor: Option<(Kp20ExpandedSupportMode, f64)>,
) -> Kp20ExpandedSupportPatchwork {
    let sector_count = macro_patches.len();
    let patch_count = (macro_patches.len() * 4) + 3 * targets.len();
    let mut patches = Vec::with_capacity(patch_count);
    let mut provenance = Vec::with_capacity(patch_count);
    for (sector, macro_patch) in macro_patches.iter().enumerate() {
        for piece in Kp20PieceId::all() {
            if piece == support_piece
                && let Some(target) = targets.iter().find(|target| target.sector == sector)
            {
                let base = *macro_patch.piece(support_piece);
                let lanes = match reauthor {
                    Some((mode, scale)) => {
                        reauthored_support_lanes(&base, *target, split_parameter, mode, scale)
                    }
                    None => split_bi6_patch_u(&base, split_parameter),
                };
                push_support_lane(
                    &mut patches,
                    &mut provenance,
                    lanes.0,
                    *target,
                    sector_count,
                    patch_count,
                    true,
                );
                push_support_lane(
                    &mut patches,
                    &mut provenance,
                    lanes.1,
                    *target,
                    sector_count,
                    patch_count,
                    false,
                );
            } else {
                let patch = *macro_patch.piece(piece);
                if piece == Kp20PieceId::NearOriginLeft
                    && targets.iter().any(|target| target.sector == sector)
                {
                    for fragment in macro_fragments(patch, split_parameter, fragment_rebuild) {
                        push_macro_fragment(
                            &mut patches,
                            &mut provenance,
                            fragment,
                            sector,
                            sector_count,
                            patch_count,
                            piece,
                        );
                    }
                } else {
                    patches.push(gregory_from_bi6(patch));
                    provenance.push(CornerSubPatchProvenance::macro_cap(
                        CornerSubPatchFamily::Kp20Bi6,
                        sector,
                        sector_count,
                        piece.index(),
                    ));
                }
            }
        }
    }
    Kp20ExpandedSupportPatchwork {
        patches,
        provenance,
    }
}

fn push_macro_fragment(
    patches: &mut Vec<GregoryPatch>,
    provenance: &mut Vec<CornerSubPatchProvenance>,
    patch: Kp20Bi6Patch,
    sector: usize,
    sector_count: usize,
    patch_count: usize,
    piece: Kp20PieceId,
) {
    let patch_index = patches.len();
    patches.push(gregory_from_bi6(patch));
    let role = CornerSubPatchProvenance::macro_cap(
        CornerSubPatchFamily::Kp20Bi6,
        sector,
        sector_count,
        piece.index(),
    )
    .role;
    provenance.push(CornerSubPatchProvenance::kp20_expanded_macro_fragment(
        sector,
        sector_count,
        patch_index,
        patch_count,
        role,
    ));
}

fn gregory_from_bi6(patch: Kp20Bi6Patch) -> GregoryPatch {
    GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
        grid: elevate_bi6_grid_to_bi9(patch.cps),
    }))
}

fn retained_fold_count(
    summary: &Kp20SamplerFoldSummary,
    targets: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> usize {
    summary
        .fold_sites
        .iter()
        .filter(|site| {
            site.piece != target_piece || !targets.iter().any(|target| target.sector == site.sector)
        })
        .count()
}

fn support_sites_by_sector(
    summary: &Kp20SamplerFoldSummary,
    target_piece: Kp20PieceId,
) -> Vec<Kp20FoldSite> {
    let mut sites = Vec::new();
    for site in summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
    {
        if let Some(existing) = sites
            .iter_mut()
            .find(|existing: &&mut Kp20FoldSite| existing.sector == site.sector)
        {
            if site.min_jacobian_normal_cos < existing.min_jacobian_normal_cos {
                *existing = site;
            }
        } else {
            sites.push(site);
        }
    }
    sites
}

fn worst_support_site(
    summary: &Kp20SamplerFoldSummary,
    target_piece: Kp20PieceId,
) -> Option<Kp20FoldSite> {
    summary
        .fold_sites
        .iter()
        .copied()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
        .min_by(|a, b| {
            a.min_jacobian_normal_cos
                .partial_cmp(&b.min_jacobian_normal_cos)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}
