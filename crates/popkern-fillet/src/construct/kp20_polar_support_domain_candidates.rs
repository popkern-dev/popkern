//! Bounded support-domain candidate families for KP20 sampling basis repair.

use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideRing,
    Kp20MacroPatch, Kp20PatchCorner, Kp20PieceCornerMaps, Kp20PieceId, Kp20SectorPieceMaps,
    TensorBorderBi3, assemble_kp20_bi6_macro_cap, kp20_nonoverlap_coupled_sector_piece_maps,
    sample_kp20_macro_patch_from_piece_maps,
};

use super::super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, is_better_summary, validation_is_clean,
};
use super::super::super::kp20_polar_support_reparameterization::{
    Kp20CoupledDomainReport, Kp20SplitDomainReport,
};
use super::{Kp20SupportDomainParams, TARGET_PIECE};

pub(super) fn refine_split_domain(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    accepted: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    best: &mut Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20SupportDomainParams,
    )>,
    anchor: Kp20SupportDomainParams,
) -> Result<Kp20SplitDomainReport, String> {
    let mut report = Kp20SplitDomainReport {
        guide_samples: 0,
        guide_reversals: 0,
        guide_min_normal_cos: 1.0,
        guide_worst_sector: 0,
        guide_worst_u: 0.0,
        guide_worst_v: 0.0,
        guide_classification: "guide_quality_not_measured",
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: accepted.0.fold_count,
        end_folds: accepted.0.fold_count,
        support_origin_u: anchor.origin[0],
        support_origin_v: anchor.origin[1],
        support_span_u: anchor.span_u,
        support_span_v: anchor.span_v,
        support_skew: anchor.skew,
        split_fraction: 0.5,
        interface_pull: 0.0,
        decline_reason: "split_domain_no_candidate",
        assembly_evaluated: 0,
        assembly_accepted: 0,
        assembly_rejected: 0,
        assembly_start_folds: accepted.0.fold_count,
        assembly_end_folds: accepted.0.fold_count,
        assembly_normal_scale: 0.0,
        assembly_radial_scale: 0.0,
        assembly_decline_reason: "preassembly_interior_not_run",
    };
    for params in Kp20SplitDomainParams::family(anchor) {
        report.evaluated += 1;
        let sampled = guide
            .sectors
            .iter()
            .map(|sector| {
                sample_kp20_macro_patch_from_piece_maps(sector, split_domain_piece_maps(params))
            })
            .collect::<Vec<_>>();
        let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border: tensor_border.clone(),
            sampled_macro_patches: sampled,
        })
        .map_err(|e| e.to_string())?;
        let summary = candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
        if !validation_is_clean(summary.1) {
            report.rejected += 1;
            report.decline_reason = "split_domain_validation_not_clean";
            continue;
        }
        if !is_better_summary(&summary, Some(accepted)) {
            report.rejected += 1;
            report.decline_reason = "split_domain_no_fold_reduction";
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_summary, _)| is_better_summary(&summary, Some(best_summary)))
        {
            *best = Some((cap.macro_patches, summary, params.support));
            report.accepted = 1;
            report.support_origin_u = params.support.origin[0];
            report.support_origin_v = params.support.origin[1];
            report.support_span_u = params.support.span_u;
            report.support_span_v = params.support.span_v;
            report.support_skew = params.support.skew;
            report.split_fraction = params.split_fraction;
            report.interface_pull = params.interface_pull;
            report.decline_reason = "none";
        } else {
            report.rejected += 1;
        }
    }
    report.rejected = report.evaluated - report.accepted;
    Ok(report)
}

pub(super) fn refine_coupled_domain(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    accepted: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    best: &mut Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
        Kp20SupportDomainParams,
    )>,
    anchor: Kp20SupportDomainParams,
) -> Result<Kp20CoupledDomainReport, String> {
    let mut report = Kp20CoupledDomainReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: accepted.0.fold_count,
        end_folds: accepted.0.fold_count,
        support_origin_u: anchor.origin[0],
        support_origin_v: anchor.origin[1],
        support_span_u: anchor.span_u,
        support_span_v: anchor.span_v,
        support_skew: anchor.skew,
        origin_shift: 0.0,
        origin_scale: 1.0,
        decline_reason: "coupled_domain_no_candidate",
    };
    for params in Kp20CoupledDomainParams::family(anchor) {
        report.evaluated += 1;
        let sampled = guide
            .sectors
            .iter()
            .map(|sector| {
                sample_kp20_macro_patch_from_piece_maps(sector, coupled_domain_piece_maps(params))
            })
            .collect::<Vec<_>>();
        let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border: tensor_border.clone(),
            sampled_macro_patches: sampled,
        })
        .map_err(|e| e.to_string())?;
        let summary = candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
        if !validation_is_clean(summary.1) {
            report.rejected += 1;
            report.decline_reason = "coupled_domain_validation_not_clean";
            continue;
        }
        if !is_better_summary(&summary, Some(accepted)) {
            report.rejected += 1;
            report.decline_reason = "coupled_domain_no_fold_reduction";
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_summary, _)| is_better_summary(&summary, Some(best_summary)))
        {
            *best = Some((cap.macro_patches, summary, params.support));
            report.accepted = 1;
            report.support_origin_u = params.support.origin[0];
            report.support_origin_v = params.support.origin[1];
            report.support_span_u = params.support.span_u;
            report.support_span_v = params.support.span_v;
            report.support_skew = params.support.skew;
            report.origin_shift = params.origin_shift;
            report.origin_scale = params.origin_scale;
            report.decline_reason = "none";
        } else {
            report.rejected += 1;
        }
    }
    report.rejected = report.evaluated - report.accepted;
    Ok(report)
}

#[derive(Debug, Clone, Copy)]
struct Kp20CoupledDomainParams {
    support: Kp20SupportDomainParams,
    origin_shift: f64,
    origin_scale: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20SplitDomainParams {
    support: Kp20SupportDomainParams,
    split_fraction: f64,
    interface_pull: f64,
}

impl Kp20SplitDomainParams {
    fn family(anchor: Kp20SupportDomainParams) -> [Self; 10] {
        [
            Self::new(anchor, 0.42, -0.02),
            Self::new(anchor, 0.42, 0.02),
            Self::new(anchor, 0.50, -0.02),
            Self::new(anchor, 0.50, 0.02),
            Self::new(anchor, 0.58, -0.02),
            Self::new(anchor, 0.58, 0.02),
            Self::new(anchor.shifted(0.0, -0.01, 0.0, 0.01, 0.0), 0.46, -0.03),
            Self::new(anchor.shifted(0.0, 0.01, 0.0, -0.01, 0.0), 0.54, 0.03),
            Self::new(anchor.shifted(0.01, 0.0, -0.01, 0.0, 0.01), 0.46, 0.0),
            Self::new(anchor.shifted(-0.01, 0.0, 0.01, 0.0, -0.01), 0.54, 0.0),
        ]
    }

    const fn new(
        support: Kp20SupportDomainParams,
        split_fraction: f64,
        interface_pull: f64,
    ) -> Self {
        Self {
            support,
            split_fraction,
            interface_pull,
        }
    }
}

impl Kp20CoupledDomainParams {
    fn family(anchor: Kp20SupportDomainParams) -> [Self; 8] {
        [
            Self::new(anchor, -0.015, 0.92),
            Self::new(anchor, 0.015, 0.92),
            Self::new(anchor, -0.015, 1.08),
            Self::new(anchor, 0.015, 1.08),
            Self::new(anchor.shifted(0.0, -0.01, 0.01, 0.0, 0.0), -0.02, 1.0),
            Self::new(anchor.shifted(0.0, 0.01, -0.01, 0.0, 0.0), 0.02, 1.0),
            Self::new(anchor.shifted(0.01, 0.0, 0.0, -0.01, 0.01), -0.01, 1.05),
            Self::new(anchor.shifted(-0.01, 0.0, 0.0, 0.01, -0.01), 0.01, 0.95),
        ]
    }

    const fn new(support: Kp20SupportDomainParams, origin_shift: f64, origin_scale: f64) -> Self {
        Self {
            support,
            origin_shift,
            origin_scale,
        }
    }
}

fn coupled_domain_piece_maps(params: Kp20CoupledDomainParams) -> Kp20SectorPieceMaps {
    let mut maps = kp20_nonoverlap_coupled_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        match *piece {
            TARGET_PIECE => *corner_maps = support_domain_corner_maps(params.support),
            Kp20PieceId::NearOriginLeft => {
                *corner_maps = origin_left_corner_maps(params.origin_shift, params.origin_scale);
            }
            Kp20PieceId::NearOriginRight => {
                *corner_maps = origin_right_corner_maps(params.origin_shift, params.origin_scale);
            }
            Kp20PieceId::NearMidpointRight => {}
        }
    }
    maps
}

fn split_domain_piece_maps(params: Kp20SplitDomainParams) -> Kp20SectorPieceMaps {
    let mut maps = kp20_nonoverlap_coupled_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        if *piece == TARGET_PIECE {
            *corner_maps = split_support_corner_maps(params);
        }
    }
    maps
}

fn support_domain_corner_maps(params: Kp20SupportDomainParams) -> Kp20PieceCornerMaps {
    square_maps(
        params.origin,
        [params.span_u, params.skew],
        [params.skew, params.span_v],
    )
}

fn split_support_corner_maps(params: Kp20SplitDomainParams) -> Kp20PieceCornerMaps {
    let support = params.support;
    let split = params.split_fraction.clamp(0.35, 0.65);
    let pull = params.interface_pull.clamp(-0.04, 0.04);
    let origin = support.origin;
    let du = [support.span_u, support.skew];
    let dv = [support.skew, support.span_v];
    let lower = scale2(dv, split);
    let upper = scale2(dv, 1.0 - split);
    let interface = add2(add2(origin, lower), [pull, -pull]);
    [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new(origin, du, lower),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new(add2(origin, du), scale2(du, -1.0), lower),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new(interface, du, upper),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new(add2(interface, du), scale2(du, -1.0), upper),
        ),
    ]
}

fn origin_left_corner_maps(shift: f64, scale: f64) -> Kp20PieceCornerMaps {
    let span = (0.25 * scale).clamp(0.20, 0.30);
    let origin = [
        (0.25 + shift).clamp(0.20, 0.30),
        (0.25 + shift).clamp(0.20, 0.30),
    ];
    square_maps(origin, [span, 0.0], [0.0, span])
}

fn origin_right_corner_maps(shift: f64, scale: f64) -> Kp20PieceCornerMaps {
    let span = (0.25 * scale).clamp(0.20, 0.30);
    let origin = [
        (0.50 - shift).clamp(0.45, 0.55),
        (0.25 + shift).clamp(0.20, 0.30),
    ];
    square_maps(origin, [-span, 0.0], [0.0, -span])
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
