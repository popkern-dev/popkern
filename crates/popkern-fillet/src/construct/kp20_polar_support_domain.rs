//! Bounded support-domain basis primitive for KP20 support-owned folds.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20Bi6MacroCapInput, Kp20GuideBuildInput, Kp20MacroPatch, Kp20PatchCorner,
    Kp20PieceCornerMaps, Kp20PieceId, Kp20SectorPieceMaps, TensorBorderBi3,
    assemble_kp20_bi6_macro_cap, build_kp20_guide_ring, kp20_nonoverlap_coupled_sector_piece_maps,
    sample_kp20_macro_patch_from_piece_maps,
};

use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, is_better_summary, validation_is_clean,
};
use super::super::kp20_polar_support_reparameterization::Kp20SupportDomainReport;

const TARGET_PIECE: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

pub(super) fn support_domain_basis_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<(Vec<Kp20MacroPatch>, Kp20SupportDomainReport), String> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point,
        tensor_border: tensor_border.clone(),
    })
    .map_err(|e| e.to_string())?;
    let baseline_sampled = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_nonoverlap_coupled_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let baseline_cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches: baseline_sampled,
    })
    .map_err(|e| e.to_string())?;
    let baseline = candidate_summary_from_cap(tensor_border, baseline_cap.macro_patches.clone())?;
    let mut report = Kp20SupportDomainReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_folds: baseline.0.fold_count,
        end_folds: baseline.0.fold_count,
        target_piece: TARGET_PIECE,
        origin_u: 0.0,
        origin_v: 0.0,
        span_u: 0.0,
        span_v: 0.0,
        skew: 0.0,
        decline_reason: "support_domain_no_candidate",
    };
    let mut best: Option<(
        Vec<Kp20MacroPatch>,
        (Kp20SamplerFoldSummary, _),
        Kp20SupportDomainParams,
    )> = None;
    for params in Kp20SupportDomainParams::family() {
        report.evaluated += 1;
        let sampled = guide
            .sectors
            .iter()
            .map(|sector| {
                sample_kp20_macro_patch_from_piece_maps(sector, support_domain_piece_maps(params))
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
            report.decline_reason = "support_domain_validation_not_clean";
            continue;
        }
        if !is_better_summary(&summary, Some(&baseline)) {
            report.rejected += 1;
            report.decline_reason = "support_domain_no_fold_reduction";
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_summary, _)| is_better_summary(&summary, Some(best_summary)))
        {
            best = Some((cap.macro_patches, summary, params));
        } else {
            report.rejected += 1;
        }
    }
    let Some((macro_patches, summary, params)) = best else {
        report.rejected = report.evaluated;
        return Ok((baseline_cap.macro_patches, report));
    };
    report.accepted = 1;
    report.rejected = report.evaluated - report.accepted;
    report.end_folds = summary.0.fold_count;
    report.origin_u = params.origin[0];
    report.origin_v = params.origin[1];
    report.span_u = params.span_u;
    report.span_v = params.span_v;
    report.skew = params.skew;
    report.decline_reason = "none";
    Ok((macro_patches, report))
}

#[derive(Debug, Clone, Copy)]
struct Kp20SupportDomainParams {
    origin: [f64; 2],
    span_u: f64,
    span_v: f64,
    skew: f64,
}

impl Kp20SupportDomainParams {
    const fn family() -> [Self; 8] {
        [
            Self::new([0.00, 0.20], 0.30, 0.30, 0.00),
            Self::new([0.00, 0.22], 0.30, 0.28, 0.02),
            Self::new([0.02, 0.20], 0.28, 0.30, -0.02),
            Self::new([0.00, 0.18], 0.32, 0.30, 0.00),
            Self::new([0.00, 0.24], 0.28, 0.28, 0.04),
            Self::new([0.04, 0.20], 0.26, 0.30, -0.04),
            Self::new([0.00, 0.16], 0.34, 0.30, 0.02),
            Self::new([0.02, 0.22], 0.30, 0.26, -0.02),
        ]
    }

    const fn new(origin: [f64; 2], span_u: f64, span_v: f64, skew: f64) -> Self {
        Self {
            origin,
            span_u,
            span_v,
            skew,
        }
    }
}

fn support_domain_piece_maps(params: Kp20SupportDomainParams) -> Kp20SectorPieceMaps {
    let mut maps = kp20_nonoverlap_coupled_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        if *piece == TARGET_PIECE {
            *corner_maps = support_domain_corner_maps(params);
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
