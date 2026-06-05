//! Pre-assembly support-domain construction candidates for KP20.

use popkern_core::geometry::gregory::{
    Kp20AffineGuideMap, Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideRing,
    Kp20MacroPatch, Kp20PatchCorner, Kp20PieceCornerMaps, Kp20SectorPieceMaps, TensorBorderBi3,
    assemble_kp20_bi6_macro_cap, kp20_nonoverlap_coupled_sector_piece_maps,
    sample_kp20_macro_patch_from_piece_maps,
};
use popkern_core::geometry::{Point3, Vec3};

use super::super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, is_better_summary, validation_is_clean,
};
use super::super::super::kp20_polar_support_reparameterization::Kp20SplitDomainReport;
use super::{Kp20SupportDomainParams, TARGET_PIECE};

type Kp20SupportDomainCandidate = (
    Vec<Kp20MacroPatch>,
    (Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    Kp20SupportDomainParams,
);

pub(super) fn refine_preassembly_interior(
    star_point: Point3,
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    accepted: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
    best: &mut Option<Kp20SupportDomainCandidate>,
    anchor: Kp20SupportDomainParams,
    report: &mut Kp20SplitDomainReport,
) -> Result<(), String> {
    report.assembly_start_folds = accepted.0.fold_count;
    report.assembly_end_folds = accepted.0.fold_count;
    report.assembly_decline_reason = "preassembly_interior_no_candidate";
    for params in Kp20PreassemblyInteriorParams::family() {
        let mut sampled = sample_with_support_layout(guide, anchor, params.layout);
        apply_preassembly_interior_field(&mut sampled, star_point, params);
        report.assembly_evaluated += 1;
        let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
            tensor_border: tensor_border.clone(),
            sampled_macro_patches: sampled,
        })
        .map_err(|e| e.to_string())?;
        let summary = candidate_summary_from_cap(tensor_border, cap.macro_patches.clone())?;
        if !validation_is_clean(summary.1) {
            report.assembly_rejected += 1;
            report.assembly_decline_reason = "preassembly_interior_validation_not_clean";
            continue;
        }
        if !is_better_summary(&summary, Some(accepted)) {
            report.assembly_rejected += 1;
            report.assembly_decline_reason = params.decline_reason();
            continue;
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_summary, _)| is_better_summary(&summary, Some(best_summary)))
        {
            *best = Some((cap.macro_patches, summary, anchor));
            report.assembly_accepted = 1;
            report.assembly_normal_scale = params.normal_scale;
            report.assembly_radial_scale = params.radial_scale;
            report.assembly_decline_reason = "none";
        } else {
            report.assembly_rejected += 1;
        }
    }
    report.assembly_rejected = report.assembly_evaluated - report.assembly_accepted;
    if let Some((_, summary, _)) = best.as_ref() {
        report.assembly_end_folds = summary.0.fold_count;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct Kp20PreassemblyInteriorParams {
    normal_scale: f64,
    radial_scale: f64,
    layout: Kp20SupportLayout,
}

#[derive(Debug, Clone, Copy)]
enum Kp20SupportLayout {
    Anchor,
    ShoulderWide,
    ShoulderTall,
    OriginPulled,
    MidpointPulled,
    BoundaryFixedBendLeft,
    BoundaryFixedBendRight,
    BoundaryFixedFanOut,
    BoundaryFixedFanIn,
}

impl Kp20PreassemblyInteriorParams {
    fn family() -> [Self; 24] {
        [
            Self::field(0.010, 0.000),
            Self::field(-0.010, 0.000),
            Self::field(0.020, 0.000),
            Self::field(-0.020, 0.000),
            Self::field(0.010, 0.010),
            Self::field(0.010, -0.010),
            Self::field(-0.010, 0.010),
            Self::field(-0.010, -0.010),
            Self::field(0.020, 0.015),
            Self::field(0.020, -0.015),
            Self::field(-0.020, 0.015),
            Self::field(-0.020, -0.015),
            Self::layout(Kp20SupportLayout::ShoulderWide),
            Self::layout(Kp20SupportLayout::ShoulderTall),
            Self::layout(Kp20SupportLayout::OriginPulled),
            Self::layout(Kp20SupportLayout::MidpointPulled),
            Self::new(0.010, 0.000, Kp20SupportLayout::ShoulderWide),
            Self::new(-0.010, 0.000, Kp20SupportLayout::OriginPulled),
            Self::layout(Kp20SupportLayout::BoundaryFixedBendLeft),
            Self::layout(Kp20SupportLayout::BoundaryFixedBendRight),
            Self::layout(Kp20SupportLayout::BoundaryFixedFanOut),
            Self::layout(Kp20SupportLayout::BoundaryFixedFanIn),
            Self::new(0.010, 0.000, Kp20SupportLayout::BoundaryFixedBendLeft),
            Self::new(-0.010, 0.000, Kp20SupportLayout::BoundaryFixedFanOut),
        ]
    }

    const fn field(normal_scale: f64, radial_scale: f64) -> Self {
        Self::new(normal_scale, radial_scale, Kp20SupportLayout::Anchor)
    }

    const fn layout(layout: Kp20SupportLayout) -> Self {
        Self::new(0.0, 0.0, layout)
    }

    const fn new(normal_scale: f64, radial_scale: f64, layout: Kp20SupportLayout) -> Self {
        Self {
            normal_scale,
            radial_scale,
            layout,
        }
    }

    fn decline_reason(self) -> &'static str {
        match self.layout {
            Kp20SupportLayout::Anchor => "preassembly_interior_no_fold_reduction",
            Kp20SupportLayout::BoundaryFixedBendLeft
            | Kp20SupportLayout::BoundaryFixedBendRight
            | Kp20SupportLayout::BoundaryFixedFanOut
            | Kp20SupportLayout::BoundaryFixedFanIn => "preassembly_warp_no_fold_reduction",
            _ => "preassembly_layout_no_fold_reduction",
        }
    }
}

fn sample_with_support_layout(
    guide: &Kp20GuideRing,
    anchor: Kp20SupportDomainParams,
    layout: Kp20SupportLayout,
) -> Vec<Kp20MacroPatch> {
    guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                support_layout_piece_maps(anchor, layout),
            )
        })
        .collect()
}

fn apply_preassembly_interior_field(
    macro_patches: &mut [Kp20MacroPatch],
    star_point: Point3,
    params: Kp20PreassemblyInteriorParams,
) {
    for macro_patch in macro_patches {
        let piece = macro_patch.piece_mut(TARGET_PIECE);
        let support_span = (piece.cps[6][3] - piece.cps[2][3])
            .length()
            .max(piece.cps[4][6].distance_to(&piece.cps[4][0]));
        let normal = support_normal(piece.cps).unwrap_or(Vec3::zero());
        for row in 4..=6 {
            for col in 1..=5 {
                let radial = (star_point - piece.cps[row][col])
                    .normalized()
                    .unwrap_or(Vec3::zero());
                let row_weight = match row {
                    4 => 0.45,
                    5 => 0.70,
                    _ => 1.00,
                };
                let col_weight = 1.0 - ((col as f64 - 3.0).abs() / 3.0);
                let displacement = normal * params.normal_scale + radial * params.radial_scale;
                piece.cps[row][col] =
                    piece.cps[row][col] + displacement * (support_span * row_weight * col_weight);
            }
        }
    }
}

fn support_normal(cps: [[Point3; 7]; 7]) -> Option<Vec3> {
    let along_depth = cps[6][3] - cps[2][3];
    let across_sector = cps[4][6] - cps[4][0];
    along_depth.cross(&across_sector).normalized()
}

fn support_layout_piece_maps(
    params: Kp20SupportDomainParams,
    layout: Kp20SupportLayout,
) -> Kp20SectorPieceMaps {
    let mut maps = kp20_nonoverlap_coupled_sector_piece_maps();
    for (piece, corner_maps) in &mut maps {
        if *piece == TARGET_PIECE {
            *corner_maps = support_layout_corner_maps(params, layout);
        }
    }
    maps
}

fn support_layout_corner_maps(
    params: Kp20SupportDomainParams,
    layout: Kp20SupportLayout,
) -> Kp20PieceCornerMaps {
    let origin = params.origin;
    let du = [params.span_u, params.skew];
    let dv = [params.skew, params.span_v];
    if matches!(
        layout,
        Kp20SupportLayout::BoundaryFixedBendLeft
            | Kp20SupportLayout::BoundaryFixedBendRight
            | Kp20SupportLayout::BoundaryFixedFanOut
            | Kp20SupportLayout::BoundaryFixedFanIn
    ) {
        return boundary_fixed_warp_corner_maps(origin, du, dv, layout);
    }
    let shoulder = match layout {
        Kp20SupportLayout::ShoulderWide => add2(add2(origin, dv), scale2(du, 0.16)),
        Kp20SupportLayout::ShoulderTall => add2(add2(origin, dv), scale2(dv, 0.14)),
        Kp20SupportLayout::OriginPulled => add2(add2(origin, dv), scale2(du, -0.08)),
        Kp20SupportLayout::MidpointPulled => add2(add2(origin, dv), scale2(du, 0.08)),
        Kp20SupportLayout::Anchor => add2(origin, dv),
        Kp20SupportLayout::BoundaryFixedBendLeft
        | Kp20SupportLayout::BoundaryFixedBendRight
        | Kp20SupportLayout::BoundaryFixedFanOut
        | Kp20SupportLayout::BoundaryFixedFanIn => unreachable!("warp layouts return above"),
    };
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
            Kp20AffineGuideMap::new(shoulder, du, scale2(dv, -1.0)),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new(add2(shoulder, du), scale2(du, -1.0), scale2(dv, -1.0)),
        ),
    ]
}

fn boundary_fixed_warp_corner_maps(
    origin: [f64; 2],
    du: [f64; 2],
    dv: [f64; 2],
    layout: Kp20SupportLayout,
) -> Kp20PieceCornerMaps {
    let (along, across) = match layout {
        Kp20SupportLayout::BoundaryFixedBendLeft => (0.18, -0.08),
        Kp20SupportLayout::BoundaryFixedBendRight => (-0.18, 0.08),
        Kp20SupportLayout::BoundaryFixedFanOut => (0.10, 0.14),
        Kp20SupportLayout::BoundaryFixedFanIn => (-0.10, -0.14),
        _ => (0.0, 0.0),
    };
    [
        (
            Kp20PatchCorner::U0V0,
            Kp20AffineGuideMap::new(
                origin,
                add2(du, scale2(dv, across)),
                add2(dv, scale2(du, along)),
            ),
        ),
        (
            Kp20PatchCorner::U1V0,
            Kp20AffineGuideMap::new(
                add2(origin, du),
                add2(scale2(du, -1.0), scale2(dv, -across)),
                add2(dv, scale2(du, -along)),
            ),
        ),
        (
            Kp20PatchCorner::U0V1,
            Kp20AffineGuideMap::new(
                add2(origin, dv),
                add2(du, scale2(dv, -across)),
                add2(scale2(dv, -1.0), scale2(du, -along)),
            ),
        ),
        (
            Kp20PatchCorner::U1V1,
            Kp20AffineGuideMap::new(
                add2(add2(origin, du), dv),
                add2(scale2(du, -1.0), scale2(dv, across)),
                add2(scale2(dv, -1.0), scale2(du, along)),
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
