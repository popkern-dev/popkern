//! Bounded KP20 support-strip interior candidate directions.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch};

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_multipatch_strip::Kp20MultiPatchStripSurface;
use super::kp20_polar_support_surface_split::{evaluate_bi6_normal, support_owner_piece};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kp20StripInteriorOptimizerMode {
    FairBlend,
    FoldNormalPositive,
    FoldNormalNegative,
}

impl Kp20StripInteriorOptimizerMode {
    pub(super) fn all() -> [Self; 3] {
        [
            Self::FairBlend,
            Self::FoldNormalPositive,
            Self::FoldNormalNegative,
        ]
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::FairBlend => "fair_blend",
            Self::FoldNormalPositive => "fold_normal_positive",
            Self::FoldNormalNegative => "fold_normal_negative",
        }
    }
}

pub(super) fn optimized_strip_surface(
    accepted_patches: &[Kp20MacroPatch],
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    interior_scale: f64,
) -> Kp20Bi6Patch {
    let mut surface = strip.surface;
    if let Some(accepted) = accepted_patches
        .get(strip.sector_index)
        .map(|patch| *patch.piece(strip.strip_piece))
    {
        let normal_direction =
            strip_fold_normal_direction(&accepted, accepted_summary, strip, mode);
        let span = average_strip_span(&accepted);
        for row in 2..=4 {
            for column in 0..7 {
                let base = accepted.cps[row][column];
                let target = strip.surface.cps[row][column];
                surface.cps[row][column] = match mode {
                    Kp20StripInteriorOptimizerMode::FairBlend => {
                        base + (target - base) * interior_scale
                    }
                    Kp20StripInteriorOptimizerMode::FoldNormalPositive
                    | Kp20StripInteriorOptimizerMode::FoldNormalNegative => {
                        let fair = base + (target - base) * interior_scale;
                        let weight =
                            strip_fold_interior_weight(accepted_summary, strip, row, column);
                        fair + normal_direction * (span * interior_scale * weight)
                    }
                };
            }
        }
    }
    surface
}

fn strip_fold_normal_direction(
    accepted: &Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
) -> Vec3 {
    let mut direction = Vec3::new(0.0, 0.0, 0.0);
    for site in accepted_summary.fold_sites.iter().filter(|site| {
        site.sector == strip.sector_index && support_owner_piece(site.piece) == strip.strip_piece
    }) {
        let normal =
            evaluate_bi6_normal(accepted, site.u, site.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let weight = (1.0 - site.min_jacobian_normal_cos).max(0.0);
        direction = direction + normal * weight;
    }
    let direction = direction.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    match mode {
        Kp20StripInteriorOptimizerMode::FoldNormalNegative => -direction,
        _ => direction,
    }
}

fn strip_fold_interior_weight(
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    row: usize,
    column: usize,
) -> f64 {
    let u = row as f64 / 6.0;
    let v = column as f64 / 6.0;
    let mut weight = 0.0_f64;
    for site in accepted_summary.fold_sites.iter().filter(|site| {
        site.sector == strip.sector_index && support_owner_piece(site.piece) == strip.strip_piece
    }) {
        let du = u - site.u;
        let dv = v - site.v;
        let r2 = du * du + dv * dv;
        let locality = (1.0 - r2 / (0.35_f64 * 0.35_f64)).max(0.0);
        let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
        weight = weight.max(locality * locality * severity);
    }
    0.08 * weight
}

fn average_strip_span(surface: &Kp20Bi6Patch) -> f64 {
    let support_width = surface.cps[1][0].distance_to(&surface.cps[1][6]);
    let interface_width = surface.cps[5][0].distance_to(&surface.cps[5][6]);
    let cross_depth = surface.cps[1][3].distance_to(&surface.cps[5][3]);
    ((support_width + interface_width + cross_depth) / 3.0).max(1.0e-9)
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};

    #[test]
    fn strip_damping_preserves_constrained_rows() {
        let mut accepted = Kp20MacroPatch::zero();
        let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
        piece.cps[1][2] = Point3::new(1.0, 0.0, 0.0);
        piece.cps[2][2] = Point3::new(2.0, 0.0, 0.0);
        piece.cps[5][2] = Point3::new(5.0, 0.0, 0.0);
        piece.cps[6][2] = Point3::new(6.0, 0.0, 0.0);
        let mut strip_surface = *accepted.piece(Kp20PieceId::NearMidpointLeft);
        strip_surface.cps[2][2] = Point3::new(10.0, 0.0, 0.0);
        let strip = Kp20MultiPatchStripSurface {
            patch_name: "strip".to_string(),
            sector_index: 0,
            strip_piece: Kp20PieceId::NearMidpointLeft,
            surface: strip_surface,
        };

        let scaled = optimized_strip_surface(
            &[accepted],
            &empty_fold_summary(),
            &strip,
            Kp20StripInteriorOptimizerMode::FairBlend,
            0.25,
        );

        assert_eq!(scaled.cps[1][2], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(scaled.cps[2][2], Point3::new(4.0, 0.0, 0.0));
        assert_eq!(scaled.cps[5][2], Point3::new(5.0, 0.0, 0.0));
        assert_eq!(scaled.cps[6][2], Point3::new(6.0, 0.0, 0.0));
    }

    #[test]
    fn fold_normal_mode_moves_only_unconstrained_rows() {
        let mut accepted = Kp20MacroPatch::zero();
        let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
        for row in 0..7 {
            for column in 0..7 {
                piece.cps[row][column] = Point3::new(row as f64, column as f64, 0.0);
            }
        }
        let strip_surface = *accepted.piece(Kp20PieceId::NearMidpointLeft);
        let strip = Kp20MultiPatchStripSurface {
            patch_name: "strip".to_string(),
            sector_index: 0,
            strip_piece: Kp20PieceId::NearMidpointLeft,
            surface: strip_surface,
        };

        let shifted = optimized_strip_surface(
            &[accepted],
            &single_fold_summary(),
            &strip,
            Kp20StripInteriorOptimizerMode::FoldNormalPositive,
            1.0,
        );

        assert_eq!(shifted.cps[1][3], Point3::new(1.0, 3.0, 0.0));
        assert!(shifted.cps[3][3].z > 0.0);
        assert_eq!(shifted.cps[5][3], Point3::new(5.0, 3.0, 0.0));
        assert_eq!(shifted.cps[6][3], Point3::new(6.0, 3.0, 0.0));
    }

    fn empty_fold_summary() -> Kp20SamplerFoldSummary {
        Kp20SamplerFoldSummary {
            fold_count: 0,
            min_jacobian_normal_cos: 1.0,
            piece_fold_counts: [0; 4],
            fold_details: Vec::new(),
            fold_sites: Vec::new(),
        }
    }

    fn single_fold_summary() -> Kp20SamplerFoldSummary {
        Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.5,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.5,
                v: 0.5,
                min_jacobian_ratio: -0.5,
                min_jacobian_normal_cos: -0.5,
            }],
        }
    }
}
