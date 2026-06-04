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
    WorstFoldNormalPositive,
    WorstFoldNormalNegative,
}

impl Kp20StripInteriorOptimizerMode {
    pub(super) fn all() -> [Self; 3] {
        [
            Self::FairBlend,
            Self::FoldNormalPositive,
            Self::FoldNormalNegative,
        ]
    }

    pub(super) fn exact_guided_all() -> [Self; 5] {
        [
            Self::FairBlend,
            Self::FoldNormalPositive,
            Self::FoldNormalNegative,
            Self::WorstFoldNormalPositive,
            Self::WorstFoldNormalNegative,
        ]
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::FairBlend => "fair_blend",
            Self::FoldNormalPositive => "fold_normal_positive",
            Self::FoldNormalNegative => "fold_normal_negative",
            Self::WorstFoldNormalPositive => "worst_fold_normal_positive",
            Self::WorstFoldNormalNegative => "worst_fold_normal_negative",
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
    let mut surface = accepted_patches
        .get(strip.sector_index)
        .map(|patch| *patch.piece(strip.strip_piece))
        .unwrap_or(strip.surface);
    if let Some(accepted) = accepted_patches
        .get(strip.sector_index)
        .map(|patch| *patch.piece(strip.strip_piece))
    {
        let span = average_strip_span(&accepted);
        for row in 3..=4 {
            for column in 0..7 {
                let base = accepted.cps[row][column];
                let target = strip.surface.cps[row][column];
                surface.cps[row][column] = match mode {
                    Kp20StripInteriorOptimizerMode::FairBlend => {
                        base + (target - base) * interior_scale
                    }
                    Kp20StripInteriorOptimizerMode::FoldNormalPositive
                    | Kp20StripInteriorOptimizerMode::FoldNormalNegative
                    | Kp20StripInteriorOptimizerMode::WorstFoldNormalPositive
                    | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative => {
                        let fair = base + (target - base) * interior_scale;
                        let displacement = strip_fold_displacement(
                            &accepted,
                            accepted_summary,
                            strip,
                            mode,
                            row,
                            column,
                        );
                        fair + displacement * (span * interior_scale)
                    }
                };
            }
        }
        for column in 0..7 {
            surface.cps[5][column] =
                shoulder_row_target(&accepted, &strip.surface, column, interior_scale);
        }
    }
    surface
}

fn shoulder_row_target(
    accepted: &Kp20Bi6Patch,
    constructed: &Kp20Bi6Patch,
    column: usize,
    interior_scale: f64,
) -> popkern_core::geometry::Point3 {
    let base = accepted.cps[5][column];
    let target = constructed.cps[5][column];
    base + (target - base) * interior_scale.clamp(0.0, 1.0)
}

fn strip_fold_displacement(
    accepted: &Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    row: usize,
    column: usize,
) -> Vec3 {
    let u = row as f64 / 6.0;
    let v = column as f64 / 6.0;
    let mut displacement = Vec3::zero();
    let worst_site = accepted_summary.worst_fold_site();
    for site in accepted_summary.fold_sites.iter().filter(|site| {
        let strip_owned = site.sector == strip.sector_index
            && support_owner_piece(site.piece) == strip.strip_piece;
        let mode_owned = match mode {
            Kp20StripInteriorOptimizerMode::WorstFoldNormalPositive
            | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative => worst_site
                .is_some_and(|worst| worst.sector == site.sector && worst.piece == site.piece),
            _ => true,
        };
        strip_owned && mode_owned
    }) {
        let normal =
            evaluate_bi6_normal(accepted, site.u, site.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let weight = fold_site_weight(site.u, site.v, u, v, site.min_jacobian_normal_cos);
        displacement = displacement + normal * weight;
    }
    let displacement = capped_displacement(displacement, 0.08);
    match mode {
        Kp20StripInteriorOptimizerMode::FoldNormalNegative
        | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative => -displacement,
        _ => displacement,
    }
}

fn fold_site_weight(
    site_u: f64,
    site_v: f64,
    control_u: f64,
    control_v: f64,
    normal_cos: f64,
) -> f64 {
    let du = control_u - site_u;
    let dv = control_v - site_v;
    let r2 = du * du + dv * dv;
    let locality = (1.0 - r2 / (0.35_f64 * 0.35_f64)).max(0.0);
    let severity = (1.0 - normal_cos).clamp(0.0, 1.0);
    locality * locality * severity
}

fn capped_displacement(displacement: Vec3, max_length: f64) -> Vec3 {
    let length = displacement.length();
    if length > max_length {
        if let Some(direction) = displacement.normalized() {
            return direction * max_length;
        }
    }
    displacement
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
    fn strip_damping_preserves_tensor_boundary_and_interface_g0_row() {
        let mut accepted = Kp20MacroPatch::zero();
        let piece = accepted.piece_mut(Kp20PieceId::NearMidpointLeft);
        piece.cps[0][2] = Point3::new(0.0, 0.0, 0.0);
        piece.cps[1][2] = Point3::new(1.0, 0.0, 0.0);
        piece.cps[2][2] = Point3::new(2.0, 0.0, 0.0);
        piece.cps[3][2] = Point3::new(3.0, 0.0, 0.0);
        piece.cps[5][2] = Point3::new(5.0, 0.0, 0.0);
        piece.cps[6][2] = Point3::new(6.0, 0.0, 0.0);
        let mut strip_surface = *accepted.piece(Kp20PieceId::NearMidpointLeft);
        strip_surface.cps[2][2] = Point3::new(10.0, 0.0, 0.0);
        strip_surface.cps[3][2] = Point3::new(11.0, 0.0, 0.0);
        strip_surface.cps[5][2] = Point3::new(9.0, 0.0, 0.0);
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

        assert_eq!(scaled.cps[0][2], Point3::new(0.0, 0.0, 0.0));
        assert_eq!(scaled.cps[1][2], Point3::new(1.0, 0.0, 0.0));
        assert_eq!(scaled.cps[2][2], Point3::new(2.0, 0.0, 0.0));
        assert_eq!(scaled.cps[3][2], Point3::new(5.0, 0.0, 0.0));
        assert_eq!(scaled.cps[5][2], Point3::new(6.0, 0.0, 0.0));
        assert_eq!(scaled.cps[6][2], Point3::new(6.0, 0.0, 0.0));
    }

    #[test]
    fn fold_normal_mode_moves_unconstrained_rows_and_preserves_interface_g0() {
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

        assert_eq!(shifted.cps[2][3], Point3::new(2.0, 3.0, 0.0));
        assert_eq!(shifted.cps[1][3], Point3::new(1.0, 3.0, 0.0));
        assert!(shifted.cps[3][3].z > 0.0);
        assert_eq!(shifted.cps[5][3], Point3::new(5.0, 3.0, 0.0));
        assert_eq!(shifted.cps[6][3], Point3::new(6.0, 3.0, 0.0));
    }

    #[test]
    fn shoulder_row_target_is_bounded_to_constructed_row() {
        let mut accepted = Kp20Bi6Patch::zero();
        let mut constructed = Kp20Bi6Patch::zero();
        accepted.cps[5][4] = Point3::new(1.0, 0.0, 0.0);
        constructed.cps[5][4] = Point3::new(5.0, 0.0, 0.0);

        assert_eq!(
            shoulder_row_target(&accepted, &constructed, 4, 0.25),
            Point3::new(2.0, 0.0, 0.0)
        );
        assert_eq!(
            shoulder_row_target(&accepted, &constructed, 4, 3.0),
            Point3::new(5.0, 0.0, 0.0)
        );
    }

    #[test]
    fn fold_site_weight_is_local_and_severity_bounded() {
        let center_weight = fold_site_weight(0.5, 0.5, 0.5, 0.5, -1.0);
        let adjacent_weight = fold_site_weight(0.5, 0.5, 0.5, 0.65, -1.0);
        let distant_weight = fold_site_weight(0.5, 0.5, 0.5, 1.0, -1.0);
        let mild_weight = fold_site_weight(0.5, 0.5, 0.5, 0.5, 0.75);

        assert_eq!(center_weight, 1.0);
        assert!(adjacent_weight > 0.0);
        assert!(adjacent_weight < center_weight);
        assert_eq!(distant_weight, 0.0);
        assert_eq!(mild_weight, 0.25);
    }

    #[test]
    fn fold_site_displacement_is_capped() {
        let capped = capped_displacement(Vec3::new(1.0, 0.0, 0.0), 0.08);
        let uncapped = capped_displacement(Vec3::new(0.01, 0.02, 0.0), 0.08);

        assert!((capped.length() - 0.08).abs() < 1.0e-12);
        assert_eq!(uncapped, Vec3::new(0.01, 0.02, 0.0));
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
