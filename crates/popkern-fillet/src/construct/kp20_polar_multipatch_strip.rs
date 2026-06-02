//! Candidate-only KP20 support-strip surface builder.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchStripBuildReport {
    pub(super) attempted_strip_patch_count: usize,
    pub(super) built_strip_patch_count: usize,
    pub(super) finite_control_point_count: usize,
    pub(super) invalid_control_point_count: usize,
    pub(super) strip_control_point_count: usize,
    pub(super) decline_reason: Option<&'static str>,
    pub(super) strip_names: Vec<String>,
    pub(super) strip_surfaces: Vec<Kp20MultiPatchStripSurface>,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchStripSurface {
    pub(super) patch_name: String,
    pub(super) sector_index: usize,
    pub(super) strip_piece: Kp20PieceId,
    pub(super) surface: Kp20Bi6Patch,
}

pub(super) fn build_multipatch_strip_surfaces(
    inventory: &Kp20MultiPatchInventoryReport,
    macro_patches: &[Kp20MacroPatch],
) -> Kp20MultiPatchStripBuildReport {
    let mut strip_surfaces = Vec::new();
    let mut finite_control_point_count = 0_usize;
    let mut invalid_control_point_count = 0_usize;
    for sector in &inventory.sector_reports {
        if sector.sector_index >= macro_patches.len() {
            continue;
        }
        for (idx, side) in sector.strip_sides.iter().copied().enumerate() {
            let surface = *macro_patches[sector.sector_index].piece(side);
            let (finite, invalid) = count_control_points(&surface);
            finite_control_point_count += finite;
            invalid_control_point_count += invalid;
            strip_surfaces.push(Kp20MultiPatchStripSurface {
                patch_name: sector.strip_names[idx].clone(),
                sector_index: sector.sector_index,
                strip_piece: side,
                surface,
            });
        }
    }
    let strip_names = strip_surfaces
        .iter()
        .map(|strip| strip.patch_name.clone())
        .collect::<Vec<_>>();
    let decline_reason = if strip_surfaces.len() != inventory.attempted_strip_patch_count {
        Some("missing_strip_surface")
    } else if invalid_control_point_count > 0 {
        Some("non_finite_strip_control_point")
    } else {
        None
    };
    Kp20MultiPatchStripBuildReport {
        attempted_strip_patch_count: inventory.attempted_strip_patch_count,
        built_strip_patch_count: strip_surfaces.len(),
        finite_control_point_count,
        invalid_control_point_count,
        strip_control_point_count: strip_surfaces.len() * BI6_CONTROL_POINTS,
        decline_reason,
        strip_names,
        strip_surfaces,
    }
}

fn count_control_points(surface: &Kp20Bi6Patch) -> (usize, usize) {
    let mut finite = 0_usize;
    let mut invalid = 0_usize;
    for row in surface.cps {
        for point in row {
            if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
                finite += 1;
            } else {
                invalid += 1;
            }
        }
    }
    (finite, invalid)
}

const BI6_CONTROL_POINTS: usize = 49;

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;

    #[test]
    fn strip_builder_clones_inventory_owned_support_surfaces() {
        let inventory = multipatch_inventory_from_fold_summary(&Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.5,
                v: 0.5,
                min_jacobian_ratio: -0.1,
                min_jacobian_normal_cos: -0.1,
            }],
        });
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[2][2] = Point3::new(1.0, 2.0, 3.0);

        let report = build_multipatch_strip_surfaces(&inventory, &[macro_patch]);

        assert_eq!(report.attempted_strip_patch_count, 1);
        assert_eq!(report.built_strip_patch_count, 1);
        assert_eq!(report.finite_control_point_count, 49);
        assert_eq!(report.invalid_control_point_count, 0);
        assert_eq!(report.decline_reason, None);
        assert_eq!(report.strip_surfaces[0].sector_index, 0);
        assert_eq!(
            report.strip_surfaces[0].strip_piece,
            Kp20PieceId::NearMidpointLeft
        );
        assert_eq!(
            report.strip_surfaces[0].surface.cps[2][2],
            Point3::new(1.0, 2.0, 3.0)
        );
    }
}
