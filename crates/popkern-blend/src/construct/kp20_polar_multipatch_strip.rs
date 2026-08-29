//! Candidate-only KP20 support-strip surface builder.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport;
#[path = "kp20_polar_multipatch_strip_rows.rs"]
mod rows;
#[cfg(test)]
use rows::{
    Kp20StripColumnGeometry, Kp20StripRowPlacementStrategy, MIN_PARAMETER_STEP,
    UNIFORM_ROW_PARAMETERS, strip_column_parameter_map,
};
use rows::{Kp20StripRowStrategyCounts, reauthor_strip_rows, row_strategy_assessment};

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchStripBuildReport {
    pub(super) attempted_strip_patch_count: usize,
    pub(super) built_strip_patch_count: usize,
    pub(super) reauthored_strip_patch_count: usize,
    pub(super) finite_control_point_count: usize,
    pub(super) invalid_control_point_count: usize,
    pub(super) max_reauthored_interior_displacement: f64,
    pub(super) min_parameter_step: f64,
    pub(super) max_parameter_skew: f64,
    pub(super) max_tangent_step_imbalance: f64,
    pub(super) min_tangent_alignment: f64,
    pub(super) uniform_row_strategy_columns: usize,
    pub(super) fold_site_row_strategy_columns: usize,
    pub(super) tangent_balance_row_strategy_columns: usize,
    pub(super) combined_row_strategy_columns: usize,
    pub(super) row_strategy_assessment: &'static str,
    pub(super) tangent_row_reauthored_columns: usize,
    pub(super) max_tangent_row_displacement: f64,
    pub(super) max_row_spacing_ratio: f64,
    pub(super) max_row_curvature_ratio: f64,
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
    fold_summary: &Kp20SamplerFoldSummary,
) -> Kp20MultiPatchStripBuildReport {
    let mut strip_surfaces = Vec::new();
    let mut finite_control_point_count = 0_usize;
    let mut invalid_control_point_count = 0_usize;
    let mut reauthored_strip_patch_count = 0_usize;
    let mut max_reauthored_interior_displacement = 0.0_f64;
    let mut min_parameter_step = f64::INFINITY;
    let mut max_parameter_skew = 0.0_f64;
    let mut max_tangent_step_imbalance = 0.0_f64;
    let mut min_tangent_alignment = f64::INFINITY;
    let mut row_strategy_counts = Kp20StripRowStrategyCounts::default();
    let mut tangent_row_reauthored_columns = 0_usize;
    let mut max_tangent_row_displacement = 0.0_f64;
    let mut max_row_spacing_ratio = 0.0_f64;
    let mut max_row_curvature_ratio = 0.0_f64;
    for sector in &inventory.sector_reports {
        if sector.sector_index >= macro_patches.len() {
            continue;
        }
        for (idx, side) in sector.strip_sides.iter().copied().enumerate() {
            let mut surface = *macro_patches[sector.sector_index].piece(side);
            let reauthor =
                reauthor_strip_rows(&mut surface, fold_summary, sector.sector_index, side);
            if reauthor.max_interior_displacement > 0.0 {
                reauthored_strip_patch_count += 1;
                max_reauthored_interior_displacement =
                    max_reauthored_interior_displacement.max(reauthor.max_interior_displacement);
            }
            min_parameter_step = min_parameter_step.min(reauthor.min_parameter_step);
            max_parameter_skew = max_parameter_skew.max(reauthor.max_parameter_skew);
            max_tangent_step_imbalance =
                max_tangent_step_imbalance.max(reauthor.max_tangent_step_imbalance);
            min_tangent_alignment = min_tangent_alignment.min(reauthor.min_tangent_alignment);
            row_strategy_counts.add(reauthor.row_strategy_counts);
            tangent_row_reauthored_columns += reauthor.tangent_row_reauthored_columns;
            max_tangent_row_displacement =
                max_tangent_row_displacement.max(reauthor.max_tangent_row_displacement);
            let (finite, invalid) = count_control_points(&surface);
            let row_shape = row_shape_metrics(&surface);
            max_row_spacing_ratio = max_row_spacing_ratio.max(row_shape.max_spacing_ratio);
            max_row_curvature_ratio = max_row_curvature_ratio.max(row_shape.max_curvature_ratio);
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
        reauthored_strip_patch_count,
        finite_control_point_count,
        invalid_control_point_count,
        max_reauthored_interior_displacement,
        min_parameter_step: if min_parameter_step.is_finite() {
            min_parameter_step
        } else {
            0.0
        },
        max_parameter_skew,
        max_tangent_step_imbalance,
        min_tangent_alignment: if min_tangent_alignment.is_finite() {
            min_tangent_alignment
        } else {
            1.0
        },
        uniform_row_strategy_columns: row_strategy_counts.uniform,
        fold_site_row_strategy_columns: row_strategy_counts.fold_site,
        tangent_balance_row_strategy_columns: row_strategy_counts.tangent_balance,
        combined_row_strategy_columns: row_strategy_counts.combined,
        row_strategy_assessment: row_strategy_assessment(
            row_strategy_counts,
            max_tangent_step_imbalance,
            min_tangent_alignment,
        ),
        tangent_row_reauthored_columns,
        max_tangent_row_displacement,
        max_row_spacing_ratio,
        max_row_curvature_ratio,
        strip_control_point_count: strip_surfaces.len() * BI6_CONTROL_POINTS,
        decline_reason,
        strip_names,
        strip_surfaces,
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20StripRowShapeMetrics {
    max_spacing_ratio: f64,
    max_curvature_ratio: f64,
}

fn row_shape_metrics(surface: &Kp20Bi6Patch) -> Kp20StripRowShapeMetrics {
    let mut min_spacing = f64::INFINITY;
    let mut max_spacing = 0.0_f64;
    let mut max_curvature = 0.0_f64;
    for column in 0..7 {
        for row in 1..5 {
            let spacing = surface.cps[row][column].distance_to(&surface.cps[row + 1][column]);
            min_spacing = min_spacing.min(spacing);
            max_spacing = max_spacing.max(spacing);
        }
        for row in 2..=4 {
            let mid = surface.cps[row][column];
            let prev = surface.cps[row - 1][column];
            let next = surface.cps[row + 1][column];
            let fair = prev + (next - prev) * 0.5;
            max_curvature = max_curvature.max(mid.distance_to(&fair));
        }
    }
    Kp20StripRowShapeMetrics {
        max_spacing_ratio: max_spacing / min_spacing.max(1.0e-12),
        max_curvature_ratio: max_curvature / max_spacing.max(1.0e-12),
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
#[path = "kp20_polar_multipatch_strip_tests.rs"]
mod tests;
