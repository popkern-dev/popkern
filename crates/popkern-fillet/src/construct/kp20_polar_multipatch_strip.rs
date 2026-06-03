//! Candidate-only KP20 support-strip surface builder.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId};

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport;

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
    let mut max_row_spacing_ratio = 0.0_f64;
    let mut max_row_curvature_ratio = 0.0_f64;
    for sector in &inventory.sector_reports {
        if sector.sector_index >= macro_patches.len() {
            continue;
        }
        for (idx, side) in sector.strip_sides.iter().copied().enumerate() {
            let mut surface = *macro_patches[sector.sector_index].piece(side);
            let reauthor =
                reauthor_strip_interior_rows(&mut surface, fold_summary, sector.sector_index, side);
            if reauthor.max_displacement > 0.0 {
                reauthored_strip_patch_count += 1;
                max_reauthored_interior_displacement =
                    max_reauthored_interior_displacement.max(reauthor.max_displacement);
            }
            min_parameter_step = min_parameter_step.min(reauthor.min_parameter_step);
            max_parameter_skew = max_parameter_skew.max(reauthor.max_parameter_skew);
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
        max_row_spacing_ratio,
        max_row_curvature_ratio,
        strip_control_point_count: strip_surfaces.len() * BI6_CONTROL_POINTS,
        decline_reason,
        strip_names,
        strip_surfaces,
    }
}

fn reauthor_strip_interior_rows(
    surface: &mut Kp20Bi6Patch,
    fold_summary: &Kp20SamplerFoldSummary,
    sector_index: usize,
    strip_piece: Kp20PieceId,
) -> Kp20StripReauthorReport {
    let original = *surface;
    let mut min_parameter_step = f64::INFINITY;
    let mut max_parameter_skew = 0.0_f64;
    for column in 0..7 {
        let support_adjacent = original.cps[1][column];
        let interface_adjacent = original.cps[5][column];
        let map = strip_column_parameter_map(fold_summary, sector_index, strip_piece, column);
        min_parameter_step = min_parameter_step.min(map.min_step());
        max_parameter_skew = max_parameter_skew.max(map.max_skew());
        for (idx, row) in (2..=4).enumerate() {
            let t = map.row_parameters[idx];
            surface.cps[row][column] =
                support_adjacent + (interface_adjacent - support_adjacent) * t;
        }
    }
    Kp20StripReauthorReport {
        max_displacement: max_interior_displacement(&original, surface),
        min_parameter_step,
        max_parameter_skew,
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20StripReauthorReport {
    max_displacement: f64,
    min_parameter_step: f64,
    max_parameter_skew: f64,
}

#[derive(Debug, Clone, Copy)]
struct Kp20StripParameterMap {
    row_parameters: [f64; 3],
}

fn strip_column_parameter_map(
    fold_summary: &Kp20SamplerFoldSummary,
    sector_index: usize,
    strip_piece: Kp20PieceId,
    column: usize,
) -> Kp20StripParameterMap {
    let control_v = column as f64 / 6.0;
    let mut weighted_u = 0.0_f64;
    let mut total_weight = 0.0_f64;
    for site in fold_summary.fold_sites.iter().filter(|site| {
        site.sector == sector_index
            && super::kp20_polar_support_surface_split::support_owner_piece(site.piece)
                == strip_piece
    }) {
        let distance = (control_v - site.v).abs();
        let weight = (1.0 - distance / 0.35).max(0.0).powi(2);
        let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
        weighted_u += site.u.clamp(0.0, 1.0) * weight * severity;
        total_weight += weight * severity;
    }
    if total_weight <= 1.0e-12 {
        return Kp20StripParameterMap::uniform();
    }
    let target_u = weighted_u / total_weight;
    let strength = (total_weight / (total_weight + 1.0)).clamp(0.0, 0.55);
    Kp20StripParameterMap::from_target(target_u, strength)
}

impl Kp20StripParameterMap {
    fn uniform() -> Self {
        Self {
            row_parameters: UNIFORM_ROW_PARAMETERS,
        }
    }

    fn from_target(target_u: f64, strength: f64) -> Self {
        let mut rows = UNIFORM_ROW_PARAMETERS;
        for (row, base) in rows.iter_mut().zip(UNIFORM_ROW_PARAMETERS) {
            *row = (base + (target_u - base) * strength).clamp(0.12, 0.88);
        }
        enforce_minimum_parameter_steps(&mut rows);
        Self {
            row_parameters: rows,
        }
    }

    fn min_step(self) -> f64 {
        let [a, b, c] = self.row_parameters;
        a.min(b - a).min(c - b).min(1.0 - c)
    }

    fn max_skew(self) -> f64 {
        self.row_parameters
            .into_iter()
            .zip(UNIFORM_ROW_PARAMETERS)
            .map(|(actual, uniform)| (actual - uniform).abs())
            .fold(0.0_f64, f64::max)
    }
}

fn enforce_minimum_parameter_steps(rows: &mut [f64; 3]) {
    rows[0] = rows[0].clamp(MIN_PARAMETER_STEP, 1.0 - 3.0 * MIN_PARAMETER_STEP);
    rows[1] = rows[1].clamp(rows[0] + MIN_PARAMETER_STEP, 1.0 - 2.0 * MIN_PARAMETER_STEP);
    rows[2] = rows[2].clamp(rows[1] + MIN_PARAMETER_STEP, 1.0 - MIN_PARAMETER_STEP);
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

fn max_interior_displacement(before: &Kp20Bi6Patch, after: &Kp20Bi6Patch) -> f64 {
    let mut max_displacement = 0.0_f64;
    for row in 2..=4 {
        for column in 0..7 {
            max_displacement =
                max_displacement.max(before.cps[row][column].distance_to(&after.cps[row][column]));
        }
    }
    max_displacement
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
const MIN_PARAMETER_STEP: f64 = 0.12;
const UNIFORM_ROW_PARAMETERS: [f64; 3] = [0.25, 0.5, 0.75];

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::{Kp20MacroPatch, Kp20PieceId};

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};
    use crate::construct::kp20_polar_multipatch_inventory::multipatch_inventory_from_fold_summary;

    #[test]
    fn strip_builder_reauthors_unconstrained_interior_rows() {
        let summary = Kp20SamplerFoldSummary {
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
        };
        let inventory = multipatch_inventory_from_fold_summary(&summary);
        let mut macro_patch = Kp20MacroPatch::zero();
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[2][2] = Point3::new(1.0, 2.0, 3.0);
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[1][2] = Point3::new(0.0, 0.0, 0.0);
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[5][2] = Point3::new(0.0, 0.0, 4.0);
        macro_patch.piece_mut(Kp20PieceId::NearMidpointLeft).cps[6][2] = Point3::new(0.0, 0.0, 5.0);

        let report = build_multipatch_strip_surfaces(&inventory, &[macro_patch], &summary);

        assert_eq!(report.attempted_strip_patch_count, 1);
        assert_eq!(report.built_strip_patch_count, 1);
        assert_eq!(report.reauthored_strip_patch_count, 1);
        assert!(report.max_reauthored_interior_displacement > 0.0);
        assert_eq!(report.finite_control_point_count, 49);
        assert_eq!(report.invalid_control_point_count, 0);
        assert!(report.max_row_spacing_ratio.is_finite());
        assert!(report.max_row_curvature_ratio.is_finite());
        assert_eq!(report.decline_reason, None);
        assert_eq!(report.strip_surfaces[0].sector_index, 0);
        assert_eq!(
            report.strip_surfaces[0].strip_piece,
            Kp20PieceId::NearMidpointLeft
        );
        assert!(report.strip_surfaces[0].surface.cps[2][2].z > 1.0);
        assert!(report.strip_surfaces[0].surface.cps[2][2].z < 2.0);
        assert_eq!(
            report.strip_surfaces[0].surface.cps[5][2],
            Point3::new(0.0, 0.0, 4.0)
        );
        assert_eq!(
            report.strip_surfaces[0].surface.cps[6][2],
            Point3::new(0.0, 0.0, 5.0)
        );
    }

    #[test]
    fn strip_parameter_map_is_monotone_and_fold_site_aware() {
        let summary = Kp20SamplerFoldSummary {
            fold_count: 1,
            min_jacobian_normal_cos: -0.9,
            piece_fold_counts: [1, 0, 0, 0],
            fold_details: Vec::new(),
            fold_sites: vec![Kp20FoldSite {
                sector: 0,
                piece: Kp20PieceId::NearOriginLeft,
                u: 0.85,
                v: 0.5,
                min_jacobian_ratio: -0.1,
                min_jacobian_normal_cos: -0.9,
            }],
        };

        let map = strip_column_parameter_map(&summary, 0, Kp20PieceId::NearMidpointLeft, 3);

        assert!(map.row_parameters[0] > UNIFORM_ROW_PARAMETERS[0]);
        assert!(map.min_step() >= MIN_PARAMETER_STEP - 1.0e-12);
        assert!(map.max_skew() > 0.0);
    }

    #[test]
    fn strip_parameter_map_uses_uniform_rows_without_owned_fold_sites() {
        let summary = Kp20SamplerFoldSummary {
            fold_count: 0,
            min_jacobian_normal_cos: 1.0,
            piece_fold_counts: [0; 4],
            fold_details: Vec::new(),
            fold_sites: Vec::new(),
        };

        let map = strip_column_parameter_map(&summary, 0, Kp20PieceId::NearMidpointLeft, 3);

        assert_eq!(map.row_parameters, UNIFORM_ROW_PARAMETERS);
        assert_eq!(map.max_skew(), 0.0);
    }
}
