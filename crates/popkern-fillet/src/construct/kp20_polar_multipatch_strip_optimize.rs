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
    WorstEdgeColumnPositive,
    WorstEdgeColumnNegative,
    WorstEdgeParameterExpand,
    WorstEdgeParameterContract,
}

impl Kp20StripInteriorOptimizerMode {
    pub(super) fn all() -> [Self; 3] {
        [
            Self::FairBlend,
            Self::FoldNormalPositive,
            Self::FoldNormalNegative,
        ]
    }

    pub(super) fn exact_guided_all() -> [Self; 9] {
        [
            Self::FairBlend,
            Self::FoldNormalPositive,
            Self::FoldNormalNegative,
            Self::WorstFoldNormalPositive,
            Self::WorstFoldNormalNegative,
            Self::WorstEdgeColumnPositive,
            Self::WorstEdgeColumnNegative,
            Self::WorstEdgeParameterExpand,
            Self::WorstEdgeParameterContract,
        ]
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::FairBlend => "fair_blend",
            Self::FoldNormalPositive => "fold_normal_positive",
            Self::FoldNormalNegative => "fold_normal_negative",
            Self::WorstFoldNormalPositive => "worst_fold_normal_positive",
            Self::WorstFoldNormalNegative => "worst_fold_normal_negative",
            Self::WorstEdgeColumnPositive => "worst_edge_column_positive",
            Self::WorstEdgeColumnNegative => "worst_edge_column_negative",
            Self::WorstEdgeParameterExpand => "worst_edge_parameter_expand",
            Self::WorstEdgeParameterContract => "worst_edge_parameter_contract",
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
                    | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative
                    | Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive
                    | Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative
                    | Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
                    | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract => {
                        let fair = base + (target - base) * interior_scale;
                        let displacement = strip_mode_displacement(
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
            surface.cps[5][column] = shoulder_row_target(
                &accepted,
                &strip.surface,
                accepted_summary,
                strip,
                mode,
                column,
                span,
                interior_scale,
            );
            surface.cps[6][column] = strip.surface.cps[6][column];
        }
        reauthor_worst_edge_parameters(
            &mut surface,
            accepted_summary,
            strip,
            mode,
            span,
            interior_scale,
        );
    }
    surface
}

fn shoulder_row_target(
    accepted: &Kp20Bi6Patch,
    constructed: &Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    column: usize,
    span: f64,
    interior_scale: f64,
) -> popkern_core::geometry::Point3 {
    let base = accepted.cps[5][column];
    let target = constructed.cps[5][column];
    let fair = base + (target - base) * interior_scale.clamp(0.0, 1.0);
    match mode {
        Kp20StripInteriorOptimizerMode::FairBlend => fair,
        Kp20StripInteriorOptimizerMode::FoldNormalPositive
        | Kp20StripInteriorOptimizerMode::FoldNormalNegative
        | Kp20StripInteriorOptimizerMode::WorstFoldNormalPositive
        | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative
        | Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive
        | Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative
        | Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
        | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract => {
            let displacement =
                strip_mode_displacement(accepted, accepted_summary, strip, mode, 5, column);
            fair + displacement * (0.5 * span * interior_scale)
        }
    }
}

fn strip_mode_displacement(
    accepted: &Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    row: usize,
    column: usize,
) -> Vec3 {
    match mode {
        Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive
        | Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative
        | Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
        | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract => {
            worst_edge_column_displacement(accepted, accepted_summary, strip, mode, row, column)
        }
        _ => strip_fold_displacement(accepted, accepted_summary, strip, mode, row, column),
    }
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
            | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative
            | Kp20StripInteriorOptimizerMode::WorstEdgeColumnPositive
            | Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative
            | Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
            | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract => worst_site
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
        | Kp20StripInteriorOptimizerMode::WorstFoldNormalNegative
        | Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative
        | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract => -displacement,
        _ => displacement,
    }
}

fn reauthor_worst_edge_parameters(
    surface: &mut Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    span: f64,
    interior_scale: f64,
) {
    if !matches!(
        mode,
        Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand
            | Kp20StripInteriorOptimizerMode::WorstEdgeParameterContract
    ) {
        return;
    }
    let Some(site) = accepted_summary.worst_fold_site() else {
        return;
    };
    if site.sector != strip.sector_index || support_owner_piece(site.piece) != strip.strip_piece {
        return;
    }
    let target_column = if site.v <= 0.5 { 1_usize } else { 5_usize };
    let site_t = ((site.u - (2.0 / 6.0)) / (4.0 / 6.0)).clamp(0.0, 1.0);
    let expand = mode == Kp20StripInteriorOptimizerMode::WorstEdgeParameterExpand;
    for column in 0..7 {
        let column_weight = edge_column_weight(column, target_column);
        if column_weight <= 0.0 {
            continue;
        }
        let row_parameters = edge_reauthor_row_parameters(site_t, expand, column_weight);
        let anchor = surface.cps[2][column];
        let interface = surface.cps[6][column];
        for (idx, row) in (3..=5).enumerate() {
            let target = anchor + (interface - anchor) * row_parameters[idx];
            let displacement = capped_displacement(
                (target - surface.cps[row][column]) * interior_scale.clamp(0.0, 1.0),
                span * 0.18,
            );
            surface.cps[row][column] = surface.cps[row][column] + displacement;
        }
    }
}

fn edge_reauthor_row_parameters(site_t: f64, expand: bool, column_weight: f64) -> [f64; 3] {
    let base = [0.25_f64, 0.5, 0.75];
    let mut parameters = base;
    for parameter in &mut parameters {
        let distance = (*parameter - site_t).abs();
        let local = (1.0 - distance / 0.42).max(0.0);
        let direction = if *parameter <= site_t { -1.0 } else { 1.0 };
        let sign = if expand { direction } else { -direction };
        *parameter += sign * 0.14 * local * column_weight;
    }
    enforce_monotone_row_parameters(parameters)
}

fn enforce_monotone_row_parameters(mut parameters: [f64; 3]) -> [f64; 3] {
    parameters[0] = parameters[0].clamp(0.08, 0.84);
    parameters[1] = parameters[1].clamp(parameters[0] + 0.08, 0.92);
    parameters[2] = parameters[2].clamp(parameters[1] + 0.08, 0.96);
    parameters
}

fn worst_edge_column_displacement(
    accepted: &Kp20Bi6Patch,
    accepted_summary: &Kp20SamplerFoldSummary,
    strip: &Kp20MultiPatchStripSurface,
    mode: Kp20StripInteriorOptimizerMode,
    row: usize,
    column: usize,
) -> Vec3 {
    let Some(site) = accepted_summary.worst_fold_site() else {
        return Vec3::zero();
    };
    if site.sector != strip.sector_index || support_owner_piece(site.piece) != strip.strip_piece {
        return Vec3::zero();
    }
    let target_column = if site.v <= 0.5 { 1_usize } else { 5_usize };
    if column.abs_diff(target_column) > 1 {
        return Vec3::zero();
    }
    let control_u = row as f64 / 6.0;
    let row_weight = (1.0 - ((control_u - site.u).abs() / 0.35)).max(0.0);
    let column_weight = edge_column_weight(column, target_column);
    let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
    let normal = evaluate_bi6_normal(accepted, site.u, site.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let displacement = capped_displacement(normal * row_weight * column_weight * severity, 0.12);
    match mode {
        Kp20StripInteriorOptimizerMode::WorstEdgeColumnNegative => -displacement,
        _ => displacement,
    }
}

fn edge_column_weight(column: usize, target_column: usize) -> f64 {
    match column.abs_diff(target_column) {
        0 => 1.0,
        1 => 0.45,
        _ => 0.0,
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
#[path = "kp20_polar_multipatch_strip_optimize_tests.rs"]
mod tests;
