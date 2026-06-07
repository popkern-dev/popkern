//! Exact-boundary support-cell primitive for KP20 folded support patches.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId, RectangularGregoryBiNonic,
    elevate_bi6_grid_to_bi9,
};
use popkern_core::geometry::{Point3, Vec3};

use super::super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::super::kp20_polar_bi6_split::{split_bi6_patch_u, split_bi6_patch_v};
use super::super::kp20_polar_fold_summary::Kp20FoldSite;
use super::super::kp20_polar_support_surface_split::support_owner_piece;
use crate::types::CornerSubPatchProvenance;

#[derive(Debug, Clone)]
pub(crate) struct Kp20SupportCellPrimitiveReport {
    pub(crate) evaluated: usize,
    pub(crate) target_count: usize,
    pub(crate) patch_count: usize,
    pub(crate) folded_patch_count: usize,
    pub(crate) cells_u: usize,
    pub(crate) cells_v: usize,
    pub(crate) construction_mode: &'static str,
    pub(crate) target_piece: Kp20PieceId,
    pub(crate) emitted_patches: Vec<GregoryPatch>,
    pub(crate) emitted_provenance: Vec<CornerSubPatchProvenance>,
    pub(crate) decline_reason: &'static str,
}

impl Default for Kp20SupportCellPrimitiveReport {
    fn default() -> Self {
        Self {
            evaluated: 0,
            target_count: 0,
            patch_count: 0,
            folded_patch_count: 0,
            cells_u: 0,
            cells_v: 0,
            construction_mode: "none",
            target_piece: Kp20PieceId::NearMidpointLeft,
            emitted_patches: Vec::new(),
            emitted_provenance: Vec::new(),
            decline_reason: "support_cell_primitive_not_evaluated",
        }
    }
}

pub(super) fn measure_support_cell_primitive(
    macro_patches: &[Kp20MacroPatch],
    fold_sites: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> Kp20SupportCellPrimitiveReport {
    let targets = support_targets_by_sector(fold_sites, target_piece);
    let mut report = Kp20SupportCellPrimitiveReport {
        evaluated: 0,
        target_count: targets.len(),
        patch_count: 0,
        folded_patch_count: 0,
        cells_u: 0,
        cells_v: 0,
        construction_mode: "none",
        target_piece,
        emitted_patches: Vec::new(),
        emitted_provenance: Vec::new(),
        decline_reason: "support_cell_primitive_no_candidate",
    };
    if targets.is_empty() {
        report.decline_reason = "support_cell_primitive_no_target_site";
        return report;
    }

    let mut best: Option<Kp20SupportCellPrimitiveTrial> = None;
    for cells_u in [2, 3] {
        for cells_v in [2, 3] {
            report.evaluated += 1;
            let trial = measure_cell_grid(macro_patches, &targets, target_piece, cells_u, cells_v);
            if best
                .as_ref()
                .is_none_or(|best| trial.folded_patch_count < best.folded_patch_count)
            {
                best = Some(trial);
            }
        }
    }
    let Some(best) = best else {
        return report;
    };
    report.patch_count = best.patch_count;
    report.folded_patch_count = best.folded_patch_count;
    report.cells_u = best.cells_u;
    report.cells_v = best.cells_v;
    report.construction_mode = best.construction_mode.label();
    report.emitted_patches = best.emitted_patches;
    report.emitted_provenance = best.emitted_provenance;
    report.decline_reason = if best.folded_patch_count == 0 {
        "support_cell_primitive_fold_free_not_emitted"
    } else {
        "support_cell_primitive_retains_folds"
    };
    report
}

#[derive(Debug, Clone)]
struct Kp20SupportCellPrimitiveTrial {
    patch_count: usize,
    folded_patch_count: usize,
    cells_u: usize,
    cells_v: usize,
    construction_mode: Kp20SupportCellConstructionMode,
    emitted_patches: Vec<GregoryPatch>,
    emitted_provenance: Vec<CornerSubPatchProvenance>,
}

fn measure_cell_grid(
    macro_patches: &[Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
    cells_u: usize,
    cells_v: usize,
) -> Kp20SupportCellPrimitiveTrial {
    let mut best = Kp20SupportCellPrimitiveTrial {
        patch_count: 0,
        folded_patch_count: usize::MAX,
        cells_u,
        cells_v,
        construction_mode: Kp20SupportCellConstructionMode::InheritedCellBoundary,
        emitted_patches: Vec::new(),
        emitted_provenance: Vec::new(),
    };
    for construction_mode in Kp20SupportCellConstructionMode::all() {
        let trial = measure_cell_grid_mode(
            macro_patches,
            targets,
            target_piece,
            cells_u,
            cells_v,
            construction_mode,
        );
        if trial.folded_patch_count < best.folded_patch_count {
            best = trial;
        }
    }
    best
}

#[derive(Debug, Clone, Copy)]
enum Kp20SupportCellConstructionMode {
    InheritedCellBoundary,
    FullSupportBoundary,
    SupportBoundaryRibbon,
    ParallelSupportBoundaryRibbon,
}

impl Kp20SupportCellConstructionMode {
    fn all() -> [Self; 4] {
        [
            Self::InheritedCellBoundary,
            Self::FullSupportBoundary,
            Self::SupportBoundaryRibbon,
            Self::ParallelSupportBoundaryRibbon,
        ]
    }

    fn label(self) -> &'static str {
        match self {
            Self::InheritedCellBoundary => "inherited_cell_boundary",
            Self::FullSupportBoundary => "full_support_boundary",
            Self::SupportBoundaryRibbon => "support_boundary_ribbon",
            Self::ParallelSupportBoundaryRibbon => "parallel_support_boundary_ribbon",
        }
    }
}

fn measure_cell_grid_mode(
    macro_patches: &[Kp20MacroPatch],
    targets: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
    cells_u: usize,
    cells_v: usize,
    construction_mode: Kp20SupportCellConstructionMode,
) -> Kp20SupportCellPrimitiveTrial {
    let mut patch_count = 0;
    let mut folded_patch_count = 0;
    let total_patch_count = targets.len() * cells_u * cells_v;
    let sector_count = macro_patches.len();
    let mut emitted_patches = Vec::with_capacity(total_patch_count);
    let mut emitted_provenance = Vec::with_capacity(total_patch_count);
    for target in targets {
        if target.sector >= macro_patches.len() {
            continue;
        }
        let base = *macro_patches[target.sector].piece(target_piece);
        let authored = match construction_mode {
            Kp20SupportCellConstructionMode::InheritedCellBoundary => base,
            Kp20SupportCellConstructionMode::FullSupportBoundary => {
                exact_boundary_coons_cell(&base)
            }
            Kp20SupportCellConstructionMode::SupportBoundaryRibbon => {
                support_boundary_ribbon(&base)
            }
            Kp20SupportCellConstructionMode::ParallelSupportBoundaryRibbon => {
                parallel_support_boundary_ribbon(&base)
            }
        };
        for cell in support_cells(&authored, cells_u, cells_v) {
            let patch_index = emitted_patches.len();
            emitted_patches.push(gregory_from_bi6(cell));
            emitted_provenance.push(CornerSubPatchProvenance::kp20_support_cell(
                target.sector,
                sector_count,
                patch_index,
                total_patch_count,
            ));
            patch_count += 1;
            if emitted_patches[patch_index]
                .self_intersection_report(SELF_INT_TOL)
                .intersects()
            {
                folded_patch_count += 1;
            }
        }
    }
    Kp20SupportCellPrimitiveTrial {
        patch_count,
        folded_patch_count,
        cells_u,
        cells_v,
        construction_mode,
        emitted_patches,
        emitted_provenance,
    }
}

fn support_cells(base: &Kp20Bi6Patch, cells_u: usize, cells_v: usize) -> Vec<Kp20Bi6Patch> {
    let mut cells = Vec::with_capacity(cells_u * cells_v);
    for u_band in split_uniform_u(base, cells_u) {
        cells.extend(split_uniform_v(&u_band, cells_v));
    }
    cells
}

fn split_uniform_u(base: &Kp20Bi6Patch, count: usize) -> Vec<Kp20Bi6Patch> {
    split_uniform(base, count, split_bi6_patch_u)
}

fn split_uniform_v(base: &Kp20Bi6Patch, count: usize) -> Vec<Kp20Bi6Patch> {
    split_uniform(base, count, split_bi6_patch_v)
}

fn split_uniform(
    base: &Kp20Bi6Patch,
    count: usize,
    split: fn(&Kp20Bi6Patch, f64) -> (Kp20Bi6Patch, Kp20Bi6Patch),
) -> Vec<Kp20Bi6Patch> {
    let mut patches = Vec::with_capacity(count);
    let mut remainder = *base;
    for index in 0..count {
        if index + 1 == count {
            patches.push(remainder);
        } else {
            let local_t = 1.0 / ((count - index) as f64);
            let (head, tail) = split(&remainder, local_t);
            patches.push(head);
            remainder = tail;
        }
    }
    patches
}

fn exact_boundary_coons_cell(cell: &Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut rebuilt = Kp20Bi6Patch::zero();
    for row in 0..7 {
        let u = row as f64 / 6.0;
        for column in 0..7 {
            let v = column as f64 / 6.0;
            let left_to_right = lerp(cell.cps[row][0], cell.cps[row][6], v);
            let bottom_to_top = lerp(cell.cps[0][column], cell.cps[6][column], u);
            let corner_blend = bilinear_corners(cell, u, v);
            rebuilt.cps[row][column] = left_to_right + (bottom_to_top - corner_blend);
        }
    }
    rebuilt
}

fn support_boundary_ribbon(base: &Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut ribbon = Kp20Bi6Patch::zero();
    let fallback_direction = average_interior_direction(base);
    for column in 0..7 {
        ribbon.cps[0][column] = base.cps[0][column];
        ribbon.cps[1][column] = base.cps[1][column];
        let support = base.cps[1][column];
        let depth = support
            .distance_to(&base.cps[6][column])
            .max(support_span(base) * 0.35);
        let direction = (base.cps[6][column] - support)
            .normalized()
            .unwrap_or(fallback_direction);
        for row in 2..7 {
            let t = (row - 1) as f64 / 5.0;
            ribbon.cps[row][column] = support + direction * (depth * smooth_step(t));
        }
    }
    ribbon
}

fn parallel_support_boundary_ribbon(base: &Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut ribbon = Kp20Bi6Patch::zero();
    let direction = average_interior_direction(base);
    let depth = average_interior_depth(base).max(support_span(base) * 0.35);
    for column in 0..7 {
        ribbon.cps[0][column] = base.cps[0][column];
        ribbon.cps[1][column] = base.cps[1][column];
        let support = base.cps[1][column];
        for row in 2..7 {
            let t = (row - 1) as f64 / 5.0;
            ribbon.cps[row][column] = support + direction * (depth * smooth_step(t));
        }
    }
    ribbon
}

fn average_interior_direction(base: &Kp20Bi6Patch) -> Vec3 {
    let mut direction = Vec3::new(0.0, 0.0, 0.0);
    for column in 0..7 {
        direction = direction + (base.cps[6][column] - base.cps[1][column]);
    }
    direction.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

fn average_interior_depth(base: &Kp20Bi6Patch) -> f64 {
    let mut depth = 0.0;
    for column in 0..7 {
        depth += base.cps[1][column].distance_to(&base.cps[6][column]);
    }
    depth / 7.0
}

fn support_span(base: &Kp20Bi6Patch) -> f64 {
    let support_width = base.cps[1][0].distance_to(&base.cps[1][6]);
    let interior_width = base.cps[6][0].distance_to(&base.cps[6][6]);
    let depth = base.cps[1][3].distance_to(&base.cps[6][3]);
    ((support_width + interior_width + depth) / 3.0).max(1.0e-9)
}

fn smooth_step(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

fn bilinear_corners(cell: &Kp20Bi6Patch, u: f64, v: f64) -> Point3 {
    let bottom = lerp(cell.cps[0][0], cell.cps[0][6], v);
    let top = lerp(cell.cps[6][0], cell.cps[6][6], v);
    lerp(bottom, top, u)
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}

fn gregory_from_bi6(patch: Kp20Bi6Patch) -> GregoryPatch {
    GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
        grid: elevate_bi6_grid_to_bi9(patch.cps),
    }))
}

fn support_targets_by_sector(
    fold_sites: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> Vec<Kp20FoldSite> {
    let mut targets = Vec::new();
    for site in fold_sites
        .iter()
        .copied()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
    {
        if let Some(existing) = targets
            .iter_mut()
            .find(|existing: &&mut Kp20FoldSite| existing.sector == site.sector)
        {
            if site.min_jacobian_normal_cos < existing.min_jacobian_normal_cos {
                *existing = site;
            }
        } else {
            targets.push(site);
        }
    }
    targets
}
