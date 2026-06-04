//! KP20 support-strip row placement and tangent-row reauthoring.

use popkern_core::geometry::gregory::{Kp20Bi6Patch, Kp20PieceId};
use popkern_core::geometry::{Point3, Vec3};

use super::super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20StripReauthorReport {
    pub(super) max_interior_displacement: f64,
    pub(super) min_parameter_step: f64,
    pub(super) max_parameter_skew: f64,
    pub(super) max_tangent_step_imbalance: f64,
    pub(super) min_tangent_alignment: f64,
    pub(super) row_strategy_counts: Kp20StripRowStrategyCounts,
    pub(super) tangent_row_reauthored_columns: usize,
    pub(super) max_tangent_row_displacement: f64,
}

pub(super) fn reauthor_strip_rows(
    surface: &mut Kp20Bi6Patch,
    fold_summary: &Kp20SamplerFoldSummary,
    sector_index: usize,
    strip_piece: Kp20PieceId,
) -> Kp20StripReauthorReport {
    let original = *surface;
    let tangent_report = reauthor_interface_tangent_row(surface);
    let tangent_adjusted = *surface;
    let mut min_parameter_step = f64::INFINITY;
    let mut max_parameter_skew = 0.0_f64;
    let mut max_tangent_step_imbalance = 0.0_f64;
    let mut min_tangent_alignment = f64::INFINITY;
    let mut row_strategy_counts = Kp20StripRowStrategyCounts::default();
    for column in 0..7 {
        let geometry = Kp20StripColumnGeometry::from_surface(&tangent_adjusted, column);
        let map =
            strip_column_parameter_map(fold_summary, sector_index, strip_piece, column, geometry);
        row_strategy_counts.record(map.strategy);
        min_parameter_step = min_parameter_step.min(map.min_step());
        max_parameter_skew = max_parameter_skew.max(map.max_skew());
        max_tangent_step_imbalance =
            max_tangent_step_imbalance.max(geometry.tangent_step_imbalance());
        min_tangent_alignment = min_tangent_alignment.min(geometry.tangent_alignment());
        for (idx, row) in (2..=4).enumerate() {
            let t = map.row_parameters[idx];
            surface.cps[row][column] = geometry.support_adjacent
                + (geometry.interface_adjacent - geometry.support_adjacent) * t;
        }
    }
    Kp20StripReauthorReport {
        max_interior_displacement: max_interior_displacement(&original, surface),
        min_parameter_step,
        max_parameter_skew,
        max_tangent_step_imbalance,
        min_tangent_alignment,
        row_strategy_counts,
        tangent_row_reauthored_columns: tangent_report.reauthored_columns,
        max_tangent_row_displacement: tangent_report.max_displacement,
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20StripTangentRowReport {
    reauthored_columns: usize,
    max_displacement: f64,
}

fn reauthor_interface_tangent_row(surface: &mut Kp20Bi6Patch) -> Kp20StripTangentRowReport {
    let mut reauthored_columns = 0_usize;
    let mut max_displacement = 0.0_f64;
    for column in 0..7 {
        let geometry = Kp20StripColumnGeometry::from_surface(surface, column);
        if !geometry.needs_tangent_row_reauthor() {
            continue;
        }
        let Some(target) = geometry.balanced_interface_adjacent() else {
            continue;
        };
        let current = surface.cps[5][column];
        let displacement = capped_vector((target - current) * TANGENT_ROW_DAMPING, geometry.cap());
        if displacement.length() <= 1.0e-12 {
            continue;
        }
        surface.cps[5][column] = current + displacement;
        reauthored_columns += 1;
        max_displacement = max_displacement.max(displacement.length());
    }
    Kp20StripTangentRowReport {
        reauthored_columns,
        max_displacement,
    }
}

fn capped_vector(displacement: Vec3, max_length: f64) -> Vec3 {
    let length = displacement.length();
    if length > max_length {
        if let Some(direction) = displacement.normalized() {
            return direction * max_length;
        }
    }
    displacement
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20StripParameterMap {
    pub(super) row_parameters: [f64; 3],
    pub(super) strategy: Kp20StripRowPlacementStrategy,
}

#[derive(Debug, Default, Clone, Copy)]
pub(super) struct Kp20StripRowStrategyCounts {
    pub(super) uniform: usize,
    pub(super) fold_site: usize,
    pub(super) tangent_balance: usize,
    pub(super) combined: usize,
}

impl Kp20StripRowStrategyCounts {
    pub(super) fn record(&mut self, strategy: Kp20StripRowPlacementStrategy) {
        match strategy {
            Kp20StripRowPlacementStrategy::UniformBlend => self.uniform += 1,
            Kp20StripRowPlacementStrategy::FoldSiteGuided => self.fold_site += 1,
            Kp20StripRowPlacementStrategy::TangentBalanced => self.tangent_balance += 1,
            Kp20StripRowPlacementStrategy::FoldSiteTangentBalanced => self.combined += 1,
        }
    }

    pub(super) fn add(&mut self, other: Self) {
        self.uniform += other.uniform;
        self.fold_site += other.fold_site;
        self.tangent_balance += other.tangent_balance;
        self.combined += other.combined;
    }

    fn total(self) -> usize {
        self.uniform + self.fold_site + self.tangent_balance + self.combined
    }

    fn tangent_aware(self) -> usize {
        self.tangent_balance + self.combined
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kp20StripRowPlacementStrategy {
    UniformBlend,
    FoldSiteGuided,
    TangentBalanced,
    FoldSiteTangentBalanced,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Kp20StripColumnGeometry {
    pub(super) support_boundary: Point3,
    pub(super) support_adjacent: Point3,
    pub(super) interface_adjacent: Point3,
    pub(super) interface_boundary: Point3,
}

impl Kp20StripColumnGeometry {
    pub(super) fn from_surface(surface: &Kp20Bi6Patch, column: usize) -> Self {
        Self {
            support_boundary: surface.cps[0][column],
            support_adjacent: surface.cps[1][column],
            interface_adjacent: surface.cps[5][column],
            interface_boundary: surface.cps[6][column],
        }
    }

    pub(super) fn tangent_step_imbalance(self) -> f64 {
        let support_step = self.support_step();
        let interface_step = self.interface_step();
        (support_step - interface_step).abs() / support_step.max(interface_step)
    }

    pub(super) fn tangent_alignment(self) -> f64 {
        let support = self.support_adjacent - self.support_boundary;
        let interface = self.interface_adjacent - self.interface_boundary;
        let Some(support) = support.normalized() else {
            return 1.0;
        };
        let Some(interface) = interface.normalized() else {
            return 1.0;
        };
        support.dot(&interface).abs().clamp(0.0, 1.0)
    }

    fn tangent_parameter_bias(self) -> f64 {
        let support_step = self.support_step();
        let interface_step = self.interface_step();
        let balance = (support_step - interface_step) / (support_step + interface_step);
        (balance * self.tangent_alignment() * 0.18).clamp(-0.12, 0.12)
    }

    fn needs_tangent_balanced_strategy(self) -> bool {
        self.tangent_step_imbalance() > 0.15 && self.tangent_alignment() > 0.35
    }

    fn needs_tangent_row_reauthor(self) -> bool {
        self.tangent_step_imbalance() > 0.25 || self.tangent_alignment() < 0.35
    }

    fn balanced_interface_adjacent(self) -> Option<Point3> {
        let support = self.support_adjacent - self.support_boundary;
        let interface = self.interface_adjacent - self.interface_boundary;
        let support_length = support.length();
        if support_length <= 1.0e-12 {
            return None;
        }
        let support_direction = support.normalized()?;
        let direction = interface
            .normalized()
            .map_or(support_direction, |interface_direction| {
                if support_direction.dot(&interface_direction) < 0.0 {
                    -support_direction
                } else {
                    support_direction
                }
            });
        Some(self.interface_boundary + direction * support_length)
    }

    fn support_step(self) -> f64 {
        self.support_adjacent
            .distance_to(&self.support_boundary)
            .max(1.0e-12)
    }

    fn interface_step(self) -> f64 {
        self.interface_adjacent
            .distance_to(&self.interface_boundary)
            .max(1.0e-12)
    }

    fn cap(self) -> f64 {
        self.support_step().max(self.interface_step()) * TANGENT_ROW_DISPLACEMENT_CAP
    }
}

pub(super) fn strip_column_parameter_map(
    fold_summary: &Kp20SamplerFoldSummary,
    sector_index: usize,
    strip_piece: Kp20PieceId,
    column: usize,
    geometry: Kp20StripColumnGeometry,
) -> Kp20StripParameterMap {
    let control_v = column as f64 / 6.0;
    let mut weighted_u = 0.0_f64;
    let mut total_weight = 0.0_f64;
    for site in fold_summary.fold_sites.iter().filter(|site| {
        site.sector == sector_index
            && super::super::kp20_polar_support_surface_split::support_owner_piece(site.piece)
                == strip_piece
    }) {
        let distance = (control_v - site.v).abs();
        let weight = (1.0 - distance / 0.35).max(0.0).powi(2);
        let severity = (1.0_f64 - site.min_jacobian_normal_cos).clamp(0.0, 1.0);
        weighted_u += site.u.clamp(0.0, 1.0) * weight * severity;
        total_weight += weight * severity;
    }
    let has_fold_site_target = total_weight > 1.0e-12;
    let needs_tangent_balance = geometry.needs_tangent_balanced_strategy();
    let strategy =
        Kp20StripRowPlacementStrategy::select(has_fold_site_target, needs_tangent_balance);
    let map = if has_fold_site_target {
        let target_u = weighted_u / total_weight;
        let strength = (total_weight / (total_weight + 1.0)).clamp(0.0, 0.55);
        Kp20StripParameterMap::from_target(target_u, strength, strategy)
    } else {
        Kp20StripParameterMap::uniform(strategy)
    };
    if needs_tangent_balance {
        map.with_tangent_bias(geometry.tangent_parameter_bias())
    } else {
        map
    }
}

impl Kp20StripParameterMap {
    fn uniform(strategy: Kp20StripRowPlacementStrategy) -> Self {
        Self {
            row_parameters: UNIFORM_ROW_PARAMETERS,
            strategy,
        }
    }

    fn from_target(target_u: f64, strength: f64, strategy: Kp20StripRowPlacementStrategy) -> Self {
        let mut rows = UNIFORM_ROW_PARAMETERS;
        for (row, base) in rows.iter_mut().zip(UNIFORM_ROW_PARAMETERS) {
            *row = (base + (target_u - base) * strength).clamp(0.12, 0.88);
        }
        enforce_minimum_parameter_steps(&mut rows);
        Self {
            row_parameters: rows,
            strategy,
        }
    }

    fn with_tangent_bias(mut self, bias: f64) -> Self {
        for row in &mut self.row_parameters {
            *row = (*row + bias * *row * (1.0 - *row)).clamp(0.12, 0.88);
        }
        enforce_minimum_parameter_steps(&mut self.row_parameters);
        self
    }

    pub(super) fn min_step(self) -> f64 {
        let [a, b, c] = self.row_parameters;
        a.min(b - a).min(c - b).min(1.0 - c)
    }

    pub(super) fn max_skew(self) -> f64 {
        self.row_parameters
            .into_iter()
            .zip(UNIFORM_ROW_PARAMETERS)
            .map(|(actual, uniform)| (actual - uniform).abs())
            .fold(0.0_f64, f64::max)
    }
}

impl Kp20StripRowPlacementStrategy {
    fn select(has_fold_site_target: bool, needs_tangent_balance: bool) -> Self {
        match (has_fold_site_target, needs_tangent_balance) {
            (false, false) => Self::UniformBlend,
            (true, false) => Self::FoldSiteGuided,
            (false, true) => Self::TangentBalanced,
            (true, true) => Self::FoldSiteTangentBalanced,
        }
    }
}

pub(super) fn row_strategy_assessment(
    counts: Kp20StripRowStrategyCounts,
    max_tangent_step_imbalance: f64,
    min_tangent_alignment: f64,
) -> &'static str {
    let total = counts.total();
    if total == 0 {
        return "no_strip_columns";
    }
    if max_tangent_step_imbalance > 0.5
        && min_tangent_alignment < 0.2
        && counts.tangent_aware() * 2 >= total
    {
        "support_tangent_alignment_limited"
    } else if max_tangent_step_imbalance > 0.5 && counts.tangent_aware() * 2 < total {
        "tangent_strategy_underactive"
    } else if counts.combined > 0 {
        "combined_strategy_active"
    } else if counts.fold_site > 0 {
        "fold_site_strategy_active"
    } else {
        "uniform_strategy_active"
    }
}

fn enforce_minimum_parameter_steps(rows: &mut [f64; 3]) {
    rows[0] = rows[0].clamp(MIN_PARAMETER_STEP, 1.0 - 3.0 * MIN_PARAMETER_STEP);
    rows[1] = rows[1].clamp(rows[0] + MIN_PARAMETER_STEP, 1.0 - 2.0 * MIN_PARAMETER_STEP);
    rows[2] = rows[2].clamp(rows[1] + MIN_PARAMETER_STEP, 1.0 - MIN_PARAMETER_STEP);
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

pub(super) const MIN_PARAMETER_STEP: f64 = 0.12;
pub(super) const UNIFORM_ROW_PARAMETERS: [f64; 3] = [0.25, 0.5, 0.75];
const TANGENT_ROW_DAMPING: f64 = 0.65;
const TANGENT_ROW_DISPLACEMENT_CAP: f64 = 0.5;
