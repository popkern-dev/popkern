//! Curvature-weighted normal solve for constructive row bands.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::DerivedRowBand;
const WEIGHTED_NORMAL_SAMPLE_LIMIT: usize = 32;
const INTERIOR_CONTROL_COUNT: usize = 6;
const WEIGHTED_NORMAL_REGULARIZATION: f64 = 1.0e-8;
const MIN_WEIGHTED_NORMAL_IMPACT: f64 = 5.0e-2;
const COMMERCIAL_G1_TARGET: f64 = 1.25e-3;
const MAX_WEIGHTED_NORMAL_CONTROL_FRACTION: f64 = 0.03;
const WEIGHTED_NORMAL_MAX_PASSES: usize = 4;
const WEIGHTED_NORMAL_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];
const WEIGHTED_NORMAL_IMPROVEMENT_EPS: f64 = 1.0e-12;
pub(super) fn repair_curvature_weighted_normal_solve(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> String {
    if bands.len() < 2 {
        return "curvature_weighted_normal_solve_status=insufficient_bands;curvature_weighted_normal_sample_count=0;curvature_weighted_normal_pass_count=0;curvature_weighted_normal_constraint_count=0;curvature_weighted_normal_repaired_row_count=0;curvature_weighted_normal_limited_row_count=0;max_curvature_weighted_normal_control_delta=0.000000e0;max_curvature_weighted_normal_impact_before=0.000000000000000e0;max_curvature_weighted_normal_impact_after=0.000000000000000e0;max_curvature_weighted_normal_g1_before=0.000000000000000e0;max_curvature_weighted_normal_g1_after=0.000000000000000e0;worst_curvature_weighted_normal_detail=[none]".to_string();
    }
    let mut before_impact = max_weighted_normal_impact(bands, chain_trims);
    let mut before_g1 = max_trim_disagreement(bands, chain_trims);
    let initial_impact = before_impact;
    let initial_g1 = before_g1;
    let mut pass_count = 0usize;
    let mut total_samples = 0usize;
    let mut total_constraints = 0usize;
    let mut total_repaired_rows = 0usize;
    let mut total_limited_rows = 0usize;
    let mut max_control_delta = 0.0_f64;
    let mut best_detail = top_weighted_normal_samples(bands, chain_trims)
        .first()
        .map(|sample| sample.detail.clone())
        .unwrap_or_else(|| "none".to_string());
    let mut decline_reason = "no_active_samples";

    for _ in 0..WEIGHTED_NORMAL_MAX_PASSES {
        let active_samples = top_weighted_normal_samples(bands, chain_trims);
        total_samples += active_samples.len();
        total_constraints += active_samples.len() * 3;
        if active_samples.is_empty() {
            if pass_count == 0 {
                decline_reason = "no_active_samples";
            }
            break;
        }
        if let Some(sample) = active_samples.first() {
            best_detail = sample.detail.clone();
        }
        let Some(solution) = solve_weighted_normal_system(bands, &active_samples) else {
            if pass_count == 0 {
                decline_reason = "weighted_normal_system_singular";
            }
            break;
        };
        let (deltas, limited_rows, max_requested_delta) = limited_row_deltas(bands, &solution);
        let mut best: Option<(f64, f64, f64, Vec<DerivedRowBand>)> = None;
        for damping in WEIGHTED_NORMAL_DAMPING_FACTORS {
            let mut candidate = bands.to_vec();
            apply_deltas(&mut candidate, &deltas, damping);
            let after_impact = max_weighted_normal_impact(&candidate, chain_trims);
            if after_impact + WEIGHTED_NORMAL_IMPROVEMENT_EPS >= before_impact {
                continue;
            }
            let after_g1 = max_trim_disagreement(&candidate, chain_trims);
            if after_g1 > before_g1.max(COMMERCIAL_G1_TARGET) + 1.0e-10 {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_impact, _, _)| after_impact < *best_impact)
            {
                best = Some((damping, after_impact, after_g1, candidate));
            }
        }
        let Some((damping, after_impact, after_g1, candidate)) = best else {
            if pass_count == 0 {
                decline_reason = "no_quality_safe_improvement";
            }
            total_limited_rows += limited_rows;
            max_control_delta = max_control_delta.max(max_requested_delta);
            break;
        };
        total_repaired_rows += write_repaired_rows(bands, candidate);
        total_limited_rows += limited_rows;
        max_control_delta = max_control_delta.max(max_requested_delta * damping);
        before_impact = after_impact;
        before_g1 = after_g1;
        pass_count += 1;
        decline_reason = "kept";
    }

    format!(
        "curvature_weighted_normal_solve_status={};curvature_weighted_normal_sample_count={};curvature_weighted_normal_pass_count={};curvature_weighted_normal_constraint_count={};curvature_weighted_normal_repaired_row_count={};curvature_weighted_normal_limited_row_count={};max_curvature_weighted_normal_control_delta={:.6e};max_curvature_weighted_normal_impact_before={:.15e};max_curvature_weighted_normal_impact_after={:.15e};max_curvature_weighted_normal_g1_before={:.15e};max_curvature_weighted_normal_g1_after={:.15e};worst_curvature_weighted_normal_detail=[{}]",
        decline_reason,
        total_samples,
        pass_count,
        total_constraints,
        total_repaired_rows,
        total_limited_rows,
        max_control_delta,
        initial_impact,
        before_impact,
        initial_g1,
        before_g1,
        best_detail
    )
}
fn top_weighted_normal_samples(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> Vec<WeightedNormalSample> {
    let band_indices = band_indices(bands);
    let mut samples = Vec::new();
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 1..8 {
            let fraction = step as f64 / 8.0;
            let chain_param = paired_chain_param(&bands[first], &bands[second], trim, fraction);
            let Some(sample_a) = row_weighted_normal_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_weighted_normal_sample(&bands[second], chain_param) else {
                continue;
            };
            let impact = weighted_normal_impact(&sample_a, &sample_b);
            if impact <= MIN_WEIGHTED_NORMAL_IMPACT {
                continue;
            }
            samples.push(WeightedNormalSample {
                first,
                second,
                fraction,
                chain_param,
                impact,
                detail: weighted_normal_detail(
                    &bands[first],
                    &bands[second],
                    fraction,
                    chain_param,
                    &sample_a,
                    &sample_b,
                    impact,
                ),
            });
        }
    }
    samples.sort_by(|a, b| {
        b.impact
            .total_cmp(&a.impact)
            .then_with(|| a.first.cmp(&b.first))
            .then_with(|| a.second.cmp(&b.second))
            .then_with(|| a.fraction.total_cmp(&b.fraction))
    });
    samples.truncate(WEIGHTED_NORMAL_SAMPLE_LIMIT);
    samples
}
fn solve_weighted_normal_system(
    bands: &[DerivedRowBand],
    samples: &[WeightedNormalSample],
) -> Option<Vec<f64>> {
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    let mut matrix = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut rhs = vec![0.0_f64; variable_count];
    for sample in samples {
        let sample_a = row_weighted_normal_sample(&bands[sample.first], sample.chain_param)?;
        let sample_b = row_weighted_normal_sample(&bands[sample.second], sample.chain_param)?;
        let orientation = normal_alignment(sample_a.normal, sample_b.normal);
        let diff = sample_a.normal - sample_b.normal * orientation;
        let weight_sqrt = sample.impact.sqrt();
        for axis in 0..3 {
            let mut row = vec![0.0_f64; variable_count];
            fill_coefficients(&mut row, sample.first, weight_sqrt, &sample_a, axis);
            fill_coefficients(
                &mut row,
                sample.second,
                -orientation * weight_sqrt,
                &sample_b,
                axis,
            );
            add_normal_equation(
                &mut matrix,
                &mut rhs,
                &row,
                -component(diff, axis) * weight_sqrt,
            );
        }
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += WEIGHTED_NORMAL_REGULARIZATION;
    }
    solve_dense(matrix, rhs)
}
fn max_weighted_normal_impact(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> f64 {
    top_weighted_normal_samples(bands, chain_trims)
        .first()
        .map(|sample| sample.impact)
        .unwrap_or(0.0)
}

fn max_trim_disagreement(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> f64 {
    let band_indices = band_indices(bands);
    let mut worst = 0.0_f64;
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 0..9 {
            let fraction = step as f64 / 8.0;
            let chain_param = paired_chain_param(&bands[first], &bands[second], trim, fraction);
            let Some(sample_a) = row_weighted_normal_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_weighted_normal_sample(&bands[second], chain_param) else {
                continue;
            };
            worst = worst.max(1.0 - sample_a.normal.dot(&sample_b.normal).abs());
        }
    }
    worst
}

fn row_weighted_normal_sample(
    band: &DerivedRowBand,
    chain_param: f64,
) -> Option<RowWeightedNormalSample> {
    let t = ((chain_param - band.interval.0) / (band.interval.1 - band.interval.0)).clamp(0.0, 1.0);
    let basis = bernstein9(t);
    let derivative = row_derivative_coefficients(&band.row0);
    let tangent = eval_ribbon_field(&derivative, t);
    let speed_sq = tangent.dot(&tangent);
    if speed_sq <= 1.0e-300 {
        return None;
    }
    let mut cross = Vec3::new(0.0, 0.0, 0.0);
    for (k, weight) in basis.iter().enumerate() {
        cross = cross + (band.row1[k] - band.row0[k]) * (9.0 * *weight);
    }
    let raw_normal = tangent.cross(&cross);
    let raw_normal_magnitude = raw_normal.dot(&raw_normal).sqrt();
    let normal = raw_normal.normalized()?;
    let second = eval_row_second_derivative(&band.row0, t);
    let directional_second = second * (1.0 / speed_sq);
    let curvature = normal.dot(&directional_second);
    let mut normal_sensitivities = [Vec3::new(0.0, 0.0, 0.0); INTERIOR_CONTROL_COUNT];
    for (column, sensitivity) in normal_sensitivities.iter_mut().enumerate() {
        let k = column + 2;
        let Some(direction) = control_normal_direction(band, k) else {
            continue;
        };
        let cross_delta = direction * (9.0 * basis[k]);
        let normal_delta = tangent.cross(&cross_delta);
        *sensitivity =
            (normal_delta - normal * normal_delta.dot(&normal)) * (1.0 / raw_normal_magnitude);
    }
    Some(RowWeightedNormalSample {
        normal,
        curvature,
        q_magnitude: directional_second.dot(&directional_second).sqrt(),
        normal_sensitivities,
    })
}

fn weighted_normal_impact(
    first: &RowWeightedNormalSample,
    second: &RowWeightedNormalSample,
) -> f64 {
    let alignment = normal_alignment(first.normal, second.normal);
    let disagreement = 1.0 - first.normal.dot(&second.normal).abs();
    let angle_estimate = (2.0 * disagreement.max(0.0)).sqrt();
    let aligned_second = second.curvature * alignment;
    let scale = 1.0 + first.curvature.abs().max(aligned_second.abs());
    first.q_magnitude.max(second.q_magnitude) * angle_estimate / scale
}

fn weighted_normal_detail(
    first: &DerivedRowBand,
    second: &DerivedRowBand,
    fraction: f64,
    chain_param: f64,
    sample_a: &RowWeightedNormalSample,
    sample_b: &RowWeightedNormalSample,
    impact: f64,
) -> String {
    let g1 = 1.0 - sample_a.normal.dot(&sample_b.normal).abs();
    format!(
        "{:?}:{:?}<->{:?}:{:?}@f={:.4},s={:.4};impact={:.6e};q={:.6e};g1={:.6e}",
        first.key.face,
        first.key.side,
        second.key.face,
        second.key.side,
        fraction,
        chain_param,
        impact,
        sample_a.q_magnitude.max(sample_b.q_magnitude),
        g1
    )
}

fn eval_row_second_derivative(row: &[popkern_core::geometry::Point3; 10], t: f64) -> Vec3 {
    let mut controls = [Vec3::new(0.0, 0.0, 0.0); 8];
    for index in 0..8 {
        let a = row[index];
        let b = row[index + 1];
        let c = row[index + 2];
        controls[index] = Vec3::new(
            c.x - 2.0 * b.x + a.x,
            c.y - 2.0 * b.y + a.y,
            c.z - 2.0 * b.z + a.z,
        ) * 72.0;
    }
    let mut working = controls;
    for level in 1..8 {
        for index in 0..(8 - level) {
            working[index] = working[index] + (working[index + 1] - working[index]) * t;
        }
    }
    working[0]
}

fn band_indices(bands: &[DerivedRowBand]) -> HashMap<ChainRowKey, usize> {
    bands
        .iter()
        .enumerate()
        .map(|(index, band)| (band.key, index))
        .collect()
}

fn paired_chain_param(
    first: &DerivedRowBand,
    second: &DerivedRowBand,
    trim: &ChainTrim<ChainRowKey, VertexHandle>,
    fraction: f64,
) -> f64 {
    let first_param = trim_chain_param(first, trim.interval_a, fraction);
    let second_param = trim_chain_param(second, trim.interval_b, 1.0 - fraction);
    0.5 * (first_param + second_param)
}

fn trim_chain_param(band: &DerivedRowBand, trim_interval: (f64, f64), fraction: f64) -> f64 {
    let row_param = trim_interval.0 + (trim_interval.1 - trim_interval.0) * fraction;
    band.interval.0 + (band.interval.1 - band.interval.0) * row_param
}

fn control_normal_direction(band: &DerivedRowBand, k: usize) -> Option<Vec3> {
    let t = k as f64 / 9.0;
    let derivative = row_derivative_coefficients(&band.row0);
    let tangent = eval_ribbon_field(&derivative, t);
    let field = eval_ribbon_field(&band.field_segment, t);
    tangent.cross(&field).normalized()
}

fn normal_alignment(first: Vec3, second: Vec3) -> f64 {
    if first.dot(&second) < 0.0 { -1.0 } else { 1.0 }
}

fn limited_row_deltas(bands: &[DerivedRowBand], solution: &[f64]) -> (Vec<[f64; 6]>, usize, f64) {
    let mut deltas = vec![[0.0_f64; 6]; bands.len()];
    let mut limited_rows = 0usize;
    let mut max_delta = 0.0_f64;
    for (row_index, band) in bands.iter().enumerate() {
        let limit = row_control_limit(band);
        let offset = row_index * INTERIOR_CONTROL_COUNT;
        let max_requested = solution[offset..offset + INTERIOR_CONTROL_COUNT]
            .iter()
            .map(|value| value.abs())
            .fold(0.0_f64, f64::max);
        let scale = if max_requested > limit {
            limited_rows += 1;
            limit / max_requested
        } else {
            1.0
        };
        for column in 0..INTERIOR_CONTROL_COUNT {
            deltas[row_index][column] = solution[offset + column] * scale;
            max_delta = max_delta.max(deltas[row_index][column].abs());
        }
    }
    (deltas, limited_rows, max_delta)
}

fn apply_deltas(bands: &mut [DerivedRowBand], deltas: &[[f64; 6]], damping: f64) {
    for (band, row_deltas) in bands.iter_mut().zip(deltas) {
        for (column, value) in row_deltas.iter().enumerate() {
            let k = column + 2;
            let Some(direction) = control_normal_direction(band, k) else {
                continue;
            };
            band.row1[k] = band.row1[k] + direction * (*value * damping);
        }
    }
}

fn write_repaired_rows(bands: &mut [DerivedRowBand], candidate: Vec<DerivedRowBand>) -> usize {
    let mut repaired_rows = 0usize;
    for (target, source) in bands.iter_mut().zip(candidate) {
        let moved = target
            .row1
            .iter()
            .zip(source.row1.iter())
            .any(|(a, b)| a.distance_to(b) > 1.0e-14);
        repaired_rows += usize::from(moved);
        target.row1 = source.row1;
    }
    repaired_rows
}

fn row_control_limit(band: &DerivedRowBand) -> f64 {
    let boundary_length: f64 = (1..10)
        .map(|k| band.row0[k].distance_to(&band.row0[k - 1]))
        .sum();
    let band_depth = band
        .row1
        .iter()
        .zip(band.row0.iter())
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    MAX_WEIGHTED_NORMAL_CONTROL_FRACTION * boundary_length.max(band_depth).max(1.0e-9)
}

fn fill_coefficients(
    row: &mut [f64],
    row_index: usize,
    weight: f64,
    sample: &RowWeightedNormalSample,
    axis: usize,
) {
    let offset = row_index * INTERIOR_CONTROL_COUNT;
    for column in 0..INTERIOR_CONTROL_COUNT {
        row[offset + column] = weight * component(sample.normal_sensitivities[column], axis);
    }
}

fn add_normal_equation(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    coefficients: &[f64],
    target: f64,
) {
    for i in 0..coefficients.len() {
        if coefficients[i].abs() <= 1.0e-18 {
            continue;
        }
        rhs[i] += coefficients[i] * target;
        for j in 0..coefficients.len() {
            if coefficients[j].abs() > 1.0e-18 {
                matrix[i][j] += coefficients[i] * coefficients[j];
            }
        }
    }
}

fn component(vector: Vec3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        _ => vector.z,
    }
}

struct WeightedNormalSample {
    first: usize,
    second: usize,
    fraction: f64,
    chain_param: f64,
    impact: f64,
    detail: String,
}

struct RowWeightedNormalSample {
    normal: Vec3,
    curvature: f64,
    q_magnitude: f64,
    normal_sensitivities: [Vec3; INTERIOR_CONTROL_COUNT],
}
