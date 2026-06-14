//! Active-set normal-curvature solve for constructive row bands.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::DerivedRowBand;

const CURVATURE_SAMPLE_LIMIT: usize = 24;
const INTERIOR_CONTROL_COUNT: usize = 6;
const CURVATURE_REGULARIZATION: f64 = 1.0e-8;
const MIN_NORMAL_ALIGNED_G2_RESIDUAL: f64 = 5.0e-2;
const COMMERCIAL_G1_TARGET: f64 = 1.25e-3;
const MAX_CURVATURE_CONTROL_FRACTION: f64 = 0.03;
const CURVATURE_MAX_PASSES: usize = 4;
const CURVATURE_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];
const CURVATURE_IMPROVEMENT_EPS: f64 = 1.0e-12;

pub(super) fn repair_curvature_solve(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> String {
    if bands.len() < 2 {
        return declined_detail("insufficient_bands", 0.0, 0.0);
    }
    let mut before_g2 = max_normal_aligned_g2_residual(bands, chain_trims);
    let mut before_g1 = max_trim_disagreement(bands, chain_trims);
    let initial_g2 = before_g2;
    let initial_g1 = before_g1;
    let mut pass_count = 0usize;
    let mut total_samples = 0usize;
    let mut total_constraints = 0usize;
    let mut total_repaired_rows = 0usize;
    let mut total_limited_rows = 0usize;
    let mut max_control_delta = 0.0_f64;
    let mut best_detail = top_curvature_samples(bands, chain_trims)
        .first()
        .map(|sample| sample.detail.clone())
        .unwrap_or_else(|| "none".to_string());
    let mut decline_reason = "no_active_samples";

    for _ in 0..CURVATURE_MAX_PASSES {
        let active_samples = top_curvature_samples(bands, chain_trims);
        total_samples += active_samples.len();
        total_constraints += active_samples.len();
        if active_samples.is_empty() {
            decline_reason = if pass_count == 0 {
                "no_active_samples"
            } else {
                "converged_no_active_samples"
            };
            break;
        }
        if let Some(sample) = active_samples.first() {
            best_detail = sample.detail.clone();
        }
        let Some(solution) = solve_curvature_system(bands, &active_samples) else {
            decline_reason = "curvature_system_singular";
            break;
        };
        let (deltas, limited_rows, max_requested_delta) = limited_row_deltas(bands, &solution);
        let mut best: Option<(f64, f64, f64, Vec<DerivedRowBand>)> = None;
        for damping in CURVATURE_DAMPING_FACTORS {
            let mut candidate = bands.to_vec();
            apply_deltas(&mut candidate, &deltas, damping);
            let after_g2 = max_normal_aligned_g2_residual(&candidate, chain_trims);
            if after_g2 + CURVATURE_IMPROVEMENT_EPS >= before_g2 {
                continue;
            }
            let after_g1 = max_trim_disagreement(&candidate, chain_trims);
            if after_g1 > before_g1.max(COMMERCIAL_G1_TARGET) + 1.0e-10 {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_g2, _, _)| after_g2 < *best_g2)
            {
                best = Some((damping, after_g2, after_g1, candidate));
            }
        }
        let Some((damping, after_g2, after_g1, candidate)) = best else {
            decline_reason = "no_quality_safe_improvement";
            total_limited_rows += limited_rows;
            max_control_delta = max_control_delta.max(max_requested_delta);
            break;
        };
        total_repaired_rows += write_repaired_rows(bands, candidate);
        total_limited_rows += limited_rows;
        max_control_delta = max_control_delta.max(max_requested_delta * damping);
        before_g2 = after_g2;
        before_g1 = after_g1;
        pass_count += 1;
        decline_reason = "kept";
    }

    format!(
        "curvature_solve_status={};curvature_sample_count={};curvature_pass_count={};curvature_constraint_count={};curvature_repaired_row_count={};curvature_limited_row_count={};max_curvature_control_delta={:.6e};max_curvature_residual_before={:.15e};max_curvature_residual_after={:.15e};max_curvature_g1_before={:.15e};max_curvature_g1_after={:.15e};worst_curvature_detail=[{}]",
        decline_reason,
        total_samples,
        pass_count,
        total_constraints,
        total_repaired_rows,
        total_limited_rows,
        max_control_delta,
        initial_g2,
        before_g2,
        initial_g1,
        before_g1,
        best_detail
    )
}

fn declined_detail(reason: &'static str, initial_g2: f64, initial_g1: f64) -> String {
    format!(
        "curvature_solve_status={reason};curvature_sample_count=0;curvature_pass_count=0;curvature_constraint_count=0;curvature_repaired_row_count=0;curvature_limited_row_count=0;max_curvature_control_delta=0.000000e0;max_curvature_residual_before={initial_g2:.15e};max_curvature_residual_after={initial_g2:.15e};max_curvature_g1_before={initial_g1:.15e};max_curvature_g1_after={initial_g1:.15e};worst_curvature_detail=[none]"
    )
}

fn top_curvature_samples(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> Vec<CurvatureSample> {
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
            let Some(sample_a) = row_curvature_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_curvature_sample(&bands[second], chain_param) else {
                continue;
            };
            let normal_alignment = normal_alignment(sample_a.normal, sample_b.normal);
            let normal_disagreement = 1.0 - sample_a.normal.dot(&sample_b.normal).abs();
            let residual = normal_aligned_g2_residual(
                sample_a.curvature,
                sample_b.curvature,
                normal_alignment,
            );
            if residual <= MIN_NORMAL_ALIGNED_G2_RESIDUAL {
                continue;
            }
            samples.push(CurvatureSample {
                first,
                second,
                fraction,
                chain_param,
                normal_alignment,
                residual,
                detail: format!(
                    "{:?}:{:?}<->{:?}:{:?}@f={:.4},s={:.4};g2_aligned={:.6e};ka={:.6e};kb={:.6e};normal_alignment={:.1};g1={:.6e}",
                    bands[first].key.face,
                    bands[first].key.side,
                    bands[second].key.face,
                    bands[second].key.side,
                    fraction,
                    chain_param,
                    residual,
                    sample_a.curvature,
                    sample_b.curvature,
                    normal_alignment,
                    normal_disagreement
                ),
            });
        }
    }
    samples.sort_by(|a, b| {
        b.residual
            .total_cmp(&a.residual)
            .then_with(|| a.first.cmp(&b.first))
            .then_with(|| a.second.cmp(&b.second))
            .then_with(|| a.fraction.total_cmp(&b.fraction))
    });
    samples.truncate(CURVATURE_SAMPLE_LIMIT);
    samples
}

fn solve_curvature_system(
    bands: &[DerivedRowBand],
    samples: &[CurvatureSample],
) -> Option<Vec<f64>> {
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    let mut matrix = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut rhs = vec![0.0_f64; variable_count];
    for sample in samples {
        let sample_a = row_curvature_sample(&bands[sample.first], sample.chain_param)?;
        let sample_b = row_curvature_sample(&bands[sample.second], sample.chain_param)?;
        let aligned_b = sample_b.curvature * sample.normal_alignment;
        let scale = 1.0 + sample_a.curvature.abs().max(aligned_b.abs());
        let diff = (sample_a.curvature - aligned_b) / scale;
        let mut row = vec![0.0_f64; variable_count];
        fill_coefficients(&mut row, sample.first, 1.0 / scale, &sample_a);
        fill_coefficients(
            &mut row,
            sample.second,
            -sample.normal_alignment / scale,
            &sample_b,
        );
        add_normal_equation(&mut matrix, &mut rhs, &row, -diff);
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += CURVATURE_REGULARIZATION;
    }
    solve_dense(matrix, rhs)
}

fn max_normal_aligned_g2_residual(
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
        for step in 1..8 {
            let fraction = step as f64 / 8.0;
            let chain_param = paired_chain_param(&bands[first], &bands[second], trim, fraction);
            let Some(sample_a) = row_curvature_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_curvature_sample(&bands[second], chain_param) else {
                continue;
            };
            let alignment = normal_alignment(sample_a.normal, sample_b.normal);
            worst = worst.max(normal_aligned_g2_residual(
                sample_a.curvature,
                sample_b.curvature,
                alignment,
            ));
        }
    }
    worst
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
            let Some(sample_a) = row_curvature_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_curvature_sample(&bands[second], chain_param) else {
                continue;
            };
            worst = worst.max(1.0 - sample_a.normal.dot(&sample_b.normal).abs());
        }
    }
    worst
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

fn row_curvature_sample(band: &DerivedRowBand, chain_param: f64) -> Option<RowCurvatureSample> {
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
    let curvature = normal.dot(&second) / speed_sq;
    let mut curvature_sensitivities = [0.0_f64; INTERIOR_CONTROL_COUNT];
    for column in 0..INTERIOR_CONTROL_COUNT {
        let k = column + 2;
        let Some(direction) = control_normal_direction(band, k) else {
            continue;
        };
        let cross_delta = direction * (9.0 * basis[k]);
        let normal_delta = tangent.cross(&cross_delta);
        let normal_sensitivity =
            (normal_delta - normal * normal_delta.dot(&normal)) * (1.0 / raw_normal_magnitude);
        curvature_sensitivities[column] = normal_sensitivity.dot(&second) / speed_sq;
    }
    Some(RowCurvatureSample {
        normal,
        curvature,
        curvature_sensitivities,
    })
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

fn normal_aligned_g2_residual(first: f64, second: f64, alignment: f64) -> f64 {
    let aligned_second = second * alignment;
    let scale = first.abs().max(aligned_second.abs());
    (first - aligned_second).abs() / (1.0 + scale)
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
    MAX_CURVATURE_CONTROL_FRACTION * boundary_length.max(band_depth).max(1.0e-9)
}

fn fill_coefficients(row: &mut [f64], row_index: usize, weight: f64, sample: &RowCurvatureSample) {
    let offset = row_index * INTERIOR_CONTROL_COUNT;
    for column in 0..INTERIOR_CONTROL_COUNT {
        row[offset + column] = weight * sample.curvature_sensitivities[column];
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

#[derive(Debug, Clone)]
struct CurvatureSample {
    first: usize,
    second: usize,
    fraction: f64,
    chain_param: f64,
    normal_alignment: f64,
    residual: f64,
    detail: String,
}

#[derive(Debug, Clone, Copy)]
struct RowCurvatureSample {
    normal: Vec3,
    curvature: f64,
    curvature_sensitivities: [f64; INTERIOR_CONTROL_COUNT],
}
