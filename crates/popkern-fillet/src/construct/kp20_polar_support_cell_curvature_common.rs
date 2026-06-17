//! Shared row-band curvature sampling and delta application.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::DerivedRowBand;

pub(super) const INTERIOR_CONTROL_COUNT: usize = 6;
const CURVATURE_REGULARIZATION: f64 = 1.0e-8;

pub(super) fn top_curvature_samples(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    active_residual: f64,
    sample_limit: usize,
    subdivisions: usize,
) -> Vec<CurvatureSample> {
    let band_indices = band_indices(bands);
    let mut samples = Vec::new();
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 1..subdivisions {
            let fraction = step as f64 / subdivisions as f64;
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
            if residual <= active_residual {
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
    samples.truncate(sample_limit);
    samples
}

pub(super) fn solve_curvature_system(
    bands: &[DerivedRowBand],
    samples: &[CurvatureSample],
) -> Option<Vec<f64>> {
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    let mut matrix = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut rhs = vec![0.0_f64; variable_count];
    add_curvature_equations(bands, samples, &mut matrix, &mut rhs)?;
    regularize(&mut matrix);
    solve_dense(matrix, rhs)
}

pub(super) fn add_curvature_equations(
    bands: &[DerivedRowBand],
    samples: &[CurvatureSample],
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
) -> Option<()> {
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    for sample in samples {
        let sample_a = row_curvature_sample(&bands[sample.first], sample.chain_param)?;
        let sample_b = row_curvature_sample(&bands[sample.second], sample.chain_param)?;
        let aligned_b = sample_b.curvature * sample.normal_alignment;
        let scale = 1.0 + sample_a.curvature.abs().max(aligned_b.abs());
        let diff = (sample_a.curvature - aligned_b) / scale;
        let mut row = vec![0.0_f64; variable_count];
        fill_curvature_coefficients(&mut row, sample.first, 1.0 / scale, &sample_a);
        fill_curvature_coefficients(
            &mut row,
            sample.second,
            -sample.normal_alignment / scale,
            &sample_b,
        );
        add_normal_equation(matrix, rhs, &row, -diff);
    }
    Some(())
}

pub(super) fn regularize(matrix: &mut [Vec<f64>]) {
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += CURVATURE_REGULARIZATION;
    }
}

pub(super) fn max_normal_aligned_g2_residual(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    subdivisions: usize,
) -> f64 {
    let band_indices = band_indices(bands);
    let mut worst = 0.0_f64;
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 1..subdivisions {
            let fraction = step as f64 / subdivisions as f64;
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

pub(super) fn max_trim_disagreement(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    subdivisions: usize,
) -> f64 {
    let band_indices = band_indices(bands);
    let mut worst = 0.0_f64;
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 0..=subdivisions {
            let fraction = step as f64 / subdivisions as f64;
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

pub(super) fn row_curvature_sample(
    band: &DerivedRowBand,
    chain_param: f64,
) -> Option<RowCurvatureSample> {
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
    for (column, sensitivity) in curvature_sensitivities.iter_mut().enumerate() {
        let k = column + 2;
        let Some(direction) = control_normal_direction(band, k) else {
            continue;
        };
        let cross_delta = direction * (9.0 * basis[k]);
        let normal_delta = tangent.cross(&cross_delta);
        let normal_sensitivity =
            (normal_delta - normal * normal_delta.dot(&normal)) * (1.0 / raw_normal_magnitude);
        *sensitivity = normal_sensitivity.dot(&second) / speed_sq;
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

pub(super) fn normal_alignment(first: Vec3, second: Vec3) -> f64 {
    if first.dot(&second) < 0.0 { -1.0 } else { 1.0 }
}

fn normal_aligned_g2_residual(first: f64, second: f64, alignment: f64) -> f64 {
    let aligned_second = second * alignment;
    let scale = first.abs().max(aligned_second.abs());
    (first - aligned_second).abs() / (1.0 + scale)
}

pub(super) fn limited_row_deltas(
    bands: &[DerivedRowBand],
    solution: &[f64],
    max_control_fraction: f64,
) -> (Vec<[f64; 6]>, usize, f64) {
    let mut deltas = vec![[0.0_f64; 6]; bands.len()];
    let mut limited_rows = 0usize;
    let mut max_delta = 0.0_f64;
    for (row_index, band) in bands.iter().enumerate() {
        let limit = row_control_limit(band, max_control_fraction);
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

pub(super) fn apply_deltas(bands: &mut [DerivedRowBand], deltas: &[[f64; 6]], damping: f64) {
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

pub(super) fn write_repaired_rows(
    bands: &mut [DerivedRowBand],
    candidate: Vec<DerivedRowBand>,
) -> usize {
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

fn row_control_limit(band: &DerivedRowBand, max_control_fraction: f64) -> f64 {
    let boundary_length: f64 = (1..10)
        .map(|k| band.row0[k].distance_to(&band.row0[k - 1]))
        .sum();
    let band_depth = band
        .row1
        .iter()
        .zip(band.row0.iter())
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max);
    max_control_fraction * boundary_length.max(band_depth).max(1.0e-9)
}

pub(super) fn fill_curvature_coefficients(
    row: &mut [f64],
    row_index: usize,
    weight: f64,
    sample: &RowCurvatureSample,
) {
    let offset = row_index * INTERIOR_CONTROL_COUNT;
    for column in 0..INTERIOR_CONTROL_COUNT {
        row[offset + column] = weight * sample.curvature_sensitivities[column];
    }
}

pub(super) fn add_normal_equation(
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

#[derive(Debug, Clone)]
pub(super) struct CurvatureSample {
    pub(super) first: usize,
    pub(super) second: usize,
    pub(super) fraction: f64,
    pub(super) chain_param: f64,
    pub(super) normal_alignment: f64,
    pub(super) residual: f64,
    pub(super) detail: String,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct RowCurvatureSample {
    pub(super) normal: Vec3,
    pub(super) curvature: f64,
    pub(super) curvature_sensitivities: [f64; INTERIOR_CONTROL_COUNT],
}
