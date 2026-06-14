//! Active-set normal solve for constructive support-cell row bands.
//!
//! Unlike the pair-normal solve, this first ranks actual sampled
//! shared trims in one sub-chain, then solves only the high-G1
//! active set while rejecting updates that move the sub-chain's
//! trim-sampled worst disagreement upward.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::topology::VertexHandle;

use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::{CellReauthorOutcome, DerivedRowBand};

const ACTIVE_NORMAL_SAMPLE_LIMIT: usize = 32;
const INTERIOR_CONTROL_COUNT: usize = 6;
const ACTIVE_NORMAL_REGULARIZATION: f64 = 1.0e-8;
const MIN_ACTIVE_NORMAL_DISAGREEMENT: f64 = 1.0e-5;
const MAX_ACTIVE_NORMAL_CONTROL_FRACTION: f64 = 0.04;
const ACTIVE_NORMAL_MAX_PASSES: usize = 3;
const ACTIVE_NORMAL_DAMPING_FACTORS: [f64; 4] = [1.0, 0.5, 0.25, 0.125];
const ACTIVE_NORMAL_IMPROVEMENT_EPS: f64 = 1.0e-12;

pub(super) fn repair_active_normal_solve(
    bands: &mut [DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
    outcome: &mut CellReauthorOutcome,
) {
    if bands.len() < 2 {
        record_decline(outcome, "insufficient_bands");
        return;
    }
    let mut before_worst = max_trim_disagreement(bands, chain_trims);
    let updates_global_worst = before_worst > outcome.max_active_normal_disagreement_before;
    outcome.max_active_normal_disagreement_before = outcome
        .max_active_normal_disagreement_before
        .max(before_worst);
    let mut kept_passes = 0usize;
    for pass_index in 0..ACTIVE_NORMAL_MAX_PASSES {
        let active_samples = top_active_samples(bands, chain_trims);
        outcome.active_normal_sample_count += active_samples.len();
        outcome.active_normal_constraint_count += active_samples.len() * 3;
        if active_samples.is_empty() {
            if kept_passes == 0 {
                record_decline(outcome, "no_active_samples");
            }
            break;
        }
        if pass_index == 0
            && updates_global_worst
            && let Some(sample) = active_samples.first()
        {
            outcome.worst_active_normal_detail = Some(format!(
                "{:?}:{:?}<->{:?}:{:?}@f={:.4},s={:.4};before={:.6e}",
                bands[sample.first].key.face,
                bands[sample.first].key.side,
                bands[sample.second].key.face,
                bands[sample.second].key.side,
                sample.fraction,
                sample.chain_param,
                sample.disagreement
            ));
        }
        let Some(solution) = solve_active_normal_system(bands, &active_samples) else {
            if kept_passes == 0 {
                record_decline(outcome, "active_normal_system_singular");
            }
            break;
        };
        let (deltas, limited_rows, max_requested_delta) = limited_row_deltas(bands, &solution);
        let mut best: Option<(f64, f64, Vec<DerivedRowBand>)> = None;
        for damping in ACTIVE_NORMAL_DAMPING_FACTORS {
            let mut candidate = bands.to_vec();
            apply_deltas(&mut candidate, &deltas, damping);
            let after_worst = max_trim_disagreement(&candidate, chain_trims);
            if after_worst + ACTIVE_NORMAL_IMPROVEMENT_EPS >= before_worst {
                continue;
            }
            if best
                .as_ref()
                .is_none_or(|(_, best_after, _)| after_worst < *best_after)
            {
                best = Some((damping, after_worst, candidate));
            }
        }
        let Some((damping, after_worst, candidate)) = best else {
            if kept_passes == 0 {
                record_decline(outcome, "no_predicted_improvement");
                outcome.active_normal_limited_row_count += limited_rows;
                outcome.max_active_normal_control_delta = outcome
                    .max_active_normal_control_delta
                    .max(max_requested_delta);
            }
            break;
        };
        let repaired_rows = write_repaired_rows(bands, candidate);
        kept_passes += 1;
        before_worst = after_worst;
        outcome.active_normal_decline_reason = Some("kept");
        outcome.active_normal_pass_count += 1;
        outcome.active_normal_repaired_row_count += repaired_rows;
        outcome.active_normal_limited_row_count += limited_rows;
        outcome.max_active_normal_control_delta = outcome
            .max_active_normal_control_delta
            .max(max_requested_delta * damping);
    }
    outcome.max_active_normal_disagreement_after = outcome
        .max_active_normal_disagreement_after
        .max(before_worst);
}

fn record_decline(outcome: &mut CellReauthorOutcome, reason: &'static str) {
    if outcome.active_normal_decline_reason != Some("kept") {
        outcome.active_normal_decline_reason = Some(reason);
    }
}

fn top_active_samples(
    bands: &[DerivedRowBand],
    chain_trims: &[ChainTrim<ChainRowKey, VertexHandle>],
) -> Vec<ActiveNormalSample> {
    let band_indices = band_indices(bands);
    let mut samples = Vec::new();
    for trim in chain_trims {
        let (Some(&first), Some(&second)) =
            (band_indices.get(&trim.row_a), band_indices.get(&trim.row_b))
        else {
            continue;
        };
        for step in 0..9 {
            let fraction = step as f64 / 8.0;
            let chain_param = trim_chain_param(&bands[first], trim.interval_a, fraction);
            let peer_chain_param =
                trim_chain_param(&bands[second], trim.interval_b, 1.0 - fraction);
            let chain_param = 0.5 * (chain_param + peer_chain_param);
            let Some(sample_a) = row_normal_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_normal_sample(&bands[second], chain_param) else {
                continue;
            };
            let disagreement = 1.0 - sample_a.normal.dot(&sample_b.normal).abs();
            if disagreement <= MIN_ACTIVE_NORMAL_DISAGREEMENT {
                continue;
            }
            samples.push(ActiveNormalSample {
                first,
                second,
                fraction,
                chain_param,
                disagreement,
            });
        }
    }
    samples.sort_by(|a, b| {
        b.disagreement
            .total_cmp(&a.disagreement)
            .then_with(|| a.first.cmp(&b.first))
            .then_with(|| a.second.cmp(&b.second))
            .then_with(|| a.fraction.total_cmp(&b.fraction))
    });
    samples.truncate(ACTIVE_NORMAL_SAMPLE_LIMIT);
    samples
}

fn solve_active_normal_system(
    bands: &[DerivedRowBand],
    samples: &[ActiveNormalSample],
) -> Option<Vec<f64>> {
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    let mut matrix = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut rhs = vec![0.0_f64; variable_count];
    for sample in samples {
        let sample_a = row_normal_sample(&bands[sample.first], sample.chain_param)?;
        let sample_b = row_normal_sample(&bands[sample.second], sample.chain_param)?;
        let orientation = if sample_a.normal.dot(&sample_b.normal) < 0.0 {
            -1.0
        } else {
            1.0
        };
        let diff = sample_a.normal - sample_b.normal * orientation;
        for axis in 0..3 {
            let mut row = vec![0.0_f64; variable_count];
            fill_coefficients(&mut row, sample.first, 1.0, &sample_a, axis);
            fill_coefficients(&mut row, sample.second, -orientation, &sample_b, axis);
            add_normal_equation(&mut matrix, &mut rhs, &row, -component(diff, axis));
        }
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += ACTIVE_NORMAL_REGULARIZATION;
    }
    solve_dense(matrix, rhs)
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
            let chain_param = trim_chain_param(&bands[first], trim.interval_a, fraction);
            let peer_chain_param =
                trim_chain_param(&bands[second], trim.interval_b, 1.0 - fraction);
            let chain_param = 0.5 * (chain_param + peer_chain_param);
            let Some(sample_a) = row_normal_sample(&bands[first], chain_param) else {
                continue;
            };
            let Some(sample_b) = row_normal_sample(&bands[second], chain_param) else {
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

fn trim_chain_param(band: &DerivedRowBand, trim_interval: (f64, f64), fraction: f64) -> f64 {
    let row_param = trim_interval.0 + (trim_interval.1 - trim_interval.0) * fraction;
    band.interval.0 + (band.interval.1 - band.interval.0) * row_param
}

fn row_normal_sample(band: &DerivedRowBand, chain_param: f64) -> Option<RowNormalSample> {
    let t = ((chain_param - band.interval.0) / (band.interval.1 - band.interval.0)).clamp(0.0, 1.0);
    let basis = bernstein9(t);
    let derivative = row_derivative_coefficients(&band.row0);
    let tangent = eval_ribbon_field(&derivative, t);
    let mut cross = Vec3::new(0.0, 0.0, 0.0);
    for (k, weight) in basis.iter().enumerate() {
        cross = cross + (band.row1[k] - band.row0[k]) * (9.0 * *weight);
    }
    let raw_normal = tangent.cross(&cross);
    let raw_normal_magnitude = raw_normal.dot(&raw_normal).sqrt();
    let normal = raw_normal.normalized()?;
    let mut sensitivities = [Vec3::new(0.0, 0.0, 0.0); INTERIOR_CONTROL_COUNT];
    for column in 0..INTERIOR_CONTROL_COUNT {
        let k = column + 2;
        let Some(direction) = control_normal_direction(band, k) else {
            continue;
        };
        let cross_delta = direction * (9.0 * basis[k]);
        let normal_delta = tangent.cross(&cross_delta);
        sensitivities[column] =
            (normal_delta - normal * normal_delta.dot(&normal)) * (1.0 / raw_normal_magnitude);
    }
    Some(RowNormalSample {
        normal,
        sensitivities,
    })
}

fn control_normal_direction(band: &DerivedRowBand, k: usize) -> Option<Vec3> {
    let t = k as f64 / 9.0;
    let derivative = row_derivative_coefficients(&band.row0);
    let tangent = eval_ribbon_field(&derivative, t);
    let field = eval_ribbon_field(&band.field_segment, t);
    tangent.cross(&field).normalized()
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
    MAX_ACTIVE_NORMAL_CONTROL_FRACTION * boundary_length.max(band_depth).max(1.0e-9)
}

fn fill_coefficients(
    row: &mut [f64],
    row_index: usize,
    orientation: f64,
    sample: &RowNormalSample,
    axis: usize,
) {
    let offset = row_index * INTERIOR_CONTROL_COUNT;
    for column in 0..INTERIOR_CONTROL_COUNT {
        row[offset + column] = orientation * component(sample.sensitivities[column], axis);
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

#[derive(Debug, Clone)]
struct ActiveNormalSample {
    first: usize,
    second: usize,
    fraction: f64,
    chain_param: f64,
    disagreement: f64,
}

#[derive(Debug, Clone, Copy)]
struct RowNormalSample {
    normal: Vec3,
    sensitivities: [Vec3; INTERIOR_CONTROL_COUNT],
}
