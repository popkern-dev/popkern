//! Paired partial-trim normal solve for constructive support-cell
//! rows. The solve is local to one authored sub-chain and may
//! move only interior row-1 controls `k = 2..=7`; boundary
//! curves, frozen legs, and dual-twist-adjacent controls remain
//! owned by their existing writers.

use popkern_core::geometry::Vec3;

use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::{CellReauthorOutcome, DerivedRowBand};

const INTERIOR_CONTROL_COUNT: usize = 6;
const NORMAL_SOLVE_REGULARIZATION: f64 = 1.0e-8;
const MIN_NORMAL_DISAGREEMENT: f64 = 1.0e-5;
const MAX_PAIR_NORMAL_CONTROL_FRACTION: f64 = 0.08;

pub(super) fn repair_pair_normal_solve(
    bands: &mut [DerivedRowBand],
    outcome: &mut CellReauthorOutcome,
) {
    if bands.len() < 2 {
        return;
    }
    let variable_count = bands.len() * INTERIOR_CONTROL_COUNT;
    let mut matrix = vec![vec![0.0_f64; variable_count]; variable_count];
    let mut rhs = vec![0.0_f64; variable_count];
    let mut constraint_count = 0usize;
    for first in 0..bands.len() {
        for second in first + 1..bands.len() {
            let Some(overlap) = partial_overlap(&bands[first], &bands[second]) else {
                continue;
            };
            for step in 1..16 {
                let chain_param = overlap.0 + (overlap.1 - overlap.0) * (step as f64 / 16.0);
                let Some(sample_a) = row_normal_sample(&bands[first], chain_param) else {
                    continue;
                };
                let Some(sample_b) = row_normal_sample(&bands[second], chain_param) else {
                    continue;
                };
                let orientation = if sample_a.normal.dot(&sample_b.normal) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let disagreement = 1.0 - sample_a.normal.dot(&sample_b.normal).abs();
                record_before(
                    outcome,
                    &bands[first],
                    &bands[second],
                    chain_param,
                    disagreement,
                );
                if disagreement <= MIN_NORMAL_DISAGREEMENT {
                    continue;
                }
                let diff = sample_a.normal - sample_b.normal * orientation;
                for axis in 0..3 {
                    let mut row = vec![0.0_f64; variable_count];
                    fill_sample_coefficients(&mut row, first, 1.0, &sample_a, axis);
                    fill_sample_coefficients(&mut row, second, -orientation, &sample_b, axis);
                    add_normal_equation(&mut matrix, &mut rhs, &row, -component(diff, axis));
                    constraint_count += 1;
                }
            }
        }
    }
    outcome.pair_normal_constraint_count += constraint_count;
    if constraint_count == 0 {
        record_after(bands, outcome);
        return;
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[index] += NORMAL_SOLVE_REGULARIZATION;
    }
    let Some(solution) = solve_dense(matrix, rhs) else {
        record_after(bands, outcome);
        return;
    };
    for (row_index, band) in bands.iter_mut().enumerate() {
        let limit = row_control_limit(band);
        let mut max_requested = 0.0_f64;
        for column in 0..INTERIOR_CONTROL_COUNT {
            max_requested =
                max_requested.max(solution[row_index * INTERIOR_CONTROL_COUNT + column].abs());
        }
        let scale = if max_requested > limit {
            outcome.pair_normal_limited_row_count += 1;
            limit / max_requested
        } else {
            1.0
        };
        let mut max_delta = 0.0_f64;
        for column in 0..INTERIOR_CONTROL_COUNT {
            let k = column + 2;
            let Some(direction) = control_normal_direction(band, k) else {
                continue;
            };
            let delta = direction * (solution[row_index * INTERIOR_CONTROL_COUNT + column] * scale);
            max_delta = max_delta.max(delta.dot(&delta).sqrt());
            band.row1[k] = band.row1[k] + delta;
        }
        if max_delta > 1.0e-14 {
            outcome.pair_normal_repaired_row_count += 1;
            outcome.max_pair_normal_control_delta =
                outcome.max_pair_normal_control_delta.max(max_delta);
        }
    }
    record_after(bands, outcome);
}

fn partial_overlap(a: &DerivedRowBand, b: &DerivedRowBand) -> Option<(f64, f64)> {
    let (lo_a, hi_a) = ordered(a.interval);
    let (lo_b, hi_b) = ordered(b.interval);
    let lo = lo_a.max(lo_b);
    let hi = hi_a.min(hi_b);
    if hi - lo <= 1.0e-6 {
        return None;
    }
    let span = hi - lo;
    let partial_a = span < hi_a - lo_a - 1.0e-6;
    let partial_b = span < hi_b - lo_b - 1.0e-6;
    (partial_a || partial_b).then_some((lo, hi))
}

fn row_normal_sample(band: &DerivedRowBand, chain_param: f64) -> Option<RowNormalSample> {
    let t = ((chain_param - band.interval.0) / (band.interval.1 - band.interval.0)).clamp(0.0, 1.0);
    let basis = bernstein9(t);
    let mut cross = Vec3::new(0.0, 0.0, 0.0);
    for (k, weight) in basis.iter().enumerate() {
        cross = cross + (band.row1[k] - band.row0[k]) * (9.0 * *weight);
    }
    let derivative = row_derivative_coefficients(&band.row0);
    let tangent = eval_ribbon_field(&derivative, t);
    let normal_vector = tangent.cross(&cross);
    let normal_magnitude = normal_vector.dot(&normal_vector).sqrt();
    let normal = normal_vector.normalized()?;
    let mut sensitivities = [Vec3::new(0.0, 0.0, 0.0); INTERIOR_CONTROL_COUNT];
    for column in 0..INTERIOR_CONTROL_COUNT {
        let k = column + 2;
        let Some(direction) = control_normal_direction(band, k) else {
            continue;
        };
        let cross_delta = direction * (9.0 * basis[k]);
        let normal_delta = tangent.cross(&cross_delta);
        sensitivities[column] =
            (normal_delta - normal * normal_delta.dot(&normal)) * (1.0 / normal_magnitude);
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

fn fill_sample_coefficients(
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

fn record_before(
    outcome: &mut CellReauthorOutcome,
    a: &DerivedRowBand,
    b: &DerivedRowBand,
    chain_param: f64,
    disagreement: f64,
) {
    outcome.pair_normal_sample_count += 1;
    if disagreement > outcome.max_pair_normal_disagreement_before {
        outcome.max_pair_normal_disagreement_before = disagreement;
        outcome.worst_pair_normal_detail = Some(format!(
            "{:?}:{:?}<->{:?}:{:?}@s={chain_param:.4};before={disagreement:.6e}",
            a.key.face, a.key.side, b.key.face, b.key.side
        ));
    }
}

fn record_after(bands: &[DerivedRowBand], outcome: &mut CellReauthorOutcome) {
    for first in 0..bands.len() {
        for second in first + 1..bands.len() {
            let Some(overlap) = partial_overlap(&bands[first], &bands[second]) else {
                continue;
            };
            for step in 1..16 {
                let chain_param = overlap.0 + (overlap.1 - overlap.0) * (step as f64 / 16.0);
                let Some(sample_a) = row_normal_sample(&bands[first], chain_param) else {
                    continue;
                };
                let Some(sample_b) = row_normal_sample(&bands[second], chain_param) else {
                    continue;
                };
                let disagreement = 1.0 - sample_a.normal.dot(&sample_b.normal).abs();
                outcome.max_pair_normal_disagreement_after =
                    outcome.max_pair_normal_disagreement_after.max(disagreement);
            }
        }
    }
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
    MAX_PAIR_NORMAL_CONTROL_FRACTION * boundary_length.max(band_depth).max(1.0e-9)
}

fn ordered(interval: (f64, f64)) -> (f64, f64) {
    (interval.0.min(interval.1), interval.0.max(interval.1))
}

fn component(vector: Vec3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        _ => vector.z,
    }
}

#[derive(Debug, Clone, Copy)]
struct RowNormalSample {
    normal: Vec3,
    sensitivities: [Vec3; INTERIOR_CONTROL_COUNT],
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;

    use super::super::super::BoundarySide;
    use super::super::super::chain_fit::ChainRowKey;
    use super::*;

    #[test]
    fn paired_solve_reduces_partial_overlap_normal_gap_without_moving_owned_controls() {
        let mut bands = [test_band(0.0, 1.0, 0.08), test_band(0.5, 1.0, 0.0)];
        let owned_a = [
            bands[0].row1[0],
            bands[0].row1[1],
            bands[0].row1[8],
            bands[0].row1[9],
        ];
        let owned_b = [
            bands[1].row1[0],
            bands[1].row1[1],
            bands[1].row1[8],
            bands[1].row1[9],
        ];
        let before = max_disagreement(&bands);
        let mut outcome = CellReauthorOutcome::default();

        repair_pair_normal_solve(&mut bands, &mut outcome);

        let after = max_disagreement(&bands);
        assert!(
            after < before,
            "normal gap must shrink: {before} -> {after}"
        );
        for (actual, expected) in [
            bands[0].row1[0],
            bands[0].row1[1],
            bands[0].row1[8],
            bands[0].row1[9],
        ]
        .iter()
        .zip(owned_a)
        {
            assert!(actual.distance_to(&expected) < 1.0e-15);
        }
        for (actual, expected) in [
            bands[1].row1[0],
            bands[1].row1[1],
            bands[1].row1[8],
            bands[1].row1[9],
        ]
        .iter()
        .zip(owned_b)
        {
            assert!(actual.distance_to(&expected) < 1.0e-15);
        }
        assert!(outcome.pair_normal_repaired_row_count > 0);
    }

    fn test_band(start: f64, end: f64, z: f64) -> DerivedRowBand {
        let mut row0 = [Point3::new(0.0, 0.0, 0.0); 10];
        let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
        for k in 0..10 {
            row0[k] = Point3::new(k as f64 / 9.0, 0.0, 0.0);
            row1[k] = row0[k] + Vec3::new(0.0, 0.2, z) * (1.0 / 9.0);
        }
        DerivedRowBand {
            key: ChainRowKey {
                face: popkern_core::topology::FaceHandle::default(),
                side: BoundarySide::UMin,
            },
            interval: (start, end),
            row0,
            row1,
            field_segment: [Vec3::new(0.0, 1.0, 0.0); 9],
        }
    }

    fn max_disagreement(bands: &[DerivedRowBand; 2]) -> f64 {
        let mut worst = 0.0_f64;
        let overlap = partial_overlap(&bands[0], &bands[1]).expect("partial overlap");
        for step in 1..16 {
            let chain_param = overlap.0 + (overlap.1 - overlap.0) * (step as f64 / 16.0);
            let a = row_normal_sample(&bands[0], chain_param).expect("sample a");
            let b = row_normal_sample(&bands[1], chain_param).expect("sample b");
            worst = worst.max(1.0 - a.normal.dot(&b.normal).abs());
        }
        worst
    }
}
