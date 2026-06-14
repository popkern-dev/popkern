//! Interior beta-envelope repair for constructive support-cell
//! row-1 bands. The solve is deliberately local: it may move
//! only controls `k = 2..=7`, leaving boundary curves, frozen end
//! legs, and dual-twist-adjacent ownership untouched.

use popkern_core::geometry::{Point3, Vec3};

use super::super::band_derive::solve_leg;
use super::super::chain_ribbon_field::{eval_ribbon_field, row_derivative_coefficients};
use super::super::chain_solver::solve_dense;
use super::super::row_fit::bernstein9;
use super::CellReauthorOutcome;

const MIN_BETA_MAGNITUDE_FRACTION: f64 = 0.12;
const BETA_ENVELOPE_REGULARIZATION: f64 = 1.0e-8;

pub(super) fn repair_beta_envelope(
    row0: &[Point3; 10],
    row1: &mut [Point3; 10],
    field_segment: &[Vec3; 9],
    outcome: &mut CellReauthorOutcome,
) {
    let before = row_beta_samples(row0, row1, field_segment);
    for sample in &before {
        record_before(outcome, sample.beta_magnitude_fraction);
    }
    let dominant_sign = dominant_beta_sign(&before);
    let constraints = beta_constraints(&before, dominant_sign);
    outcome.beta_envelope_constraint_count += constraints.len();
    if constraints.is_empty() {
        for sample in &before {
            record_after(outcome, sample.beta_magnitude_fraction);
        }
        return;
    }
    let derivative = row_derivative_coefficients(row0);
    let control_fields = interior_control_fields(field_segment);
    let mut matrix = vec![vec![0.0_f64; 6]; 6];
    let mut rhs = vec![0.0_f64; 6];
    for constraint in constraints {
        let tangent = eval_ribbon_field(&derivative, constraint.t);
        let field = eval_ribbon_field(field_segment, constraint.t);
        let basis = bernstein9(constraint.t);
        let mut row = [0.0_f64; 6];
        for (column, control_field) in control_fields.iter().enumerate() {
            let vector = *control_field * basis[column + 2];
            row[column] = solve_leg(vector, tangent, field)
                .map(|(_, beta)| beta)
                .unwrap_or(0.0);
        }
        for i in 0..6 {
            rhs[i] += row[i] * constraint.delta_beta;
            for j in 0..6 {
                matrix[i][j] += row[i] * row[j];
            }
        }
    }
    for (i, row) in matrix.iter_mut().enumerate() {
        row[i] += BETA_ENVELOPE_REGULARIZATION;
    }
    let Some(solution) = solve_dense(matrix, rhs) else {
        for sample in &before {
            record_after(outcome, sample.beta_magnitude_fraction);
        }
        return;
    };
    let mut max_delta = 0.0_f64;
    for (column, value) in solution.iter().enumerate() {
        let k = column + 2;
        let delta = control_fields[column] * (*value / 9.0);
        max_delta = max_delta.max(delta.dot(&delta).sqrt());
        row1[k] = row1[k] + delta;
    }
    if max_delta > 0.0 {
        outcome.beta_envelope_repaired_row_count += 1;
        outcome.max_beta_envelope_control_delta =
            outcome.max_beta_envelope_control_delta.max(max_delta);
    }
    for sample in row_beta_samples(row0, row1, field_segment) {
        record_after(outcome, sample.beta_magnitude_fraction);
    }
}

fn row_beta_samples(
    row0: &[Point3; 10],
    row1: &[Point3; 10],
    field_segment: &[Vec3; 9],
) -> Vec<BetaEnvelopeSample> {
    let derivative = row_derivative_coefficients(row0);
    let mut samples = Vec::new();
    for step in 1..16 {
        let t = step as f64 / 16.0;
        let basis = bernstein9(t);
        let mut cross = Vec3::new(0.0, 0.0, 0.0);
        for k in 0..10 {
            cross = cross + (row1[k] - row0[k]) * (9.0 * basis[k]);
        }
        let tangent = eval_ribbon_field(&derivative, t);
        let field = eval_ribbon_field(field_segment, t);
        let cross_magnitude = cross.dot(&cross).sqrt();
        let field_magnitude = field.dot(&field).sqrt();
        if cross_magnitude <= 1.0e-300 || field_magnitude <= 1.0e-300 {
            continue;
        }
        let Some((_, beta)) = solve_leg(cross, tangent, field) else {
            continue;
        };
        samples.push(BetaEnvelopeSample {
            t,
            beta,
            beta_magnitude_fraction: beta.abs() * field_magnitude / cross_magnitude,
            target_beta_magnitude: MIN_BETA_MAGNITUDE_FRACTION * cross_magnitude / field_magnitude,
        });
    }
    samples
}

fn dominant_beta_sign(samples: &[BetaEnvelopeSample]) -> f64 {
    let signed: f64 = samples
        .iter()
        .map(|sample| sample.beta.signum() * sample.beta_magnitude_fraction)
        .sum();
    if signed < 0.0 { -1.0 } else { 1.0 }
}

fn beta_constraints(
    samples: &[BetaEnvelopeSample],
    dominant_sign: f64,
) -> Vec<BetaEnvelopeConstraint> {
    samples
        .iter()
        .filter_map(|sample| {
            let desired_beta = dominant_sign * sample.target_beta_magnitude;
            let needs_lift = sample.beta * dominant_sign < sample.target_beta_magnitude
                || sample.beta_magnitude_fraction < MIN_BETA_MAGNITUDE_FRACTION;
            needs_lift.then_some(BetaEnvelopeConstraint {
                t: sample.t,
                delta_beta: desired_beta - sample.beta,
            })
        })
        .collect()
}

fn interior_control_fields(field_segment: &[Vec3; 9]) -> [Vec3; 6] {
    [
        eval_ribbon_field(field_segment, 2.0 / 9.0),
        eval_ribbon_field(field_segment, 3.0 / 9.0),
        eval_ribbon_field(field_segment, 4.0 / 9.0),
        eval_ribbon_field(field_segment, 5.0 / 9.0),
        eval_ribbon_field(field_segment, 6.0 / 9.0),
        eval_ribbon_field(field_segment, 7.0 / 9.0),
    ]
}

fn record_before(outcome: &mut CellReauthorOutcome, fraction: f64) {
    if outcome.beta_envelope_sample_count == 0
        || fraction < outcome.min_beta_envelope_fraction_before
    {
        outcome.min_beta_envelope_fraction_before = fraction;
    }
    outcome.beta_envelope_sample_count += 1;
}

fn record_after(outcome: &mut CellReauthorOutcome, fraction: f64) {
    if outcome.beta_envelope_after_sample_count == 0
        || fraction < outcome.min_beta_envelope_fraction_after
    {
        outcome.min_beta_envelope_fraction_after = fraction;
    }
    outcome.beta_envelope_after_sample_count += 1;
}

#[derive(Debug, Clone, Copy)]
struct BetaEnvelopeSample {
    t: f64,
    beta: f64,
    beta_magnitude_fraction: f64,
    target_beta_magnitude: f64,
}

#[derive(Debug, Clone, Copy)]
struct BetaEnvelopeConstraint {
    t: f64,
    delta_beta: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repair_lifts_interior_beta_without_moving_owned_controls() {
        let mut row0 = [Point3::new(0.0, 0.0, 0.0); 10];
        let mut row1 = [Point3::new(0.0, 0.0, 0.0); 10];
        for k in 0..10 {
            row0[k] = Point3::new(k as f64 / 9.0, 0.0, 0.0);
            row1[k] = row0[k] + Vec3::new(1.0, 0.01, 0.0) * (1.0 / 9.0);
        }
        let field = [Vec3::new(0.0, 1.0, 0.0); 9];
        let original = row1;
        let before = row_beta_samples(&row0, &row1, &field)
            .into_iter()
            .map(|sample| sample.beta_magnitude_fraction)
            .fold(f64::INFINITY, f64::min);
        let mut outcome = CellReauthorOutcome::default();

        repair_beta_envelope(&row0, &mut row1, &field, &mut outcome);

        let after = row_beta_samples(&row0, &row1, &field)
            .into_iter()
            .map(|sample| sample.beta_magnitude_fraction)
            .fold(f64::INFINITY, f64::min);
        assert!(
            after > before,
            "repair should lift beta: {before} -> {after}"
        );
        for k in [0usize, 1, 8, 9] {
            assert!(
                row1[k].distance_to(&original[k]) < 1.0e-15,
                "owned control {k} must not move"
            );
        }
        assert_eq!(outcome.beta_envelope_repaired_row_count, 1);
    }
}
