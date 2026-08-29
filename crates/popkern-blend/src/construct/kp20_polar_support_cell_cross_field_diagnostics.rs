//! Cross-field beta diagnostics for constructive support-cell
//! rows. These helpers measure the realized row-1 cross
//! derivative against the row-local `(C', R)` span without
//! changing geometry.

use popkern_core::geometry::Vec3;

use super::super::band_derive::solve_leg;
use super::super::cell_reauthor::CellReauthorOutcome;

const NEAR_TANGENT_FRACTION: f64 = 0.98;
const BETA_SIGN_FRACTION_TOL: f64 = 1.0e-6;

#[derive(Debug, Clone, Copy)]
pub(super) struct CrossFieldSample {
    pub(super) beta: f64,
    pub(super) tangential_fraction: f64,
    pub(super) beta_magnitude_fraction: f64,
    pub(super) cross_magnitude: f64,
    pub(super) field_magnitude: f64,
}

pub(super) fn sample_cross_field(
    cross: Vec3,
    tangent: Vec3,
    field: Vec3,
) -> Option<CrossFieldSample> {
    let cross_magnitude = cross.dot(&cross).sqrt();
    let tangent_magnitude = tangent.dot(&tangent).sqrt();
    let field_magnitude = field.dot(&field).sqrt();
    if cross_magnitude <= 1.0e-300 || tangent_magnitude <= 1.0e-300 || field_magnitude <= 1.0e-300 {
        return None;
    }
    let (_, beta) = solve_leg(cross, tangent, field)?;
    Some(CrossFieldSample {
        beta,
        tangential_fraction: cross.dot(&tangent).abs() / (tangent_magnitude * cross_magnitude),
        beta_magnitude_fraction: beta.abs() * field_magnitude / cross_magnitude,
        cross_magnitude,
        field_magnitude,
    })
}

pub(super) fn record_cross_field_sample(
    outcome: &mut CellReauthorOutcome,
    sample: CrossFieldSample,
    detail: impl FnOnce() -> String,
) {
    let first = outcome.cross_field_sample_count == 0;
    outcome.cross_field_sample_count += 1;
    if sample.tangential_fraction >= NEAR_TANGENT_FRACTION {
        outcome.near_tangent_cross_count += 1;
    }
    if first || sample.beta_magnitude_fraction < outcome.min_beta_magnitude_fraction {
        outcome.min_beta_magnitude_fraction = sample.beta_magnitude_fraction;
    }
    if first || sample.tangential_fraction > outcome.max_tangential_fraction {
        outcome.max_tangential_fraction = sample.tangential_fraction;
        outcome.worst_cross_field_detail = Some(detail());
    }
}

pub(super) fn record_beta_sign_flip(
    outcome: &mut CellReauthorOutcome,
    previous_sign: &mut i8,
    sample: CrossFieldSample,
) {
    let sign = beta_sign(sample);
    if *previous_sign != 0 && sign != 0 && sign != *previous_sign {
        outcome.beta_sign_flip_count += 1;
    }
    if sign != 0 {
        *previous_sign = sign;
    }
}

fn beta_sign(sample: CrossFieldSample) -> i8 {
    if sample.beta_magnitude_fraction <= BETA_SIGN_FRACTION_TOL {
        0
    } else if sample.beta < 0.0 {
        -1
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_exposes_tangent_dominated_beta() {
        let tangent = Vec3::new(1.0, 0.0, 0.0);
        let field = Vec3::new(0.0, 1.0, 0.0);
        let cross = Vec3::new(10.0, 0.1, 0.0);

        let sample = sample_cross_field(cross, tangent, field).expect("basis is non-degenerate");

        assert!(
            sample.tangential_fraction > 0.999,
            "cross should be tangent dominated: {sample:?}"
        );
        assert!(
            (sample.beta - 0.1).abs() < 1.0e-12,
            "beta must recover the transverse coefficient"
        );
        assert!(
            sample.beta_magnitude_fraction < 0.02,
            "beta fraction should expose transverse collapse: {sample:?}"
        );
        assert_eq!(beta_sign(sample), 1);
    }
}
