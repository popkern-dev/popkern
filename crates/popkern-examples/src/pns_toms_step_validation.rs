use std::path::Path;

use popkern_examples::step_validation::{
    StepRoundTripReport, StepRoundTripValidationOptions, validate_step_bytes,
};

pub(super) fn validate_step_file(
    label: &str,
    path: &Path,
    expected_bodies: usize,
) -> Result<StepRoundTripReport, String> {
    let bytes = std::fs::read(path).map_err(|err| err.to_string())?;
    validate_step_bytes(
        label,
        &bytes,
        expected_bodies,
        &StepRoundTripValidationOptions::default(),
    )
    .map_err(|err| err.to_string())
}
