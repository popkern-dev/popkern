//! The build tolerance is a length in the model's unit (PK-462): a mm
//! model's default tolerance is 1e-6 mm, so a part means the same
//! precision relative to its features whatever unit it is authored in.

use popkern_history::unit::Unit;
use popkern_script::Model;

#[test]
fn a_millimetre_model_builds_at_a_millimetre_scaled_tolerance() {
    let si = Model::new();
    let mm = Model::new().unit(Unit::Millimeters);
    assert_eq!(si.tolerance().linear(), 1.0e-6);
    assert!((mm.tolerance().linear() - 1.0e-9).abs() < 1e-24);
}

#[test]
fn an_explicit_tolerance_is_read_in_the_model_unit() {
    use popkern_core::tolerance::Tolerance;
    // Authored as 1e-4 in a mm document: 1e-4 mm.
    let model =
        Model::with_tolerance(Tolerance::default().with_linear(1.0e-4)).unit(Unit::Millimeters);
    assert!((model.tolerance().linear() - 1.0e-7).abs() < 1e-22);
    // Changing unit again rescales from the current unit.
    let inches = model.unit(Unit::Inches);
    assert!((inches.tolerance().linear() - 1.0e-4 * 0.0254).abs() < 1e-18);
}

#[test]
fn non_length_units_leave_the_tolerance_alone() {
    let model = Model::new().unit(Unit::Degrees);
    assert_eq!(model.tolerance().linear(), 1.0e-6);
}
