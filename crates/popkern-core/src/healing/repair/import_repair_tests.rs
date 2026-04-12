use super::*;
use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Helper: check if a step label is a specific standard op.
fn is_standard(label: &StepLabel, expected: ImportRepairOp) -> bool {
    matches!(label, StepLabel::Standard(op) if *op == expected)
}

/// Helper to extract a StepLabel from a config's operations.
fn step_label(config: &ImportRepairConfig, idx: usize) -> StepLabel {
    match &config.operations[idx] {
        RepairOperation::Standard(op) => StepLabel::Standard(*op),
        RepairOperation::Custom(op) => StepLabel::Custom(op.name().to_string()),
    }
}

#[test]
fn default_config_on_valid_cube_does_not_increase_defects() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
    let config = ImportRepairConfig::default_config(tol);
    let report = run_import_repair(&mut brep, &config);

    // The pipeline must not make things worse.
    assert!(
        report.total_defects_after <= report.total_defects_before,
        "import repair should not increase defects: before={}, after={}",
        report.total_defects_before,
        report.total_defects_after,
    );
    // Every step should also not increase defects.
    for (label, before, after) in &report.steps {
        assert!(
            *after <= *before,
            "{label:?} increased defects: {before} -> {after}"
        );
    }
}

#[test]
fn empty_ops_returns_immediately() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
    let config = ImportRepairConfig {
        operations: vec![],
        max_displacement: tol.linear() * 10.0,
        tolerance: tol,
    };
    let report = run_import_repair(&mut brep, &config);

    assert!(
        report.steps.is_empty(),
        "empty operations should produce no steps"
    );
    assert_eq!(report.total_defects_before, report.total_defects_after);
}

#[test]
fn default_config_has_correct_operation_order() {
    let tol = Tolerance::default();
    let config = ImportRepairConfig::default_config(tol);

    assert_eq!(config.operations.len(), 6);
    assert!(is_standard(
        &step_label(&config, 0),
        ImportRepairOp::OrientationFix
    ));
    assert!(is_standard(
        &step_label(&config, 1),
        ImportRepairOp::PeriodicTopologyRestore
    ));
    assert!(is_standard(
        &step_label(&config, 2),
        ImportRepairOp::ToleranceNormalization
    ));
    assert!(is_standard(
        &step_label(&config, 3),
        ImportRepairOp::GapClosure
    ));
    assert!(is_standard(
        &step_label(&config, 4),
        ImportRepairOp::SurfaceSimplification
    ));
    assert!(is_standard(
        &step_label(&config, 5),
        ImportRepairOp::ToleranceEnforcement
    ));
}

#[test]
fn max_displacement_scales_with_tolerance() {
    let tol = Tolerance::default();
    let config = ImportRepairConfig::default_config(tol);
    let expected = tol.linear() * 10.0;
    assert!(
        (config.max_displacement - expected).abs() < f64::EPSILON,
        "max_displacement should be 10x linear tolerance"
    );
}

#[test]
fn single_op_pipeline_runs() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
    let config = ImportRepairConfig {
        operations: vec![RepairOperation::Standard(
            ImportRepairOp::ToleranceEnforcement,
        )],
        max_displacement: tol.linear() * 10.0,
        tolerance: tol,
    };
    let report = run_import_repair(&mut brep, &config);
    assert_eq!(report.steps.len(), 1);
    assert!(is_standard(
        &report.steps[0].0,
        ImportRepairOp::ToleranceEnforcement
    ));
}

#[test]
fn empty_brep_returns_no_defects() {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    let config = ImportRepairConfig::default_config(tol);
    let report = run_import_repair(&mut brep, &config);

    assert_eq!(report.total_defects_before, 0);
    assert_eq!(report.total_defects_after, 0);
}

// --- Custom repair op tests ---

/// A test custom op that does nothing (succeeds).
struct NoopCustomOp {
    name: String,
}

impl CustomRepairOp for NoopCustomOp {
    fn name(&self) -> &str {
        &self.name
    }

    fn execute(&self, _brep: &mut BRep, _tolerance: &Tolerance) -> Result<(), CoreError> {
        Ok(())
    }
}

#[test]
fn custom_op_can_be_inserted_and_executed() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
    let mut config = ImportRepairConfig::default_config(tol);

    // Insert a custom op at position 2 (after orientation + periodic).
    config.insert_custom_op(
        2,
        Box::new(NoopCustomOp {
            name: "my-custom-op".into(),
        }),
    );

    assert_eq!(config.operations.len(), 7);
    let report = run_import_repair(&mut brep, &config);
    assert_eq!(report.steps.len(), 7);

    // The third step should be our custom op.
    assert!(
        matches!(&report.steps[2].0, StepLabel::Custom(name) if name == "my-custom-op"),
        "expected custom op at position 2, got {:?}",
        report.steps[2].0,
    );
}

#[test]
fn custom_op_push_appends_to_end() {
    let tol = Tolerance::default();
    let mut config = ImportRepairConfig::default_config(tol);
    let original_len = config.operations.len();

    config.push_custom_op(Box::new(NoopCustomOp {
        name: "tail-op".into(),
    }));

    assert_eq!(config.operations.len(), original_len + 1);
    assert!(matches!(
        config.operations.last(),
        Some(RepairOperation::Custom(_))
    ));
}

#[test]
fn standard_op_count_excludes_custom() {
    let tol = Tolerance::default();
    let mut config = ImportRepairConfig::default_config(tol);
    let original_count = config.standard_op_count();

    config.push_custom_op(Box::new(NoopCustomOp {
        name: "extra".into(),
    }));

    assert_eq!(config.standard_op_count(), original_count);
}

// --- Format-specific profile tests ---

#[test]
fn step_ap203_profile_has_non_empty_operations() {
    let tol = Tolerance::default();
    let config = ImportRepairConfig::for_format(ImportFormat::StepAP203, tol);
    assert!(
        !config.operations.is_empty(),
        "AP203 config should have operations"
    );
    assert!(config.standard_op_count() >= 2);
}

#[test]
fn step_ap214_profile_has_non_empty_operations() {
    let tol = Tolerance::default();
    let config = ImportRepairConfig::for_format(ImportFormat::StepAP214, tol);
    assert!(
        !config.operations.is_empty(),
        "AP214 config should have operations"
    );
    assert!(config.standard_op_count() >= 3);
}

#[test]
fn iges_profile_has_non_empty_operations() {
    let tol = Tolerance::default();
    let config = ImportRepairConfig::for_format(ImportFormat::Iges, tol);
    assert!(
        !config.operations.is_empty(),
        "IGES config should have operations"
    );
    assert!(config.standard_op_count() >= 3);
}

#[test]
fn iges_has_larger_displacement_budget() {
    let tol = Tolerance::default();
    let default_config = ImportRepairConfig::default_config(tol);
    let iges_config = ImportRepairConfig::for_format(ImportFormat::Iges, tol);

    assert!(
        iges_config.max_displacement > default_config.max_displacement,
        "IGES should have larger displacement budget"
    );
}

#[test]
fn format_profiles_run_on_valid_cube() {
    let tol = Tolerance::default();
    for format in [
        ImportFormat::StepAP203,
        ImportFormat::StepAP214,
        ImportFormat::Iges,
    ] {
        let mut brep = BRep::unit_cube(tol).expect("unit cube creation failed");
        let config = ImportRepairConfig::for_format(format, tol);
        let report = run_import_repair(&mut brep, &config);
        assert!(
            report.total_defects_after <= report.total_defects_before,
            "{format:?} repair should not increase defects"
        );
    }
}
