//! PK-639: the file's geometric uncertainty is read (writers type it,
//! `LENGTH_MEASURE(1.E-06)`, which once never parsed, so every import
//! ran at the default) and is the trim-recovery allowance when looser
//! than the kernel default (reported); the imported body keeps the default
//! tolerance with per-edge tolerances where imprecision is measured.
//! Angular uncertainties are skipped and a value in another length
//! unit is converted.

use popkern_io::{StepImportOptions, import_step};

const CUBE: &str = include_str!("fixtures/cube.step");

/// The mapping allowance the import reported (mm), or `None` when it
/// kept the default; the body tolerance is always the default.
fn linear_tolerance(text: &str) -> Option<f64> {
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    solid.validate().expect("valid");
    let default = popkern_core::tolerance::Tolerance::default().linear();
    assert!(
        (solid.tolerance.linear() - default).abs() < 1e-18,
        "the body keeps the default"
    );
    imported.diagnostics.entries().iter().find_map(|e| {
        e.message
            .strip_prefix("trim allowance widened to the file's uncertainty ")?
            .trim_end_matches(" mm")
            .parse::<f64>()
            .ok()
    })
}

fn with_uncertainty(entity: &str) -> String {
    let old = "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #401, 'tolerance');";
    assert!(CUBE.contains(old));
    CUBE.replacen(old, entity, 1)
}

#[test]
fn a_looser_file_uncertainty_is_the_import_tolerance() {
    let text = with_uncertainty(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-03), #401, 'tolerance');",
    );
    assert!((linear_tolerance(&text).expect("loosened") - 1e-3).abs() < 1e-15);
}

#[test]
fn a_tighter_file_uncertainty_keeps_the_kernel_default() {
    let text = with_uncertainty(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-12), #401, 'tolerance');",
    );
    assert_eq!(linear_tolerance(&text), None);
}

#[test]
fn an_angular_uncertainty_is_not_a_length() {
    let text = CUBE
        .replacen(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#402))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#404, #402))",
            1,
        )
        .replacen(
            "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #401, 'tolerance');",
            "#404 = UNCERTAINTY_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.5), #405, 'angle');\n\
             #405 = (NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($, .RADIAN.));\n\
             #402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(2.0E-03), #401, 'tolerance');",
            1,
        );
    assert!((linear_tolerance(&text).expect("loosened") - 2e-3).abs() < 1e-15);
}

#[test]
fn an_uncertainty_in_another_length_unit_is_converted() {
    // 1e-6 metre = 1e-3 mm.
    let text = with_uncertainty(
        "#406 = (LENGTH_UNIT() NAMED_UNIT(#400) SI_UNIT($, .METRE.));\n\
         #402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #406, 'tolerance');",
    );
    assert!((linear_tolerance(&text).expect("loosened") - 1e-3).abs() < 1e-12);
}

/// The cube with edge E0's line moved 1e-5 off both its vertices (and
/// off one of its faces by as much).
fn cube_with_an_edge_off_its_vertices(uncertainty: &str) -> String {
    let text = with_uncertainty(&format!(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty}), #401, 'tolerance');"
    ));
    let text = text.replacen(
        "#40 = LINE('', #1, #30);",
        "#409 = CARTESIAN_POINT('', (0.0, 1.0E-05, 0.0));\n#40 = LINE('', #409, #30);",
        1,
    );
    assert!(text.contains("#409"));
    text
}

#[test]
fn an_edge_off_its_vertices_within_the_file_precision_imports_watertight() {
    let text = cube_with_an_edge_off_its_vertices("1.0E-03");
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    solid
        .validate()
        .expect("valid within the recorded edge tolerance");
    let mesh = popkern_mesh::tessellate(solid, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "corners meet at the vertices"
    );
}

#[test]
fn the_same_gap_beyond_the_file_precision_is_not_absorbed() {
    // Written to 1e-6, a 1e-5 gap is a defect: nothing records it.
    let text = cube_with_an_edge_off_its_vertices("1.0E-06");
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    assert!(
        solid.validate().is_err(),
        "a gap beyond the declared precision must fail validation"
    );
}
