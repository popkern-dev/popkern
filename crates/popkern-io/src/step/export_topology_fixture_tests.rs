use super::*;

use crate::step::import::{StepImportOptions, import_step};

#[test]
fn imported_fixture_cube_exports_shared_manifold_topology() {
    let imported = import_step(
        include_bytes!("../../tests/fixtures/cube.step"),
        &StepImportOptions {
            healing: crate::step::healing::HealingLevel::None,
            ..Default::default()
        },
    )
    .expect("fixture cube imports");
    let brep = imported.document.products[0].shapes[0].solids[0].clone();
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(&brep).unwrap();

    let count =
        |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

    assert_eq!(
        count("VERTEX_POINT"),
        8,
        "fixture cube must keep shared vertices"
    );
    assert_eq!(
        count("EDGE_CURVE"),
        12,
        "fixture cube must keep shared edges"
    );
    assert_eq!(
        count("CLOSED_SHELL"),
        1,
        "fixture cube must export as closed shell"
    );
    assert_eq!(
        count("MANIFOLD_SOLID_BREP"),
        1,
        "fixture cube must export as a manifold solid"
    );
    assert_eq!(
        count("SHELL_BASED_SURFACE_MODEL"),
        0,
        "coherent fixture cube must not need face-local open-shell fallback"
    );
}
