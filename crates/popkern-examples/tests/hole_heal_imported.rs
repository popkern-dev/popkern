//! PK-195 acceptance 299 (unblocked by PK-199): filled holes are
//! watertight and OCC-clean on an **imported** damaged model.
//!
//! The damaged model arrives through the STEP import path: an open
//! shell (cube with a face removed) exported as
//! `SHELL_BASED_SURFACE_MODEL`/`OPEN_SHELL` and re-imported with
//! healing disabled, so the hole and its null-sentinel boundary
//! edges are the *importer's* product, not a hand-built fixture.
//! `heal_holes` then detects, fills, and commits the cap through the
//! Euler fan operator, and the closed body survives the serving
//! gates, the welded-mesh watertight count, and the OCC oracle.

use std::sync::Arc;

use popkern_core::healing::repair::heal_holes;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    HealingLevel, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};
use popkern_part::Part;

#[test]
fn imported_damaged_model_heals_watertight_and_occ_clean() {
    let tolerance = Tolerance::default();

    // Damage travels through STEP: export the open shell, re-import
    // it with healing off.
    let damaged = BRep::open_box(tolerance).expect("open box fixture");
    let document = StepDocument {
        products: vec![StepProductNode {
            name: "damaged".into(),
            id: "damaged".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![damaged],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let exported = export_step(&document, &StepExportOptions::default()).expect("export");
    let imported = import_step(
        &exported,
        &StepImportOptions {
            healing: HealingLevel::None,
            ..Default::default()
        },
    )
    .expect("import");
    assert!(
        !imported.diagnostics.has_errors(),
        "import diagnostics: {:?}",
        imported.diagnostics.entries()
    );
    let mut body = imported
        .document
        .products
        .iter()
        .flat_map(|product| &product.shapes)
        .flat_map(|shape| &shape.solids)
        .next()
        .expect("one imported body")
        .clone();
    assert_eq!(
        body.boundary_edges.len(),
        4,
        "the hole must arrive through the import"
    );

    let report = heal_holes(&mut body, tolerance.linear());
    assert_eq!(report.holes_detected, 1);
    assert_eq!(report.holes_committed, 1);
    assert!(body.boundary_edges.is_empty());
    body.validate_euler_poincare()
        .expect("strict genus-0 Euler");

    let part = Part::from_brep(body, tolerance, "healed import").expect("serving gates");
    assert_eq!(
        part.open_edge_count().expect("welded mesh"),
        0,
        "healed imported body is watertight"
    );

    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("fixtures");
    std::fs::create_dir_all(&out_dir).expect("fixture dir");
    let path = out_dir.join("healed_imported_open_box.step");
    part.export_step(&path).expect("STEP export");

    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[hole-heal-import-occ] OCP oracle unavailable — leg skipped");
        return;
    }
    let output = std::process::Command::new(&oracle_python)
        .arg(oracle_script)
        .arg(&path)
        .output()
        .expect("oracle invocation");
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("\"passed\": true"),
        "OCC rejected the healed imported STEP: {}",
        report.chars().take(400).collect::<String>()
    );
}
