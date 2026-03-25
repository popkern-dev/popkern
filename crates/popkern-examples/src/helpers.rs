//! Shared helpers for popkern showcase examples.

use std::path::Path;
use std::sync::Arc;

use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

/// Export a BRep solid as a STEP file with the given product name.
pub fn export(brep: &BRep, name: &str, filename: &str) {
    let output_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&output_dir).expect("create output dir");
    let path = output_dir.join(filename);

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: format!("{}-001", name.to_uppercase()),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    match export_step_file(&doc, &path, &StepExportOptions::default()) {
        Ok(()) => {
            println!("  Exported: {}", path.display());
            // Validate by re-importing the STEP file.
            validate_step_file(&path);
        }
        Err(e) => eprintln!("  STEP export failed: {e}"),
    }
}

/// Validate a STEP file by re-importing it and checking for errors.
fn validate_step_file(path: &Path) {
    use popkern_io::{StepImportOptions, import_step};

    match std::fs::read(path) {
        Ok(data) => match import_step(&data, &StepImportOptions::default()) {
            Ok(result) => {
                let total_solids: usize = result
                    .document
                    .products
                    .iter()
                    .flat_map(|p| &p.shapes)
                    .map(|s| s.solids.len())
                    .sum();
                println!("  STEP re-import OK: {total_solids} solid(s)");
            }
            Err(e) => eprintln!("  STEP re-import FAILED: {e}"),
        },
        Err(e) => eprintln!("  Could not read STEP file for validation: {e}"),
    }
}

/// Try to render a BRep to PNG via the headless renderer.
/// Returns Ok if successful, or prints a warning and returns Err if GPU is unavailable.
pub fn try_render(brep: &BRep, filename: &str) {
    let output_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&output_dir).expect("create output dir");
    let path = output_dir.join(filename);

    let camera = popkern_viz::camera::Camera::default_iso();
    match popkern_viz::headless::render_to_png_gpu(brep, &camera, 800, 600, &path, true) {
        Ok(()) => println!("  Rendered: {}", path.display()),
        Err(e) => eprintln!("  PNG render skipped (no GPU?): {e}"),
    }
}
