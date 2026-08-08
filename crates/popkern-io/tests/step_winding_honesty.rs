//! PK-245: winding-flag honesty at the STEP import serving boundary.
//!
//! Import maps `same_sense` from the file into the stored orientation
//! flag, so a dishonest exporter (or a mapping defect) would ship
//! dishonest flags into every downstream operator. Import is a
//! serving boundary for UNTRUSTED external data: the flag must hold
//! the interior-left invariant (outer loop CCW in UV ⇔ flag true) on
//! the imported body regardless of what the file claimed.

use std::sync::Arc;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    HealingLevel, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};

fn dishonest_report(brep: &BRep) -> Vec<String> {
    brep.audit_winding_honesty()
        .into_iter()
        .filter(|v| v.honest == Some(false))
        .map(|v| {
            let surface = brep
                .faces
                .get(v.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                .unwrap_or_else(|| "<no-surface>".into());
            format!(
                "{:?} [{surface}] area {:.3e} flag {}",
                v.face, v.signed_uv_area, v.orientation
            )
        })
        .collect()
}

fn assert_honest(brep: &BRep, fixture: &str) {
    let report = dishonest_report(brep);
    assert!(
        report.is_empty(),
        "{fixture}: {} dishonest face(s):\n  {}",
        report.len(),
        report.join("\n  ")
    );
}

fn document_with_body(brep: BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "body".into(),
            id: "body".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

fn import_options() -> StepImportOptions {
    StepImportOptions {
        healing: HealingLevel::None,
        ..Default::default()
    }
}

fn first_body(result: &popkern_io::StepImportResult) -> &BRep {
    result
        .document
        .products
        .iter()
        .flat_map(|product| &product.shapes)
        .flat_map(|shape| &shape.solids)
        .next()
        .expect("imported body")
}

/// Honest body roundtrip: mapping must not INTRODUCE dishonesty.
#[test]
fn honest_cube_roundtrip_imports_honest() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    assert_honest(&brep, "unit_cube pre-export");
    let exported = export_step(&document_with_body(brep), &StepExportOptions::default()).unwrap();
    let imported = import_step(&exported, &import_options()).unwrap();
    assert!(
        !imported.diagnostics.has_errors(),
        "import diagnostics: {:?}",
        imported.diagnostics.entries()
    );
    assert_honest(first_body(&imported), "unit_cube roundtrip");
}

/// Dishonest-exporter simulation: flip one stored flag before export
/// so the file carries a `same_sense` that contradicts the loop
/// winding. The imported body must still serve honest flags — import
/// may not trust the file's claim over the measured invariant.
#[test]
fn dishonest_flag_body_imports_honest() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let victim = brep.faces.keys().next().unwrap();
    {
        let sof = brep.faces[victim]
            .surface
            .as_mut()
            .expect("cube face has surface");
        sof.orientation = !sof.orientation;
    }
    assert_eq!(
        dishonest_report(&brep).len(),
        1,
        "fixture precondition: exactly the flipped face is dishonest"
    );
    let exported = export_step(&document_with_body(brep), &StepExportOptions::default()).unwrap();
    let imported = import_step(&exported, &import_options()).unwrap();
    assert!(
        !imported.diagnostics.has_errors(),
        "import diagnostics: {:?}",
        imported.diagnostics.entries()
    );
    assert_honest(first_body(&imported), "dishonest-flag roundtrip");
}
