//! PK-143 — blend identity survives a STEP round-trip.
//!
//! Identify is curvature/tangency based (no metadata channel yet), so
//! a fillet exported to STEP and re-imported through the full
//! parse → resolve → map → heal pipeline must still classify: same
//! kind, same radius, live support pair. This is the
//! provenance-through-STEP gate: downstream editing (modify /
//! unfillet) starts from identify, so identification loss across
//! serialization silently strands imported models.

use std::sync::Arc;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_fillet::{
    BlendKind, FilletOptions, FilletRadius, FilletSpec, fillet_edges, identify_blend,
};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn step_document(brep: &BRep, name: &str) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: name.into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
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

fn roundtrip(brep: &BRep, name: &str) -> BRep {
    let bytes = export_step(&step_document(brep, name), &StepExportOptions::default())
        .expect("STEP export");
    let imported = import_step(&bytes, &StepImportOptions::default()).expect("STEP import");
    imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("round-trip returns a solid")
        .clone()
}

fn cylinder_faces(brep: &BRep) -> Vec<FaceHandle> {
    brep.faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(h, _)| h)
        .collect()
}

#[test]
fn cylinder_fillet_identity_survives_step_roundtrip() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let edge = cube.edges.keys().next().expect("cube edge");
    let filleted = fillet_edges(
        &cube,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        }],
        &FilletOptions::default(),
        tol(),
        None,
    )
    .expect("baseline fillet")
    .brep;

    let imported = roundtrip(&filleted, "PK143-CYL-FILLET");
    imported
        .validate()
        .expect("round-tripped body must validate");

    let cyls = cylinder_faces(&imported);
    assert_eq!(cyls.len(), 1, "one fillet band expected after round-trip");
    let info = identify_blend(&imported, cyls[0])
        .expect("identify must recognize the fillet after STEP round-trip");
    assert_eq!(info.kind, BlendKind::Fillet);
    let radius = info.radius.expect("round-tripped fillet has a radius");
    assert!(
        (radius - 0.1).abs() < 1e-6,
        "identified radius must survive serialization: {radius}"
    );
    let (sa, sb) = info
        .support_faces
        .expect("support pair must resolve on the imported body");
    for support in [sa, sb] {
        assert!(
            imported.faces.contains_key(support),
            "support face {support:?} must be live in the imported body"
        );
    }
}

#[test]
fn chamfer_identity_survives_step_roundtrip() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let edge = cube.edges.keys().next().expect("cube edge");
    let chamfered = popkern_fillet::chamfer_edges(
        &cube,
        &[popkern_fillet::ChamferSpec {
            edge,
            kind: popkern_fillet::ChamferKind::Equal(0.1),
        }],
        &popkern_fillet::ChamferOptions::default(),
        tol(),
        None,
    )
    .expect("baseline chamfer")
    .brep;

    let imported = roundtrip(&chamfered, "PK143-CHAMFER");
    imported
        .validate()
        .expect("round-tripped chamfer body must validate");

    // The chamfer band is planar; find it as the face identify
    // classifies as Chamfer.
    let mut found = None;
    for (fh, _) in imported.faces.iter() {
        if let Some(info) = identify_blend(&imported, fh)
            && info.kind == BlendKind::Chamfer
        {
            found = Some(info);
            break;
        }
    }
    let info = found.expect("identify must recognize the chamfer after STEP round-trip");
    let half_width = info.radius.expect("chamfer half-width");
    assert!(
        half_width.is_finite() && half_width > 0.0,
        "chamfer width must be a finite positive measurement: {half_width}"
    );
}
