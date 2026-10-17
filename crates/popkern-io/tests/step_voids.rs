//! Solids with cavities through STEP (PK-880).
//!
//! A popkern solid owns its outer shell and any cavity shells. The
//! writer exported one `MANIFOLD_SOLID_BREP` per shell, so a cavity
//! became a second, positive solid (OpenCASCADE read a hollowed dome
//! as two bodies), and the reader ignored `BREP_WITH_VOIDS`, so a
//! hollow part from another system imported as nothing. A solid with
//! cavities now travels as one `BREP_WITH_VOIDS`.

use std::sync::Arc;

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    HealingLevel, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};

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

/// A 3-cube with a unit-cube cavity at (1, 1, 1).
fn hollow_cube() -> BRep {
    let tol = Tolerance::default();
    let outer = BRep::unit_cube(tol)
        .unwrap()
        .transform(&Transform3::scale_uniform(3.0))
        .unwrap();
    let inner = BRep::unit_cube(tol)
        .unwrap()
        .transform(&Transform3::translation(1.0, 1.0, 1.0))
        .unwrap();
    popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        outer,
        inner,
        tol,
        None,
    )
    .expect("hollow cube")
}

#[test]
fn a_hollow_solid_exports_as_one_brep_with_voids() {
    let hollow = hollow_cube();
    assert_eq!((hollow.solids.len(), hollow.shells.len()), (1, 2));
    let bytes = export_step(&document_with_body(hollow), &StepExportOptions::default()).unwrap();
    let text = String::from_utf8_lossy(&bytes);
    let count = |name: &str| text.matches(&format!("={name}(")).count();
    assert_eq!(count("BREP_WITH_VOIDS"), 1, "one solid with a void");
    assert_eq!(count("MANIFOLD_SOLID_BREP"), 0, "the cavity is not a solid");
    assert_eq!(count("ORIENTED_CLOSED_SHELL"), 1);
    assert_eq!(
        count("ORIENTED_FACE"),
        0,
        "void faces are written reversed inline"
    );
    assert_eq!(count("CLOSED_SHELL"), 2, "outer + void shells");
    assert!(text.contains(",.F.);"), "voids are used reversed");
}

#[test]
fn a_hollow_solid_round_trips_with_its_cavity_facing_the_void() {
    let bytes = export_step(
        &document_with_body(hollow_cube()),
        &StepExportOptions::default(),
    )
    .unwrap();
    let imported = import_step(
        &bytes,
        &StepImportOptions {
            healing: HealingLevel::None,
            ..Default::default()
        },
    )
    .unwrap();
    let bodies: Vec<&BRep> = imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| &s.solids)
        .collect();
    assert_eq!(bodies.len(), 1, "one body");
    let body = bodies[0];
    assert_eq!(
        (body.solids.len(), body.shells.len(), body.faces.len()),
        (1, 2, 12)
    );
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let volume = welded_mesh_volume(body);
    assert!((volume - 26.0).abs() < 1e-9, "volume {volume}");
}

/// The served mesh's signed volume (exact for a planar body): the
/// cavity counts negative only if its faces face the void.
fn welded_mesh_volume(body: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v = 0.0;
    for fm in &mesh.face_meshes {
        for t in &fm.triangles {
            let [a, b, c] = t.map(|i| fm.vertices[i].position);
            v += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x);
        }
    }
    v / 6.0
}
