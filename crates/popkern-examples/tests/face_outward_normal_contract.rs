//! PK-772: the documented meaning of `SurfaceOnFace::orientation`
//! matches the bodies. `orientation` says the outer loop winds
//! counter-clockwise about the CHART normal; the out-of-material normal
//! is `folded_normal × component_material_sense().outward_multiplier()`.
//! Checked face by face against the tessellation's winding (the mesher
//! orients the whole mesh outward) on a primitive, a STEP-imported box
//! and a boolean result.

use std::sync::Arc;

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

fn reimported(brep: &BRep) -> BRep {
    let document = StepDocument {
        products: vec![StepProductNode {
            name: "box".into(),
            id: "box".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&document, &StepExportOptions::default()).expect("exports");
    let imported = import_step(&bytes, &StepImportOptions::default()).expect("imports");
    imported.document.products[0].shapes[0].solids[0].clone()
}

/// Every face's out-of-material normal agrees with its mesh winding.
fn assert_outward_contract(name: &str, brep: &BRep) {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellates");
    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let p = |i: usize| mesh.vertices[t[i]].position;
            let (a, b, c) = (p(0), p(1), p(2));
            a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x)
        })
        .sum();
    assert!(volume > 0.0, "{name}: the mesh is outward");
    let seed = brep.faces.keys().next().expect("a face");
    let outward = brep
        .component_material_sense(seed)
        .expect("material sense")
        .outward_multiplier();
    for fm in &mesh.face_meshes {
        let surface = brep.faces[fm.face].surface.as_ref().expect("surface");
        let mut agree = 0.0;
        for t in &fm.triangles {
            let [a, b, c] = t.map(|i| &fm.vertices[i]);
            let facet = (b.position - a.position).cross(&(c.position - a.position));
            let uv = (
                (a.uv.0 + b.uv.0 + c.uv.0) / 3.0,
                (a.uv.1 + b.uv.1 + c.uv.1) / 3.0,
            );
            if let Some(n) = surface.folded_normal(uv.0, uv.1) {
                agree += facet.dot(&(n * outward));
            }
        }
        assert!(
            agree > 0.0,
            "{name}: face {:?} ({}) out-of-material normal opposes its mesh",
            fm.face,
            surface.geometry.kind_name()
        );
    }
}

#[test]
fn outward_normals_agree_with_the_mesh_on_primitive_imported_and_boolean_bodies() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    assert_outward_contract("unit cube", &cube);
    assert_outward_contract("imported box", &reimported(&cube));
    let bore = BRep::cylinder(0.2, 3.0, tol)
        .expect("cylinder")
        .transform(&Transform3::translation(0.5, 0.5, -1.0))
        .expect("moved");
    let bored = popkern_bool::subtraction(cube, bore, tol, None).expect("bored cube");
    assert_outward_contract("bored cube", &bored);
}
