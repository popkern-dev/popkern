//! PK-629: assembly product structure round-trips through STEP. The
//! writer places every child through a context-dependent shape
//! representation; the reader resolves that placement (it returned
//! identity for every occurrence before) and the instance name.

use std::sync::Arc;

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_io::{
    AssemblyChild, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};

fn cube() -> popkern_core::topology::BRep {
    let bytes = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/cube.step"
    ))
    .expect("fixture");
    let imported = import_step(&bytes, &StepImportOptions::default()).expect("import");
    imported.document.products[0].shapes[0].solids[0].clone()
}

fn product(
    name: &str,
    children: Vec<AssemblyChild>,
    solids: Vec<popkern_core::topology::BRep>,
) -> StepProductNode {
    StepProductNode {
        name: name.into(),
        id: name.to_uppercase(),
        children,
        shapes: if solids.is_empty() {
            vec![]
        } else {
            vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids,
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })]
        },
        metadata: ProductMetadata::default(),
    }
}

fn placements() -> Vec<Transform3> {
    let turn =
        Transform3::rotation_around_axis(Point3::new(1.0, 2.0, 0.0), Vec3::new(1.0, 1.0, 1.0), 0.7)
            .expect("rotation");
    vec![
        Transform3::identity(),
        Transform3::translation(10.0, -3.0, 0.5),
        Transform3::translation(0.0, 4.0, 2.0).compose(&turn),
    ]
}

#[test]
fn an_assembly_round_trips_its_occurrences_placements_and_names() {
    let children: Vec<AssemblyChild> = placements()
        .into_iter()
        .enumerate()
        .map(|(i, transform)| AssemblyChild {
            product_index: 2,
            transform,
            instance_name: Some(format!("block_{i}")),
        })
        .collect();
    // root -> sub (one occurrence, translated) -> block (three).
    let doc = StepDocument {
        products: vec![
            product(
                "root",
                vec![AssemblyChild {
                    product_index: 1,
                    transform: Transform3::translation(0.0, 0.0, 5.0),
                    instance_name: Some("stage".into()),
                }],
                vec![],
            ),
            product("sub", children, vec![]),
            product("block", vec![], vec![cube()]),
        ],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default()).expect("export");
    let text = String::from_utf8(bytes.clone()).expect("ascii");
    let count = |e: &str| text.matches(&format!("={e}(")).count();
    assert_eq!(count("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), 4);
    assert_eq!(count("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"), 4);
    assert_eq!(count("ITEM_DEFINED_TRANSFORMATION"), 4);
    assert_eq!(
        count("MANIFOLD_SOLID_BREP"),
        1,
        "the block's solid is written once"
    );

    let back = import_step(&bytes, &StepImportOptions::default()).expect("re-import");
    let products = &back.document.products;
    let by_name = |name: &str| products.iter().position(|p| p.name == name).expect(name);
    let (root, sub, block) = (by_name("root"), by_name("sub"), by_name("block"));
    assert_eq!(back.document.roots, vec![root]);
    let stage = &products[root].children;
    assert_eq!(stage.len(), 1);
    assert_eq!(stage[0].product_index, sub);
    assert_eq!(stage[0].instance_name.as_deref(), Some("stage"));
    assert!(
        stage[0]
            .transform
            .approx_eq(&Transform3::translation(0.0, 0.0, 5.0), 1e-12)
    );
    let blocks = &products[sub].children;
    assert_eq!(blocks.len(), 3);
    for (i, (child, expected)) in blocks.iter().zip(placements()).enumerate() {
        assert_eq!(child.product_index, block);
        assert_eq!(
            child.instance_name.as_deref(),
            Some(format!("block_{i}").as_str())
        );
        assert!(
            child.transform.approx_eq(&expected, 1e-12),
            "block_{i}: {:?} vs {:?}",
            child.transform,
            expected
        );
    }
    assert_eq!(products[block].shapes[0].solids.len(), 1);
}

#[test]
fn a_reader_relating_parent_first_is_inverted() {
    // Same file with rep_1/rep_2 swapped and the IDT items swapped:
    // the operator now maps the parent frame onto the child's, and the
    // reader must still recover the placement.
    let doc = StepDocument {
        products: vec![
            product(
                "root",
                vec![AssemblyChild {
                    product_index: 1,
                    transform: placements()[2],
                    instance_name: Some("b".into()),
                }],
                vec![],
            ),
            product("block", vec![], vec![cube()]),
        ],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let text = String::from_utf8(export_step(&doc, &StepExportOptions::default()).expect("export"))
        .expect("ascii");
    let swapped = swap_relationship(&text);
    assert_ne!(swapped, text);
    let back = import_step(swapped.as_bytes(), &StepImportOptions::default()).expect("import");
    let root = back
        .document
        .products
        .iter()
        .find(|p| p.name == "root")
        .expect("root");
    assert!(
        root.children[0]
            .transform
            .approx_eq(&placements()[2], 1e-12)
    );
}

/// Swap rep_1/rep_2 in the relationship and item_1/item_2 in the
/// transformation (the parent-first spelling some writers use).
fn swap_relationship(text: &str) -> String {
    let mut out = String::new();
    for line in text.split_inclusive(';') {
        let swapped = if let Some(at) = line.find("REPRESENTATION_RELATIONSHIP('','',") {
            swap_last_two_refs(line, at + "REPRESENTATION_RELATIONSHIP('','',".len())
        } else if let Some(at) = line.find("ITEM_DEFINED_TRANSFORMATION('','',") {
            swap_last_two_refs(line, at + "ITEM_DEFINED_TRANSFORMATION('','',".len())
        } else {
            line.to_string()
        };
        out.push_str(&swapped);
    }
    out
}

fn swap_last_two_refs(line: &str, start: usize) -> String {
    let rest = &line[start..];
    let close = rest.find(')').expect("close");
    let args: Vec<&str> = rest[..close].split(',').map(str::trim).collect();
    assert_eq!(args.len(), 2, "{line}");
    format!(
        "{}{},{}{}",
        &line[..start],
        args[1],
        args[0],
        &rest[close..]
    )
}
