//! PK-629 substrate: definitions built once, named rigid instances,
//! nesting, flattening, the `.pkn` v5 round trip, and STEP product
//! structure checked by re-import and (where installed) OpenCASCADE.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_assembly::{Assembly, AssemblyError, Document, InstancePath, Parent, Placement};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn plate(w: f64, h: f64, t: f64) -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(w, h)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(&mut m, End::Blind(t.into()))
        .expect("extrude");
    m
}

fn turned() -> Placement {
    Placement::rotation(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), 30.0)
        .expect("rotation")
        .then(&Placement::translation(5.0, 0.0, 0.0))
}

/// root: plate_1, plate_2 (turned), stack_a, stack_b (raised)
/// stack: plate (low), plate (high)
fn sample() -> Assembly {
    let mut asm = Assembly::new("sample");
    let plate_def = asm.add_part("plate", plate(2.0, 1.0, 0.25)).expect("part");
    asm.place(Parent::Root, plate_def, None, Placement::identity())
        .expect("p1");
    asm.place(Parent::Root, plate_def, None, turned())
        .expect("p2");
    let stack = asm.add_subassembly("stack").expect("sub");
    asm.place(
        Parent::Definition(stack),
        plate_def,
        Some("low"),
        Placement::identity(),
    )
    .expect("low");
    asm.place(
        Parent::Definition(stack),
        plate_def,
        Some("high"),
        Placement::translation(0.0, 0.0, 1.0),
    )
    .expect("high");
    asm.place(
        Parent::Root,
        stack,
        Some("stack_a"),
        Placement::translation(0.0, 5.0, 0.0),
    )
    .expect("a");
    asm.place(
        Parent::Root,
        stack,
        Some("stack_b"),
        Placement::translation(0.0, -5.0, 2.0),
    )
    .expect("b");
    asm
}

fn bounds(brep: &BRep) -> [[f64; 3]; 2] {
    let mut b = [[f64::MAX; 3], [f64::MIN; 3]];
    for v in brep.vertices.values() {
        for (k, x) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            b[0][k] = b[0][k].min(x);
            b[1][k] = b[1][k].max(x);
        }
    }
    b
}

#[test]
fn instances_share_one_definition_and_flatten_to_world_bodies() {
    let mut asm = sample();
    assert_eq!(asm.definitions().len(), 2);
    let rows = asm.tree();
    let paths: Vec<String> = rows.iter().map(|r| r.path.to_string()).collect();
    assert_eq!(
        paths,
        [
            "plate_1",
            "plate_2",
            "stack_a",
            "stack_a/low",
            "stack_a/high",
            "stack_b",
            "stack_b/low",
            "stack_b/high"
        ]
    );
    let flat = asm.flatten().expect("flatten");
    assert_eq!(flat.len(), 6, "six plate occurrences, one body each");
    let local = bounds(&flat[0].brep);
    let high_b = flat
        .iter()
        .find(|f| f.path.to_string() == "stack_b/high")
        .expect("high");
    let b = bounds(&high_b.brep);
    for k in 0..3 {
        let offset = [0.0, -5.0, 3.0][k];
        assert!((b[0][k] - (local[0][k] + offset)).abs() < 1e-12, "{k}");
    }
    // The turned instance is the definition's body under the turn.
    let turned_body = &flat[1].brep;
    for (_, v) in &flat[0].brep.vertices {
        let want = turned().apply_point(&v.position);
        assert!(
            turned_body
                .vertices
                .values()
                .any(|w| (w.position - want).length() < 1e-12),
            "{want:?}"
        );
    }
    // Hiding a sub-assembly hides what it holds.
    let a = asm
        .resolve(&InstancePath::parse("stack_a"))
        .expect("resolve");
    asm.set_visible(a, false).expect("hide");
    let hidden: Vec<String> = asm
        .tree()
        .into_iter()
        .filter(|r| !r.visible)
        .map(|r| r.path.to_string())
        .collect();
    assert_eq!(hidden, ["stack_a", "stack_a/low", "stack_a/high"]);
}

#[test]
fn names_definitions_and_cycles_are_checked() {
    let mut asm = sample();
    let plate_def = asm.definition_named("plate").expect("plate");
    let stack = asm.definition_named("stack").expect("stack");
    assert!(matches!(
        asm.place(
            Parent::Root,
            plate_def,
            Some("plate_1"),
            Placement::identity()
        ),
        Err(AssemblyError::DuplicateInstanceName { .. })
    ));
    assert!(matches!(
        asm.place(Parent::Root, plate_def, Some("a/b"), Placement::identity()),
        Err(AssemblyError::InvalidName { .. })
    ));
    assert!(matches!(
        asm.add_part("plate", plate(1.0, 1.0, 1.0)),
        Err(AssemblyError::DuplicateDefinitionName { .. })
    ));
    assert!(matches!(
        asm.add_part("nothing", Model::new()),
        Err(AssemblyError::EmptyDefinition { .. })
    ));
    // stack inside stack, and stack inside outer inside stack.
    match asm.place(
        Parent::Definition(stack),
        stack,
        None,
        Placement::identity(),
    ) {
        Err(AssemblyError::Cycle { chain, .. }) => assert_eq!(chain, "stack > stack"),
        other => panic!("{other:?}"),
    }
    let outer = asm.add_subassembly("outer").expect("outer");
    asm.place(
        Parent::Definition(outer),
        stack,
        None,
        Placement::identity(),
    )
    .expect("stack in outer");
    match asm.place(
        Parent::Definition(stack),
        outer,
        None,
        Placement::identity(),
    ) {
        Err(AssemblyError::Cycle { chain, .. }) => assert_eq!(chain, "stack > outer > stack"),
        other => panic!("{other:?}"),
    }
    // Placing into a part is refused.
    assert!(matches!(
        asm.place(
            Parent::Definition(plate_def),
            plate_def,
            None,
            Placement::identity()
        ),
        Err(AssemblyError::Cycle { .. } | AssemblyError::WrongDefinitionKind { .. })
    ));
    // Auto names skip taken names.
    let third = asm
        .place(Parent::Root, plate_def, None, Placement::identity())
        .expect("auto");
    assert_eq!(asm.instance(third).expect("instance").name, "plate_3");
}

#[test]
fn the_document_round_trips_through_pkn_version_5() {
    let mut asm = sample();
    let a = asm
        .resolve(&InstancePath::parse("stack_a"))
        .expect("resolve");
    asm.set_visible(a, false).expect("hide");
    let bytes = asm.to_pkn_bytes().expect("save");
    assert!(
        String::from_utf8_lossy(&bytes)
            .starts_with("{\"version\":5,\"naming_scheme\":3,\"assembly\"")
    );
    assert!(matches!(
        popkern_history::FeatureTree::load(&bytes),
        Err(popkern_history::HistoryError::AssemblyDocument)
    ));
    let Document::Assembly(mut back) = Document::from_pkn_bytes(&bytes).expect("open") else {
        panic!("an assembly document");
    };
    assert_eq!(
        back.tree(),
        asm.tree(),
        "paths, placements (bitwise), visibility"
    );
    let (before, after) = (asm.flatten().expect("a"), back.flatten().expect("b"));
    assert_eq!(before.len(), after.len());
    for (x, y) in before.iter().zip(&after) {
        assert_eq!(bounds(&x.brep), bounds(&y.brep), "{}", x.path);
    }
    assert_eq!(back.to_pkn_bytes().expect("again"), bytes, "stable bytes");

    // A part document opens through the same door.
    let part = plate(1.0, 1.0, 1.0).to_pkn_bytes().expect("part");
    assert!(matches!(
        Document::from_pkn_bytes(&part),
        Ok(Document::Part(_))
    ));
}

#[test]
fn step_export_writes_the_product_structure() {
    let mut asm = sample();
    let bytes = asm.to_step_bytes("sample").expect("step");
    let text = String::from_utf8(bytes.clone()).expect("ascii");
    let count = |e: &str| text.matches(&format!("={e}(")).count();
    assert_eq!(count("PRODUCT"), 3, "root + plate + stack");
    assert_eq!(
        count("NEXT_ASSEMBLY_USAGE_OCCURRENCE"),
        6,
        "4 at the root + 2 in stack"
    );
    assert_eq!(count("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"), 6);
    assert_eq!(count("MANIFOLD_SOLID_BREP"), 1, "the plate is written once");

    let back =
        popkern_io::import_step(&bytes, &popkern_io::StepImportOptions::default()).expect("import");
    let products = &back.document.products;
    let root = &products[back.document.roots[0]];
    assert_eq!(root.name, "sample");
    let names: Vec<_> = root
        .children
        .iter()
        .map(|c| c.instance_name.clone().unwrap_or_default())
        .collect();
    assert_eq!(names, ["plate_1", "plate_2", "stack_a", "stack_b"]);
    assert!(
        root.children[1]
            .transform
            .approx_eq(turned().transform(), 1e-12)
    );

    // An independent reader places every occurrence where we do.
    let flat = asm.flatten().expect("flatten");
    if let Some(report) = occ_assembly(&text) {
        let instances = report["instances"].as_array().expect("instances");
        assert_eq!(instances.len(), flat.len());
        assert_eq!(report["parts"], 1);
        for (row, body) in instances.iter().zip(&flat) {
            assert_eq!(row["path"], body.path.to_string());
            assert_eq!(row["part"], "plate");
            let ours = bounds(&body.brep);
            for (k, corner) in ours.iter().enumerate() {
                for (axis, value) in corner.iter().enumerate() {
                    let theirs = row["bbox"][k][axis].as_f64().expect("bbox");
                    assert!(
                        (theirs - value).abs() < 1e-7,
                        "{}: {theirs} vs {value}",
                        body.path
                    );
                }
            }
        }
    }
}

/// Run the OCC XCAF oracle on `text` when the OCP venv is installed.
pub fn occ_assembly(text: &str) -> Option<serde_json::Value> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let python =
        PathBuf::from(std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let script = root.join("tools/occ-oracle/occ_assembly.py");
    if !python.exists() || !script.exists() {
        eprintln!("[pk629] OCC oracle absent; checked by re-import only");
        return None;
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk629");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{:x}.step", fnv(text)));
    std::fs::write(&path, text).expect("write");
    let output = Command::new(python)
        .arg(script)
        .arg(&path)
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&output.stdout).expect("oracle json"))
}

fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf29ce484222325u64, |h, b| {
        (h ^ b as u64).wrapping_mul(0x100000001b3)
    })
}
