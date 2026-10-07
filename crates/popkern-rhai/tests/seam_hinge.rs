//! PK-702 Rhai surface: `m.hinge_bodies` grows a hinge on a box and its
//! lid; labelled faces of both bodies still resolve after it; the
//! document emits as code (`seam_hinge`, `seam_hinge_component`) that
//! rebuilds the same bodies.

use popkern_rhai::{ScriptHost, to_rhai};
use popkern_script::prelude::*;

fn run(source: &str) -> Model {
    ScriptHost::new()
        .run(source, "seam.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()))
        .model
        .expect("model")
}

fn geometry(m: &mut Model) -> Vec<Vec<[f64; 3]>> {
    let mut all: Vec<Vec<[f64; 3]>> = m
        .bodies()
        .into_iter()
        .map(|b| {
            let mut v: Vec<[f64; 3]> = m
                .brep_of(b)
                .expect("brep")
                .vertices
                .values()
                .map(|v| [v.position.x, v.position.y, v.position.z])
                .collect();
            v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
            v
        })
        .collect();
    all.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    all
}

const BOX: &str = include_str!("../examples/hinged_box.rhai");

#[test]
fn the_hinged_box_builds_and_keeps_its_face_labels() {
    let source = BOX.replace(
        "let h = m.hinge_bodies(",
        "m.set_active(base);\nm.label_face(\"bottom\", faces().facing_neg(\"z\"));\nm.set_active(lid);\nm.label_face(\"lid_top\", faces().facing(\"z\").at_height(33.0));\nlet h = m.hinge_bodies(",
    );
    let mut m = run(&source);
    // The box and its lid, both knuckled; print-in-place: no pin body.
    assert_eq!(m.bodies().len(), 2);
    // The box's bottom face, labelled before the hinge, still resolves.
    m.entity_by_label("bottom")
        .expect("the box's label survives the hinge");
    m.entity_by_label("lid_top")
        .expect("the lid's label survives the hinge");
    m.assert_watertight().expect("watertight");
}

#[test]
fn a_seam_hinge_emits_as_code_that_rebuilds_the_same_bodies() {
    let source = r#"let m = model();
let base = m.profile(xy()).rect(60.0, 40.0).sketch("base").extrude(blind(20.0));
let lid = m.profile(xy_at(20.0)).rect(60.0, 40.0).sketch("lid").extrude(blind(5.0));
let h = m.hinge_bodies("butt", base, lid, axis(edges().straight().direction("x").at_height(20.0).at_vertex([30.0, 20.0, 20.0])),
    #{ knuckles: 5, pin: 2.0, wall: 1.0, clearance: 0.3, label: "hinge" });
m
"#;
    let mut original = run(source);
    assert_eq!(original.bodies().len(), 3, "box, lid, pin");
    let emitted = to_rhai(&mut original).expect("emit");
    assert!(emitted.contains(".seam_hinge(\"butt\""), "{emitted}");
    assert!(
        emitted.contains("seam_hinge_component(\"hinge\", \"leaf_b\""),
        "{emitted}"
    );
    assert!(
        emitted.contains("m.seam_hinge_component(\"hinge\", \"pin\""),
        "{emitted}"
    );
    let mut again = run(&emitted);
    assert_eq!(geometry(&mut original), geometry(&mut again), "{emitted}");
    assert_eq!(to_rhai(&mut again).expect("emit again"), emitted);
}
