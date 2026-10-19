//! A body touching another face to face, with a small real overlap
//! beside the contact (PK-850): a box sits on a block (opposite
//! coplanar contact at z = 10) and carries a peg dipping 0.2 into the
//! block. The peg's wall meets the contact plane along the box's hole
//! rim, which reads like a contact ending there; the faces beyond that
//! rim decide. The touching shortcut once took it for a pure contact:
//! the intersection served empty and the subtraction served the block
//! untouched, the peg's overlap lost. The pure contact (no peg) still
//! reads touching.

use popkern_rhai::ScriptHost;

fn volumes(op: &str, with_peg: bool) -> Vec<f64> {
    let peg = if with_peg {
        "let peg = m.tool(|t| t.profile(xy_at(9.8)).circle(1.0).at(2.0, 1.0).sketch(\"p\").extrude(blind(1.0)));\n\
         let b = b.union(peg, #{ label: \"bp\" });"
    } else {
        ""
    };
    let src = format!(
        "let m = model();\n\
         let a = m.profile(xy()).rect(10.0, 10.0).sketch(\"a\").extrude(blind(10.0));\n\
         let b = m.tool(|t| t.profile(xy_at(10.0)).rect(10.0, 10.0).sketch(\"b\").extrude(blind(10.0)));\n\
         {peg}\n\
         let x = a.{op}(b, #{{ label: \"x\" }});\nm"
    );
    let mut model = ScriptHost::new()
        .run(&src, "contact.rhai")
        .unwrap_or_else(|e| panic!("{op}: {}", e.render()))
        .model
        .expect("model");
    let Ok(parts) = model.to_parts_all() else {
        return Vec::new();
    };
    parts
        .iter()
        .map(|p| {
            popkern_measure::body_volume(p.brep(), &popkern_core::tolerance::Tolerance::default())
                .expect("volume")
        })
        .collect()
}

#[test]
fn a_peg_beside_a_face_contact_overlaps() {
    let peg = std::f64::consts::PI * 0.2;
    let inter = volumes("intersect", true);
    assert_eq!(inter.len(), 1, "the overlap serves: {inter:?}");
    assert!(
        (inter[0] - peg).abs() < 1e-6,
        "intersection {inter:?} vs {peg}"
    );
    let sub = volumes("subtract", true);
    assert_eq!(sub.len(), 1);
    assert!(
        (sub[0] - (1000.0 - peg)).abs() < 1e-6,
        "subtraction {sub:?}"
    );
    // Control: the pure face contact is touching (empty intersection,
    // the block untouched by subtraction).
    assert!(volumes("intersect", false).is_empty());
    let pure = volumes("subtract", false);
    assert!((pure[0] - 1000.0).abs() < 1e-9, "{pure:?}");
}
