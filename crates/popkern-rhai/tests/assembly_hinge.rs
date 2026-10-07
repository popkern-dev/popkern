//! PK-701 Rhai surface: `a.hinge`, `a.limits`, `a.swing`, and angle /
//! distance mates driven by `p("assembly parameter")`.

use popkern_assembly::{Assembly, Document};
use popkern_rhai::ScriptHost;

fn run(source: &str) -> (Assembly, Vec<String>) {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "asm.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()));
    let printed = outcome.printed.clone();
    match outcome.take_document() {
        Some(Document::Assembly(a)) => (*a, printed),
        _ => panic!("the script returns an assembly"),
    }
}

#[test]
fn an_assembly_hinge_opens_by_its_angle_parameter_and_swings() {
    let (mut asm, printed) = run(r#"let a = assembly("box");
let h = a.hinge("lid", "print_in_place", 60.0, #{ angle: 110.0 });
print(h.angle);
print(h.pin == ());
a.limits("lid_angle", 0.0, 300.0);
let s = a.swing("lid_angle", #{ from: 180.0, to: 300.0, step: 10.0 });
print(s.first_contact);
a
"#);
    assert_eq!(printed[0], "lid_angle");
    assert_eq!(printed[1], "true");
    let contact: f64 = printed[2].parse().expect("a contact angle");
    assert!((contact - 273.86).abs() < 0.5, "{contact}");
    let p = asm
        .params()
        .iter()
        .find(|p| p.name == "lid_angle")
        .expect("param");
    assert!(
        (p.display_value() - 110.0).abs() < 1e-9,
        "the swing restores the start"
    );
    assert_eq!(p.bounds, Some([0.0, 300.0]));
    assert!(asm.set_assembly_param("lid_angle", 301.0).is_err());
    assert!(
        asm.set_assembly_param("lid_angle", 45.0)
            .expect("45")
            .is_empty()
    );
}

#[test]
fn mates_take_an_assembly_parameter() {
    let (mut asm, _) = run(r#"let a = assembly("slide");
a.param("gap", 0.5, "m");
a.param("turn", 30.0, "deg");
let block = a.part("block", |m| { m.profile(xy()).rect(1.0, 1.0).sketch("s").extrude(blind(1.0)); });
let base = a.place(block, "base");
let top = a.place(block, "top");
a.ground(base);
a.distance(top.datum("xy"), base.datum("xy"), p("gap"), #{ aligned: true });
a.concentric(top.datum("z"), base.datum("z"));
a.angle(base.datum("x"), top.datum("x"), p("turn"), #{ about: base.datum("z") });
a
"#);
    assert!(asm.drives_mates("gap") && asm.drives_mates("turn"));
    let top =
        |asm: &Assembly| asm.instances(popkern_assembly::Parent::Root).expect("root")[1].placement;
    let z = |asm: &Assembly| {
        top(asm)
            .apply_point(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0))
            .z
    };
    let z0 = z(&asm);
    asm.set_assembly_param("gap", 1.5).expect("gap");
    assert!(
        ((z(&asm) - z0).abs() - 1.0).abs() < 1e-9,
        "the gap moved the top block by 1"
    );
    let x =
        |asm: &Assembly| top(asm).apply_point(&popkern_core::geometry::Point3::new(1.0, 0.0, 0.0));
    asm.set_assembly_param("turn", 90.0).expect("turn");
    let tip = x(&asm);
    let origin = top(&asm).apply_point(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
    assert!(
        (tip.x - origin.x).abs() < 1e-9 && ((tip.y - origin.y).abs() - 1.0).abs() < 1e-9,
        "{tip:?}"
    );
}
