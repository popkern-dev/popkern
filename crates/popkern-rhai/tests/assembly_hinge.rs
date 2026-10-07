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

/// PK-711: show-as-code keeps the ground, every mate (driven ones with
/// their parameter, scale and offset) and parameter limits, so the
/// emitted assembly's sliders still move it; emission is a fixed point.
#[test]
fn assembly_show_as_code_keeps_mates_drives_and_limits() {
    for source in [
        r#"let a = assembly("box");
let h = a.hinge("lid", "print_in_place", 60.0, #{ angle: 110.0 });
a.limits("lid_angle", 0.0, 300.0);
a
"#,
        r#"let a = assembly("slide");
a.param("gap", 0.5, "m");
a.param("turn", 30.0, "deg");
let block = a.part("block", |m| { m.profile(xy()).rect(1.0, 1.0).sketch("s").extrude(blind(1.0)); });
let base = a.place(block, "base");
let top = a.place(block, "top");
a.ground(base);
a.distance(top.datum("xy"), base.datum("xy"), p("gap"), #{ aligned: true });
a.concentric(top.datum("z"), base.datum("z"));
a.angle(base.datum("x"), top.datum("x"), p("turn"), #{ about: base.datum("z"), scale: 2.0, offset: 10.0 });
a
"#,
    ] {
        let (mut original, _) = run(source);
        let emitted = popkern_rhai::assembly_to_rhai(&mut original).expect("emit");
        let (mut again, _) = run(&emitted);
        assert_eq!(original.mates(), again.mates(), "{emitted}");
        assert_eq!(original.params(), again.params(), "{emitted}");
        assert_eq!(original.ground(), again.ground(), "{emitted}");
        let world = |asm: &Assembly| -> Vec<String> {
            asm.instances(popkern_assembly::Parent::Root)
                .expect("root")
                .iter()
                .map(|i| format!("{:?}", i.placement))
                .collect()
        };
        assert_eq!(world(&original), world(&again), "{emitted}");
        assert_eq!(
            popkern_rhai::assembly_to_rhai(&mut again).expect("emit again"),
            emitted
        );
        // The emitted assembly's first driven parameter still moves it.
        let driven = again
            .mates()
            .iter()
            .find_map(|m| m.drive.as_ref().map(|d| d.param.clone()))
            .expect("a driven mate");
        let before = world(&again);
        let p = again
            .params()
            .iter()
            .find(|p| p.name == driven)
            .expect("p")
            .display_value();
        again
            .set_assembly_param(&driven, p * 0.5 + 1.0)
            .expect("slider");
        assert_ne!(world(&again), before, "{driven} moves the emitted assembly");
    }
}

/// PK-712: a box and its lid as separate parts, hinged along their seam
/// in the assembly — the lid opens by `case_angle` without any part
/// replay, the swing finds where it meets the box (the cross-section
/// value from popkern-hinge's tests/seam.rs: 221.02°), and show-as-code
/// keeps the literal-axis mates and the drive.
#[test]
fn a_box_and_lid_hinged_in_an_assembly_open_by_a_slider() {
    let source = r#"let a = assembly("case");
let boxp = a.part("box", |m| { m.profile(xy()).rect(80.0, 50.0).at(40.0, 25.0).sketch("b").extrude(blind(30.0)); });
let lidp = a.part("lid", |m| { m.profile(xy_at(30.0)).rect(80.0, 50.0).at(40.0, 25.0).sketch("l").extrude(blind(10.0)); });
let bx = a.place(boxp, "box");
let ld = a.place(lidp, "lid");
a.ground(bx);
let h = a.hinge_bodies("case", "print_in_place", bx, ld, axis_at([8.0, 50.0, 30.0], [1.0, 0.0, 0.0]), #{ length: 64.0, pin: 3.0, wall: 1.6 });
print(h.angle);
let s = a.swing("case_angle", #{ from: 1.0, to: 300.0, step: 5.0 });
print(s.first_contact);
a
"#;
    let (mut asm, printed) = run(source);
    assert_eq!(printed[0], "case_angle");
    let contact: f64 = printed[1].parse().expect("a contact angle");
    assert!((contact - 221.02).abs() < 0.5, "first contact {contact}");
    // Open to 90°: no part replays, the lid's front edge rises.
    let lid_front = |asm: &Assembly| {
        let i = &asm.instances(popkern_assembly::Parent::Root).expect("root")[1];
        i.placement
            .apply_point(&popkern_core::geometry::Point3::new(40.0, 0.0, 40.0))
    };
    let closed = lid_front(&asm);
    assert!(
        asm.set_assembly_param("case_angle", 90.0)
            .expect("open")
            .is_empty()
    );
    let open = lid_front(&asm);
    assert!(open.z > closed.z + 40.0, "{closed:?} -> {open:?}");
    let emitted = popkern_rhai::assembly_to_rhai(&mut asm).expect("emit");
    let (again, _) = run(&emitted);
    assert_eq!(asm.mates(), again.mates(), "{emitted}");
}
