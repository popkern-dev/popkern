//! Revolute joints in Rhai (six-axis-arm §4.3): `a.revolute(...)`
//! makes a slider (or follows one at a ratio), and show-as-code emits
//! the joint as one verb — a fixed point that still moves.

use popkern_assembly::{Assembly, Document, JointDrive};
use popkern_rhai::{ScriptHost, assembly_to_rhai};

fn run(source: &str) -> Assembly {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "joints.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()));
    match outcome.take_document() {
        Some(Document::Assembly(a)) => *a,
        _ => panic!("the script returns an assembly"),
    }
}

const ARM: &str = r#"let a = assembly("arm");
let block = a.part("block", |m| { m.profile(xy()).rect(4.0, 1.0).sketch("s").extrude(blind(1.0)); });
let base = a.place(block, "base");
let link = a.place(block, "link").translate(1.0, 2.0, 0.0);
let rotor = a.place(block, "rotor").translate(-3.0, 0.0, 0.0);
a.ground(base);
let j = a.revolute("j1", base.axis([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]), link.datum("z"),
    #{ zero: [base.datum("x"), link.datum("x")], angle: 20.0, limits: [-170.0, 170.0], offset: 1.0 });
let f = a.revolute("j1_rotor", base.datum("z"), rotor.datum("z"),
    #{ zero: [base.datum("x"), rotor.datum("x")], follow: "j1", ratio: -21.0 });
print(j.param);
print(f.param == ());
a
"#;

fn turn(asm: &Assembly, name: &str) -> (f64, [f64; 3]) {
    let row = asm
        .leaves()
        .into_iter()
        .find(|r| r.path.to_string() == name)
        .expect("occurrence");
    let x = row.world.x_axis();
    let o = row.world.origin();
    (x.y.atan2(x.x).to_degrees(), [o.x, o.y, o.z])
}

fn same_angle(a: f64, b: f64) -> bool {
    ((a - b + 540.0).rem_euclid(360.0) - 180.0).abs() < 1e-7
}

#[test]
fn revolute_makes_a_slider_and_a_follower() {
    let mut asm = run(ARM);
    let (t, o) = turn(&asm, "link");
    assert!(same_angle(t, 20.0), "link at 20: {t}");
    // The seat: base axis origin z = 1, axial offset 1.
    assert!(
        (o[2] - 2.0).abs() < 1e-9 && o[0].abs() < 1e-9 && o[1].abs() < 1e-9,
        "seated: {o:?}"
    );
    let rotor_start = turn(&asm, "rotor").0;
    asm.set_assembly_param("j1", 30.0).expect("slider");
    assert!(same_angle(turn(&asm, "link").0, 30.0));
    let rotor = turn(&asm, "rotor").0;
    assert!(
        same_angle(rotor - rotor_start, -210.0),
        "rotor −21 × 10°: {rotor_start} -> {rotor}"
    );
    assert!(asm.set_assembly_param("j1", 171.0).is_err(), "limits hold");
    assert_eq!(asm.joints().len(), 2);
    assert!(matches!(asm.joints()[1].drive, JointDrive::Follow { ratio, .. } if ratio == -21.0));
}

#[test]
fn show_as_code_emits_joints_as_one_verb_and_is_a_fixed_point() {
    let mut asm = run(ARM);
    asm.set_assembly_param("j1", -45.0).expect("pose");
    let code = assembly_to_rhai(&mut asm).expect("emit");
    assert_eq!(code.matches("a.revolute(").count(), 2, "{code}");
    assert!(
        !code.contains("a.concentric("),
        "joint mates are not emitted raw:\n{code}"
    );
    assert!(
        !code.contains("a.param(\"j1\""),
        "the joint makes its own slider:\n{code}"
    );
    let mut again = run(&code);
    for name in ["base", "link", "rotor"] {
        let (a, b) = (turn(&asm, name), turn(&again, name));
        assert!(same_angle(a.0, b.0), "{name}: {} vs {}", a.0, b.0);
        for k in 0..3 {
            assert!((a.1[k] - b.1[k]).abs() < 1e-9, "{name} origin");
        }
    }
    assert_eq!(
        assembly_to_rhai(&mut again).expect("emit again"),
        code,
        "fixed point"
    );
    // The emitted joints still move.
    again.set_assembly_param("j1", 60.0).expect("slider");
    assert!(same_angle(turn(&again, "link").0, 60.0));
}

#[test]
fn revolute_refusals_are_typed() {
    let mut host = ScriptHost::new();
    let err = host
        .run(
            r#"let a = assembly("x");
let block = a.part("block", |m| { m.profile(xy()).rect(1.0, 1.0).sketch("s").extrude(blind(1.0)); });
let base = a.place(block, "base");
let link = a.place(block, "link");
a.revolute("j", base.datum("xy"), link.datum("z"));
a"#,
            "bad.rhai",
        )
        .expect_err("a plane is no axis");
    assert!(
        err.payload.to_string().contains("joint_entity"),
        "{}",
        err.payload
    );
}
