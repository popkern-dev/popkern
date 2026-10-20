//! Released joints within their limits (PK-808): a released joint
//! turns with the other mates and reads back its angle; a solve that
//! would take it past its limits is refused with nothing moved; driven
//! again it resumes from the angle it has.

use popkern_assembly::{
    Assembly, AssemblyError, Datum, Mate, MateEntity, MateKind, MateRef, Parent, Placement,
    RevoluteOptions,
};
use popkern_core::geometry::{PortableMath, Vec3};
use popkern_script::Plane;
use popkern_script::prelude::*;

fn block() -> Model {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .rect(1.0, 0.5)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.2.into()))
        .expect("extrude");
    m
}

fn datum(instance: &str, datum: Datum) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum },
    }
}

fn value(asm: &Assembly, name: &str) -> f64 {
    asm.params()
        .iter()
        .find(|p| p.name == name)
        .expect("param")
        .display_value()
}

fn world(asm: &Assembly, name: &str) -> Placement {
    asm.leaves()
        .into_iter()
        .find(|r| r.path.to_string() == name)
        .expect("occurrence")
        .world
}

/// A base with an arm on joint `j` (limits -90..90) and a rotor on
/// joint `t`, both about the base's z axis; `j` released, the arm held
/// parallel to the rotor.
fn following() -> Assembly {
    let mut asm = Assembly::new("release");
    let d = asm.add_part("block", block()).expect("part");
    for name in ["base", "arm", "rotor"] {
        asm.place(Parent::Root, d, Some(name), Placement::identity())
            .expect("place");
    }
    // The first root, `base`, is the ground.
    for (name, limits, offset) in [("j", Some([-90.0, 90.0]), 0.5), ("t", None, 1.0)] {
        let target = if name == "j" { "arm" } else { "rotor" };
        asm.add_revolute(
            name,
            datum("base", Datum::Z),
            datum(target, Datum::Z),
            &RevoluteOptions {
                zero: Some((datum("base", Datum::X), datum(target, Datum::X))),
                angle: Some(if name == "j" { 0.0 } else { 45.0 }),
                limits,
                offset,
                ..RevoluteOptions::default()
            },
        )
        .expect("revolute");
    }
    asm.release_joint("j").expect("release");
    asm.add_mate(Mate {
        kind: MateKind::Angle {
            degrees: 0.0,
            about: Some(datum("base", Datum::Z)),
        },
        a: datum("arm", Datum::X),
        b: datum("rotor", Datum::X),
        drive: None,
    })
    .expect("the arm follows the rotor");
    asm
}

#[test]
fn a_released_joint_follows_the_solve_and_reads_back_its_angle() {
    let mut asm = following();
    assert_eq!(asm.released_joints(), ["j"]);
    assert!(
        (value(&asm, "j") - 45.0).abs() < 1e-6,
        "{}",
        value(&asm, "j")
    );
    asm.set_assembly_param("t", 60.0)
        .expect("within the limits");
    assert!(
        (value(&asm, "j") - 60.0).abs() < 1e-6,
        "{}",
        value(&asm, "j")
    );
    // Its parameter reads back the solve; it cannot be set.
    assert!(asm.set_assembly_param("j", 10.0).is_err());
}

#[test]
fn a_solve_past_a_released_joints_limits_is_refused_and_nothing_moves() {
    let mut asm = following();
    asm.set_assembly_param("t", 60.0).expect("within");
    let (arm, rotor) = (world(&asm, "arm"), world(&asm, "rotor"));
    match asm.set_assembly_param("t", 120.0) {
        Err(AssemblyError::JointLimit {
            joint,
            angle,
            limits,
        }) => {
            assert_eq!(joint, "j");
            assert!((angle - 120.0).abs() < 1e-6, "{angle}");
            assert_eq!(limits, [-90.0, 90.0]);
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(world(&asm, "arm"), arm);
    assert_eq!(world(&asm, "rotor"), rotor);
    assert!((value(&asm, "t") - 60.0).abs() < 1e-9);
    assert!((value(&asm, "j") - 60.0).abs() < 1e-6);
}

#[test]
fn driven_again_a_joint_resumes_from_the_angle_it_has() {
    let mut asm = following();
    asm.set_assembly_param("t", 30.0).expect("within");
    asm.drive_joint("j").expect("drive");
    assert!(asm.released_joints().is_empty());
    assert!((value(&asm, "j") - 30.0).abs() < 1e-6);
    // Now held at 30 by its own mate, it cannot follow the rotor.
    assert!(asm.set_assembly_param("t", 50.0).is_err());
    assert!((value(&asm, "t") - 30.0).abs() < 1e-9);
}

#[test]
fn a_released_joints_follower_turns_with_its_solved_angle() {
    // A reducer stage on the released joint: `cage` follows `j` at half
    // its angle (from the pose it has when added). When the solve turns
    // `j`, the cage must turn by half of what `j` reads back, not stay at
    // the angle its parameter had before.
    let mut asm = following();
    let d = asm.add_part("cage_block", block()).expect("part");
    asm.place(Parent::Root, d, Some("cage"), Placement::identity())
        .expect("place");
    asm.add_revolute(
        "j_cage",
        datum("base", Datum::Z),
        datum("cage", Datum::Z),
        &RevoluteOptions {
            zero: Some((datum("base", Datum::X), datum("cage", Datum::X))),
            follow: Some(("j".into(), 0.5)),
            offset: 1.5,
            ..RevoluteOptions::default()
        },
    )
    .expect("follower");
    let mut seen = Vec::new();
    for t in [60.0, -30.0] {
        asm.set_assembly_param("t", t).expect("within the limits");
        let x = world(&asm, "cage")
            .transform()
            .apply_direction(&Vec3::new(1.0, 0.0, 0.0));
        seen.push((x.y.patan2(x.x).to_degrees(), value(&asm, "j")));
    }
    let (turned, j) = (seen[0].0 - seen[1].0, seen[0].1 - seen[1].1);
    assert!((j - 90.0).abs() < 1e-6, "j turned {j}°");
    assert!(
        (turned - 0.5 * j).abs() < 1e-6,
        "the cage turned {turned}° as j turned {j}°"
    );
}
