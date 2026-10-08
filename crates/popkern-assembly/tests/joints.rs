//! Revolute joints (six-axis-arm §4.3): one verb for concentric +
//! axial seat + driven angle; followers at a ratio; transactional
//! refusals; `.pkn` round trip.

use popkern_assembly::{
    Assembly, AssemblyError, Datum, Document, JointDrive, MateEntity, MateRef, Parent, Placement,
    RevoluteOptions,
};
use popkern_core::geometry::{Point3, Vec3};
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

/// `base` at the origin; `arm` and `rotor` placed off-axis and turned.
fn three() -> Assembly {
    let mut asm = Assembly::new("joints");
    let d = asm.add_part("block", block()).expect("part");
    asm.place(Parent::Root, d, Some("base"), Placement::identity())
        .expect("base");
    let turned = |x: f64, deg: f64| {
        let r = Placement::rotation(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), deg)
            .expect("rotation");
        Placement::translation(x, 0.3, 0.5).then(&r)
    };
    asm.place(Parent::Root, d, Some("arm"), turned(2.0, 17.0))
        .expect("arm");
    asm.place(Parent::Root, d, Some("rotor"), turned(-2.0, -40.0))
        .expect("rotor");
    asm
}

fn z_axis(instance: &str) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum: Datum::Z },
    }
}

fn x_axis(instance: &str) -> MateRef {
    MateRef {
        instance: instance.into(),
        entity: MateEntity::Datum { datum: Datum::X },
    }
}

fn world(asm: &Assembly, name: &str) -> Placement {
    asm.leaves()
        .into_iter()
        .find(|r| r.path.to_string() == name)
        .expect("occurrence")
        .world
}

/// The occurrence's turn about +z (degrees) and its origin.
fn pose(asm: &Assembly, name: &str) -> (f64, Point3) {
    let p = world(asm, name);
    let x = p.x_axis();
    (x.y.atan2(x.x).to_degrees(), p.origin())
}

fn close(a: f64, b: f64) -> bool {
    let d = (a - b + 540.0).rem_euclid(360.0) - 180.0;
    d.abs() < 1e-9
}

fn joint_on_base(asm: &mut Assembly, angle: Option<f64>) {
    asm.add_revolute(
        "j",
        z_axis("base"),
        z_axis("arm"),
        &RevoluteOptions {
            zero: Some((x_axis("base"), x_axis("arm"))),
            angle,
            limits: Some([-180.0, 180.0]),
            ..RevoluteOptions::default()
        },
    )
    .expect("revolute");
}

#[test]
fn revolute_seats_the_axis_and_drives_the_angle() {
    let mut asm = three();
    joint_on_base(&mut asm, Some(30.0));
    let (turn, origin) = pose(&asm, "arm");
    assert!(close(turn, 30.0), "arm turned {turn}");
    assert!(
        origin.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-9,
        "seated at {origin:?}"
    );
    asm.set_assembly_param("j", -75.0).expect("slider");
    let (turn, origin) = pose(&asm, "arm");
    assert!(close(turn, -75.0), "arm turned {turn}");
    assert!(origin.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-9);
    // The joint names its slider and owns three mates.
    assert_eq!(asm.joints().len(), 1);
    assert_eq!(asm.joints()[0].drive, JointDrive::Param);
    assert_eq!(asm.mates().len(), 3);
    assert_eq!(asm.params()[0].bounds, Some([-180.0, 180.0]));
}

#[test]
fn revolute_without_an_angle_keeps_the_current_pose() {
    let mut asm = three();
    joint_on_base(&mut asm, None);
    let value = asm.params()[0].display_value();
    let (turn, _) = pose(&asm, "arm");
    assert!(
        close(value, 17.0),
        "the slider starts at the current 17°: {value}"
    );
    assert!(close(turn, 17.0), "nothing turned: {turn}");
}

#[test]
fn follower_turns_at_its_ratio_from_its_phase() {
    let mut asm = three();
    joint_on_base(&mut asm, Some(0.0));
    let report = asm
        .add_revolute(
            "rotor_follow",
            z_axis("base"),
            z_axis("rotor"),
            &RevoluteOptions {
                zero: Some((x_axis("base"), x_axis("rotor"))),
                follow: Some(("j".into(), -29.0)),
                ..RevoluteOptions::default()
            },
        )
        .expect("follower");
    assert!(report.param.is_none(), "a follower creates no slider");
    // Phase: the rotor's current −40° at j = 0.
    let (turn, _) = pose(&asm, "rotor");
    assert!(close(turn, -40.0), "rotor kept its phase: {turn}");
    asm.set_assembly_param("j", 10.0).expect("slider");
    let (turn, _) = pose(&asm, "rotor");
    assert!(close(turn, -40.0 - 290.0), "rotor at −40 − 29·10: {turn}");
    let (arm, _) = pose(&asm, "arm");
    assert!(close(arm, 10.0), "arm at 10: {arm}");
}

#[test]
fn default_zero_is_perpendicular_to_each_axis() {
    let mut asm = three();
    asm.add_revolute(
        "j",
        z_axis("base"),
        z_axis("arm"),
        &RevoluteOptions::default(),
    )
    .expect("revolute");
    let before = pose(&asm, "arm").0;
    let start = asm.params()[0].display_value();
    asm.set_assembly_param("j", start + 25.0).expect("slider");
    let after = pose(&asm, "arm").0;
    assert!(
        close(after - before, 25.0),
        "the slider turns the arm 1:1: {before} -> {after}"
    );
}

#[test]
fn refusals_leave_the_assembly_untouched() {
    let mut asm = three();
    joint_on_base(&mut asm, Some(0.0));
    let mates = asm.mates().len();
    let params = asm.params().len();
    let arm = world(&asm, "arm");
    let refuse = |asm: &mut Assembly, name: &str, a: MateRef, b: MateRef, o: RevoluteOptions| {
        asm.add_revolute(name, a, b, &o).expect_err("must refuse")
    };
    // A plane is no axis.
    let plane = MateRef {
        instance: "rotor".into(),
        entity: MateEntity::Datum { datum: Datum::Xy },
    };
    let e = refuse(
        &mut asm,
        "k",
        z_axis("base"),
        plane,
        RevoluteOptions::default(),
    );
    assert!(
        matches!(e, AssemblyError::JointEntity { role: "axis", .. }),
        "{e}"
    );
    // A zero direction along the axis.
    let e = refuse(
        &mut asm,
        "k",
        z_axis("base"),
        z_axis("rotor"),
        RevoluteOptions {
            zero: Some((x_axis("base"), z_axis("rotor"))),
            ..RevoluteOptions::default()
        },
    );
    assert!(
        matches!(e, AssemblyError::JointEntity { role: "zero", .. }),
        "{e}"
    );
    // Zero ratio, unknown leader, duplicate name.
    for (name, follow) in [
        ("k", Some(("j".to_string(), 0.0))),
        ("k", Some(("nope".into(), 2.0))),
        ("j", None),
    ] {
        let e = refuse(
            &mut asm,
            name,
            z_axis("base"),
            z_axis("rotor"),
            RevoluteOptions {
                follow,
                ..RevoluteOptions::default()
            },
        );
        assert!(matches!(e, AssemblyError::Param { .. }), "{e}");
    }
    // Limits that exclude the starting angle: refused after the first
    // two mates were added — all of them come back out.
    let e = refuse(
        &mut asm,
        "k",
        z_axis("base"),
        z_axis("rotor"),
        RevoluteOptions {
            angle: Some(50.0),
            limits: Some([0.0, 10.0]),
            ..RevoluteOptions::default()
        },
    );
    assert!(matches!(e, AssemblyError::Param { .. }), "{e}");
    assert_eq!(asm.mates().len(), mates);
    assert_eq!(asm.params().len(), params);
    assert_eq!(asm.joints().len(), 1);
    assert_eq!(world(&asm, "arm"), arm);
}

#[test]
fn joints_round_trip_through_pkn() {
    let mut asm = three();
    joint_on_base(&mut asm, Some(12.5));
    let bytes = asm.to_pkn_bytes().expect("save");
    let Document::Assembly(mut back) = Document::from_pkn_bytes(&bytes).expect("open") else {
        panic!("an assembly");
    };
    assert_eq!(back.joints(), asm.joints());
    back.set_assembly_param("j", 40.0)
        .expect("slider after reopen");
    assert!(close(pose(&back, "arm").0, 40.0));
}
